function out = simple_mag_nexrad_asr11_rev1(opts)
%SIMPLE_MAG_NEXRAD_ASR11_REV1 Simple mutual antenna gain between a WSR-88D and an ASR-11.
%   out = simple_mag_nexrad_asr11_rev1() computes the mutual antenna gain
%   (MAG) between a WSR-88D (NEXRAD) and an ASR-11 in the simplest useful
%   form, for explaining the method before the full model
%   (nexrad_asr11_mag_rev1):
%
%     MAG = G_WSR-88D(toward the ASR-11) + G_ASR-11(toward the WSR-88D)  (dBi)
%
%   Both antennas rotate at constant rates on a horizontal path. The WSR-88D
%   scans one cut at a fixed tilt; the ASR-11 beam nose is on the horizon,
%   so only azimuth matters for it. Each radar uses one pattern to transmit
%   and receive, so the MAG is the same in both directions.
%
%   Patterns, written out so the slide equation is the code:
%     WSR-88D (ROC specification envelope, as the full model's default):
%       G = 45.5 - 12 (theta/0.925)^2 dBi until it meets the sidelobe
%       envelope: 29 dB below peak (16.5 dBi) out to 10 deg off axis,
%       40 dB below peak (5.5 dBi) beyond. theta is the off-axis angle from
%       the beam axis (at the tilt) to the path.
%     ASR-11 (M.1464-2 Radar C levels, azimuth):
%       G = 34 - 12 (phi/1.45)^2 dBi until it meets the first sidelobe,
%       +9.5 dBi, held out to 3.55 deg; +3.5 dBi (remote sidelobe) beyond.
%       3.55 deg is where the M.1851-2 pattern of the full model reaches the
%       remote level. Compared with that pattern, this one is up to 2.4 dB
%       lower from 1.4 to 2.5 deg and up to 3 dB higher from 2.5 to 3.55 deg.
%
%   The distribution of the MAG takes both azimuths as uniform and
%   independent, which holds when the two rotation rates drift apart (with
%   exactly constant rates the pointing pair repeats; see out.ts). It is
%   computed exactly, without sampling: each pattern's gain distribution is
%   built on a fine azimuth grid (0.0001 deg within 6 deg of the beam,
%   0.01 deg elsewhere), binned on a 0.001 dB grid, and the two are
%   convolved. Percentiles out to 99.9999% are read from it.
%
%   opts (name-value or struct; defaults in parentheses):
%     nexrad_rate      "paper" (3 rpm, 18 deg/s: the 2.7 GHz link budget
%                      working paper, Table 4.4-1) | "vcp212" (21.149 deg/s:
%                      the 0.5 deg surveillance cut of VCP 212,
%                      nexrad_vcp_rev1)
%     nexrad_tilt_deg  WSR-88D beam tilt (0.5, the lowest VCP cut)
%     duration_s       time series length (60)
%     dt_s             time series step (0.001)
%     align_s          time at which both beams point along the path, so
%                      the window shows a main-beam crossing (30); NaN for
%                      random start phases drawn with seed
%     seed             RNG seed for random start phases (1)
%     fig_file         "" (no figure) or a file name for the three-panel
%                      figure (patterns, time series, exceedance curve)
%
%   out fields:
%     params     table of the inputs and their sources
%     pattern    table: az_off_deg, g_nexrad_dBi, g_asr11_dBi
%     peak_dBi   both main beams on the path (G_WSR-88D at the tilt + 34)
%     pct        table: pct, mag_dBi, for 50, 90, 99, 99.9, 99.99,
%                99.999, and 99.9999
%     exceed     table: thr_dBi, time_pct (share of time at or above),
%                for 40, 50, 60, 70 dBi
%     dist       table: mag_dBi (0.001 dB grid), prob, exceed (P(MAG >= x))
%     ts         table: t_s, az_nexrad_deg, az_asr11_deg, mag_dBi
%     repeat_s   period after which the pointing pair repeats at the
%                constant rates (Inf if the rate ratio isn't rational to
%                1e-9)
%
%   Example:
%     out = simple_mag_nexrad_asr11_rev1(fig_file="simple_mag.png");
%     out.pct          % 50%: 9.0 dBi (5.5 + 3.5), 99.9999%: near 76 dBi
%
%   DRAFT//PRE-DECISIONAL
%
%   See also nexrad_asr11_mag_rev1, test_simple_mag_rev1.

arguments
    opts.nexrad_rate (1,1) string {mustBeMember(opts.nexrad_rate, ["paper", "vcp212"])} = "paper"
    opts.nexrad_tilt_deg (1,1) double {mustBeInRange(opts.nexrad_tilt_deg, -1, 20)} = 0.5
    opts.duration_s (1,1) double {mustBePositive} = 60
    opts.dt_s (1,1) double {mustBePositive} = 1e-3
    opts.align_s (1,1) double = 30
    opts.seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.fig_file (1,1) string = ""
end

% ---- Inputs and sources
P.nex_g_max = 45.5;     % dBi, ROC_TECH (WSR-88D technical description), 2850 MHz incl. radome
P.nex_bw = 0.925;       % deg, ROC_TECH
P.nex_sl1 = -29;        % dB below peak, out to 10 deg, ROC_TECH specification
P.nex_sl2 = -40;        % dB below peak, beyond 10 deg, ROC_TECH specification
P.nex_sl_break = 10;    % deg
P.asr_g_max = 34;       % dBi, M.1464-2 Radar C
P.asr_bw = 1.45;        % deg, M.1464-2 Radar C azimuth beamwidth
P.asr_sl1 = 9.5;        % dBi, M.1464-2 Radar C first sidelobe
P.asr_remote = 3.5;     % dBi, M.1464-2 Radar C remote sidelobe
P.asr_remote_from = 3.55;   % deg, where the full model's M.1851-2 pattern reaches the remote level
P.asr_rate = 75;        % deg/s, M.1464-2 Radar C (12.5 rpm)
switch opts.nexrad_rate
    case "paper"
        P.nex_rate = 18;        % deg/s, 3 rpm, 2.7 GHz link budget working paper Table 4.4-1
        rate_src = "WP_REV3 Table 4.4-1 (3 rpm)";
    case "vcp212"
        P.nex_rate = 21.149;    % deg/s, VCP 212 cut 1 (0.5 deg, contiguous surveillance), nexrad_vcp_rev1
        rate_src = "VCP 212 cut 1 (nexrad_vcp_rev1)";
end
params = table(["WSR-88D peak gain"; "WSR-88D beamwidth"; "WSR-88D sidelobes to 10 deg"; ...
    "WSR-88D sidelobes beyond 10 deg"; "WSR-88D tilt"; "WSR-88D rotation rate"; "ASR-11 peak gain"; ...
    "ASR-11 azimuth beamwidth"; "ASR-11 first sidelobe"; "ASR-11 remote sidelobe"; "ASR-11 remote from"; ...
    "ASR-11 rotation rate"; "Path elevation"], ...
    [P.nex_g_max; P.nex_bw; P.nex_sl1; P.nex_sl2; opts.nexrad_tilt_deg; P.nex_rate; P.asr_g_max; P.asr_bw; ...
    P.asr_sl1; P.asr_remote; P.asr_remote_from; P.asr_rate; 0], ...
    ["dBi"; "deg"; "dB"; "dB"; "deg"; "deg/s"; "dBi"; "deg"; "dBi"; "dBi"; "deg"; "deg/s"; "deg"], ...
    ["ROC_TECH"; "ROC_TECH"; "ROC_TECH"; "ROC_TECH"; "lowest VCP cut"; rate_src; ...
    "M.1464-2 Radar C"; "M.1464-2 Radar C"; "M.1464-2 Radar C"; "M.1464-2 Radar C"; "M.1851-2 (full model)"; ...
    "M.1464-2 Radar C"; "horizontal path"], 'VariableNames', ["parameter", "value", "unit", "source"]);

gn = @(az) nexrad_gain(az, opts.nexrad_tilt_deg, P);
ga = @(az) asr_gain(az, P);

% ---- Pattern table for display
az = unique([(-180:0.05:180)'; (-6:0.001:6)']);
out.params = params;
out.pattern = table(az, gn(az), ga(az), 'VariableNames', ["az_off_deg", "g_nexrad_dBi", "g_asr11_dBi"]);
out.peak_dBi = gn(0) + ga(0);

% ---- Exact distribution: fine azimuth grid with cell weights, gains
% binned on a 0.001 dB grid, then the two distributions convolved.
% (A histogram convolution, conv, is exact for independent uniform
% azimuths up to the binning; sampling would need 1e9 draws for 99.9999%.)
dg = 1e-3;
[an, wn] = az_grid();
[aa, wa] = az_grid();
lo = floor(min([gn(an); ga(aa)])/dg)*dg;
hn = accum_hist(gn(an), wn, lo, dg);
ha = accum_hist(ga(aa), wa, lo, dg);
pm = conv(hn, ha);                                  % MAG bin k is 2*lo + (k - 1)*dg
pm = max(pm, 0)/sum(pm);
x = 2*lo + (0:numel(pm) - 1)'*dg;
ex = flip(cumsum(flip(pm)));                        % P(MAG >= x)
keep = pm > 0 | ex > 0;
k1 = find(keep, 1, "last");
out.dist = table(x(1:k1), pm(1:k1), ex(1:k1), 'VariableNames', ["mag_dBi", "prob", "exceed"]);
pcts = [50 90 99 99.9 99.99 99.999 99.9999]';
cdf = cumsum(pm);
mag_p = zeros(size(pcts));
for i = 1:numel(pcts)
    mag_p(i) = x(find(cdf >= pcts(i)/100 - 1e-12, 1));
end
out.pct = table(pcts, mag_p, 'VariableNames', ["pct", "mag_dBi"]);
thr = [40 50 60 70]';
tp = zeros(size(thr));
for i = 1:numel(thr)
    tp(i) = 100*sum(pm(x >= thr(i) - dg/2));
end
out.exceed = table(thr, tp, 'VariableNames', ["thr_dBi", "time_pct"]);

% ---- Time series: both beams rotate at constant rates. Azimuths are
% measured from the direction to the other radar, so 0 is on the path.
t = (0:opts.dt_s:opts.duration_s)';
if isnan(opts.align_s)
    rs = RandStream("mrg32k3a", "Seed", opts.seed);
    ph = 360*rand(rs, 1, 2);
else
    ph = -[P.nex_rate P.asr_rate]*opts.align_s;     % both on the path at align_s
end
az_n = mod(ph(1) + P.nex_rate*t + 180, 360) - 180;
az_a = mod(ph(2) + P.asr_rate*t + 180, 360) - 180;
out.ts = table(t, az_n, az_a, gn(az_n) + ga(az_a), 'VariableNames', ["t_s", "az_nexrad_deg", "az_asr11_deg", "mag_dBi"]);
[nr, dr] = rat(P.asr_rate/P.nex_rate, 1e-9);
out.repeat_s = Inf;
if dr <= 1e4
    out.repeat_s = 360*nr/P.asr_rate;               % the ASR-11 makes nr turns while the WSR-88D makes dr
end

if opts.fig_file ~= ""
    draw(out, opts, P);
end
end

% =========================================================================
function g = nexrad_gain(az, tilt, P)
% Off-axis angle from the beam axis (azimuth az from the path, elevation
% tilt) to the path (azimuth 0, elevation 0): spherical law of cosines.
th = acosd(min(max(cosd(tilt).*cosd(az), -1), 1));
env = P.nex_sl1*(th <= P.nex_sl_break) + P.nex_sl2*(th > P.nex_sl_break);
g = P.nex_g_max + max(-12*(th/P.nex_bw).^2, env);
end

function g = asr_gain(az, P)
a = abs(az);
main = P.asr_g_max - 12*(a/P.asr_bw).^2;
sl = P.asr_sl1*(a <= P.asr_remote_from) + P.asr_remote*(a > P.asr_remote_from);
g = max(main, sl);
end

function [a, w] = az_grid()
% Azimuth cell centers and weights (cell width / 360): fine near the beam.
e = unique([(-180:0.01:-6)'; (-6:1e-4:6)'; (6:0.01:180)']);
a = (e(1:end - 1) + e(2:end))/2;
w = diff(e)/360;
end

function h = accum_hist(g, w, lo, dg)
% Probability in each 0.001 dB bin; bin k covers lo + (k - 1)*dg, rounded.
k = round((g - lo)/dg) + 1;
h = accumarray(k, w);
end

function draw(out, opts, P)
c = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, "gray", [138 141 145]/255, ...
    "purple", [112 48 160]/255, "grid", [0.9 0.9 0.9]);
f = figure("Color", "w", "Units", "inches", "Visible", "off");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(f, "light");
end
f.Position(3:4) = [6.5 7.2];
tl = tiledlayout(f, 3, 1, "TileSpacing", "compact", "Padding", "compact");

% (a) Patterns
ax = nexttile(tl);
p = out.pattern;
k = abs(p.az_off_deg) <= 20;
plot(ax, p.az_off_deg(k), p.g_nexrad_dBi(k), "Color", c.blue, "LineWidth", 1.4); hold(ax, "on");
plot(ax, p.az_off_deg(k), p.g_asr11_dBi(k), "Color", c.gray, "LineWidth", 1.4);
text(ax, 1.2, P.nex_g_max - 2, sprintf("WSR-88D, tilt %.1f%s", opts.nexrad_tilt_deg, char(176)), ...
    "Color", c.blue, "FontWeight", "bold", "FontName", "Times New Roman", "FontSize", 8);
text(ax, 2.6, 26, "ASR-11", "Color", c.gray, "FontWeight", "bold", "FontName", "Times New Roman", "FontSize", 8);
style(ax, c);
xlim(ax, [-20 20]); ylim(ax, [0 50]);
xlabel(ax, "Beam pointing from the path (degrees)");
ylabel(ax, "Gain toward the other radar (dBi)");
title(ax, "(a) Each antenna's gain toward the other", "FontName", "Times New Roman", "FontSize", 9, ...
    "FontWeight", "normal", "Color", c.dark);

% (b) Time series
ax = nexttile(tl);
ts = out.ts;
plot(ax, ts.t_s, ts.mag_dBi, "Color", c.dark, "LineWidth", 0.6); hold(ax, "on");
[mx, im] = max(ts.mag_dBi);
plot(ax, ts.t_s(im), mx, "o", "Color", c.purple, "MarkerFaceColor", c.purple, "MarkerSize", 4);
text(ax, ts.t_s(im) + 0.8, mx - 1, sprintf("Both main beams on the path: %.1f dBi", mx), ...
    "Color", c.purple, "FontWeight", "bold", "FontName", "Times New Roman", "FontSize", 8);
style(ax, c);
xlim(ax, [0 opts.duration_s]); ylim(ax, [0 85]);
xlabel(ax, "Time (s)");
ylabel(ax, "MAG (dBi)");
title(ax, sprintf("(b) Two beams turning: WSR-88D at %.3g and ASR-11 at %g degrees per second", ...
    P.nex_rate, P.asr_rate), "FontName", "Times New Roman", "FontSize", 9, "FontWeight", "normal", "Color", c.dark);

% (c) Exceedance
ax = nexttile(tl);
d = out.dist;
k = d.exceed > 1e-9;
semilogy(ax, d.mag_dBi(k), d.exceed(k), "Color", c.blue, "LineWidth", 1.4); hold(ax, "on");
q = out.pct(out.pct.pct >= 99, :);
plot(ax, q.mag_dBi, 1 - q.pct/100, "o", "Color", c.dark, "MarkerFaceColor", "w", "MarkerSize", 4);
for i = 1:height(q)
    text(ax, q.mag_dBi(i) - 1, (1 - q.pct(i)/100)*0.45, sprintf("%g%%: %.1f", q.pct(i), q.mag_dBi(i)), ...
        "Color", c.dark, "FontName", "Times New Roman", "FontSize", 7, "HorizontalAlignment", "right");
end
style(ax, c);
ax.YScale = "log";
xlim(ax, [0 85]); ylim(ax, [1e-8 1]);
xlabel(ax, "MAG (dBi)");
ylabel(ax, "Share of time at or above");
title(ax, "(c) How often each MAG is reached, with random, independent pointing", ...
    "FontName", "Times New Roman", "FontSize", 9, "FontWeight", "normal", "Color", c.dark);
exportgraphics(f, opts.fig_file, "Resolution", 300);
close(f);
end

function style(ax, c)
box(ax, "off");
ax.FontName = "Times New Roman";
ax.FontSize = 8;
ax.TickDir = "out";
ax.XColor = c.dark;
ax.YColor = c.dark;
ax.YGrid = "on";
ax.GridColor = c.grid;
ax.GridAlpha = 1;
end
