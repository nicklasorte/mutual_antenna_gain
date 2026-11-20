%REPORT_FIGURES_MACD_REV2 Figures and tables for the MACD walkthrough working paper.
%   Runs the MACD steps (macd_rev1) with steps 100 times finer than the C#
%   code, checks them with a Monte Carlo of up to 1e10 random pointings per
%   case (macd_monte_carlo_rev1), and compares both with the C# settings.
%   Writes to out_dir:
%
%     fig_macd_sampling_rev2.png   Step 1: fine pointings and the C# pointings
%     fig_macd_hist_rev2.png       Step 2: gain distributions, 0.01 dB bins
%     fig_macd_modes_rev2.png      Step 3: VCP 212 cuts averaged by time
%     fig_macd_convgrid_rev2.png   Step 4: every pair of bins (shown at
%                                  0.25 dB) and the MAG distribution
%     fig_macd_lookup_rev2.png     Step 5: reading the 99% and 99.999% values
%     fig_macd_mc_rev2.png         Monte Carlo against the convolution
%     fig_macd_mc_conv_rev2.png    Monte Carlo 99.999% value against draws
%     fig_macd_csharp_rev2.png     C# settings against the fine result
%     macd_pct_rev2.csv            percentiles: case, method, settings
%     macd_mc_convergence_rev2.csv Monte Carlo percentiles at each draw count
%
%   Cases: "one cut" (WSR-88D at 0.5 deg only) and "VCP 212" (VCP 212
%   SAILSx3 AVSET 6.4, cuts weighted by volume time), horizontal path,
%   Section 1.1 patterns of the methodology report.
%
%   Set macd2_overrides (struct: out_dir, walk_step_div, walk_bin_dB,
%   ref_step_div, ref_bin_dB, mc_draws, mc_chunk, parallel_flag) for a
%   quick run.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also macd_rev1, macd_monte_carlo_rev1, report_figures_macd_rev1.

if ~exist("macd2_overrides", "var"), macd2_overrides = struct(); end
here = fileparts(which("report_figures_macd_rev2"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "walk_step_div", 100, "walk_bin_dB", 0.01, ...
    "ref_step_div", 400, "ref_bin_dB", 0.001, "mc_draws", 1e10, "mc_chunk", 1e7, "parallel_flag", 1);
f = fieldnames(macd2_overrides);
for i = 1:numel(f)
    ov.(f{i}) = macd2_overrides.(f{i});
end
if ~isfolder(ov.out_dir), mkdir(ov.out_dir); end
if ov.parallel_flag == 1 && isempty(gcp("nocreate"))
    parpool("Processes");
end

nex_gain = @(az, tilt) 45.5 + max(-12*(acosd(min(cosd(tilt).*cosd(az), 1))/0.925).^2, ...
    -29*(acosd(min(cosd(tilt).*cosd(az), 1)) <= 10) - 40*(acosd(min(cosd(tilt).*cosd(az), 1)) > 10));
asr_gain = @(az, ~) max(34 - 12*(abs(az)/1.45).^2, 9.5*(abs(az) <= 3.55) + 3.5*(abs(az) > 3.55));
nex = struct("name", "WSR-88D", "g_max_dBi", 45.5, "bw_deg", 0.925, "gain_fn", nex_gain, ...
    "backlobe_dBi", 5.5, "modes", table(0.5, 1, 'VariableNames', ["tilt_deg", "prob"]));
asr = struct("name", "ASR-11", "g_max_dBi", 34, "bw_deg", 1.45, "gain_fn", asr_gain, "backlobe_dBi", 3.5);
v = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
nexv = nex;
nexv.modes = table(v.cuts.el_deg, v.cuts.dur_s, 'VariableNames', ["tilt_deg", "prob"]);
cases = {"one cut", nex; "VCP 212", nexv};
pcts = [50 90 99 99.9 99.99 99.999 99.9999]';
mcp = [99 99.9 99.99 99.999 99.9999]';

R = {}; CV = {};
res = struct();
for c = 1:2
    a1 = cases{c, 2};
    res(c).cs = macd_rev1(a1, asr, pcts=pcts);
    res(c).wk = macd_rev1(a1, asr, pcts=pcts, step_div=ov.walk_step_div, bin_dB=ov.walk_bin_dB);
    res(c).rf = macd_rev1(a1, asr, pcts=pcts, step_div=ov.ref_step_div, bin_dB=ov.ref_bin_dB);
    cps = 10.^(7:log10(ov.mc_draws));
    res(c).mc = macd_monte_carlo_rev1(a1, asr, n_draws=ov.mc_draws, chunk=min(ov.mc_chunk, ov.mc_draws), ...
        checkpoints=cps, pcts=mcp, parallel_flag=ov.parallel_flag, seed=20261120 + c);
    small = macd_monte_carlo_rev1(a1, asr, n_draws=1e6, chunk=1e4, checkpoints=[1e4 1e5 1e6], pcts=mcp, ...
        seed=20261130 + c);
    res(c).conv = [small.pct; res(c).mc.pct];
    names = ["C# (step BW/4, 1 dB)", "walkthrough (step BW/" + ov.walk_step_div + ", " + ov.walk_bin_dB + " dB)", ...
        "reference (step BW/" + ov.ref_step_div + ", " + ov.ref_bin_dB + " dB)"];
    oo = {res(c).cs, res(c).wk, res(c).rf};
    for m = 1:3
        for i = 1:numel(pcts)
            R(end + 1, :) = {cases{c, 1}, names(m), oo{m}.pct.pct(i), oo{m}.pct.interp_dBi(i), oo{m}.pct.upper_edge_dBi(i), NaN, NaN}; %#ok<SAGROW>
        end
    end
    M = res(c).mc.pct(res(c).mc.pct.n_draws == res(c).mc.n_draws, :);
    for i = 1:height(M)
        R(end + 1, :) = {cases{c, 1}, "Monte Carlo (" + sprintf("%.0e", res(c).mc.n_draws) + " draws)", M.pct(i), ...
            M.mag_dBi(i), NaN, M.ci_lo_dBi(i), M.ci_hi_dBi(i)}; %#ok<SAGROW>
    end
    T = res(c).conv;
    T.case = repmat(string(cases{c, 1}), height(T), 1);
    CV{end + 1} = T; %#ok<SAGROW>
end
writetable(cell2table(R, 'VariableNames', ["case", "method", "pct", "mag_dBi", "upper_edge_dBi", "ci_lo_dBi", "ci_hi_dBi"]), ...
    fullfile(ov.out_dir, "macd_pct_rev2.csv"));
writetable(vertcat(CV{:}), fullfile(ov.out_dir, "macd_mc_convergence_rev2.csv"));

% ---- Figures
C = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, "gray", [138 141 145]/255, ...
    "purple", [112 48 160]/255, "grid", [0.9 0.9 0.9]);
fs = {"FontName", "Times New Roman", "FontSize", 8};
wk = res(1).wk; cs = res(1).cs;

% Step 1: fine pointings with the C# pointings
fh = new_fig(6.5, 3.6);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
A2 = {nex, asr};
for a = 1:2
    ax = nexttile(tl);
    s = wk.samples{a}; k = abs(s.rel_az_deg) <= 2.5;
    plot(ax, s.rel_az_deg(k), s.gain_dBi(k), ".", "Color", C.blue, "MarkerSize", 3); hold(ax, "on");
    s2 = cs.samples{a}; k2 = abs(s2.rel_az_deg) <= 2.5;
    plot(ax, s2.rel_az_deg(k2), s2.gain_dBi(k2), "o", "Color", C.dark, "MarkerSize", 4, "LineWidth", 0.8);
    style_ax(ax, C); xlim(ax, [-2.5 2.5]); ylim(ax, [0 50]);
    xlabel(ax, "Pointing from the path (degrees)", fs{:});
    if a == 1
        ylabel(ax, "Gain toward the other radar (dBi)", fs{:});
        text(ax, -2.4, 47, "Dots: walkthrough pointings", fs{:}, "Color", C.blue, "FontWeight", "bold");
        text(ax, -2.4, 44, "Circles: C# pointings", fs{:}, "Color", C.dark, "FontWeight", "bold");
    end
    title(ax, sprintf("%s: a sample every %.4f%s (C#: %.3f%s)", A2{a}.name, A2{a}.bw_deg/ov.walk_step_div, char(176), ...
        A2{a}.bw_deg/4, char(176)), fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_sampling_rev2.png"), "Resolution", 300); close(fh);

% Step 2: gain distributions as share per dB
fh = new_fig(6.5, 3.4);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for a = 1:2
    ax = nexttile(tl);
    h = cs.hist{a};
    bar(ax, h.bin_lo_dBi + 0.5, h.prob, 1, "FaceColor", C.grid, "EdgeColor", "w"); hold(ax, "on");
    [xb, pb] = rebin(wk.hist{a}.bin_lo_dBi, wk.hist{a}.prob, 0.25);
    stairs(ax, xb, pb/0.25, "Color", C.blue, "LineWidth", 1.2);
    ax.YScale = "log"; style_ax(ax, C);
    ylim(ax, [1e-5 1e2]); xlim(ax, [0 48]);
    xlabel(ax, "Gain toward the other radar (dBi)", fs{:});
    if a == 1
        ylabel(ax, "Share of time per dB", fs{:});
        text(ax, 9, 30, sprintf("Blue: walkthrough, %.2f dB bins", ov.walk_bin_dB), fs{:}, "Color", C.blue, "FontWeight", "bold");
        text(ax, 9, 8, "Gray bars: C# 1 dB bins", fs{:}, "Color", C.gray, "FontWeight", "bold");
    end
    title(ax, sprintf("%s: %s pointings", A2{a}.name, num2str(wk.n_samples(a))), fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_hist_rev2.png"), "Resolution", 300); close(fh);

% Step 3: modes
fh = new_fig(6.5, 3.2);
ax = axes(fh);
wv = res(2).wk; H = wv.mode_hist{1}; h = wv.hist{1};
[ue, ia] = unique(nexv.modes.tilt_deg, "stable");
for e = [0.5 0.9 1.3]
    j = ia(ue == e);
    [xb, pb] = rebin(h.bin_lo_dBi, H(:, j), 0.25);
    stairs(ax, xb, pb/0.25, "Color", C.gray, "LineWidth", 0.7); hold(ax, "on");
    kmax = find(H(:, j) > 0, 1, "last");
    text(ax, h.bin_lo_dBi(kmax) + 0.6, 3e-3, sprintf("%.1f%s cut", e, char(176)), fs{:}, "Color", C.gray, "FontWeight", "bold");
end
[xb, pb] = rebin(h.bin_lo_dBi, h.prob, 0.25);
stairs(ax, xb, pb/0.25, "Color", C.blue, "LineWidth", 1.4);
text(ax, 7, 30, "Time-weighted average of all cuts", fs{:}, "Color", C.blue, "FontWeight", "bold");
ax.YScale = "log"; style_ax(ax, C);
ylim(ax, [1e-5 1e2]); xlim(ax, [0 48]);
xlabel(ax, "WSR-88D gain toward the ASR-11 (dBi)", fs{:}); ylabel(ax, "Share of time per dB", fs{:});
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_modes_rev2.png"), "Resolution", 300); close(fh);

% Step 4: pairs of bins (0.25 dB for display) and the MAG distribution
og = macd_rev1(nex, asr, step_div=ov.walk_step_div, bin_dB=0.25);
h1 = og.hist{1}; h2 = og.hist{2};
J = h2.prob*h1.prob.';
fh = new_fig(6.5, 3.4);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
ax = nexttile(tl);
Jl = log10(J); Jl(J == 0) = NaN;
imagesc(ax, h1.bin_lo_dBi + 0.125, h2.bin_lo_dBi + 0.125, Jl, "AlphaData", ~isnan(Jl));
ax.YDir = "normal"; hold(ax, "on");
colormap(ax, flipud(gray(256))*0.85 + 0.15*C.blue);
cb = colorbar(ax); cb.Label.String = "log_{10} share of time"; cb.FontName = "Times New Roman"; cb.FontSize = 7;
for g = [50 60 70]
    x = 5:46; plot(ax, x, g - x, "Color", C.purple, "LineWidth", 0.8);
    text(ax, g - 35.5, 35.5, sprintf("%d dBi", g), fs{:}, "Color", C.purple, "FontWeight", "bold", "HorizontalAlignment", "center");
end
text(ax, 6.5, 5.2, sprintf("Both far sidelobes: %.1f%%", 100*J(1, 1)), fs{:}, "Color", C.dark, "FontWeight", "bold");
style_ax(ax, C); ax.YGrid = "off"; xlim(ax, [4 46]); ylim(ax, [2 37]);
xlabel(ax, "WSR-88D gain bin (dBi)", fs{:}); ylabel(ax, "ASR-11 gain bin (dBi)", fs{:});
title(ax, "(a) Each square: one bin of each antenna", fs{:}, "FontWeight", "normal", "Color", C.dark);
ax = nexttile(tl);
m = wk.mag;
[xb, pb] = rebin(m.bin_lo_dBi, m.prob, 0.25);
stairs(ax, xb, pb/0.25, "Color", C.blue, "LineWidth", 0.8);
ax.YScale = "log"; style_ax(ax, C); ylim(ax, [1e-7 1e2]); xlim(ax, [0 80]);
xlabel(ax, "MAG (dBi): sum of the two gains", fs{:}); ylabel(ax, "Share of time per dB", fs{:});
title(ax, "(b) Add up the squares on each diagonal", fs{:}, "FontWeight", "normal", "Color", C.dark);
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_convgrid_rev2.png"), "Resolution", 300); close(fh);

% Step 5: reading the percentiles (linear scale, percent)
fh = new_fig(6.5, 3.2);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for q = [99 99.999]
    ax = nexttile(tl);
    r = wk.pct(wk.pct.pct == q, :);
    x0 = r.interp_dBi; k = m.bin_lo_dBi >= x0 - 0.15 & m.bin_lo_dBi <= x0 + 0.15;
    stairs(ax, m.bin_lo_dBi(k), 100*m.exceed(k), "Color", C.blue, "LineWidth", 1.2); hold(ax, "on");
    yline(ax, 100 - q, "Color", C.purple, "LineWidth", 0.9);
    plot(ax, r.interp_dBi, 100 - q, "o", "Color", C.dark, "MarkerFaceColor", C.dark, "MarkerSize", 4);
    yy = 100*m.exceed(k); yl = [min(yy) max(yy)];
    text(ax, r.interp_dBi + 0.01, 100 - q + 0.15*diff(yl), sprintf("%g%% value: %.2f dBi", q, r.interp_dBi), fs{:}, ...
        "Color", C.dark, "FontWeight", "bold");
    style_ax(ax, C); xlim(ax, [x0 - 0.15 x0 + 0.15]); ylim(ax, yl + [-0.02 0.02]*diff(yl));
    ax.YAxis.Exponent = 0; ytickformat(ax, "%.5g");
    xlabel(ax, "MAG (dBi)", fs{:});
    if q == 99, ylabel(ax, "Share of time at or above (%)", fs{:}); end
    title(ax, sprintf("Line at %g%%: crosses one %.2f dB bin", 100 - q, ov.walk_bin_dB), fs{:}, ...
        "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_lookup_rev2.png"), "Resolution", 300); close(fh);

% Monte Carlo against the convolution
fh = new_fig(6.5, 3.6);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
ttl = ["WSR-88D at 0.5° only", "VCP 212, cuts weighted by time"];
for c = 1:2
    ax = nexttile(tl);
    e = res(c).mc.exceed; k = e.exceed > 0;
    semilogy(ax, e.mag_dBi(k), e.exceed(k), "Color", C.gray, "LineWidth", 2.4); hold(ax, "on");
    mm = res(c).rf.mag; k = mm.exceed > 1e-10;
    semilogy(ax, mm.bin_lo_dBi(k), mm.exceed(k), "Color", C.blue, "LineWidth", 0.9);
    style_ax(ax, C); ax.YScale = "log"; xlim(ax, [0 85]); ylim(ax, [1e-9 1]);
    xlabel(ax, "MAG (dBi)", fs{:});
    if c == 1
        ylabel(ax, "Share of time at or above", fs{:});
        text(ax, 40, 2e-1, sprintf("Gray: Monte Carlo, %.0e draws", res(c).mc.n_draws), fs{:}, "Color", C.gray, "FontWeight", "bold");
        text(ax, 40, 5e-2, "Blue: convolution", fs{:}, "Color", C.blue, "FontWeight", "bold");
    end
    title(ax, ttl(c), fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_mc_rev2.png"), "Resolution", 300); close(fh);

% Monte Carlo 99.999% value against the number of draws
fh = new_fig(6.5, 3.4);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for c = 1:2
    ax = nexttile(tl);
    T = res(c).conv; T = T(T.pct == 99.999 & T.n_draws >= 1e5, :);
    patch(ax, [T.n_draws; flip(T.n_draws)], [T.ci_lo_dBi; flip(T.ci_hi_dBi)], C.blue, "FaceAlpha", 0.2, "EdgeColor", "none"); hold(ax, "on");
    semilogx(ax, T.n_draws, T.mag_dBi, "o-", "Color", C.blue, "MarkerFaceColor", C.blue, "MarkerSize", 3);
    yref = res(c).rf.pct.interp_dBi(res(c).rf.pct.pct == 99.999);
    ycs = res(c).cs.pct.upper_edge_dBi(res(c).cs.pct.pct == 99.999);
    yline(ax, yref, "Color", C.dark, "LineWidth", 0.9);
    yline(ax, ycs, "--", "Color", C.purple, "LineWidth", 0.9);
    text(ax, 1e8, yref - 1.2, sprintf("Convolution: %.1f dBi", yref), fs{:}, "Color", C.dark, "FontWeight", "bold");
    text(ax, 1e8, ycs + 1.2, sprintf("C# settings: %.0f dBi", ycs), fs{:}, "Color", C.purple, "FontWeight", "bold");
    ax.XScale = "log"; style_ax(ax, C); xlim(ax, [5e4 2*max(T.n_draws)]); ylim(ax, [60 80]);
    xlabel(ax, "Monte Carlo draws", fs{:});
    if c == 1, ylabel(ax, "99.999% MAG (dBi), with 95% interval", fs{:}); end
    title(ax, ttl(c), fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_mc_conv_rev2.png"), "Resolution", 300); close(fh);

% C# settings against the fine result
fh = new_fig(6.5, 3.6);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for c = 1:2
    ax = nexttile(tl);
    mf = res(c).rf.mag; mc = res(c).cs.mag;
    k = mf.exceed > 1e-9; semilogy(ax, mf.bin_lo_dBi(k), mf.exceed(k), "Color", C.gray, "LineWidth", 1.2); hold(ax, "on");
    k = mc.exceed > 1e-9; stairs(ax, mc.bin_lo_dBi(k), mc.exceed(k), "Color", C.blue, "LineWidth", 1.4);
    q = res(c).cs.pct(res(c).cs.pct.pct >= 99 & res(c).cs.pct.pct <= 99.999, :);
    plot(ax, q.upper_edge_dBi, 1 - q.pct/100, "o", "Color", C.dark, "MarkerFaceColor", "w", "MarkerSize", 4);
    for i = 1:height(q)
        text(ax, q.upper_edge_dBi(i) - 1, (1 - q.pct(i)/100)*0.4, sprintf("%g%%: %.0f", q.pct(i), q.upper_edge_dBi(i)), ...
            fs{:}, "Color", C.dark, "HorizontalAlignment", "right", "FontSize", 7);
    end
    style_ax(ax, C); ax.YScale = "log"; xlim(ax, [0 85]); ylim(ax, [1e-7 1]);
    xlabel(ax, "MAG (dBi)", fs{:});
    if c == 1
        ylabel(ax, "Share of time at or above", fs{:});
        text(ax, 50, 3e-2, "C# settings", fs{:}, "Color", C.blue, "FontWeight", "bold");
        text(ax, 50, 1e-2, "Fine reference", fs{:}, "Color", C.gray, "FontWeight", "bold");
    end
    title(ax, ttl(c), fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_csharp_rev2.png"), "Resolution", 300); close(fh);

% -------------------------------------------------------------------------
function fh = new_fig(w, h)
fh = figure("Color", "w", "Units", "inches", "Visible", "off");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(fh, "light");
end
fh.Position(3:4) = [w h];
end

function [xb, pb] = rebin(x, p, w)
% Sum fine bins into display bins of width w (lower edges).
k = floor((x - x(1))/w + 1e-9) + 1;
pb = accumarray(k, p);
xb = x(1) + (0:numel(pb) - 1)'*w;
pb(pb == 0) = 1e-12;                                % so the log axis draws empty bins at the floor
end

function style_ax(ax, C)
box(ax, "off");
ax.FontName = "Times New Roman";
ax.FontSize = 8;
ax.TickDir = "out";
ax.XColor = C.dark;
ax.YColor = C.dark;
ax.YGrid = "on";
ax.GridColor = C.grid;
ax.GridAlpha = 1;
end
