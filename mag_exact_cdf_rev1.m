function out = mag_exact_cdf_rev1(sys1, sys2, opts)
%MAG_EXACT_CDF_REV1 Exact random-mode MAG distribution of two rotating radars.
%   out = mag_exact_cdf_rev1(sys1, sys2, opts) returns the time-weighted
%   mutual antenna gain (MAG) distribution of two mechanically rotating
%   radars in the Rec. ITU-R M.2069-0 Annex 1 random mode, without
%   sampling. It is the exact counterpart of radar_mutual_gain_cdf_rev1
%   with opts.mode = "random" and returns the same output fields, plus the
%   distribution itself:
%       g12 = G_tx,1(toward 2) + G_rx,2(toward 1)
%       g21 = G_tx,2(toward 1) + G_rx,1(toward 2)       (dBi)
%
%   Method: in random mode the two pointing states are independent, so each
%   MAG distribution is the convolution of two single-antenna
%   distributions. Each single-antenna distribution is a time-weighted mix
%   of exact gain slices (gain_slice_rev1), one per elevation step, with
%   weights n_rev*sector_deg/|rate| (the step_dur_s rule of
%   radar_mutual_gain_cdf_rev1). Sector scans (sector_deg < 360) restrict
%   each slice to the sector's azimuths, uniform in time, as the random
%   mode draws them. Gains are rounded up to the next bin edge, so every
%   percentile is conservative by at most 4*bin_dB: per antenna, one bin
%   from rounding up and, for closed-form slices, one bin from taking each
%   cell's larger end gain (red team F1 found the earlier 2*bin_dB
%   statement exceeded by up to 0.05 dB on coarse tables).
%
%   sys1, sys2  radar structs as for radar_mutual_gain_cdf_rev1: .tx
%               (antenna), optional .rx (default .tx), .scan of type
%               "mech_rotating" (rate_deg_s or rate_rpm, optional tilt_deg,
%               n_rev, sector_deg, sector_center_deg, sweep_mode,
%               turnaround_s). Arrays and beam steering are not supported
%               (use the Monte Carlo engine).
%   opts        struct, all fields optional:
%     .el12_deg   elevation of radar 2 seen from radar 1 (takeoff angle at
%                 radar 1, deg; default 0)
%     .el21_deg   elevation of radar 1 seen from radar 2 (default 0)
%     .pct        percentiles (default [50 90 99 99.9 99.99 99.999])
%     .cdf_pct    percent grid of cdf_array12/21 (default 0:1:100)
%     .bin_dB     gain bin width (default 0.05)
%     .az_step_deg  azimuth step for separable and grid slices (0.001)
%
%   out struct (fields of radar_mutual_gain_cdf_rev1, random mode):
%     .g12_dBi, .g21_dBi   [] (no samples)
%     .cdf12, .cdf21       tables: gain_dBi (0.1 dB grid), cdf = P(MAG <=
%                          gain), exceed_pct
%     .pct12, .pct21       tables: pct, gain_dBi, ci95_lo_dBi and
%                          ci95_hi_dBi (equal to gain_dBi: no sampling
%                          error), n_tail (NaN)
%     .cdf_array12/21      [percent, MAG dBi]
%     .dist12, .dist21     structs: gain_dBi (bin values, upper edges),
%                          pmf, exceed_ge = P(MAG >= gain_dBi)
%     .geom, .opts, .method ("exact"), .bias_dB
%     .events12, .ts, .sync  [] (random mode)
%
%   Example:
%     o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, ...
%         "duration_s", 60, "n_trials", 1));
%     r = mag_exact_cdf_rev1(o.sys_nexrad, o.sys_asr11, ...
%         struct("el12_deg", 1.7, "el21_deg", 2.0));
%     disp(r.pct21)
%
%   See also radar_mutual_gain_cdf_rev1, gain_slice_rev1, mag_pair_cdf_rev1.

arguments
    sys1 (1, 1) struct
    sys2 (1, 1) struct
    opts (1, 1) struct = struct()
end
opts = fill_opts(opts);
[sys1, sc1] = check_system(sys1, "sys1");
[sys2, sc2] = check_system(sys2, "sys2");
check_rates(sc1, sc2);
so = struct("bin_dB", opts.bin_dB, "az_step_deg", opts.az_step_deg);
S1t = gain_slice_prep_rev1(sys1.tx, so);
S1r = S1t;
if ~isequal(sys1.rx, sys1.tx), S1r = gain_slice_prep_rev1(sys1.rx, so); end
S2t = gain_slice_prep_rev1(sys2.tx, so);
S2r = S2t;
if ~isequal(sys2.rx, sys2.tx), S2r = gain_slice_prep_rev1(sys2.rx, so); end

% Location-agnostic geometry, as radar_mutual_gain_cdf_rev1: radar 2 at
% azimuth 0 from radar 1, radar 1 at azimuth 180 from radar 2.
[k1t, p1t] = side_dist(S1t, sc1, opts.el12_deg, 0);
[k2r, p2r] = side_dist(S2r, sc2, opts.el21_deg, 180);
[k2t, p2t] = side_dist(S2t, sc2, opts.el21_deg, 180);
[k1r, p1r] = side_dist(S1r, sc1, opts.el12_deg, 0);
d12 = combine(k1t, p1t, k2r, p2r, opts.bin_dB);
d21 = combine(k2t, p2t, k1r, p1r, opts.bin_dB);

out = struct();
out.g12_dBi = [];
out.g21_dBi = [];
out.cdf12 = cdf_table(d12);
out.cdf21 = cdf_table(d21);
out.pct12 = pct_table(d12, opts.pct);
out.pct21 = pct_table(d21, opts.pct);
out.cdf_array12 = [opts.cdf_pct(:), quantile_of(d12, opts.cdf_pct(:))];
out.cdf_array21 = [opts.cdf_pct(:), quantile_of(d21, opts.cdf_pct(:))];
out.dist12 = d12;
out.dist21 = d21;
out.geom = struct("az12_deg", 0, "el12_deg", opts.el12_deg, "az21_deg", 180, "el21_deg", opts.el21_deg);
out.events12 = [];
out.ts = [];
out.sync = [];
out.opts = opts;
out.method = "exact";
out.bias_dB = 4*opts.bin_dB;                    % 2 bins per antenna (see header)
end

% =========================================================================
function [k0, p] = side_dist(S, sc, el, az_to)
% Time-weighted mix of the gain slices of every elevation step.
lim = [];
if sc.sector_deg < 360
    c = az_to - sc.sector_center_deg;           % az_off = az_to - beam azimuth
    lim = [c - sc.sector_deg/2, c + sc.sector_deg/2];
end
n = numel(sc.tilt_u);
hs = cell(n, 1);
for i = 1:n
    if isempty(lim)
        hs{i} = gain_slice_rev1(S, sc.tilt_u(i), el);
    else
        hs{i} = gain_slice_rev1(S, sc.tilt_u(i), el, lim, sc.edge_u(i));
    end
end
k0 = min(cellfun(@(h) h.k_lo, hs));
nb = max(cellfun(@(h) h.k_lo + numel(h.p), hs)) - k0;
p = zeros(nb, 1);
for i = 1:n
    r = hs{i}.k_lo - k0 + (1:numel(hs{i}.p));
    p(r) = p(r) + sc.w_u(i)*hs{i}.p;
end
end

% -------------------------------------------------------------------------
function d = combine(ka, pa, kb, pb, dg)
pmf = max(conv(pa, pb), 0);
pmf = pmf/sum(pmf);
gain = (ka + kb + (0:numel(pmf) - 1)')*dg;
d = struct("gain_dBi", gain, "pmf", pmf, "exceed_ge", cumsum(pmf, "reverse"));
end

% -------------------------------------------------------------------------
function q = quantile_of(d, pct)
% Smallest gain with P(MAG > gain) <= 1 - pct/100; pct 0 gives the
% smallest gain with mass, pct 100 the largest.
E = d.exceed_ge - d.pmf;                         % P(MAG > gain), tail-accurate
q = zeros(size(pct));
first = find(d.pmf > 0, 1);
for j = 1:numel(pct)
    if pct(j) <= 0
        q(j) = d.gain_dBi(first);
    else
        q(j) = d.gain_dBi(find(E <= 1 - pct(j)/100 + 1e-13, 1));
    end
end
end

% -------------------------------------------------------------------------
function T = pct_table(d, pct)
q = quantile_of(d, pct(:));
T = table(pct(:), q, q, q, NaN(numel(pct), 1), 'VariableNames', ...
    ["pct", "gain_dBi", "ci95_lo_dBi", "ci95_hi_dBi", "n_tail"]);
end

% -------------------------------------------------------------------------
function T = cdf_table(d)
% P(MAG <= g) on a 0.1 dB grid, the grid of radar_mutual_gain_cdf_rev1.
nz = find(d.pmf > 0);
g = (floor(d.gain_dBi(nz(1))*10)/10:0.1:ceil(d.gain_dBi(nz(end))*10)/10)';
F = zeros(size(g));
C = cumsum(d.pmf);
for i = 1:numel(g)
    j = find(d.gain_dBi <= g(i) + 1e-9, 1, "last");
    if ~isempty(j), F(i) = C(j); end
end
F = min(F, 1);
T = table(g, F, 100*(1 - F), 'VariableNames', ["gain_dBi", "cdf", "exceed_pct"]);
end

% -------------------------------------------------------------------------
function [s, sc] = check_system(s, label)
if ~isfield(s, "tx") || ~isfield(s, "scan")
    error("mag_exact:missingField", "%s needs fields tx and scan.", label);
end
if ~isfield(s, "rx") || isempty(s.rx)
    s.rx = s.tx;
end
scan = s.scan;
if ~isfield(scan, "type") || string(scan.type) ~= "mech_rotating"
    error("mag_exact:scanType", "%s: only ""mech_rotating"" scans have an exact form here; use radar_mutual_gain_cdf_rev1.", label);
end
if isfield(scan, "beam_mode") && ~isempty(scan.beam_mode) && string(scan.beam_mode) ~= "none"
    error("mag_exact:steering", "%s: beam steering is not supported; use radar_mutual_gain_cdf_rev1.", label);
end
if isfield(scan, "rate_rpm") && ~isempty(scan.rate_rpm)
    rate = 6*scan.rate_rpm(:);
elseif isfield(scan, "rate_deg_s")
    rate = scan.rate_deg_s(:);
else
    error("mag_exact:missingRate", "%s: scan needs rate_deg_s or rate_rpm.", label);
end
if any(rate == 0)
    error("mag_exact:zeroRate", "%s: rotation rate must be nonzero.", label);
end
n = numel(rate);
tilt = field_or(scan, "tilt_deg", 0);
n_rev = field_or(scan, "n_rev", 1);
if isscalar(tilt), tilt = tilt + zeros(n, 1); end
if isscalar(n_rev), n_rev = n_rev + zeros(n, 1); end
sec = field_or(scan, "sector_deg", 360);
if ~isscalar(sec) || sec <= 0 || sec > 360
    error("mag_exact:badSector", "%s: sector_deg must be a scalar in (0, 360].", label);
end
mode = string(field_or(scan, "sweep_mode", "restart"));
turn = field_or(scan, "turnaround_s", 0);
if ~ismember(mode, ["restart", "bidirectional"]) || (mode == "bidirectional" && sec >= 360) || ...
        turn < 0 || (turn > 0 && mode ~= "bidirectional")
    error("mag_exact:badSweep", "%s: sweep_mode must be ""restart"" or ""bidirectional"" (sector below 360), turnaround_s >= 0 only for bidirectional.", label);
end
w = n_rev(:).*(sec./abs(rate) + turn);           % step durations (radar_mutual_gain_cdf_rev1)
w_edge = n_rev(:)*turn;                          % time parked at the edges per step
% uniquetol: tilts from VCP tables can differ by round-off.
[tilt_u, ~, iu] = uniquetol(tilt(:), 1e-6, "DataScale", 1);
w_u = accumarray(iu, w);
sc = struct("rate_deg_s", rate, "tilt_u", tilt_u, "w_u", w_u/sum(w), ...
    "edge_u", accumarray(iu, w_edge)./w_u, ...
    "sector_deg", sec, "sector_center_deg", field_or(scan, "sector_center_deg", 0));
end

% -------------------------------------------------------------------------
function check_rates(sc1, sc2)
% The random-mode validity rules of radar_mutual_gain_cdf_rev1
% (check_random_mode_valid), for single-rate scans.
if ~isscalar(sc1.rate_deg_s) || ~isscalar(sc2.rate_deg_s)
    return
end
r1 = abs(sc1.rate_deg_s);
r2 = abs(sc2.rate_deg_s);
if ismembertol(r1, r2, 1e-9)
    warning("radar_coupling:equalRates", join(["Equal rotation rates keep a fixed phase offset; " ...
        "random mode averages over an unknown offset. For a known offset use " ...
        "radar_mutual_gain_cdf_rev1 with opts.scan_sync = ""synchronized""."], ""));
    return
end
[p, q] = rat(r1/r2, 1e-9*r1/r2);
if p <= 100 && q <= 100
    warning("radar_coupling:commensurateRates", join(["Rotation rates are in a %d:%d ratio; " ...
        "random mode assumes the rates drift apart. Check with radar_mutual_gain_cdf_rev1, " ...
        "opts.mode = ""time_series""."], ""), p, q);
end
end

% -------------------------------------------------------------------------
function v = field_or(s, f, d)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end

% -------------------------------------------------------------------------
function opts = fill_opts(opts)
d = struct("el12_deg", 0, "el21_deg", 0, "pct", [50 90 99 99.9 99.99 99.999], ...
    "cdf_pct", 0:1:100, "bin_dB", 0.05, "az_step_deg", 0.001);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("mag_exact:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
if abs(opts.el12_deg) > 90 || abs(opts.el21_deg) > 90
    error("mag_exact:badElevation", "el12_deg and el21_deg must lie in [-90, 90] deg.");
end
if any(opts.pct < 0 | opts.pct >= 100) || any(opts.cdf_pct < 0 | opts.cdf_pct > 100)
    error("mag_exact:badPct", "pct must lie in [0, 100) and cdf_pct in [0, 100].");
end
end
