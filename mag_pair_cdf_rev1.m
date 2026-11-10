function R = mag_pair_cdf_rev1(T_tx, T_rx, el_tx_deg, el_rx_deg, opts)
%MAG_PAIR_CDF_REV1 Mutual antenna gain CDF of radar pairs from two gain tables.
%   R = mag_pair_cdf_rev1(T_tx, T_rx, el_tx_deg, el_rx_deg, opts) returns
%   the mutual antenna gain (MAG) distribution of each transmitter-receiver
%   pair, path tx -> rx:
%       MAG = G_tx(toward rx, el_tx) + G_rx(toward tx, el_rx)   (dBi)
%   from two tables made once by gain_dist_table_rev1. No Monte Carlo run
%   is needed per pair: in the M.2069-0 Annex 1 random mode the two
%   pointing states are independent, so the MAG distribution is the
%   convolution of the two single-antenna gain distributions (see
%   gain_dist_table_rev1).
%
%   T_tx, T_rx   tables from gain_dist_table_rev1 (same bin_dB). The
%                transmit pattern of T_tx and the receive pattern of T_rx
%                are used (T.P_rx when present, otherwise T.P_tx).
%   el_tx_deg    elevation of the receiver seen from the transmitter
%                (takeoff angle at the transmitter, deg), one per pair
%   el_rx_deg    elevation of the transmitter seen from the receiver
%                (takeoff angle at the receiver, deg), one per pair. Both
%                map to the "Takeoff Angle (deg)" columns of the 2.7 GHz
%                RFS tool and to el12_deg / el21_deg of
%                radar_mutual_gain_cdf_rev1. Values must lie in each
%                table's elevation grid.
%   opts         struct, all fields optional:
%     .pct          percentiles (default [50 90 99 99.9 99.99 99.999])
%     .brg_tx_deg   for a sector-scanning transmitter: bearing of the
%                   receiver measured clockwise from the center of the
%                   transmitter's sector, deg (scalar or one per pair;
%                   default 0, the path at the sector center, the worst
%                   case for a sector aimed at the other station)
%     .brg_rx_deg   the same for a sector-scanning receiver (default 0).
%                   Full-rotation tables ignore both. Sector slices are
%                   made for each pair with gain_slice_rev1 (exact).
%     .tx_modes     T_tx mode names to include (default all)
%     .rx_modes     T_rx mode names to include (default all)
%     .per_mode     true also returns every mode pair's percentiles
%                   (default false)
%     .return_curve true also returns the envelope exceedance curve of
%                   each pair (default false)
%     .interp       "quantile" (default) interpolates between the two
%                   neighboring grid elevations in quantile space (no false
%                   double peaks); "nearest" uses the nearest grid
%                   elevation; "exact" computes every slice at the pair's
%                   own elevations with gain_slice_rev1 (no interpolation
%                   error; slower for many tilts; tables made after
%                   gain_slice_rev1 only)
%     .parallel_flag 1 runs the pair groups in parfor (default 0)
%     .conv_method  "fft" (default: one FFT per slice column, one inverse
%                   FFT for every tilt pair) | "direct" (conv; same
%                   results, kept for checks)
%     .el_round_deg  round takeoff angles to this step before grouping
%                   (default 0: exact angles). Pairs that share the
%                   heavier side's angle share its slices; rounding to,
%                   for example, 0.01 deg makes more pairs share them
%
%   Scan modes: every pair of included modes (tx mode, rx mode) is
%   convolved, and the result is the worst-mode envelope: at each gain the
%   largest exceedance over the mode pairs. Its p percentile is the largest
%   mode-pair p percentile, which bounds any mix of the modes (the bound of
%   mag_mode_envelope_rev1). For one mode, pass tx_modes and rx_modes.
%
%   R struct:
%     .pct          table, one row per pair: el_tx_deg, el_rx_deg, one
%                   column per percentile (dBi, envelope), max_on_path_dBi
%                   (largest MAG along the path: both antennas' best
%                   pointing at their takeoff angles; between grid
%                   elevations the per-tilt maxima are interpolated
%                   linearly)
%     .worst_mode   string, pairs by percentiles: "tx mode | rx mode" that
%                   sets each envelope percentile
%     .mode_pairs   string, the mode pairs in the envelope
%     .mode_pct     pairs by mode pairs by percentiles (dBi), per_mode only
%     .curve        pairs-by-1 cell of tables (gain_dBi, exceed_pct), when
%                   return_curve; exceed_pct = 100*P(MAG > gain_dBi)
%     .bias_dB      upper bound of the conservative grid bias of every
%                   percentile (dB): one bin per table from rounding gains
%                   up, plus one bin per table from quantile interpolation
%     .notes        method notes
%
%   Example:
%     R = mag_pair_cdf_rev1(T_asr11, T_nexrad, 2.0, 1.7);   % envelope, all VCP modes
%     disp(R.pct); disp(R.worst_mode)
%
%   See also gain_dist_table_rev1, nexrad_mode_scans_rev1,
%   radar_mutual_gain_cdf_rev1, mag_mode_envelope_rev1.

arguments
    T_tx (1, 1) struct
    T_rx (1, 1) struct
    el_tx_deg {mustBeNumeric, mustBeReal, mustBeFinite}
    el_rx_deg {mustBeNumeric, mustBeReal, mustBeFinite}
    opts (1, 1) struct = struct()
end

opts = fill_opts(opts);
el_tx_deg = el_tx_deg(:);
el_rx_deg = el_rx_deg(:);
if numel(el_tx_deg) ~= numel(el_rx_deg)
    error("mag_pair:sizeMismatch", "el_tx_deg and el_rx_deg must have one value per pair.");
end
if abs(T_tx.bin_dB - T_rx.bin_dB) > 1e-12
    error("mag_pair:binMismatch", "T_tx and T_rx must use the same bin_dB (%g and %g).", T_tx.bin_dB, T_rx.bin_dB);
end
n_pair = numel(el_tx_deg);
brg_t = expand_brg(opts.brg_tx_deg, n_pair, "brg_tx_deg");
brg_r = expand_brg(opts.brg_rx_deg, n_pair, "brg_rx_deg");
is_sec_t = is_sector(T_tx);
is_sec_r = is_sector(T_rx);
if opts.interp == "exact" || is_sec_t || is_sec_r
    check_range(el_tx_deg, [-90; 90], "el_tx_deg", T_tx.name);
    check_range(el_rx_deg, [-90; 90], "el_rx_deg", T_rx.name);
else
    check_range(el_tx_deg, T_tx.el_grid_deg, "el_tx_deg", T_tx.name);
    check_range(el_rx_deg, T_rx.el_grid_deg, "el_rx_deg", T_rx.name);
end

% ---- Mode selection and weights. Only tilts with weight are kept.
mt = select_modes(T_tx, opts.tx_modes, "tx_modes");
mr = select_modes(T_rx, opts.rx_modes, "rx_modes");
check_sync(T_tx.mode_rates(mt), T_rx.mode_rates(mr));
kt = find(any(T_tx.W(mt, :) > 0, 1));
kr = find(any(T_rx.W(mr, :) > 0, 1));
sd_t = side(T_tx, kt, mt, false, opts.interp);
sd_r = side(T_rx, kr, mr, true, opts.interp);
Wk = kron(sd_r.Wc, sd_t.Wc);                     % mode-pair weights on slice-column pairs, row (n-1)*Mt + m
[nm_t, nm_r] = ndgrid(T_tx.mode_names(mt), T_rx.mode_names(mr));
mode_pairs = nm_t(:) + " | " + nm_r(:);
dg = T_tx.bin_dB;

% ---- Pairs
n_pct = numel(opts.pct);
n_mp = size(Wk, 1);
pct_env = NaN(n_pair, n_pct);
worst = zeros(n_pair, n_pct);
mode_pct = NaN(n_pair, n_mp, n_pct);
max_path = NaN(n_pair, 1);
curve = cell(n_pair, 1);
o = opts;
% Batching (plan T5): group the pairs by the elevation (and bearing) of the
% side with more slice columns, usually the WSR-88D's 21 tilts, and make
% that side's slices and their FFT once per group; the other side is made
% once per distinct elevation (and bearing) and reused.
et = el_tx_deg;
er = el_rx_deg;
if opts.el_round_deg > 0
    et = round(et/opts.el_round_deg)*opts.el_round_deg;
    er = round(er/opts.el_round_deg)*opts.el_round_deg;
end
heavy_r = n_cols(sd_r) >= n_cols(sd_t);
if heavy_r
    [kh, kl] = deal([er, brg_r.*is_sec_r], [et, brg_t.*is_sec_t]);
else
    [kh, kl] = deal([et, brg_t.*is_sec_t], [er, brg_r.*is_sec_r]);
end
[uh, ~, gh] = unique(kh, "rows");
[ul, ~, gl] = unique(kl, "rows");
light = cell(size(ul, 1), 1);
sd_l = sd_t;
sd_h = sd_r;
if ~heavy_r
    [sd_l, sd_h] = deal(sd_r, sd_t);
end
pr = [];                                          % keys are sorted, so one interval's keys are adjacent
for i = 1:size(ul, 1)
    [light{i}.p, light{i}.k0, light{i}.g, pr] = slices_at(sd_l, ul(i, 1), ul(i, 2), opts.interp, pr);
end
L_light = max(cellfun(@(x) size(x.p, 1), light));
members = accumarray(gh, (1:n_pair)', [size(uh, 1) 1], @(x) {x});
% Groups whose heavy elevations fall in the same table interval share the
% part of the quantile interpolation that does not depend on the position
% in the interval (one sort per interval, not per elevation).
n_g = size(uh, 1);
if opts.interp == "quantile" && sd_h.sector_deg == 360
    ch = discretize(uh(:, 1), sd_h.grid);
else
    ch = (1:n_g)';
end
[~, ~, ch] = unique(ch);
chunks = accumarray(ch, (1:n_g)', [], @(x) {x});
keys_c = cellfun(@(x) uh(x, :), chunks, "UniformOutput", false);   % sliced for parfor
mem_c = cellfun(@(x) members(x), chunks, "UniformOutput", false);
G = cell(numel(chunks), 1);
if opts.parallel_flag
    parfor ic = 1:numel(chunks)
        G{ic} = run_chunk(sd_h, keys_c{ic}, mem_c{ic}, light, gl, heavy_r, L_light, Wk, dg, o);
    end
else
    for ic = 1:numel(chunks)
        G{ic} = run_chunk(sd_h, keys_c{ic}, mem_c{ic}, light, gl, heavy_r, L_light, Wk, dg, o);
    end
end
for ic = 1:numel(G)
    for j = 1:numel(G{ic})
        ip = members{chunks{ic}(j)};
        g = G{ic}{j};
        pct_env(ip, :) = g.q_env;
        worst(ip, :) = g.worst;
        mode_pct(ip, :, :) = g.q_mode;
        max_path(ip) = g.max_path;
        curve(ip) = g.curve;
    end
end

vn = "p" + replace(compose("%.10g", opts.pct), ".", "_");
R = struct();
R.pct = [table(el_tx_deg, el_rx_deg), array2table(pct_env, "VariableNames", vn), ...
    table(max_path, 'VariableNames', "max_on_path_dBi")];
R.worst_mode = reshape(mode_pairs(worst), n_pair, n_pct);
R.mode_pairs = mode_pairs;
if opts.per_mode
    R.mode_pct = mode_pct;
end
if opts.return_curve
    R.curve = curve;
end
n_interp = 2*(opts.interp == "quantile");
R.bias_dB = (2 + n_interp)*dg;
R.pct_list = opts.pct;
R.notes = [
    "Path " + T_tx.name + " (tx) -> " + T_rx.name + " (rx); M.2069-0 Annex 1 random mode (independent pointing)"
    "Envelope over " + n_mp + " mode pairs: at each gain the largest exceedance (bounds any mode mix)"
    "Percentiles are conservative: at most " + compose("%.3g", R.bias_dB) + " dB above the exact value (grid rounding)"
    "Not covered: persistence, volume locks and per-volume statistics (use nexrad_asr11_mag_rev1 time series)"];
end

% =========================================================================
function out = run_chunk(sd_h, keys, mem, light, gl, heavy_r, L_light, Wk, dg, opts)
% Groups of one heavy-side table interval (or a single group).
pr = [];
out = cell(size(keys, 1), 1);
for j = 1:size(keys, 1)
    [ph, k0h, gmh, pr] = slices_at(sd_h, keys(j, 1), keys(j, 2), opts.interp, pr);
    out{j} = run_group(ph, k0h, gmh, mem{j}, light, gl, heavy_r, L_light, Wk, dg, opts);
end
end

% -------------------------------------------------------------------------
function out = run_group(ph, k0h, gmh, ip, light, gl, heavy_r, L_light, Wk, dg, opts)
% Every pair of one group: the heavy side's slices (ph) and FFT once.
Lh = size(ph, 1);
nfft = 2^nextpow2(Lh + L_light - 1);
Fh = [];
if opts.conv_method == "fft"
    Fh = fft(ph, nfft, 1);                          % dim 1: a slice can be one bin
end
sh = find(any(ph > 0, 2));                         % rows of the heavy side's support
n = numel(ip);
n_pct = numel(opts.pct);
n_mp = size(Wk, 1);
out = struct("q_env", NaN(n, n_pct), "worst", zeros(n, n_pct), "q_mode", NaN(n, n_mp, n_pct), ...
    "max_path", NaN(n, 1), "curve", {cell(n, 1)});
for i = 1:n
    li = light{gl(ip(i))};
    L = Lh + size(li.p, 1) - 1;
    if heavy_r                                     % tx = light, rx = heavy
        pt = li.p;
        pr = ph;
    else
        pt = ph;
        pr = li.p;
    end
    nt = size(pt, 2);
    nr = size(pr, 2);
    if opts.conv_method == "fft"
        % Plan T4: one FFT per light column, one inverse FFT for every
        % tilt pair, columns in tx-fastest order as Wk expects.
        Fl = fft(li.p, nfft, 1);
        if heavy_r
            Ft = Fl;
            Fr = Fh;
        else
            Ft = Fh;
            Fr = Fl;
        end
        Gp = reshape(Ft.*reshape(Fr, nfft, 1, nr), nfft, nt*nr);
        C = real(ifft(Gp, [], 1));
        % Rows outside the exact support hold only FFT round-off; zero them
        % so the envelope curve ends where the direct convolution's does.
        sl = find(any(li.p > 0, 2));
        C([1:sl(1) + sh(1) - 2, sl(end) + sh(end):end], :) = 0;
        C = C(1:L, :);
    else
        C = zeros(L, nt*nr);
        for s = 1:nr
            for t = 1:nt
                C(:, (s - 1)*nt + t) = conv(pt(:, t), pr(:, s));
            end
        end
    end
    [q_env, worst, q_mode, curve] = finish_pair(C, li.k0 + k0h, dg, Wk, opts);
    out.q_env(i, :) = q_env;
    out.worst(i, :) = worst;
    out.q_mode(i, :, :) = q_mode;
    out.curve{i} = curve;
    out.max_path(i) = max(li.g) + max(gmh);
end
end

% -------------------------------------------------------------------------
function [q_env, worst, q_mode, curve] = finish_pair(C, k_off, dg, Wk, opts)
% Mode-pair distributions from the tilt-pair convolutions; percentiles,
% worst mode and the envelope curve. Conv bin c holds (k_off + c - 1)*bin_dB.
L = size(C, 1);
PM = max(C*Wk', 0);                                   % bins by mode pairs; clip round-off
E = cumsum(PM, 1, "reverse") - PM;                    % P(MAG > gain of bin), tail-accurate
gain = (k_off + (0:L - 1)')*dg;
n_pct = numel(opts.pct);
n_mp = size(PM, 2);
q_mode = NaN(1, n_mp, n_pct);
for j = 1:n_pct
    lim = 1 - opts.pct(j)/100;
    % Smallest gain with exceedance <= 1 - p, every mode pair at once
    % (max of a logical returns the first true row).
    [~, b] = max(E <= lim + 1e-13, [], 1);
    q_mode(1, :, j) = gain(b);
end
[q_env, worst] = max(q_mode, [], 2);
q_env = reshape(q_env, 1, n_pct);
worst = reshape(worst, 1, n_pct);
curve = [];
if opts.return_curve
    e_env = max(E, [], 2);
    keep = find(e_env > 0, 1, "last");
    keep = min(L, keep + 1);
    curve = table(gain(1:keep), 100*e_env(1:keep), 'VariableNames', ["gain_dBi", "exceed_pct"]);
end
end

% -------------------------------------------------------------------------
function n = n_cols(sd)
n = numel(sd.tilt)*(1 + sd.has_edge);
end

% -------------------------------------------------------------------------
function [p, k0, gmax, pr] = slices_at(sd, el, brg, how, pr)
% Distributions of every slice column at elevation el (bins by columns,
% double), the bin index k0 of row 1, and each column's largest gain
% toward the station. Columns are the tilts; for a sector with turnaround,
% the tilts' uniform parts, then their edge parts. brg is the bearing of
% the other station from the sector center (sector sides only). pr, when
% given, is the interpolation setup of el's table interval (reused when
% the next call falls in the same interval).
if nargin < 5
    pr = [];
end
if sd.sector_deg < 360
    n_tilt = numel(sd.tilt);
    lim = [brg - sd.sector_deg/2, brg + sd.sector_deg/2];   % az_off = bearing - beam azimuth
    hs = cell(n_tilt*(1 + sd.has_edge), 1);
    for t = 1:n_tilt
        hs{t} = gain_slice_rev1(sd.S, sd.tilt(t), el, lim);
        if sd.has_edge
            hs{n_tilt + t} = gain_slice_rev1(sd.S, sd.tilt(t), el, lim, 1);   % all time at the edges
        end
    end
    gmax = cellfun(@(x) x.gmax, hs)';
    k0 = min(cellfun(@(x) x.k_lo, hs));
    n = max(cellfun(@(x) x.k_lo + numel(x.p), hs)) - k0;
    p = zeros(n, numel(hs));
    for t = 1:numel(hs)
        p(hs{t}.k_lo - k0 + (1:numel(hs{t}.p)), t) = hs{t}.p;
    end
    return
end
if how == "exact"
    n_tilt = numel(sd.tilt);
    hs = cell(n_tilt, 1);
    gmax = zeros(1, n_tilt);
    for t = 1:n_tilt
        hs{t} = gain_slice_rev1(sd.S, sd.tilt(t), el);
        gmax(t) = hs{t}.gmax;
    end
    k0 = min(cellfun(@(x) x.k_lo, hs));
    n = max(cellfun(@(x) x.k_lo + numel(x.p), hs)) - k0;
    p = zeros(n, n_tilt);
    for t = 1:n_tilt
        p(hs{t}.k_lo - k0 + (1:numel(hs{t}.p)), t) = hs{t}.p;
    end
    return
end
k0 = sd.k0;
gmax = gmax_at(sd.G, sd.grid, el);
[i, a] = bracket(sd.grid, el);
n_tilt = size(sd.P, 3);
if how == "nearest" || a == 0 || a == 1
    if how == "nearest"
        i = i + (a >= 0.5);
    else
        i = i + a;
    end
    p = double(reshape(sd.P(:, i, :), [], n_tilt));
    return
end
if isempty(pr) || pr.i ~= i
    p0 = double(reshape(sd.P(:, i, :), [], n_tilt));
    p1 = double(reshape(sd.P(:, i + 1, :), [], n_tilt));
    pr = interp_prep(p0, p1);
    pr.i = i;
end
p = interp_apply(pr, a);
end

% -------------------------------------------------------------------------
function pr = interp_prep(p0, p1)
% Displacement interpolation between two distributions on the same bins,
% every column (tilt) at once: Q(u) = (1 - a)*Q0(u) + a*Q1(u) for every
% probability level u, rounded up to a bin (conservative). Considered:
% linear mixing of the two pmfs ((1 - a)*p0 + a*p1), which leaves two
% peaks where the gain should move; three sorts per column (rev 2).
% This is the part that does not depend on a (interp_apply). One sort of both
% CDFs gives every breakpoint u. For each u, the smallest bin k with
% F(k) >= u is 1 + (number of F values strictly below u): each CDF's count
% of values before the start of u's run of equal values in the sorted
% order (cummax of the run starts). Considered: a merge of the two sorted
% CDFs on their support only (find, two 1-D sorts): equal, but slower.
[n, K] = size(p0);
F0 = cumsum(p0, 1);
F0 = F0./F0(end, :);
F1 = cumsum(p1, 1);
F1 = F1./F1(end, :);
[u, I] = sort([F0; F1], 1);
m2 = 2*n;
from0 = I <= n;
c0 = [zeros(1, K); cumsum(from0, 1)];
c1 = [zeros(1, K); cumsum(~from0, 1)];
pos = repmat((1:m2)', 1, K);
run0 = cummax([true(1, K); diff(u, 1, 1) > 0].*pos, 1);   % first position of each tie run
lin = run0 + (m2 + 1)*(0:K - 1);                % row just before the run in c0, c1
pr = struct("n", n, "K", K, "Q0", c0(lin) + 1, "Q1", c1(lin) + 1, ...
    "m", diff([zeros(1, K); u], 1, 1), "col", repmat(1:K, m2, 1), "i", 0);
end

% -------------------------------------------------------------------------
function p = interp_apply(pr, a)
% Interpolated distribution at fraction a of the interval (interp_prep).
q = (1 - a)*pr.Q0 + a*pr.Q1;
b = min(max(ceil(q - 1e-9), 1), pr.n);
p = accumarray([b(:), pr.col(:)], pr.m(:), [pr.n pr.K]);
end

% -------------------------------------------------------------------------
function g = gmax_at(G, grid, el)
% Per-tilt maximum gain toward the station at el (linear in elevation).
[i, a] = bracket(grid, el);
if a == 0
    g = G(i, :);
else
    g = (1 - a)*G(i, :) + a*G(i + 1, :);
end
end

% -------------------------------------------------------------------------
function sd = side(T, k, m, is_rx, how)
% One end of the path: its stored slices (full rotation), or a slice cache
% for "exact" and for sectors; Wc holds the selected modes' weights on the
% slice columns.
sec = 360;
if isfield(T, "sector")
    sec = T.sector.sector_deg;
end
sd = struct("grid", T.el_grid_deg, "k0", T.k0, "tilt", T.tilt_deg(k), "S", [], ...
    "sector_deg", sec, "has_edge", false, "P", [], "G", []);
W = T.W(m, k);
sd.Wc = W;
if sec < 360
    We = T.W_edge(m, k);
    sd.has_edge = any(We(:) > 0);
    if sd.has_edge
        sd.Wc = [W - We, We];                    % uniform parts, then edge parts
    end
end
use_rx = is_rx && ~isempty(T.ant_rx);
if sec >= 360 && how ~= "exact"
    if is_rx && ~isempty(T.P_rx)
        sd.P = T.P_rx(:, :, k);
        sd.G = T.gmax_rx(:, k);
    else
        sd.P = T.P_tx(:, :, k);
        sd.G = T.gmax_tx(:, k);
    end
    return
end
if ~isfield(T, "slice_opts")
    error("mag_pair:oldTable", "Table %s predates exact slices; rebuild it with gain_dist_table_rev1.", T.name);
end
ant = T.ant;
if use_rx
    ant = T.ant_rx;
end
sd.S = gain_slice_prep_rev1(ant, T.slice_opts);
end

% -------------------------------------------------------------------------
function tf = is_sector(T)
tf = isfield(T, "sector") && T.sector.sector_deg < 360;
end

% -------------------------------------------------------------------------
function b = expand_brg(b, n, label)
if isscalar(b)
    b = b + zeros(n, 1);
end
b = b(:);
if numel(b) ~= n || any(~isfinite(b))
    error("mag_pair:badBearing", "opts.%s must be a scalar or one finite value per pair.", label);
end
end

% -------------------------------------------------------------------------
function [i, a] = bracket(grid, el)
% Grid interval [grid(i), grid(i+1)] that holds el, and the fraction a.
i = discretize(el, grid);
if isnan(i)
    error("mag_pair:outOfGrid", "Elevation %g deg is outside the table grid.", el);
end
a = (el - grid(i))/(grid(i + 1) - grid(i));
if a < 1e-9
    a = 0;
elseif a > 1 - 1e-9
    a = 1;
end
end

% -------------------------------------------------------------------------
function m = select_modes(T, names, label)
if isempty(names)
    m = 1:numel(T.mode_names);
    return
end
names = string(names);
[ok, m] = ismember(names, T.mode_names);
if ~all(ok)
    error("mag_pair:unknownMode", "opts.%s: mode(s) not in table %s: %s.", label, T.name, strjoin(names(~ok), ", "));
end
end

% -------------------------------------------------------------------------
function check_sync(rt, rr)
% Random mode needs unequal rates not in a small-integer ratio (the
% check_random_mode_valid / joint_period rules of radar_mutual_gain_cdf_rev1,
% repeated here because those are local functions). Volume scans with
% several rates are not checked, as in that function.
for i = 1:numel(rt)
    for j = 1:numel(rr)
        if ~isscalar(rt{i}) || ~isscalar(rr{j})
            continue
        end
        r1 = abs(rt{i});
        r2 = abs(rr{j});
        if ismembertol(r1, r2, 1e-9)
            warning("mag_pair:equalRates", join(["Equal rotation rates keep a fixed phase offset. " ...
                "The convolution averages over an unknown offset (MACD). For a known offset use " ...
                "radar_mutual_gain_cdf_rev1 with opts.scan_sync = ""synchronized""."], ""));
            return
        end
        [p, q] = rat(r1/r2, 1e-9*r1/r2);
        if p <= 100 && q <= 100
            warning("mag_pair:commensurateRates", join(["Rotation rates are in a %d:%d ratio; the pointing " ...
                "pair repeats every %.1f s. The convolution assumes the rates drift apart. Check with " ...
                "radar_mutual_gain_cdf_rev1, opts.mode = ""time_series""."], ""), p, q, 360*p/r1);
            return
        end
    end
end
end

% -------------------------------------------------------------------------
function check_range(el, grid, label, name)
if any(el < grid(1) - 1e-12 | el > grid(end) + 1e-12)
    error("mag_pair:outOfGrid", "%s must lie in the %s table grid [%g, %g] deg.", label, name, grid(1), grid(end));
end
end

% -------------------------------------------------------------------------
function opts = fill_opts(opts)
d = struct("pct", [50 90 99 99.9 99.99 99.999], "brg_tx_deg", 0, "brg_rx_deg", 0, "tx_modes", strings(1, 0), ...
    "rx_modes", strings(1, 0), "per_mode", false, "return_curve", false, ...
    "interp", "quantile", "parallel_flag", 0, "conv_method", "fft", "el_round_deg", 0);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i})
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("mag_pair:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
opts.interp = string(opts.interp);
if ~ismember(opts.interp, ["quantile", "nearest", "exact"])
    error("mag_pair:badInterp", "opts.interp must be ""quantile"", ""nearest"" or ""exact"".");
end
opts.conv_method = string(opts.conv_method);
if ~ismember(opts.conv_method, ["fft", "direct"])
    error("mag_pair:badConv", "opts.conv_method must be ""fft"" or ""direct"".");
end
if ~isscalar(opts.el_round_deg) || opts.el_round_deg < 0
    error("mag_pair:badRound", "opts.el_round_deg must be a scalar >= 0.");
end
if any(opts.pct < 0 | opts.pct >= 100)
    error("mag_pair:badPct", "opts.pct values must lie in [0, 100).");
end
end
