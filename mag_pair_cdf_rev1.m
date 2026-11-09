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
%     .tx_modes     T_tx mode names to include (default all)
%     .rx_modes     T_rx mode names to include (default all)
%     .per_mode     true also returns every mode pair's percentiles
%                   (default false)
%     .return_curve true also returns the envelope exceedance curve of
%                   each pair (default false)
%     .interp       "quantile" (default) interpolates between the two
%                   neighboring grid elevations in quantile space (no false
%                   double peaks); "nearest" uses the nearest grid
%                   elevation
%     .parallel_flag 1 runs the pairs in parfor (default 0)
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
check_range(el_tx_deg, T_tx.el_grid_deg, "el_tx_deg", T_tx.name);
check_range(el_rx_deg, T_rx.el_grid_deg, "el_rx_deg", T_rx.name);

% ---- Mode selection and weights. Only tilts with weight are kept.
mt = select_modes(T_tx, opts.tx_modes, "tx_modes");
mr = select_modes(T_rx, opts.rx_modes, "rx_modes");
check_sync(T_tx.mode_rates(mt), T_rx.mode_rates(mr));
Wt = T_tx.W(mt, :);
Wr = T_rx.W(mr, :);
kt = find(any(Wt > 0, 1));
kr = find(any(Wr > 0, 1));
Wt = Wt(:, kt);
Wr = Wr(:, kr);
Wk = kron(Wr, Wt);                               % mode-pair weights on tilt pairs, row (n-1)*Mt + m
[nm_t, nm_r] = ndgrid(T_tx.mode_names(mt), T_rx.mode_names(mr));
mode_pairs = nm_t(:) + " | " + nm_r(:);

P_tx = T_tx.P_tx(:, :, kt);
G_tx = T_tx.gmax_tx(:, kt);
if isempty(T_rx.P_rx)
    P_rx = T_rx.P_tx(:, :, kr);
    G_rx = T_rx.gmax_tx(:, kr);
else
    P_rx = T_rx.P_rx(:, :, kr);
    G_rx = T_rx.gmax_rx(:, kr);
end
k_off = T_tx.k0 + T_rx.k0;                       % conv bin c holds (k_off + c - 1)*bin_dB
dg = T_tx.bin_dB;

% ---- Pairs
n_pair = numel(el_tx_deg);
n_pct = numel(opts.pct);
n_mp = size(Wk, 1);
pct_env = NaN(n_pair, n_pct);
worst = zeros(n_pair, n_pct);
mode_pct = NaN(n_pair, n_mp, n_pct);
max_path = NaN(n_pair, 1);
curve = cell(n_pair, 1);
grid_t = T_tx.el_grid_deg;
grid_r = T_rx.el_grid_deg;
o = opts;
if opts.parallel_flag
    parfor ip = 1:n_pair
        [pct_env(ip, :), worst(ip, :), mode_pct(ip, :, :), max_path(ip), curve{ip}] = ...
            one_pair(P_tx, G_tx, grid_t, el_tx_deg(ip), P_rx, G_rx, grid_r, el_rx_deg(ip), ...
            Wk, k_off, dg, o);
    end
else
    for ip = 1:n_pair
        [pct_env(ip, :), worst(ip, :), mode_pct(ip, :, :), max_path(ip), curve{ip}] = ...
            one_pair(P_tx, G_tx, grid_t, el_tx_deg(ip), P_rx, G_rx, grid_r, el_rx_deg(ip), ...
            Wk, k_off, dg, o);
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
function [q_env, worst, q_mode, max_path, curve] = one_pair(P_tx, G_tx, grid_t, et, ...
    P_rx, G_rx, grid_r, er, Wk, k_off, dg, opts)
pt = slices_at(P_tx, grid_t, et, opts.interp);        % bins by tx tilts
pr = slices_at(P_rx, grid_r, er, opts.interp);        % bins by rx tilts
nt = size(pt, 2);
nr = size(pr, 2);
L = size(pt, 1) + size(pr, 1) - 1;
C = zeros(L, nt*nr);
for s = 1:nr
    for t = 1:nt
        C(:, (s - 1)*nt + t) = conv(pt(:, t), pr(:, s));
    end
end
PM = max(C*Wk', 0);                                   % bins by mode pairs; clip conv round-off
E = cumsum(PM, 1, "reverse") - PM;                    % P(MAG > gain of bin), tail-accurate
gain = (k_off + (0:L - 1)')*dg;
n_pct = numel(opts.pct);
n_mp = size(PM, 2);
q_mode = NaN(1, n_mp, n_pct);
for j = 1:n_pct
    lim = 1 - opts.pct(j)/100;
    for m = 1:n_mp
        b = find(E(:, m) <= lim + 1e-13, 1);          % smallest gain with exceedance <= 1 - p
        q_mode(1, m, j) = gain(b);
    end
end
[q_env, worst] = max(q_mode, [], 2);
q_env = reshape(q_env, 1, n_pct);
worst = reshape(worst, 1, n_pct);
max_path = max(gmax_at(G_tx, grid_t, et)) + max(gmax_at(G_rx, grid_r, er));
curve = [];
if opts.return_curve
    e_env = max(E, [], 2);
    keep = find(e_env > 0, 1, "last");
    keep = min(L, keep + 1);
    curve = table(gain(1:keep), 100*e_env(1:keep), 'VariableNames', ["gain_dBi", "exceed_pct"]);
end
end

% -------------------------------------------------------------------------
function p = slices_at(P, grid, el, how)
% Distributions of every tilt at elevation el (bins by tilts, double).
[i, a] = bracket(grid, el);
n_tilt = size(P, 3);
if how == "nearest" || a == 0 || a == 1
    if how == "nearest"
        i = i + (a >= 0.5);
    else
        i = i + a;
    end
    p = double(reshape(P(:, i, :), [], n_tilt));
    return
end
p0 = double(reshape(P(:, i, :), [], n_tilt));
p1 = double(reshape(P(:, i + 1, :), [], n_tilt));
p = zeros(size(p0));
for t = 1:n_tilt
    p(:, t) = quantile_interp(p0(:, t), p1(:, t), a);
end
end

% -------------------------------------------------------------------------
function p = quantile_interp(p0, p1, a)
% Displacement interpolation between two distributions on the same bins:
% Q(u) = (1 - a)*Q0(u) + a*Q1(u) for every probability level u, rounded up
% to a bin (conservative). Considered: linear mixing of the two pmfs
% ((1 - a)*p0 + a*p1), which leaves two peaks where the gain should move.
F0 = cumsum(p0);
F0 = F0/F0(end);
F1 = cumsum(p1);
F1 = F1/F1(end);
u = unique([F0; F1]);
u = u(u > 0);
q = (1 - a)*q_index(F0, u) + a*q_index(F1, u);
b = ceil(q - 1e-9);
p = accumarray(b, diff([0; u]), size(p0));
end

% -------------------------------------------------------------------------
function k = q_index(F, u)
% Smallest bin k with F(k) >= u, for each u in (0, 1] (F nondecreasing,
% F(end) = 1). The leading 0 (bin 0) gives interp1 at least two points when
% all the mass is in one bin; u > 0 never maps to bin 0.
[Fu, ia] = unique([0; F], "first");
k = interp1(Fu, ia - 1, u, "next");
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
d = struct("pct", [50 90 99 99.9 99.99 99.999], "tx_modes", strings(1, 0), ...
    "rx_modes", strings(1, 0), "per_mode", false, "return_curve", false, ...
    "interp", "quantile", "parallel_flag", 0);
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
if ~ismember(opts.interp, ["quantile", "nearest"])
    error("mag_pair:badInterp", "opts.interp must be ""quantile"" or ""nearest"".");
end
if any(opts.pct < 0 | opts.pct >= 100)
    error("mag_pair:badPct", "opts.pct values must lie in [0, 100).");
end
end
