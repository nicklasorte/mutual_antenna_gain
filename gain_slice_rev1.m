function h = gain_slice_rev1(S, tilt_deg, el_deg, az_lim_deg, edge_frac)
%GAIN_SLICE_REV1 Exact time-weighted gain distribution at one tilt and elevation.
%   h = gain_slice_rev1(S, tilt_deg, el_deg) returns the distribution of the
%   gain of a mechanically rotating antenna (prepared by
%   gain_slice_prep_rev1) toward another station at elevation el_deg, while
%   the beam, at tilt tilt_deg, sweeps azimuth at a constant rate through a
%   full revolution. Gains are rounded up to the next bin edge (bin_dB), so
%   every percentile is conservative by at most one bin.
%
%   h = gain_slice_rev1(S, tilt_deg, el_deg, az_lim_deg) restricts the sweep
%   to an azimuth interval, for sector scans. az_lim_deg = [lo hi] is the
%   range of az_off = (bearing of the other station) - (beam azimuth), deg,
%   the offset antenna_gain_rev1 uses, with 0 < hi - lo <= 360. Time is
%   uniform across the interval (constant sweep rate, no turnaround).
%
%   h = gain_slice_rev1(S, tilt_deg, el_deg, az_lim_deg, edge_frac) also
%   puts edge_frac of the time at the two ends of the interval, half at
%   each (a bidirectional sweep that dwells at each edge, sector_sweep_rev1);
%   the rest stays uniform across the interval.
%
%   S          struct from gain_slice_prep_rev1
%   tilt_deg   beam elevation (deg), scalar
%   el_deg     elevation of the other station above the local horizontal
%              (takeoff angle, deg), scalar
%   az_lim_deg optional [lo hi] (default full rotation)
%   edge_frac  optional share of time at the interval ends, 0 to 1
%              (default 0; needs az_lim_deg)
%
%   h struct:
%     .k_lo   integer index of the first bin: bin b holds gain
%             (k_lo + b - 1)*bin_dB, the upper edge of its interval
%     .p      probability mass per bin (column, sums to 1)
%     .gmax   largest gain toward the station over the sweep (dBi; exact
%             for closed_form and separable, grid maximum for grid)
%
%   Example:
%     S = gain_slice_prep_rev1(ant_nexrad);
%     h = gain_slice_rev1(S, 0.5, 1.7);                % full rotation
%     h = gain_slice_rev1(S, 0.5, 1.7, [-30 30]);      % 60 deg sector on the path
%
%   See also gain_slice_prep_rev1, gain_dist_table_rev1, mag_pair_cdf_rev1.

arguments
    S (1, 1) struct
    tilt_deg (1, 1) double {mustBeReal, mustBeFinite}
    el_deg (1, 1) double {mustBeReal, mustBeFinite}
    az_lim_deg double = []
    edge_frac (1, 1) double {mustBeInRange(edge_frac, 0, 1)} = 0
end
if edge_frac > 0 && isempty(az_lim_deg)
    error("gain_slice:edgeNeedsLim", "edge_frac needs a sector (az_lim_deg).");
end
full = isempty(az_lim_deg);
if ~full
    [lo, L] = check_lim(az_lim_deg);
    full = L >= 360 - 1e-12;
end
dg = S.bin_dB;
switch S.method
    case "closed_form"
        if full
            h = closed_form(S, tilt_deg, el_deg, -180, 360, dg);
        else
            h = closed_form(S, tilt_deg, el_deg, lo, L, dg);
        end
    case "separable"
        if full
            h = separable_full(S, tilt_deg, el_deg, dg);
        else
            m = in_lim(S.az, lo, L);
            h = separable_subset(S, m, tilt_deg, el_deg, dg);
        end
    case "grid"
        m = true(size(S.az));
        if ~full
            m = in_lim(S.az, lo, L);
        end
        h = grid_slice(S, m, tilt_deg, el_deg, dg);
    otherwise
        error("gain_slice:badMethod", "Unknown method ""%s"" in S.", S.method);
end
if edge_frac > 0
    h = add_edges(S, h, tilt_deg, el_deg, az_lim_deg, edge_frac, dg);
end
end

% =========================================================================
function h = add_edges(S, h, tilt, el, lim, f, dg)
% Mix in point masses f/2 at each end of the interval (beam parked there).
a = lim(:);
z = zeros(2, 1);
g = antenna_gain_rev1(S.pa, 0, el, struct("mech_az", -a, "tilt", z + tilt, "steer_az", z, "steer_el", z));
k = ceil(g/dg - 1e-9);
k_lo = min([h.k_lo; k]);
n = max([h.k_lo + numel(h.p) - 1; k]) - k_lo + 1;
p = zeros(n, 1);
p(h.k_lo - k_lo + (1:numel(h.p))) = (1 - f)*h.p;
p = p + accumarray(k - k_lo + 1, f/2, [n 1]);
h = struct("k_lo", k_lo, "p", p, "gmax", max([h.gmax; g]));
end

% =========================================================================
function h = closed_form(S, tilt, el, lo, L, dg)
A = sind(el)*sind(tilt);
B = cosd(el)*cosd(tilt);
d_min = dist_min(lo, L);                         % smallest |az_off| in the interval
d_max = dist_max(lo, L);                         % largest |az_off|
gam_lo = acosd(min(max(A + B*cosd(d_min), -1), 1));
gam_hi = acosd(min(max(A + B*cosd(d_max), -1), 1));
h_gmax = gain_at(S, gam_lo);
if B < 1e-12 || gam_hi - gam_lo < 1e-12
    % Beam or station at the zenith (or a point sector): one gain.
    k = ceil(h_gmax/dg - 1e-9);
    h = struct("k_lo", k, "p", 1, "gmax", h_gmax);
    return
end
in = S.gam_b > gam_lo & S.gam_b < gam_hi;
gam = [gam_lo; S.gam_b(in); gam_hi];
g = [h_gmax; S.g_b(in); gain_at(S, gam_hi)];
[gam, g] = refine_cells(gam, g, dg);
F = share(gam, A, B, lo, L);
F(1) = 0;
F(end) = 1;
m = diff(F);
gu = max(g(1:end - 1), g(2:end));               % cell maximum (pattern monotone in each cell)
keep = m > 0;
k = ceil(gu(keep)/dg - 1e-9);
k_lo = min(k);
p = accumarray(k - k_lo + 1, m(keep));
h = struct("k_lo", k_lo, "p", p/sum(p), "gmax", h_gmax);
end

% -------------------------------------------------------------------------
function [gam, g] = refine_cells(gam, g, dg)
% Split every cell whose gain changes by more than one bin into equal
% sub-cells in gamma, with gains interpolated linearly (exact for a table,
% which is linear between nodes). Each sub-cell then spans at most one bin,
% so taking its larger end gain moves mass up by at most one bin.
n_sub = max(1, ceil(abs(diff(g))/dg));
if all(n_sub == 1)
    return
end
n_c = numel(n_sub);
idx = repelem((1:n_c)', n_sub);                  % parent cell of each sub-cell
first = cumsum([1; n_sub(1:end - 1)]);
j = (1:sum(n_sub))' - first(idx);                % 0 .. n_sub - 1 within the cell
f = j./n_sub(idx);
gam = [gam(idx) + f.*(gam(idx + 1) - gam(idx)); gam(end)];
g = [g(idx) + f.*(g(idx + 1) - g(idx)); g(end)];
end

% -------------------------------------------------------------------------
function F = share(gam, A, B, lo, L)
% Share of the sweep interval with off-axis angle <= gam: the interval's
% overlap with |az_off| <= phi(gam), counting the copy at +360 deg (lo is
% in [-180, 180), so the interval ends before 540 deg).
phi = acosd(min(max((cosd(gam) - A)/B, -1), 1));
hi = lo + L;
ov0 = max(0, min(hi, phi) - max(lo, -phi));
ov1 = max(0, min(hi, 360 + phi) - max(lo, 360 - phi));
F = (ov0 + ov1)/L;
end

% -------------------------------------------------------------------------
function g = gain_at(S, gam)
z = zeros(size(gam));
g = antenna_gain_rev1(S.pa, 0, 0, struct("mech_az", -gam, "tilt", z, "steer_az", z, "steer_el", z));
end

% -------------------------------------------------------------------------
function h = separable_full(S, tilt, el, dg)
% Every slice is the sorted azimuth cut shifted by the elevation term and
% clamped at the floor; bins are counted by binary search (discretize).
s = elev_term(S, tilt, el);
v_min = S.g_max + max(S.g_az_sorted(1) + s, S.floor_dB);
v_max = S.g_max + max(S.g_az_sorted(end) + s, S.floor_dB);
k_lo = ceil(v_min/dg - 1e-9);
k_hi = ceil(v_max/dg - 1e-9);
k = (k_lo:k_hi)';
% Bin k holds v in ((k - 1 + 1e-9)*dg, (k + 1e-9)*dg]; above the floor
% v <= x  <=>  g_az <= x - g_max - s. Clamped samples all land in bin k_lo
% or above, so counts below k_lo are zero.
x = (k + 1e-9)*dg - S.g_max - s;
c = count_le(S.u_az, S.cnt_az, x);
c(end) = numel(S.g_az_sorted);
p = diff([0; c])/numel(S.g_az_sorted);
% Exact peak: the gain at az_off = 0, g_max + s (the grid excludes 0).
h = struct("k_lo", k_lo, "p", p, "gmax", S.g_max + s);
end

% -------------------------------------------------------------------------
function h = separable_subset(S, m, tilt, el, dg)
s = elev_term(S, tilt, el);
v = S.g_max + max(S.g_az(m) + s, S.floor_dB);
h = bin_samples(v, dg);
end

% -------------------------------------------------------------------------
function s = elev_term(S, tilt, el)
% g_el(el, tilt) clamped at the floor, from the gain at az_off = 0 (where
% g_az = 0); clamping does not change max(g_az + g_el, floor).
g0 = antenna_gain_rev1(S.pa, 0, el, struct("mech_az", 0, "tilt", tilt, "steer_az", 0, "steer_el", 0));
s = g0 - S.g_max;
end

% -------------------------------------------------------------------------
function h = grid_slice(S, m, tilt, el, dg)
az = S.az(m);
z = zeros(size(az));
g = antenna_gain_rev1(S.pa, 0, el, struct("mech_az", -az, "tilt", z + tilt, "steer_az", z, "steer_el", z));
h = bin_samples(g, dg);
end

% -------------------------------------------------------------------------
function h = bin_samples(v, dg)
if isempty(v)
    error("gain_slice:emptySector", "The azimuth interval holds no grid point; widen it or refine az_step_deg.");
end
k = ceil(v/dg - 1e-9);
k_lo = min(k);
p = accumarray(k - k_lo + 1, 1)/numel(v);
h = struct("k_lo", k_lo, "p", p, "gmax", max(v));
end

% -------------------------------------------------------------------------
function c = count_le(u, cnt, x)
% Number of samples <= x, for each x, from the unique sorted values u and
% the count of samples <= each (binary search with discretize).
% Considered: histcounts (one pass over all samples per slice).
b = discretize(x, [-Inf; u(:); Inf]);           % e(b) <= x < e(b + 1)
c_ext = [0; cnt(:)];
c = c_ext(b);
end

% -------------------------------------------------------------------------
function m = in_lim(az, lo, L)
m = mod(az - lo, 360) <= L;
end

% -------------------------------------------------------------------------
function d = dist_min(lo, L)
hi = lo + L;
if (lo <= 0 && hi >= 0) || (lo <= 360 && hi >= 360)
    d = 0;
else
    d = min(abs([lo, hi, lo - 360, hi - 360]));
end
end

% -------------------------------------------------------------------------
function d = dist_max(lo, L)
hi = lo + L;
if (lo <= 180 && hi >= 180) || lo <= -180
    d = 180;
else
    d = max(abs(wrapTo180([lo, hi])));
end
end

% -------------------------------------------------------------------------
function [lo, L] = check_lim(lim)
if numel(lim) ~= 2 || ~all(isfinite(lim))
    error("gain_slice:badLim", "az_lim_deg must be [lo hi] (deg).");
end
L = lim(2) - lim(1);
if L <= 0 || L > 360 + 1e-9
    error("gain_slice:badLim", "az_lim_deg must satisfy 0 < hi - lo <= 360 (got %g).", L);
end
L = min(L, 360);
lo = mod(lim(1) + 180, 360) - 180;              % [-180, 180)
end
