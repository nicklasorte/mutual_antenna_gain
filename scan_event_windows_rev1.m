function [idx, info] = scan_event_windows_rev1(C, duration_s, dt_s, opts)
%SCAN_EVENT_WINDOWS_REV1 Time steps where either beam is near the path.
%   [idx, info] = scan_event_windows_rev1(C, duration_s, dt_s, opts)
%   returns the indices idx into the fixed time grid t = (0:dt_s:duration_s)'
%   of every step at which the WSR-88D beam azimuth is within
%   opts.delta_n_deg of the bearing to the ATC radar, or the ATC beam
%   azimuth is within opts.delta_a_deg of the bearing to the NEXRAD, for
%   the scan clock C (scan_clock_rev1). Away from these windows both beams
%   are off the path, so the mutual antenna gain (MAG) is at most the sum
%   of the two antennas' gains beyond info.delta_eff_n_deg and
%   info.delta_eff_a_deg (ts_floor_rev1). Every MAG maximum above that
%   floor therefore lies in a window, and the windows are a subset of the
%   fixed grid: maxima taken over idx equal fixed-grid maxima exactly.
%
%   Method: evaluate the clock on a coarse grid (opts.dt_coarse_s), flag
%   coarse samples within delta of the path, and keep every fine step in
%   a coarse interval that touches a flagged sample. Between two unflagged
%   coarse samples a beam moving at rate r stays at least delta - r*dt_c
%   from the path (it cannot cross the path in one interval while
%   r*dt_c < 2*delta), which gives delta_eff. Works for rotation, sectors,
%   turnaround dwells and rate jitter alike.
%
%   opts struct, all fields optional:
%     .delta_n_deg  NEXRAD window half-width in azimuth (default 4)
%     .delta_a_deg  ATC window half-width in azimuth (default 4)
%     .dt_coarse_s  coarse step (default 0.02)
%
%   idx   column vector of grid indices (sorted)
%   info  struct: delta_eff_n_deg, delta_eff_a_deg, n_grid, share (of the
%         grid kept)
%
%   See also scan_clock_rev1, scan_event_library_rev1, ts_floor_rev1,
%   mag_pair_timeseries_rev1.

arguments
    C (1, 1) struct
    duration_s (1, 1) double {mustBePositive}
    dt_s (1, 1) double {mustBePositive}
    opts (1, 1) struct = struct()
end
opts = fill_opts(opts);
r_n = max(abs(C.rate_n));
r_a = C.rate_a + C.jit_a;
dtc = opts.dt_coarse_s;
if r_n*dtc >= 2*opts.delta_n_deg || r_a*dtc >= 2*opts.delta_a_deg
    error("scan_events:coarseStep", "dt_coarse_s too large: a beam would move over 2*delta in one step.");
end
tc = (0:dtc:duration_s)';
if tc(end) < duration_s
    tc(end + 1) = duration_s;
end
[st_n, st_a] = scan_clock_state_rev1(C, tc);
flag = abs(wrapTo180(0 - st_n.mech_az)) <= opts.delta_n_deg | ...
    abs(wrapTo180(180 - st_a.mech_az)) <= opts.delta_a_deg;
% Coarse intervals [tc(m), tc(m+1)] that touch a flagged sample.
keep = flag(1:end - 1) | flag(2:end);
n_grid = numel(0:dt_s:duration_s);              % the grid of nexrad_asr11_mag_rev1
mask = false(n_grid, 1);
m = find(keep);
k0 = max(1, floor(tc(m)/dt_s) + 1);
k1 = min(n_grid, ceil(tc(m + 1)/dt_s) + 1);
% Mark each fine run with a difference array (no loop over intervals).
d = zeros(n_grid + 1, 1);
d = d + accumarray(k0, 1, [n_grid + 1, 1]) - accumarray(k1 + 1, 1, [n_grid + 1, 1]);
mask(cumsum(d(1:n_grid)) > 0) = true;
idx = find(mask);
info = struct("delta_eff_n_deg", opts.delta_n_deg - r_n*dtc, ...
    "delta_eff_a_deg", opts.delta_a_deg - r_a*dtc, "n_grid", n_grid, "share", numel(idx)/n_grid);
end

% -------------------------------------------------------------------------
function opts = fill_opts(opts)
d = struct("delta_n_deg", 4, "delta_a_deg", 4, "dt_coarse_s", 0.02);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("scan_events:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
end
