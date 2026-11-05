function out = radar_mutual_gain_cdf_rev1(sys1, sys2, opts)
%RADAR_MUTUAL_GAIN_CDF_REV1 Mutual antenna gain CDF between two scanning radars.
%   out = radar_mutual_gain_cdf_rev1(sys1, sys2, opts) returns the
%   distribution of mutual antenna gain (antenna coupling) between two
%   radars whose antennas rotate or scan:
%       g12 = G_tx,1(toward 2) + G_rx,2(toward 1)   (dBi), path 1 -> 2
%       g21 = G_tx,2(toward 1) + G_rx,1(toward 2)   (dBi), path 2 -> 1
%   g12 equals g21 when each radar uses one pattern for transmit and receive.
%
%   Method: Rec. ITU-R M.2069-0 Annex 1. Two modes:
%     "random"       Independent, time-weighted pointing states for the two
%                    radars (the M.2069-0 Fig. 10 method, 5 million
%                    samples). Valid when the rotation rates are unequal and
%                    not in a simple ratio, or for equal rates averaged
%                    over an unknown phase (MACD, see opts.scan_sync).
%     "time_series"  Each antenna rotates at its own rate from a random
%                    initial phase (M.2069-0 Figs 1-8). Captures equal or
%                    nearly equal rates, and returns how long and how often
%                    high coupling lasts.
%   Scan synchronization (opts.scan_sync), for radars that rotate at the
%   same rate:
%     "unsynchronized"  The relative rotation phase is unknown and uniform
%                       (default). In random mode this is the Mutual Antenna
%                       Coupling Distribution (MACD) of the 2.7 GHz link
%                       budget working paper, Appendix B: each antenna's
%                       pointing is independent and uniform.
%     "synchronized"    Both antennas rotate at one rate with a fixed phase
%                       offset, az2(t) = az1(t) + opts.sync_offset_deg. The
%                       default offset 0 starts both antennas at the same
%                       angle (the working paper's same-rate method, Section
%                       3.1); out.sync holds the maximum coupling over one
%                       revolution, the MAG that method uses.
%   The default antenna pattern is Statgain (NTIA TM-13-489 Table 4-1);
%   Rec. ITU-R F.1245, Rec. ITU-R M.1851-2 patterns and planar arrays are
%   also available.
%   See antenna_prepare_rev1.
%
%   Geometry is location-agnostic. Azimuths are measured from the line
%   joining the radars (radar 2 lies at 0 deg from radar 1; radar 1 at
%   180 deg from radar 2). Because each antenna sweeps all azimuths, the
%   link bearing only shifts the rotation phase and does not change the
%   CDF. The path elevation at each end is an input (default 0 deg, a
%   horizontal line-of-sight path between ground radars). Polarization
%   mismatch is not modeled.
%
%   sys1, sys2  structs:
%     .name       char or string
%     .tx         antenna model struct (antenna_prepare_rev1)
%     .rx         antenna model struct; defaults to .tx when absent or empty
%     .scan       scan model struct:
%        .type        "mech_rotating"  reflector or array that rotates in
%                                      azimuth
%                     "array_multiface" fixed faces, beam steered over 360
%        .rate_deg_s  azimuth rate, positive clockwise. A vector sets one
%                     rate per elevation step of a volume scan. .rate_rpm
%                     can be given instead (rate_deg_s = 6*rate_rpm).
%        .tilt_deg    reflector beam elevation per step (default 0)
%        .n_rev       revolutions per step (default 1, integers). With a
%                     sector, sweeps per step.
%        .sector_deg  azimuth sector width (default 360, full rotation).
%                     Below 360, each step sweeps the sector once in the
%                     direction of rotation and restarts at its edge, so a
%                     step lasts n_rev*sector_deg/|rate|. Elevation-step
%                     time shares are unchanged by the sector.
%        .sector_center_deg  sector center azimuth (default 0)
%        .rate_jitter_deg_s  time_series only: each trial offsets the rate
%                     by U(-j, +j) deg/s (default 0). M.2069-0 Annex 1
%                     quotes a +/-0.5 deg/s drive accuracy.
%        Arrays only (defaults in parentheses):
%        .beam_mode   "fixed" (default) | "random" (mech_rotating)
%                     "random" (default) | "sequential" (array_multiface)
%        .steer_az_deg, .steer_el_deg  fixed steering (0, 0: broadside on
%                                      the horizon)
%        .beam_az_lim_deg [min max]    random steering azimuth relative to
%                                      the face normal ([0 0])
%        .beam_el_lim_deg [min max]    random steering elevation ([0 0],
%                                      horizon only)
%        .dwell_s     beam dwell time (0.05 s, within the "tens to hundreds
%                     of ms" of M.1464-2 Section 2.3)
%        .n_faces     array_multiface faces (4, M.1464-2 Radar F1)
%        .face_az0_deg  normal of face 1 relative to the link (0)
%        .rate_deg_s  (array_multiface "sequential") effective beam sweep
%                     rate across the faces
%
%   opts struct:
%     .seed           required RNG seed
%     .mode           "random" (default) | "time_series"
%     .n_samples      random mode samples (default 5e6, as in M.2069-0)
%     .duration_s     time_series duration per trial (default 3600)
%     .dt_s           time_series step (default: one tenth of the shortest
%                     main-beam transit or dwell time)
%     .n_trials       time_series trials, each with new phases (default 4)
%     .el12_deg       elevation of radar 2 seen from radar 1 (default 0)
%     .el21_deg       elevation of radar 1 seen from radar 2 (default 0)
%     .parallel_flag  0 (default) or 1 (parfor over chunks or trials)
%     .chunk_size     random mode samples per chunk (default 2.5e5)
%     .pct            percentiles to report (default [50 90 99 99.9 99.99 99.999])
%     .cdf_pct        percent grid of the CDF arrays, 0 to 100 (default 0:1:100
%                     plus 99.5, 99.9, 99.95, 99.99, 99.995 and 99.999)
%     .event_thresh_dBi  coupling thresholds for event statistics
%                     (time_series; default [20 30 40])
%     .scan_sync      "unsynchronized" (default) | "synchronized". The
%                     synchronized option needs two mech_rotating radars
%                     with one equal rate, no random beam steering and no
%                     rate jitter on sys2 (sys2 follows sys1's drive).
%     .sync_offset_deg  synchronized phase offset az2 - az1 (deg, default
%                     0). Both azimuths use the common frame in which radar
%                     2 lies at 0 deg from radar 1: offset 0 points both
%                     antennas the same way, so a main beam always meets the
%                     other antenna's back lobe; offset 180 gives main beam
%                     to main beam once per revolution.
%
%   out struct:
%     .g12_dBi, .g21_dBi  coupling samples, column vectors
%     .cdf12, .cdf21      tables: gain_dBi, cdf, exceed_pct (0.1 dB grid)
%     .cdf_array12, .cdf_array21  N-by-2 double: column 1 is the percent
%                         of time (0 to 100, opts.cdf_pct), column 2 the MAG
%                         (dBi) not exceeded for that percent of time
%                         (prctile; 0 gives the minimum, 100 the maximum)
%     .pct12, .pct21      tables: pct, gain_dBi, ci95_lo_dBi, ci95_hi_dBi
%                         (95% confidence interval from order statistics),
%                         n_tail (samples above the
%                         percentile; < 10 means the tail is poorly sampled)
%     .geom               az/el each way, beat period
%                         360/|rate1-rate2|, joint_period_s (repeat period
%                         of the pointing pair when the rates are in a
%                         small-integer ratio, else Inf) and ratio [p q]
%     .events12           (time_series) per-threshold event statistics
%     .ts                 (time_series) trial 1: t_s, g12_dBi, g21_dBi
%     .sync               (synchronized) deterministic sweep of one
%                         revolution: max_g12_dBi, max_g21_dBi,
%                         az1_at_max12_deg, az1_at_max21_deg, offset_deg,
%                         step_deg; [] when unsynchronized
%     .opts               options used
%
%   See also antenna_prepare_rev1, antenna_gain_rev1, ecdf, prctile.

arguments
    sys1 (1,1) struct
    sys2 (1,1) struct
    opts (1,1) struct
end

opts = fill_opts(opts);
sys1 = prep_system(sys1, "sys1");
sys2 = prep_system(sys2, "sys2");

% Location-agnostic geometry: bearings along the line joining the radars,
% path elevations from opts.
geom = struct("az12_deg", 0, "el12_deg", opts.el12_deg, "az21_deg", 180, ...
    "el21_deg", opts.el21_deg, "beat_period_s", beat_period(sys1.scan, sys2.scan));
[geom.joint_period_s, geom.ratio] = joint_period(sys1.scan, sys2.scan);
is_sync = opts.scan_sync == "synchronized";
if is_sync
    check_sync_valid(sys1.scan, sys2.scan);
end

if opts.mode == "random"
    if ~is_sync
        check_random_mode_valid(sys1.scan, sys2.scan, geom);
    end
    [g12, g21] = run_random(sys1, sys2, geom, opts);
    events12 = [];
    ts = [];
else
    [g12, g21, events12, ts, dt_s] = run_time_series(sys1, sys2, geom, opts);
    opts.dt_s = dt_s;
end

out = struct();
out.g12_dBi = g12;
out.g21_dBi = g21;
out.cdf12 = cdf_table(g12);
out.cdf21 = cdf_table(g21);
out.pct12 = pct_table(g12, opts.pct);
out.pct21 = pct_table(g21, opts.pct);
out.cdf_array12 = cdf_array(g12, opts.cdf_pct);
out.cdf_array21 = cdf_array(g21, opts.cdf_pct);
out.geom = geom;
out.events12 = events12;
out.ts = ts;
out.sync = [];
if is_sync
    out.sync = sync_sweep(sys1, sys2, geom, opts);
end
out.opts = opts;
end

% =========================================================================
function s = sync_sweep(sys1, sys2, geom, opts)
% Deterministic sweep of one revolution of synchronized antennas: radar 1
% at azimuth x, radar 2 at x + offset. The step is 1/100 of the narrowest
% azimuth beamwidth, and the grid includes x = 0 (radar 1 pointing at
% radar 2) and x = 180 - offset (radar 2 pointing at radar 1).
step = min(sys1.min_az_bw_deg, sys2.min_az_bw_deg)/100;
x = unique([(0:step:360-step)'; 0; wrapTo360(180 - opts.sync_offset_deg)]);
st1 = sync_state(sys1.scan, x);
st2 = sync_state(sys2.scan, mod(x + opts.sync_offset_deg, 360));
[g12, g21] = coupling(sys1, sys2, geom, st1, st2);
[m12, i12] = max(g12);
[m21, i21] = max(g21);
s = struct("max_g12_dBi", m12, "max_g21_dBi", m21, "az1_at_max12_deg", x(i12), ...
    "az1_at_max21_deg", x(i21), "offset_deg", opts.sync_offset_deg, "step_deg", step);
end

% -------------------------------------------------------------------------
function st = sync_state(scan, az)
z = zeros(numel(az), 1);
st = struct("mech_az", az, "tilt", scan.tilt_deg + z, "steer_az", z, "steer_el", z);
if scan.beam_mode == "fixed"
    st.steer_az = scan.steer_az_deg + z;
    st.steer_el = scan.steer_el_deg + z;
end
end

% =========================================================================
function [g12, g21] = run_random(sys1, sys2, geom, opts)
n_chunks = ceil(opts.n_samples/opts.chunk_size);
c12 = cell(n_chunks, 1);
c21 = cell(n_chunks, 1);
if opts.parallel_flag == 1
    parfor k = 1:n_chunks
        [c12{k}, c21{k}] = random_chunk(k, sys1, sys2, geom, opts);
    end
else
    for k = 1:n_chunks
        [c12{k}, c21{k}] = random_chunk(k, sys1, sys2, geom, opts);
    end
end
g12 = vertcat(c12{:});
g21 = vertcat(c21{:});
end

% -------------------------------------------------------------------------
function [g12, g21] = random_chunk(k, sys1, sys2, geom, opts)
% One independent RandStream substream per chunk, so serial and parfor
% runs give identical samples (no rng(seed + k)).
rs = RandStream("mrg32k3a", "Seed", opts.seed);
rs.Substream = k;
n = min(opts.chunk_size, opts.n_samples - (k - 1)*opts.chunk_size);
st1 = draw_state_random(sys1.scan, n, rs);
st2 = draw_state_random(sys2.scan, n, rs);
if opts.scan_sync == "synchronized"
    st2.mech_az = mod(st1.mech_az + opts.sync_offset_deg, 360);
end
[g12, g21] = coupling(sys1, sys2, geom, st1, st2);
end

% -------------------------------------------------------------------------
function [g12, g21] = coupling(sys1, sys2, geom, st1, st2)
gt1 = antenna_gain_rev1(sys1.tx_p, geom.az12_deg, geom.el12_deg, st1);
gr2 = antenna_gain_rev1(sys2.rx_p, geom.az21_deg, geom.el21_deg, st2);
gt2 = antenna_gain_rev1(sys2.tx_p, geom.az21_deg, geom.el21_deg, st2);
gr1 = antenna_gain_rev1(sys1.rx_p, geom.az12_deg, geom.el12_deg, st1);
g12 = gt1 + gr2;
g21 = gt2 + gr1;
end

% -------------------------------------------------------------------------
function st = draw_state_random(scan, n, rs)
% Long-run (time-weighted) pointing state for one radar.
z = zeros(n, 1);
st = struct("mech_az", z, "tilt", z, "steer_az", z, "steer_el", z);
switch scan.type
    case "mech_rotating"
        st.mech_az = sector_az(scan, scan.sector_deg*rand(rs, n, 1));
        w = scan.step_dur_s/sum(scan.step_dur_s);    % time fraction per step
        if numel(w) > 1
            step = discretize(rand(rs, n, 1), [0; cumsum(w(:))]);
        else
            step = ones(n, 1);
        end
        st.tilt = scan.tilt_deg(step);
        st.tilt = st.tilt(:);
        [st.steer_az, st.steer_el] = draw_steer(scan, n, rs);
    case "array_multiface"
        beam_az = 360*rand(rs, n, 1);
        [st.mech_az, st.steer_az] = face_for_beam(scan, beam_az);
        st.steer_el = uniform(rs, n, scan.beam_el_lim_deg);
end
end

% -------------------------------------------------------------------------
function [s_az, s_el] = draw_steer(scan, n, rs)
switch scan.beam_mode
    case "none"
        s_az = zeros(n, 1);
        s_el = zeros(n, 1);
    case "fixed"
        s_az = scan.steer_az_deg + zeros(n, 1);
        s_el = scan.steer_el_deg + zeros(n, 1);
    case "random"
        s_az = uniform(rs, n, scan.beam_az_lim_deg);
        s_el = uniform(rs, n, scan.beam_el_lim_deg);
end
end

% -------------------------------------------------------------------------
function [face_az, steer_az] = face_for_beam(scan, beam_az)
% Active face is the one whose normal is nearest the beam azimuth.
spacing = 360/scan.n_faces;
idx = mod(round(wrapTo360(beam_az - scan.face_az0_deg)/spacing), scan.n_faces);
face_az = wrapTo360(scan.face_az0_deg + idx*spacing);
steer_az = wrapTo180(beam_az - face_az);
end

% =========================================================================
function [g12, g21, events12, ts, dt_s] = run_time_series(sys1, sys2, geom, opts)
dt_s = opts.dt_s;
if isnan(dt_s)
    dt_s = min(transit_time(sys1), transit_time(sys2))/10;
end
n_t = floor(opts.duration_s/dt_s) + 1;
if n_t*opts.n_trials > 1e8
    error("radar_coupling:tooManySteps", ...
        "time_series needs %.3g samples; shorten duration_s, raise dt_s or lower n_trials.", n_t*opts.n_trials);
end
t = (0:n_t-1)'*dt_s;
c12 = cell(opts.n_trials, 1);
c21 = cell(opts.n_trials, 1);
if opts.parallel_flag == 1
    parfor j = 1:opts.n_trials
        [c12{j}, c21{j}] = trial(j, t, sys1, sys2, geom, opts);
    end
else
    for j = 1:opts.n_trials
        [c12{j}, c21{j}] = trial(j, t, sys1, sys2, geom, opts);
    end
end
events12 = event_stats(c12, dt_s, opts.event_thresh_dBi);
ts = struct("t_s", t, "g12_dBi", c12{1}, "g21_dBi", c21{1});
g12 = vertcat(c12{:});
g21 = vertcat(c21{:});
end

% -------------------------------------------------------------------------
function [g12, g21] = trial(j, t, sys1, sys2, geom, opts)
rs = RandStream("mrg32k3a", "Seed", opts.seed);
rs.Substream = j;
st1 = state_time_series(sys1.scan, t, rs);
st2 = state_time_series(sys2.scan, t, rs);
if opts.scan_sync == "synchronized"
    % One drive: radar 2 follows radar 1 (including any sys1 rate jitter).
    st2.mech_az = mod(st1.mech_az + opts.sync_offset_deg, 360);
end
[g12, g21] = coupling(sys1, sys2, geom, st1, st2);
end

% -------------------------------------------------------------------------
function st = state_time_series(scan, t, rs)
n = numel(t);
z = zeros(n, 1);
st = struct("mech_az", z, "tilt", z, "steer_az", z, "steer_el", z);
switch scan.type
    case "mech_rotating"
        az0 = 360*rand(rs);
        if scan.rate_jitter_deg_s > 0
            % Per-trial rate offset within the drive tolerance, uniform in
            % +/- jitter (M.2069-0 Annex 1 cites +/-0.5 deg/s accuracy).
            dr = scan.rate_jitter_deg_s*(2*rand(rs) - 1);
            scan.rate_deg_s = sign(scan.rate_deg_s).*(abs(scan.rate_deg_s) + dr);
            scan.step_dur_s = scan.n_rev*scan.sector_deg./abs(scan.rate_deg_s);
        end
        t_cycle = sum(scan.step_dur_s);
        t0 = t_cycle*rand(rs);
        tau = mod(t + t0, t_cycle);
        edges = [0; cumsum(scan.step_dur_s(:))];
        step = discretize(tau, edges);
        step(isnan(step)) = numel(scan.step_dur_s);   % tau == t_cycle edge
        % Integer revolutions per step: azimuth restarts each step mod 360.
        if scan.sector_deg < 360
            % Sector: position within the sweep, from the leading edge.
            u = mod(az0/360*scan.sector_deg + abs(scan.rate_deg_s(step)).*(tau - edges(step)), scan.sector_deg);
            st.mech_az = sector_az(scan, u, sign(scan.rate_deg_s(step)));
        else
            st.mech_az = mod(az0 + scan.rate_deg_s(step).*(tau - edges(step)), 360);
        end
        st.tilt = scan.tilt_deg(step);
        st.tilt = st.tilt(:);
        [st.steer_az, st.steer_el] = steer_time_series(scan, t, rs);
    case "array_multiface"
        [didx, n_dw] = dwell_index(scan, t, rs);
        switch scan.beam_mode
            case "random"
                beam_az = 360*rand(rs, n_dw, 1);
            case "sequential"
                beam_az = mod(360*rand(rs) + scan.rate_deg_s*scan.dwell_s*(0:n_dw-1)', 360);
        end
        el = uniform(rs, n_dw, scan.beam_el_lim_deg);
        [st.mech_az, st.steer_az] = face_for_beam(scan, beam_az(didx));
        st.steer_el = el(didx);
end
end

% -------------------------------------------------------------------------
function [s_az, s_el] = steer_time_series(scan, t, rs)
n = numel(t);
switch scan.beam_mode
    case "none"
        s_az = zeros(n, 1);
        s_el = zeros(n, 1);
    case "fixed"
        s_az = scan.steer_az_deg + zeros(n, 1);
        s_el = scan.steer_el_deg + zeros(n, 1);
    case "random"
        [didx, n_dw] = dwell_index(scan, t, rs);
        a = uniform(rs, n_dw, scan.beam_az_lim_deg);
        e = uniform(rs, n_dw, scan.beam_el_lim_deg);
        s_az = a(didx);
        s_el = e(didx);
end
end

% -------------------------------------------------------------------------
function [didx, n_dw] = dwell_index(scan, t, rs)
off = scan.dwell_s*rand(rs);
didx = floor((t + off)/scan.dwell_s) + 1;
n_dw = didx(end);
end

% -------------------------------------------------------------------------
function ev = event_stats(c, dt_s, thresh)
% Runs of consecutive samples at or above each threshold. Runs cut by the
% trial boundaries count as events.
n_th = numel(thresh);
[n_events, per_hour, mean_dur, max_dur, frac_time] = deal(zeros(n_th, 1));
total_s = sum(cellfun(@numel, c))*dt_s;
for i = 1:n_th
    dur = [];
    above = 0;
    for j = 1:numel(c)
        x = c{j} >= thresh(i);
        d = diff([false; x; false]);
        len = find(d == -1) - find(d == 1);
        dur = [dur; len*dt_s]; %#ok<AGROW> one short vector per trial
        above = above + sum(x);
    end
    n_events(i) = numel(dur);
    per_hour(i) = n_events(i)/(total_s/3600);
    if ~isempty(dur)
        mean_dur(i) = mean(dur);
        max_dur(i) = max(dur);
    end
    frac_time(i) = above*dt_s/total_s;
end
ev = table(thresh(:), n_events, per_hour, mean_dur, max_dur, 100*frac_time, ...
    'VariableNames', ["thresh_dBi", "n_events", "events_per_hour", ...
    "mean_duration_s", "max_duration_s", "time_above_pct"]);
end

% =========================================================================
function T = cdf_table(g)
% histcounts on a 0.1 dB grid instead of ecdf, which returns one row per
% unique sample (millions of rows); the grid keeps the output compact.
lo = floor(min(g)*10)/10 - 0.05;
hi = ceil(max(g)*10)/10 + 0.05;
edges = lo:0.1:(hi + 0.1);
counts = histcounts(g, edges);
F = cumsum(counts(:))/numel(g);
T = table(edges(2:end)', F, 100*(1 - F), ...
    'VariableNames', ["gain_dBi", "cdf", "exceed_pct"]);
end

% -------------------------------------------------------------------------
function A = cdf_array(g, pct)
% CDF as a plain numeric array: [percent, MAG dBi], percent 0 to 100.
A = [pct(:), prctile(g, pct(:))];
end

% -------------------------------------------------------------------------
function T = pct_table(g, pct)
% Percentiles with distribution-free 95% confidence intervals from order
% statistics: the p-quantile lies between order statistics n*p -/+ 1.96 *
% sqrt(n*p*(1-p)) (normal approximation to the binomial count).
q = prctile(g, pct(:));
n = numel(g);
n_tail = n*(1 - pct(:)/100);
gs = sort(g(:));
p = pct(:)/100;
half = 1.96*sqrt(n*p.*(1 - p));
lo = gs(min(max(floor(n*p - half), 1), n));
hi = gs(min(max(ceil(n*p + half), 1), n));
T = table(pct(:), q(:), lo(:), hi(:), n_tail, 'VariableNames', ...
    ["pct", "gain_dBi", "ci95_lo_dBi", "ci95_hi_dBi", "n_tail"]);
end

% =========================================================================
function s = prep_system(s, label)
need = ["tx", "scan"];
missing = need(~isfield(s, need));
if ~isempty(missing)
    error("radar_coupling:missingField", "%s is missing field(s): %s.", label, strjoin(missing, ", "));
end
if ~isfield(s, "name") || isempty(s.name)
    s.name = label;
end
if ~isfield(s, "rx") || isempty(s.rx)
    s.rx = s.tx;
end
s.tx_p = antenna_prepare_rev1(s.tx);
s.rx_p = antenna_prepare_rev1(s.rx);
is_array = s.tx_p.type == "array_planar";
if is_array ~= (s.rx_p.type == "array_planar")
    error("radar_coupling:txRxMismatch", "%s: tx and rx must both be arrays or both apertures.", label);
end
s.scan = prep_scan(s.scan, is_array, label);
s.min_az_bw_deg = min(beamwidth_az(s.tx_p), beamwidth_az(s.rx_p));
end

% -------------------------------------------------------------------------
function scan = prep_scan(scan, is_array, label)
if ~isfield(scan, "type")
    error("radar_coupling:missingField", "%s.scan.type is required.", label);
end
scan.type = string(scan.type);
if isfield(scan, "rate_rpm") && ~isempty(scan.rate_rpm)
    scan.rate_deg_s = 6*scan.rate_rpm;
end
switch scan.type
    case "mech_rotating"
        req(scan, "rate_deg_s", label);
        scan.rate_deg_s = scan.rate_deg_s(:);
        if any(scan.rate_deg_s == 0)
            error("radar_coupling:zeroRate", "%s: rotation rate must be nonzero.", label);
        end
        n_step = numel(scan.rate_deg_s);
        scan.tilt_deg = expand(field_or(scan, "tilt_deg", 0), n_step, "tilt_deg", label);
        scan.n_rev = expand(field_or(scan, "n_rev", 1), n_step, "n_rev", label);
        if any(scan.n_rev ~= round(scan.n_rev) | scan.n_rev < 1)
            error("radar_coupling:badNRev", "%s: n_rev must be positive integers.", label);
        end
        scan.sector_deg = field_or(scan, "sector_deg", 360);
        scan.sector_center_deg = field_or(scan, "sector_center_deg", 0);
        if ~isscalar(scan.sector_deg) || scan.sector_deg <= 0 || scan.sector_deg > 360
            error("radar_coupling:badSector", "%s: sector_deg must be a scalar in (0, 360].", label);
        end
        scan.step_dur_s = scan.n_rev*scan.sector_deg./abs(scan.rate_deg_s);
        scan.rate_jitter_deg_s = field_or(scan, "rate_jitter_deg_s", 0);
        if scan.rate_jitter_deg_s < 0 || scan.rate_jitter_deg_s >= min(abs(scan.rate_deg_s))
            error("radar_coupling:badJitter", "%s: rate_jitter_deg_s must be >= 0 and below the slowest rate.", label);
        end
        if is_array
            scan.beam_mode = string(field_or(scan, "beam_mode", "fixed"));
            switch scan.beam_mode
                case "fixed"
                    scan.steer_az_deg = field_or(scan, "steer_az_deg", 0);
                    scan.steer_el_deg = field_or(scan, "steer_el_deg", 0);
                case "random"
                    scan = array_scan_defaults(scan);
                    scan.beam_az_lim_deg = field_or(scan, "beam_az_lim_deg", [0 0]);
                otherwise
                    error("radar_coupling:badBeamMode", "%s: beam_mode must be ""fixed"" or ""random"" for mech_rotating arrays.", label);
            end
        else
            scan.beam_mode = "none";
        end
    case "array_multiface"
        if ~is_array
            error("radar_coupling:scanAntennaMismatch", "%s: array_multiface needs an array_planar antenna.", label);
        end
        scan.n_faces = field_or(scan, "n_faces", 4);          % M.1464-2 Radar F1
        scan.face_az0_deg = field_or(scan, "face_az0_deg", 0);
        scan = array_scan_defaults(scan);
        scan.beam_mode = string(field_or(scan, "beam_mode", "random"));
        if scan.beam_mode == "sequential"
            req(scan, "rate_deg_s", label);
        elseif scan.beam_mode ~= "random"
            error("radar_coupling:badBeamMode", "%s: beam_mode must be ""random"" or ""sequential"" for array_multiface.", label);
        end
        scan.step_dur_s = scan.dwell_s;
    otherwise
        error("radar_coupling:badScanType", ...
            "%s: scan.type must be ""mech_rotating"" or ""array_multiface"", got ""%s"".", label, scan.type);
end
end

% -------------------------------------------------------------------------
function scan = array_scan_defaults(scan)
% Electronic beam defaults: horizon only, 50 ms dwell (within the "tens to
% hundreds of ms" of M.1464-2 Section 2.3).
scan.beam_el_lim_deg = field_or(scan, "beam_el_lim_deg", [0 0]);
scan.dwell_s = field_or(scan, "dwell_s", 0.05);
end

% -------------------------------------------------------------------------
function opts = fill_opts(opts)
if ~isfield(opts, "seed") || isempty(opts.seed)
    error("radar_coupling:missingSeed", "opts.seed is required for reproducibility.");
end
d = struct("mode", "random", "n_samples", 5e6, "duration_s", 3600, ...
    "dt_s", NaN, "n_trials", 4, "el12_deg", 0, "el21_deg", 0, "parallel_flag", 0, ...
    "chunk_size", 2.5e5, "pct", [50 90 99 99.9 99.99 99.999], ...
    "event_thresh_dBi", [20 30 40], "scan_sync", "unsynchronized", "sync_offset_deg", 0, ...
    "cdf_pct", unique([0:1:100, 99.5 99.9 99.95 99.99 99.995 99.999]));
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
        opts.(f{i}) = d.(f{i});
    end
end
opts.mode = string(opts.mode);
if ~ismember(opts.mode, ["random", "time_series"])
    error("radar_coupling:badMode", "opts.mode must be ""random"" or ""time_series"".");
end
if ~ismember(opts.parallel_flag, [0 1])
    error("radar_coupling:badParallelFlag", "opts.parallel_flag must be 0 or 1.");
end
opts.scan_sync = string(opts.scan_sync);
if ~ismember(opts.scan_sync, ["unsynchronized", "synchronized"])
    error("radar_coupling:badScanSync", "opts.scan_sync must be ""unsynchronized"" or ""synchronized"".");
end
if ~isnumeric(opts.cdf_pct) || isempty(opts.cdf_pct) || any(opts.cdf_pct < 0 | opts.cdf_pct > 100) ...
        || any(diff(opts.cdf_pct(:)) <= 0)
    error("radar_coupling:badCdfPct", "opts.cdf_pct must increase strictly within 0 to 100.");
end
if ~isscalar(opts.sync_offset_deg) || ~isfinite(opts.sync_offset_deg)
    error("radar_coupling:badSyncOffset", "opts.sync_offset_deg must be a finite scalar.");
end
end

% -------------------------------------------------------------------------
function check_sync_valid(scan1, scan2)
% Synchronized scanning: one mechanical rate shared by both antennas.
if scan1.type ~= "mech_rotating" || scan2.type ~= "mech_rotating" ...
        || ~isscalar(scan1.rate_deg_s) || ~isscalar(scan2.rate_deg_s) ...
        || ~ismembertol(scan1.rate_deg_s, scan2.rate_deg_s, 1e-9)
    error("radar_coupling:syncNeedsEqualRates", ...
        "Synchronized scanning needs two mech_rotating radars with one equal rate (same sign).");
end
if scan1.sector_deg < 360 || scan2.sector_deg < 360
    error("radar_coupling:syncSector", "Synchronized scanning needs full 360 deg rotation (no sector_deg).");
end
if scan1.beam_mode == "random" || scan2.beam_mode == "random"
    error("radar_coupling:syncRandomBeam", ...
        "Synchronized scanning does not support random beam steering; use beam_mode ""fixed"".");
end
if scan2.rate_jitter_deg_s > 0
    error("radar_coupling:syncJitter", ...
        "Synchronized radars share one drive; set rate jitter on sys1 only.");
end
end

% -------------------------------------------------------------------------
function check_random_mode_valid(scan1, scan2, geom)
% Unsynchronized random mode treats the two pointing angles as independent.
% Equal single rates keep a fixed but unknown phase offset (M.2069-0
% Annex 1); independent sampling then averages over that unknown offset,
% which is the MACD method. Any one pair of radars sees a single offset, so
% warn and point to the synchronized option. Rates in a small-integer
% ratio (for example 75 and 90 deg/s, 5:6) retrace the same closed path
% every joint period and visit only part of the azimuth-azimuth plane, so
% the random-mode CDF describes them only if the rates wander (drive
% tolerance) over the time of interest.
if scan1.type == "mech_rotating" && scan2.type == "mech_rotating" ...
        && isscalar(scan1.rate_deg_s) && isscalar(scan2.rate_deg_s) ...
        && ismembertol(scan1.rate_deg_s, scan2.rate_deg_s, 1e-9)
    warning("radar_coupling:equalRates", join(["Equal rotation rates keep a fixed phase offset. " ...
        "Unsynchronized random mode averages over an unknown offset (MACD). For a known " ...
        "offset, set opts.scan_sync = ""synchronized"" and opts.sync_offset_deg."], ""));
    return
end
if isfinite(geom.joint_period_s)
    warning("radar_coupling:commensurateRates", ...
        join(["Rotation rates are in a %d:%d ratio; the pointing pair repeats every %.1f s. " ...
        "Random mode assumes the rates drift apart. Check with opts.mode = ""time_series"" " ...
        "and scan.rate_jitter_deg_s."], ""), geom.ratio(1), geom.ratio(2), geom.joint_period_s);
end
end

% -------------------------------------------------------------------------
function az = sector_az(scan, u, sgn)
% Azimuth at distance u (0 to sector_deg) into the sector sweep. A full
% 360 deg sector is the ordinary rotation (u is the azimuth).
if nargin < 3
    sgn = 1;
end
if scan.sector_deg >= 360
    az = mod(u, 360);
else
    az = mod(scan.sector_center_deg + sgn.*(u - scan.sector_deg/2), 360);
end
end

% -------------------------------------------------------------------------
function [tj, ratio] = joint_period(scan1, scan2)
% Repeat period of the joint pointing (az1, az2) when the two single
% rates are in a ratio p:q with q <= 100; Inf otherwise.
tj = Inf;
ratio = [NaN NaN];
if scan1.type == "mech_rotating" && scan2.type == "mech_rotating" ...
        && isscalar(scan1.rate_deg_s) && isscalar(scan2.rate_deg_s)
    r1 = abs(scan1.rate_deg_s);
    r2 = abs(scan2.rate_deg_s);
    [p, q] = rat(r1/r2, 1e-9*r1/r2);        % r1/r2 = p/q
    if q <= 100 && p <= 100
        tj = 360*p/r1;                      % radar 1 makes p turns, radar 2 makes q
        ratio = [p q];
    end
end
end

% -------------------------------------------------------------------------
function tb = beat_period(scan1, scan2)
% Repeat period of the relative pointing, 360/|rate1 - rate2|
% (M.2069-0 Annex 1), for single-rate mechanical scans only.
tb = NaN;
if scan1.type == "mech_rotating" && scan2.type == "mech_rotating" ...
        && isscalar(scan1.rate_deg_s) && isscalar(scan2.rate_deg_s)
    dr = abs(scan1.rate_deg_s - scan2.rate_deg_s);
    if dr > 0
        tb = 360/dr;
    else
        tb = Inf;
    end
end
end

% -------------------------------------------------------------------------
function tt = transit_time(s)
% Shortest time the beam takes to sweep across its own 3 dB width, or the
% dwell time for electronically scanned beams.
scan = s.scan;
tt = Inf;
if scan.type == "mech_rotating"
    tt = s.min_az_bw_deg/max(abs(scan.rate_deg_s));
end
if isfield(scan, "dwell_s") && scan.beam_mode ~= "fixed" && scan.beam_mode ~= "none"
    tt = min(tt, scan.dwell_s);
end
end

% -------------------------------------------------------------------------
function bw = beamwidth_az(ant_p)
switch ant_p.type
    case "statgain"
        % Implied 3 dB width of the Table 4-1 main beam:
        % 4e-4*10^(G/10)*(bw/2)^2 = 3.
        bw = 2*sqrt(3/(4e-4*10^(ant_p.g_max_dBi/10)));
    case "f1245"
        % F.1245 main beam: 2.5e-3*(D/lambda*bw/2)^2 = 3, capped at phi_m.
        bw = 2*min(sqrt(3/2.5e-3)/ant_p.f.dl, ant_p.f.phi_m);
    case {"rect", "parabolic", "offaxis_table"}
        bw = ant_p.az_bw_deg;
    case "measured"
        bw = field_or(ant_p, "az_bw_deg", 1);
        if ~isfield(ant_p, "az_bw_deg")
            warning("radar_coupling:measuredBw", "Measured antenna has no az_bw_deg; time step assumes 1 deg.");
        end
    case "array_planar"
        % Broadside 3 dB width of a uniform line, ~0.886*lambda/(N*d) rad;
        % used only to pick the time step.
        bw = rad2deg(0.886*ant_p.lambda_m/(ant_p.nx*ant_p.dx_m));
end
end

% -------------------------------------------------------------------------
function x = uniform(rs, n, lim)
x = lim(1) + (lim(2) - lim(1))*rand(rs, n, 1);
end

function v = field_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function v = expand(v, n, name, label)
v = v(:);
if isscalar(v)
    v = repmat(v, n, 1);
elseif numel(v) ~= n
    error("radar_coupling:stepMismatch", "%s.scan.%s must be scalar or match rate_deg_s.", label, name);
end
end

function req(s, names, label)
missing = names(~isfield(s, names));
if ~isempty(missing)
    error("radar_coupling:missingField", "%s.scan is missing field(s): %s.", label, strjoin(missing, ", "));
end
end
