function T = gain_dist_table_rev1(ant, mode_scans, opts)
%GAIN_DIST_TABLE_REV1 Exact time-weighted gain distribution of one scanning antenna.
%   T = gain_dist_table_rev1(ant, mode_scans, opts) returns the distribution
%   of the gain of one rotating radar antenna toward another station, as a
%   function of the elevation of that station, for every scan mode in
%   mode_scans. Run it once per radar system and save T. Then
%   mag_pair_cdf_rev1 combines two tables into the mutual antenna gain
%   (MAG) CDF of any pair from the two takeoff angles alone, with no Monte
%   Carlo run per pair.
%
%   Why this is exact: in the M.2069-0 Annex 1 random mode (the method of
%   radar_mutual_gain_cdf_rev1 with opts.mode = "random") the two antennas
%   point independently, so MAG = G_tx + G_rx is the sum of two independent
%   gains and its distribution is the convolution of the two single-antenna
%   distributions (the MACD of the 2.7 GHz working paper, Appendix B; the
%   continuous_macd check in crosscheck_csharp_rev1). Each single-antenna
%   distribution depends only on that antenna's pattern, its scan, and the
%   elevation of the other station.
%
%   Method: for a mechanically rotating antenna the time-weighted pointing
%   is uniform in azimuth within each elevation step, with step weights
%   proportional to the step durations (draw_state_random in
%   radar_mutual_gain_cdf_rev1). This function replaces the random draw
%   with a midpoint azimuth grid, so the distribution has no sampling noise:
%     1. For each unique beam tilt and each grid elevation, evaluate
%        antenna_gain_rev1 on n = 360/az_step_deg azimuths and histogram the
%        gains on a bin_dB grid. Gains are rounded UP to the next bin edge,
%        so every percentile is conservative by at most bin_dB.
%     2. Scan modes differ only in how long they dwell at each tilt, so a
%        mode is a weight vector over the tilts (T.W, modes by tilts). Any
%        mode's distribution is T.W(m, :) applied to the tilt
%        distributions. All modes cost the same as one.
%
%   ant         antenna model struct (antenna_prepare_rev1 input), used for
%               transmit and receive unless opts.ant_rx is given
%   mode_scans  scan struct, or cell array of scan structs, one per scan
%               mode (for example nexrad_mode_scans_rev1). Each scan must
%               be "mech_rotating" over 360 deg (sector scans make the
%               distribution depend on the bearing; not supported).
%   opts        struct, all fields optional:
%     .name          system name (default "radar")
%     .mode_names    string array, one per scan (default "mode 1", ...)
%     .ant_rx        receive antenna struct when it differs from ant
%     .el_grid_deg   elevations of the other station, deg above the local
%                    horizontal (default -5:0.05:10)
%     .bin_dB        gain bin width, dB (default 0.05)
%     .az_step_deg   azimuth grid step, deg (default 0.001, about 1/1000
%                    of the WSR-88D beamwidth)
%     .parallel_flag 1 runs the elevations in parfor (default 0)
%
%   T struct:
%     .name, .mode_names (1-by-M), .tilt_deg (K-by-1 unique tilts)
%     .el_grid_deg (E-by-1), .bin_dB, .az_step_deg
%     .k0           integer index of the first gain bin: bin b holds gain
%                   (k0 + b - 1)*bin_dB, the upper edge of its interval
%     .P_tx         single, B-by-E-by-K probability mass per bin
%     .P_rx         same for the receive pattern ([] when ant_rx is absent)
%     .gmax_tx      E-by-K maximum gain toward the station (dBi, exact)
%     .gmax_rx      same for receive ([] when ant_rx is absent)
%     .W            M-by-K time weight of each tilt in each mode (rows sum
%                   to 1)
%     .mode_rates   1-by-M cell of each mode's per-step azimuth rates
%                   (deg/s), for the synchronization check in
%                   mag_pair_cdf_rev1
%     .ant, .ant_rx, .scans  the inputs, for traceability
%     .notes        string array of method notes
%     .created      creation time
%
%   Example:
%     o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, ...
%         "duration_s", 60, "n_trials", 1));      % builds the Rev 5 patterns
%     [scans, names] = nexrad_mode_scans_rev1();
%     Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scans, ...
%         struct("name", "WSR-88D", "mode_names", names));
%     Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, ...
%         struct("name", "ASR-11"));
%     R = mag_pair_cdf_rev1(Ta, Tn, 2.0, 1.7);
%
%   See also mag_pair_cdf_rev1, nexrad_mode_scans_rev1, antenna_gain_rev1,
%   radar_mutual_gain_cdf_rev1.

arguments
    ant (1, 1) struct
    mode_scans
    opts (1, 1) struct = struct()
end

opts = fill_opts(opts);
if isstruct(mode_scans)
    mode_scans = num2cell(mode_scans(:)');
end
if ~iscell(mode_scans) || isempty(mode_scans)
    error("gain_dist:badScans", "mode_scans must be a scan struct or a nonempty cell array of scan structs.");
end
n_mode = numel(mode_scans);
if isempty(opts.mode_names)
    opts.mode_names = "mode " + string(1:n_mode);
    if n_mode == 1
        opts.mode_names = string(opts.name);
    end
end
opts.mode_names = string(opts.mode_names(:)');
if numel(opts.mode_names) ~= n_mode
    error("gain_dist:badNames", "opts.mode_names must have one name per scan mode (%d).", n_mode);
end

% ---- Scan modes -> unique tilts and time weights
[steps, mode_rates] = deal(cell(1, n_mode));
for m = 1:n_mode
    [steps{m}, mode_rates{m}] = scan_steps(mode_scans{m}, opts.mode_names(m));
end
all_tilt = cell2mat(cellfun(@(s) s.tilt_deg, steps, "UniformOutput", false)');
% uniquetol: tilts decoded from different VCP tables can differ by
% round-off; 1e-6 deg is far below any physical tilt difference.
tilt = uniquetol(all_tilt, 1e-6, "DataScale", 1);
n_tilt = numel(tilt);
W = zeros(n_mode, n_tilt);
for m = 1:n_mode
    [~, it] = ismembertol(steps{m}.tilt_deg, tilt, 1e-6, "DataScale", 1);
    W(m, :) = accumarray(it(:), steps{m}.dur_s(:), [n_tilt 1])';
    W(m, :) = W(m, :)/sum(W(m, :));
end

% ---- Antennas
pa_tx = antenna_prepare_rev1(ant);
check_antenna(pa_tx, "ant");
has_rx = ~isempty(opts.ant_rx);
pa_rx = pa_tx;                                  % unused unless has_rx (defined for parfor)
if has_rx
    pa_rx = antenna_prepare_rev1(opts.ant_rx);
    check_antenna(pa_rx, "opts.ant_rx");
end

% ---- Gain histograms
el = opts.el_grid_deg(:);
n_el = numel(el);
n_az = round(360/opts.az_step_deg);
az = ((1:n_az)' - 0.5)*360/n_az;                % midpoint rule: uniform time weight
[H_tx, H_rx] = deal(cell(n_el, 1));
[gmax_tx, gmax_rx] = deal(NaN(n_el, n_tilt));
dg = opts.bin_dB;
if opts.parallel_flag
    parfor ie = 1:n_el
        [H_tx{ie}, gmax_tx(ie, :)] = el_hist(pa_tx, el(ie), tilt, az, dg);
        if has_rx
            [H_rx{ie}, gmax_rx(ie, :)] = el_hist(pa_rx, el(ie), tilt, az, dg);
        end
    end
else
    for ie = 1:n_el
        [H_tx{ie}, gmax_tx(ie, :)] = el_hist(pa_tx, el(ie), tilt, az, dg);
        if has_rx
            [H_rx{ie}, gmax_rx(ie, :)] = el_hist(pa_rx, el(ie), tilt, az, dg);
        end
    end
end
H = H_tx;
if has_rx
    H = [H_tx; H_rx];
end
k_lo = min(cellfun(@(h) h.k_lo, H));
k_hi = max(cellfun(@(h) h.k_lo + size(h.p, 1) - 1, H));
n_bin = k_hi - k_lo + 1;
P_tx = stack(H_tx, k_lo, n_bin, n_tilt);
P_rx = [];
if has_rx
    P_rx = stack(H_rx, k_lo, n_bin, n_tilt);
else
    gmax_rx = [];
end

T = struct();
T.name = string(opts.name);
T.mode_names = opts.mode_names;
T.tilt_deg = tilt;
T.el_grid_deg = el;
T.bin_dB = opts.bin_dB;
T.az_step_deg = 360/n_az;
T.k0 = k_lo;
T.P_tx = P_tx;
T.P_rx = P_rx;
T.gmax_tx = gmax_tx;
T.gmax_rx = gmax_rx;
T.W = W;
T.mode_rates = mode_rates;
T.ant = ant;
T.ant_rx = opts.ant_rx;
T.scans = mode_scans;
T.notes = [
    "M.2069-0 Annex 1 random mode: time-weighted pointing, uniform in azimuth within each elevation step"
    "Deterministic midpoint azimuth grid, " + compose("%.6g", T.az_step_deg) + " deg (no sampling noise)"
    "Gains rounded up to the next " + compose("%.6g", T.bin_dB) + " dB bin edge (conservative by at most one bin)"
    "Full 360 deg rotation only; polarization mismatch not modeled (as in radar_mutual_gain_cdf_rev1)"];
T.created = datetime("now");
end

% =========================================================================
function [h, gmax] = el_hist(pa, el, tilt, az, dg)
% Gain histograms at one elevation for every tilt (bins are columns of h.p).
n_tilt = numel(tilt);
z = zeros(size(az));
k = cell(n_tilt, 1);
gmax = zeros(1, n_tilt);
for it = 1:n_tilt
    st = struct("mech_az", az, "tilt", z + tilt(it), "steer_az", z, "steer_el", z);
    g = antenna_gain_rev1(pa, 0, el, st);
    % Round up to the bin edge; the 1e-9 bin tolerance keeps gains that
    % sit on an edge (round-off) from moving up a full bin.
    k{it} = ceil(g/dg - 1e-9);
    gmax(it) = max(g);
end
k_lo = min(cellfun(@min, k));
k_hi = max(cellfun(@max, k));
p = zeros(k_hi - k_lo + 1, n_tilt);
for it = 1:n_tilt
    p(:, it) = accumarray(k{it} - k_lo + 1, 1, [k_hi - k_lo + 1, 1])/numel(az);
end
h = struct("k_lo", k_lo, "p", p);
end

% -------------------------------------------------------------------------
function P = stack(H, k_lo, n_bin, n_tilt)
% Place each elevation's histograms on the common bin grid.
n_el = numel(H);
P = zeros(n_bin, n_el, n_tilt, "single");
for ie = 1:n_el
    r = H{ie}.k_lo - k_lo + (1:size(H{ie}.p, 1));
    P(r, ie, :) = reshape(single(H{ie}.p), [], 1, n_tilt);
end
end

% -------------------------------------------------------------------------
function [s, rates] = scan_steps(scan, name)
% Elevation steps of one rotating scan and their durations (the step_dur_s
% rule of radar_mutual_gain_cdf_rev1 prep_scan: n_rev*360/|rate|).
if ~isfield(scan, "type") || string(scan.type) ~= "mech_rotating"
    error("gain_dist:badScanType", "%s: scan.type must be ""mech_rotating"".", name);
end
if isfield(scan, "rate_rpm") && ~isempty(scan.rate_rpm)
    rates = 6*scan.rate_rpm(:);
elseif isfield(scan, "rate_deg_s")
    rates = scan.rate_deg_s(:);
else
    error("gain_dist:missingRate", "%s: scan needs rate_deg_s or rate_rpm.", name);
end
if any(rates == 0)
    error("gain_dist:zeroRate", "%s: rotation rate must be nonzero.", name);
end
n = numel(rates);
tilt = field_or(scan, "tilt_deg", 0);
n_rev = field_or(scan, "n_rev", 1);
if isscalar(tilt), tilt = tilt + zeros(n, 1); end
if isscalar(n_rev), n_rev = n_rev + zeros(n, 1); end
if numel(tilt) ~= n || numel(n_rev) ~= n
    error("gain_dist:badSteps", "%s: tilt_deg and n_rev must have one value per rate step.", name);
end
if field_or(scan, "sector_deg", 360) < 360
    error("gain_dist:sector", "%s: sector scans are not supported (the distribution depends on the bearing).", name);
end
if string(field_or(scan, "beam_mode", "none")) ~= "none"
    error("gain_dist:steering", "%s: electronic beam steering is not supported.", name);
end
s = struct("tilt_deg", tilt(:), "dur_s", n_rev(:)*360./abs(rates));
end

% -------------------------------------------------------------------------
function check_antenna(pa, label)
if pa.type == "array_planar"
    error("gain_dist:array", "%s: planar arrays are not supported (beam steering).", label);
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
d = struct("name", "radar", "mode_names", strings(1, 0), "ant_rx", [], ...
    "el_grid_deg", -5:0.05:10, "bin_dB", 0.05, "az_step_deg", 0.001, "parallel_flag", 0);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i})
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("gain_dist:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
if ~isnumeric(opts.el_grid_deg) || numel(opts.el_grid_deg) < 2 || any(diff(opts.el_grid_deg(:)) <= 0)
    error("gain_dist:badGrid", "opts.el_grid_deg must be an increasing vector with at least 2 values.");
end
if any(abs(opts.el_grid_deg) > 90)
    error("gain_dist:badGrid", "opts.el_grid_deg values must lie in [-90, 90] deg.");
end
if ~isscalar(opts.bin_dB) || opts.bin_dB <= 0 || opts.bin_dB > 1
    error("gain_dist:badBin", "opts.bin_dB must be a scalar in (0, 1] dB.");
end
if ~isscalar(opts.az_step_deg) || opts.az_step_deg <= 0 || opts.az_step_deg > 0.1
    error("gain_dist:badAzStep", "opts.az_step_deg must be a scalar in (0, 0.1] deg.");
end
end
