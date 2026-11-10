function R = mag_pair_timeseries_rev1(lib, ant_n, ant_a, el_n_deg, el_a_deg, opts)
%MAG_PAIR_TIMESERIES_REV1 Per-scan MAG maxima and harm metric for one radar pair.
%   R = mag_pair_timeseries_rev1(lib, ant_n, ant_a, el_n_deg, el_a_deg,
%   opts) returns the clock-based outputs of the WSR-88D / ATC radar time
%   series for one pair: the maximum mutual antenna gain (MAG) in each ATC
%   revolution, NEXRAD cut and NEXRAD volume, the persistence of high
%   coupling, and a harm metric: the number of NEXRAD radials at the lowest
%   cut (0.5 deg) that see the MAG at or above each threshold.
%
%   It evaluates the gains only in the event windows of lib
%   (scan_event_library_rev1), a subset of the fixed time grid of
%   nexrad_asr11_mag_rev1. Outside the windows both beams are off the path
%   and the MAG is at most R.floor21_dBi (ts_floor_rev1), so every maximum
%   above the floor equals the fixed-grid maximum exactly. A scan with no
%   window sample is reported as NaN (censored: at or below the floor).
%
%   lib         event library (scan_event_library_rev1); one per pair of
%               scan modes, shared by every radar pair with those modes
%   ant_n       WSR-88D antenna struct (antenna_prepare_rev1 input), tx = rx
%   ant_a       ATC radar antenna struct, tx = rx
%   el_n_deg    takeoff angle at the NEXRAD (elevation of the ATC radar)
%   el_a_deg    takeoff angle at the ATC radar (elevation of the NEXRAD)
%   opts        struct, all fields optional:
%     .persist_thr_dBi  thresholds ([40 45 50 55 60 65 70]); each must be
%                       above the floor
%     .radial_el_deg    cut of the harm metric (0.5)
%     .radial_deg       radial width (1, the WSR-88D legacy radial;
%                       0.5 for super-resolution)
%     .radial_thr_dBi   harm-metric thresholds ([40 50 60])
%
%   R struct:
%     .asr_rev, .nexrad_cut, .nexrad_volume  tables: trial, scan number,
%                     max_g12_dBi, max_g21_dBi (NaN = at or below the
%                     floor), complete
%     .persistence, .persistence_volume  per threshold: share of complete
%                     scans reaching it, number of runs, longest run (scans
%                     and seconds)
%     .radial_cuts    table, one row per 0.5 deg cut: trial, volume, cut,
%                     complete, and n_radials_ge_<thr> per threshold
%     .radial_summary per threshold: share of complete 0.5 deg cuts with at
%                     least one radial at or above it, mean and largest
%                     number of such radials per cut
%     .floor21_dBi, .floor12_dBi  MAG bound outside the windows
%     .share          share of the fixed grid evaluated
%     .notes
%
%   See also scan_event_library_rev1, scan_event_windows_rev1,
%   ts_floor_rev1, nexrad_asr11_mag_rev1.

arguments
    lib (1, 1) struct
    ant_n (1, 1) struct
    ant_a (1, 1) struct
    el_n_deg (1, 1) double {mustBeReal, mustBeFinite}
    el_a_deg (1, 1) double {mustBeReal, mustBeFinite}
    opts (1, 1) struct = struct()
end
opts = fill_opts(opts);
spec = lib.spec;
pn = struct("tx", antenna_prepare_rev1(ant_n), "rx", antenna_prepare_rev1(ant_n));
pa = struct("tx", antenna_prepare_rev1(ant_a), "rx", antenna_prepare_rev1(ant_a));
[f21, f12] = ts_floor_rev1(pn, pa, spec.el_n, el_n_deg, spec.tilt_a, el_a_deg, ...
    lib.delta_eff_n_deg, lib.delta_eff_a_deg);
if min([opts.persist_thr_dBi(:); opts.radial_thr_dBi(:)]) <= max(f21, f12)
    error("mag_pair_ts:thresholdBelowFloor", join(["A threshold (%.1f dBi) is at or below the event " ...
        "floor (%.1f dBi); widen the windows (delta_n_deg, delta_a_deg) or raise the threshold."], ""), ...
        min([opts.persist_thr_dBi(:); opts.radial_thr_dBi(:)]), max(f21, f12));
end
t_full = (0:lib.dt_s:lib.duration_s)';
n_cut = numel(spec.rate_n);
n_tr = numel(lib.trials);
[rev_tab, cut_tab, vol_tab, rad_tab] = deal(cell(n_tr, 1));
n_eval = 0;
for jj = 1:n_tr
    j = lib.trials(jj);
    ru = double(lib.runs{jj});
    idx = cell2mat(arrayfun(@(a, b) (a:b)', ru(:, 1), ru(:, 2), "UniformOutput", false));
    t = t_full(idx);
    n_eval = n_eval + numel(t);
    C = scan_clock_rev1(spec, lib.seed, j, lib.duration_s);
    [st_n, st_a, ck] = scan_clock_state_rev1(C, t);
    [~, ~, ends] = scan_clock_state_rev1(C, [0; t_full(end)]);
    gN = antenna_gain_rev1(pn.rx, 0, el_n_deg, st_n);
    gA = antenna_gain_rev1(pa.tx, 180, el_a_deg, st_a);
    g21 = gA + gN;
    g12 = g21;                                  % one pattern each for tx and rx

    % ---- per ATC revolution (every revolution of the trial, censored ones
    %      as NaN)
    kr = ck.rev - ends.rev(1) + 1;
    n_r = ends.rev(2) - ends.rev(1) + 1;
    rev_tab{jj} = scan_rows(j, kr, n_r, g12, g21, "rev");

    % ---- per NEXRAD cut and volume
    key = (ck.vol - ends.vol(1))*n_cut + ck.cut;
    key0 = ends.cut(1);
    key1 = (ends.vol(2) - ends.vol(1))*n_cut + ends.cut(2);
    kc = key - key0 + 1;
    n_c = key1 - key0 + 1;
    T = scan_rows(j, kc, n_c, g12, g21, "cut");
    keys = (key0:key1)';
    vrel = floor((keys - 1)/n_cut);
    T.volume = vrel + 1;
    T.cut = keys - vrel*n_cut;
    T.el_deg = spec.el_n(T.cut);
    cut_tab{jj} = movevars(T, ["volume", "cut", "el_deg"], "After", "trial");
    kv = ck.vol - ends.vol(1) + 1;
    n_v = ends.vol(2) - ends.vol(1) + 1;
    vol_tab{jj} = scan_rows(j, kv, n_v, g12, g21, "volume");

    % ---- harm metric: radials of the lowest cut at or above each threshold
    low = abs(spec.el_n(T.cut) - opts.radial_el_deg) < 1e-6;
    in_low = abs(spec.el_n(ck.cut) - opts.radial_el_deg) < 1e-6;
    radial = floor(mod(st_n.mech_az, 360)/opts.radial_deg) + 1;
    n_rad = ceil(360/opts.radial_deg);
    counts = zeros(n_c, numel(opts.radial_thr_dBi));
    for i = 1:numel(opts.radial_thr_dBi)
        hit = in_low & g21 >= opts.radial_thr_dBi(i);
        if any(hit)
            u = unique((kc(hit) - 1)*n_rad + radial(hit));     % one count per (cut, radial)
            counts(:, i) = accumarray(floor((u - 1)/n_rad) + 1, 1, [n_c 1]);
        end
    end
    Rt = T(low, ["trial", "volume", "cut", "complete"]);
    Rt = [Rt, array2table(counts(low, :), "VariableNames", "n_radials_ge_" + string(opts.radial_thr_dBi))]; %#ok<AGROW>
    rad_tab{jj} = Rt;
end
R = struct();
R.asr_rev = vertcat(rev_tab{:});
R.nexrad_cut = vertcat(cut_tab{:});
R.nexrad_volume = vertcat(vol_tab{:});
t_rev = spec.sec_a/spec.rate_a;
dur_n = spec.sec_n./spec.rate_n + spec.turn_n;
R.persistence = persistence(R.asr_rev, opts.persist_thr_dBi, t_rev, "rev");
R.persistence_volume = persistence(R.nexrad_volume, opts.persist_thr_dBi, sum(dur_n), "vol");
R.radial_cuts = vertcat(rad_tab{:});
Rc = R.radial_cuts(R.radial_cuts.complete, :);
thr = opts.radial_thr_dBi(:);
[share_hit, mean_n, max_n] = deal(zeros(numel(thr), 1));
for i = 1:numel(thr)
    v = Rc.("n_radials_ge_" + string(thr(i)));
    share_hit(i) = 100*mean(v > 0);
    mean_n(i) = mean(v);
    max_n(i) = max([v; 0]);
end
R.radial_summary = table(thr, share_hit, mean_n, max_n, 'VariableNames', ["threshold_dBi", ...
    "cuts_with_hit_pct", "mean_radials_per_cut", "max_radials_per_cut"]);
R.radial_summary.Properties.Description = sprintf("%d complete %.1f deg cuts; radial width %g deg", ...
    height(Rc), opts.radial_el_deg, opts.radial_deg);
R.floor21_dBi = f21;
R.floor12_dBi = f12;
R.share = n_eval/(numel(t_full)*n_tr);
R.notes = [
    "Gains evaluated only in event windows (" + compose("%.1f", 100*R.share) + "% of the fixed grid); maxima above the floor equal fixed-grid maxima"
    "NaN maxima are censored: at or below the floor (" + compose("%.1f", max(f21, f12)) + " dBi)"
    "Harm metric: radials of the " + compose("%g", opts.radial_el_deg) + " deg cut whose MAG reaches each threshold while the NEXRAD beam points there (g21: ATC transmits, NEXRAD receives)"];
end

% =========================================================================
function T = scan_rows(j, k, n, g12, g21, name)
m12 = accumarray(k, g12, [n 1], @max, NaN);
m21 = accumarray(k, g21, [n 1], @max, NaN);
c = true(n, 1);
c([1 n]) = false;                               % first and last scans are cut by the trial
T = table(j + zeros(n, 1), (1:n)', m12, m21, c, 'VariableNames', ...
    ["trial", name, "max_g12_dBi", "max_g21_dBi", "complete"]);
end

% -------------------------------------------------------------------------
function T = persistence(Rv, thr, t_scan, unit)
% Runs of consecutive complete scans whose maximum reaches each threshold
% (the persistence rule of nexrad_asr11_mag_rev1). Censored (NaN) maxima
% are below every threshold.
Rc = Rv(Rv.complete, :);
n_thr = numel(thr);
[pct, n_runs, max_run] = deal(zeros(n_thr, 1));
hours = height(Rc)*t_scan/3600;
for i = 1:n_thr
    len = [];
    for j = unique(Rc.trial)'
        hit = Rc.max_g21_dBi(Rc.trial == j) >= thr(i);
        b = find(diff(double(hit)) ~= 0);
        f = [1; b + 1];
        l = [b; numel(hit)];
        k = hit(f);
        len = [len; l(k) - f(k) + 1]; %#ok<AGROW>
    end
    pct(i) = 100*mean(Rc.max_g21_dBi >= thr(i));
    n_runs(i) = numel(len);
    max_run(i) = max([len; 0]);
end
T = table(thr(:), pct, n_runs, n_runs/hours, max_run, max_run*t_scan, 'VariableNames', ...
    ["threshold_dBi", unit + "_pct", "n_runs", "runs_per_hour", "max_run_" + unit + "s", "max_run_s"]);
T.Properties.Description = sprintf("%.1f h of complete scans (scan period %.3f s)", hours, t_scan);
end

% -------------------------------------------------------------------------
function opts = fill_opts(opts)
d = struct("persist_thr_dBi", [40 45 50 55 60 65 70], "radial_el_deg", 0.5, "radial_deg", 1, ...
    "radial_thr_dBi", [40 50 60]);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("mag_pair_ts:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
end
