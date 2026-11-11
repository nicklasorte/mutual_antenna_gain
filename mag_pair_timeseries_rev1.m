function R = mag_pair_timeseries_rev1(lib, ant_n, ant_a, el_n_deg, el_a_deg, opts)
%MAG_PAIR_TIMESERIES_REV1 Per-scan MAG maxima and harm metric for radar pairs.
%   R = mag_pair_timeseries_rev1(lib, ant_n, ant_a, el_n_deg, el_a_deg,
%   opts) returns the clock-based outputs of the WSR-88D / ATC radar time
%   series for one pair, or for a batch of pairs that share the two scan
%   modes: the maximum mutual antenna gain (MAG) in each ATC revolution,
%   NEXRAD cut and NEXRAD volume, the persistence of high coupling, the
%   exceedance of the per-scan maxima, and a harm metric: the number of
%   NEXRAD radials at the lowest cut (0.5 deg) that see the MAG at or above
%   each threshold.
%
%   It evaluates the gains only in the event windows of lib
%   (scan_event_library_rev1), a subset of the fixed time grid of
%   nexrad_asr11_mag_rev1. Outside the windows both beams are off the path
%   and the MAG is at most R.floor21_dBi (ts_floor_rev1), so every maximum
%   above the floor equals the fixed-grid maximum exactly. A scan with no
%   window sample is reported as NaN (censored: at or below the floor).
%
%   Batches: el_n_deg and el_a_deg can be vectors (one element per pair;
%   a scalar is repeated). The scan clock and pointing states of each
%   trial are computed once for the whole batch, and the gains of all
%   pairs in a chunk are evaluated in one call (antenna_gain_rev1 with a
%   row of takeoff angles), so a pair in a batch costs far less than a
%   single-pair call. Each pair's result equals its single-pair result.
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
%     .exceed_step_dB   gain step of the exceedance tables (0.1)
%     .keep_scans       true returns the per-scan tables (default: true
%                       for one pair, false for a batch, to save memory;
%                       the summaries are always returned)
%     .batch_size       pairs per gain call (64)
%     .parallel_flag    1 runs the chunks of pairs in parfor (0)
%
%   R struct (struct array, one element per pair, for a batch):
%     .el_n_deg, .el_a_deg  the pair's takeoff angles
%     .asr_rev, .nexrad_cut, .nexrad_volume  tables: trial, scan number,
%                     max_g12_dBi, max_g21_dBi (NaN = at or below the
%                     floor), complete (empty when keep_scans is false)
%     .persistence, .persistence_volume  per threshold: share of complete
%                     scans reaching it, number of runs, longest run (scans
%                     and seconds). A run cut by the end of a trial counts
%                     to the trial end; .censored_runs says how many runs
%                     touch a trial end, so long runs can be read as lower
%                     bounds (red-team finding P2)
%     .exceedance     table: gain_dBi, rev_pct, cut_pct, vol_pct = share
%                     of complete revolutions, cuts and volumes whose
%                     maximum is at or above gain_dBi, from the floor up
%                     (red-team finding P3: curves, not single thresholds)
%     .radial_cuts    table, one row per 0.5 deg cut: trial, volume, cut,
%                     complete, and n_radials_ge_<thr> per threshold
%                     (empty when keep_scans is false)
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
    el_n_deg (:, 1) double {mustBeReal, mustBeFinite}
    el_a_deg (:, 1) double {mustBeReal, mustBeFinite}
    opts (1, 1) struct = struct()
end
n_p = max(numel(el_n_deg), numel(el_a_deg));
if isscalar(el_n_deg), el_n_deg = repmat(el_n_deg, n_p, 1); end
if isscalar(el_a_deg), el_a_deg = repmat(el_a_deg, n_p, 1); end
if numel(el_n_deg) ~= n_p || numel(el_a_deg) ~= n_p || n_p == 0
    error("mag_pair_ts:sizeMismatch", "el_n_deg and el_a_deg must have the same number of elements (or one be scalar).");
end
opts = fill_opts(opts, n_p);
spec = lib.spec;
pn = struct("tx", antenna_prepare_rev1(ant_n), "rx", antenna_prepare_rev1(ant_n));
pa = struct("tx", antenna_prepare_rev1(ant_a), "rx", antenna_prepare_rev1(ant_a));
[f21, f12] = ts_floor_rev1(pn, pa, spec.el_n, el_n_deg, spec.tilt_a, el_a_deg, ...
    lib.delta_eff_n_deg, lib.delta_eff_a_deg);
thr_min = min([opts.persist_thr_dBi(:); opts.radial_thr_dBi(:)]);
if thr_min <= max([f21; f12])
    error("mag_pair_ts:thresholdBelowFloor", join(["A threshold (%.1f dBi) is at or below the event " ...
        "floor (%.1f dBi); widen the windows (delta_n_deg, delta_a_deg) or raise the threshold."], ""), ...
        thr_min, max([f21; f12]));
end

% ---- chunks of pairs; each chunk gets every trial's states once
starts = 1:opts.batch_size:n_p;
n_ch = numel(starts);
chunk = cell(n_ch, 1);
if opts.parallel_flag == 1
    parfor c = 1:n_ch
        % Broadcast is fine: the angle vectors and opts are small next to lib.
        ch = starts(c):min(starts(c) + opts.batch_size - 1, n_p); %#ok<PFBNS>
        chunk{c} = run_chunk(lib, pn, pa, el_n_deg(ch), el_a_deg(ch), opts); %#ok<PFBNS>
    end
else
    for c = 1:n_ch
        ch = starts(c):min(starts(c) + opts.batch_size - 1, n_p);
        chunk{c} = run_chunk(lib, pn, pa, el_n_deg(ch), el_a_deg(ch), opts);
    end
end

% ---- one result per pair
t_rev = spec.sec_a/spec.rate_a;
t_vol = sum(spec.sec_n./spec.rate_n + spec.turn_n);
R = cell(n_p, 1);
for c = 1:n_ch
    K = chunk{c};
    for q = 1:K.n_pairs
        p = starts(c) + q - 1;
        R{p} = pair_result(K, q, el_n_deg(p), el_a_deg(p), f21(p), f12(p), t_rev, t_vol, opts);
    end
end
R = vertcat(R{:});
end

% =========================================================================
function K = run_chunk(lib, pn, pa, el_n, el_a, opts)
% Per-scan maxima and harm counts of every trial for a chunk of pairs.
spec = lib.spec;
n_pc = numel(el_n);
t_full = (0:lib.dt_s:lib.duration_s)';
n_cut = numel(spec.rate_n);
n_tr = numel(lib.trials);
n_thr = numel(opts.radial_thr_dBi);
n_rad = ceil(360/opts.radial_deg);
K = struct("n_pairs", n_pc, "n_eval", 0, "n_grid", numel(t_full)*n_tr, ...
    "rev", {cell(n_tr, 1)}, "cut", {cell(n_tr, 1)}, "vol", {cell(n_tr, 1)}, "rad", {cell(n_tr, 1)});
for jj = 1:n_tr
    j = lib.trials(jj);
    idx = run_index(double(lib.runs{jj}));
    t = t_full(idx);
    K.n_eval = K.n_eval + numel(t);
    C = scan_clock_rev1(spec, lib.seed, j, lib.duration_s);
    [st_n, st_a, ck] = scan_clock_state_rev1(C, t);
    [~, ~, ends] = scan_clock_state_rev1(C, [0; t_full(end)]);
    % Gains for every pair of the chunk at once: samples by pairs.
    g21 = antenna_gain_rev1(pa.tx, 180, el_a', st_a) + antenna_gain_rev1(pn.rx, 0, el_n', st_n);
    % One pattern each for tx and rx, so g12 = g21.

    % ---- scan keys (shared by the chunk)
    kr = ck.rev - ends.rev(1) + 1;
    n_r = ends.rev(2) - ends.rev(1) + 1;
    key = (ck.vol - ends.vol(1))*n_cut + ck.cut;
    key0 = ends.cut(1);
    key1 = (ends.vol(2) - ends.vol(1))*n_cut + ends.cut(2);
    kc = key - key0 + 1;
    n_c = key1 - key0 + 1;
    keys = (key0:key1)';
    vrel = floor((keys - 1)/n_cut);
    cut_of = keys - vrel*n_cut;
    kv = ck.vol - ends.vol(1) + 1;
    n_v = ends.vol(2) - ends.vol(1) + 1;
    K.rev{jj} = struct("trial", j, "max", col_max(kr, n_r, g21));
    K.cut{jj} = struct("trial", j, "volume", vrel + 1, "cut", cut_of, "el_deg", spec.el_n(cut_of), ...
        "max", col_max(kc, n_c, g21));
    K.vol{jj} = struct("trial", j, "max", col_max(kv, n_v, g21));

    % ---- harm metric: radials of the lowest cut at or above each threshold
    low = abs(spec.el_n(cut_of) - opts.radial_el_deg) < 1e-6;
    in_low = abs(spec.el_n(ck.cut) - opts.radial_el_deg) < 1e-6;
    radial = floor(mod(st_n.mech_az, 360)/opts.radial_deg) + 1;
    counts = zeros(n_c, n_pc, n_thr);
    cell_id = (kc - 1)*n_rad + radial;                      % one id per (cut, radial)
    for i = 1:n_thr
        [r_, p_] = find(in_low & g21 >= opts.radial_thr_dBi(i));
        if ~isempty(r_)
            u = unique([p_, cell_id(r_)], "rows");             % one count per (pair, cut, radial)
            counts(:, :, i) = accumarray([floor((u(:, 2) - 1)/n_rad) + 1, u(:, 1)], 1, [n_c n_pc]);
        end
    end
    K.rad{jj} = struct("trial", j, "low", find(low), "volume", vrel(low) + 1, "cut", cut_of(low), ...
        "n_c", n_c, "counts", counts(low, :, :));
end
end

% -------------------------------------------------------------------------
function idx = run_index(ru)
% Grid indices of the runs [first last] (one row per run), in order.
% Considered: arrayfun with a colon per run (10% of the old run time).
len = ru(:, 2) - ru(:, 1) + 1;
n = sum(len);
step = ones(n, 1);
pos = cumsum([1; len(1:end - 1)]);                % first output position of each run
step(pos) = ru(:, 1) - [0; ru(1:end - 1, 2)];     % jump to the run's first index
idx = cumsum(step);
end

% -------------------------------------------------------------------------
function M = col_max(k, n, G)
% Largest value of each column of G within each scan k (1..n); NaN for a
% scan with no sample. Considered: accumarray per column (one call per pair).
n_pc = size(G, 2);
sub = [repmat(k, n_pc, 1), reshape(repmat(1:n_pc, numel(k), 1), [], 1)];   % repelem of a scalar returns a row
M = accumarray(sub, G(:), [n n_pc], @max, NaN);
end

% =========================================================================
function R = pair_result(K, q, el_n, el_a, f21, f12, t_rev, t_vol, opts)
% Tables and summaries of pair q of chunk K.
% Arrays are concatenated over the trials first and each table is built
% once (building tables per trial was 60% of a batch run).
n_tr = numel(K.rev);
[ar, av, ac, ah] = deal(cell(n_tr, 1));
for jj = 1:n_tr
    ar{jj} = scan_cols(K.rev{jj}.trial, K.rev{jj}.max(:, q));
    av{jj} = scan_cols(K.vol{jj}.trial, K.vol{jj}.max(:, q));
    a = K.cut{jj};
    x = scan_cols(a.trial, a.max(:, q));
    ac{jj} = [x(:, 1), a.volume, a.cut, a.el_deg, x(:, 3:4)];
    h = K.rad{jj};
    cmp = true(h.n_c, 1);
    cmp([1 h.n_c]) = false;
    ah{jj} = [h.trial + zeros(numel(h.low), 1), h.volume, h.cut, cmp(h.low), ...
        reshape(h.counts(:, q, :), numel(h.low), [])];
end
ar = vertcat(ar{:});
av = vertcat(av{:});
ac = vertcat(ac{:});
ah = vertcat(ah{:});
R = struct();
R.el_n_deg = el_n;
R.el_a_deg = el_a;
R.asr_rev = scan_table(ar, "rev");
R.nexrad_volume = scan_table(av, "volume");
R.nexrad_cut = table(ac(:, 1), ac(:, 2), ac(:, 3), ac(:, 4), ac(:, 5), ac(:, 5), logical(ac(:, 6)), ...
    'VariableNames', ["trial", "volume", "cut", "el_deg", "max_g12_dBi", "max_g21_dBi", "complete"]);
R.persistence = persistence(R.asr_rev, opts.persist_thr_dBi, t_rev, "rev");
R.persistence_volume = persistence(R.nexrad_volume, opts.persist_thr_dBi, t_vol, "vol");
R.exceedance = exceedance(R, max(f21, f12), opts.exceed_step_dB);
R.radial_cuts = [table(ah(:, 1), ah(:, 2), ah(:, 3), logical(ah(:, 4)), ...
    'VariableNames', ["trial", "volume", "cut", "complete"]), ...
    array2table(ah(:, 5:end), "VariableNames", "n_radials_ge_" + string(opts.radial_thr_dBi))];
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
R.share = K.n_eval/K.n_grid;
R.notes = [
    "Gains evaluated only in event windows (" + compose("%.1f", 100*R.share) + "% of the fixed grid); maxima above the floor equal fixed-grid maxima"
    "NaN maxima are censored: at or below the floor (" + compose("%.1f", max(f21, f12)) + " dBi)"
    "Harm metric: radials of the " + compose("%g", opts.radial_el_deg) + " deg cut whose MAG reaches each threshold while the NEXRAD beam points there (g21: ATC transmits, NEXRAD receives)"
    "Runs cut by a trial end are counted to the trial end (censored_runs); use trials long enough that the longest run ends inside one"];
if ~opts.keep_scans
    R.asr_rev = R.asr_rev([], :);
    R.nexrad_cut = R.nexrad_cut([], :);
    R.nexrad_volume = R.nexrad_volume([], :);
    R.radial_cuts = R.radial_cuts([], :);
end
end

% -------------------------------------------------------------------------
function X = scan_cols(j, m)
% [trial, scan, max, complete] of one trial; the first and last scans are
% cut by the trial.
n = numel(m);
c = true(n, 1);
c([1 n]) = false;
X = [j + zeros(n, 1), (1:n)', m, c];
end

function T = scan_table(X, name)
T = table(X(:, 1), X(:, 2), X(:, 3), X(:, 3), logical(X(:, 4)), 'VariableNames', ...
    ["trial", name, "max_g12_dBi", "max_g21_dBi", "complete"]);
end

% -------------------------------------------------------------------------
function T = persistence(Rv, thr, t_scan, unit)
% Runs of consecutive complete scans whose maximum reaches each threshold
% (the persistence rule of nexrad_asr11_mag_rev1). Censored (NaN) maxima
% are below every threshold. A run that touches the first or last
% complete scan of a trial is censored by the trial length.
Rc = Rv(Rv.complete, :);
n_thr = numel(thr);
[pct, n_runs, max_run, n_cens] = deal(zeros(n_thr, 1));
hours = height(Rc)*t_scan/3600;
tr = Rc.trial;
m = Rc.max_g21_dBi;
u_tr = unique(tr)';
for i = 1:n_thr
    len = [];
    cens = 0;
    for j = u_tr
        hit = m(tr == j) >= thr(i);
        b = find(diff(double(hit)) ~= 0);
        f = [1; b + 1];
        l = [b; numel(hit)];
        k = hit(f);
        len = [len; l(k) - f(k) + 1]; %#ok<AGROW>
        cens = cens + sum(k & (f == 1 | l == numel(hit)));
    end
    pct(i) = 100*mean(m >= thr(i));
    n_runs(i) = numel(len);
    max_run(i) = max([len; 0]);
    n_cens(i) = cens;
end
T = table(thr(:), pct, n_runs, n_runs/hours, max_run, max_run*t_scan, n_cens, 'VariableNames', ...
    ["threshold_dBi", unit + "_pct", "n_runs", "runs_per_hour", "max_run_" + unit + "s", "max_run_s", "censored_runs"]);
T.Properties.Description = sprintf("%.1f h of complete scans (scan period %.3f s)", hours, t_scan);
end

% -------------------------------------------------------------------------
function T = exceedance(R, fl, step)
% Share of complete revolutions, cuts and volumes whose maximum is at or
% above each gain, from the floor up to the largest maximum.
mr = R.asr_rev.max_g21_dBi(R.asr_rev.complete);
mc = R.nexrad_cut.max_g21_dBi(R.nexrad_cut.complete);
mv = R.nexrad_volume.max_g21_dBi(R.nexrad_volume.complete);
top = max([mr; mc; mv; fl]);
g = (ceil(fl/step)*step:step:ceil(top/step)*step)';
% Count of maxima at or above each gain: histcounts bins are closed on the
% left, so the reverse cumulative sum over [g_i, g_i+1) counts x >= g_i.
% NaN (censored) is not counted, so it is below every gain. Considered:
% mean(x >= g) over a samples-by-gains matrix (slow for long runs).
share = @(x) 100*flipud(cumsum(flipud(histcounts(x, [g; Inf])')))/numel(x);
if isempty(mv)
    T = table(g, NaN(size(g)), NaN(size(g)), NaN(size(g)), 'VariableNames', ["gain_dBi", "rev_pct", "cut_pct", "vol_pct"]);
    return
end
T = table(g, share(mr), share(mc), share(mv), 'VariableNames', ["gain_dBi", "rev_pct", "cut_pct", "vol_pct"]);
end

% -------------------------------------------------------------------------
function opts = fill_opts(opts, n_p)
d = struct("persist_thr_dBi", [40 45 50 55 60 65 70], "radial_el_deg", 0.5, "radial_deg", 1, ...
    "radial_thr_dBi", [40 50 60], "exceed_step_dB", 0.1, "keep_scans", n_p == 1, ...
    "batch_size", 64, "parallel_flag", 0);
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
if ~(isscalar(opts.batch_size) && opts.batch_size >= 1 && opts.batch_size == fix(opts.batch_size))
    error("mag_pair_ts:badBatch", "opts.batch_size must be a positive integer.");
end
end
