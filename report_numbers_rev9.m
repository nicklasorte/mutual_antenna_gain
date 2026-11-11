%REPORT_NUMBERS_REV9 Tables and figures for methodology report rev9.
%   Answers the independent recheck of revision 8
%   (docs/rev8_independent_review_rev1.md, findings S1 to S11) and writes
%   every evidence file that revision 9 adds to docs/figures. Runs under
%   the revision 6 defaults (ROC specification, 216 modes, ASR-11, exact
%   CDF) and reads ts_all_modes_rev8.csv, written by report_numbers_rev8.
%
%     per_hour_modes_rev9.csv     for each of the 216 modes and both paths:
%                                 high-coupling volumes, exposed 0.5 deg
%                                 cuts, and ATC revolutions per hour at 60
%                                 dBi beside the shares, with ranks (S1, S4)
%     per_hour_envelope_rev9.csv  largest value of each per-hour and share
%                                 statistic, the mode setting it, and the
%                                 rank and ratio of the CDF-worst mode
%     run_modes_rev9.csv          spread over the 216 modes of the run
%                                 statistics in ts_all_modes_rev8.csv (S3)
%     seeds_runs_rev9.csv         20 seeds x ten 24 h trials, three drift
%                                 models, two modes, with runs cut by a
%                                 trial end flagged as censored (S2)
%     seeds_summary_rev9.csv      spread over seeds, longest uncensored and
%                                 censored runs
%     lock_trials_rev9.csv        fixed-offset trials: the ATC rate drawn,
%                                 the ATC phase slip per WSR-88D volume, and
%                                 the longest run in the trial (S2)
%     lock_fit_rev9.csv           share of fixed-offset 24 h trials whose
%                                 longest run reaches N volumes, and how
%                                 many of those runs the trial end cut
%     timing_modes_rev9.csv       time series cost at six modes and the
%                                 event library cost (S10)
%     timing_rev9_machine.txt     machine and settings of the timings
%     fig13_runs_vs_drift_rev9.png, fig_lock_vs_slip_rev9.png
%
%   The seed runs use parfor when a pool is open (or
%   rev9_overrides.parallel_flag = 1 opens one). Set rev9_overrides
%   (struct: out_dir, n_seeds, drift_trials, drift_hours, n_modes_time,
%   parallel_flag) before running for a small run. reuse_seeds = true
%   and reuse_timing = true read the seed and timing tables already in
%   out_dir instead of rerunning them.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also report_numbers_rev8, report_figures_rev9, scan_clock_rev1.

if ~exist("rev9_overrides", "var"), rev9_overrides = struct(); end
here = fileparts(which("report_numbers_rev9"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "in_dir", fullfile(here, "docs", "figures"), ...
    "n_seeds", 20, "drift_trials", 10, "drift_hours", 24, "n_modes_time", 6, "time_hours", 120, ...
    "parallel_flag", 1, "reuse_seeds", false, "reuse_timing", false);
f = fieldnames(rev9_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev9_overrides.(f{i});
end
out_dir = ov.out_dir;
if ~isfolder(out_dir), mkdir(out_dir); end
tiny = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1);
if ov.parallel_flag == 1 && ~(ov.reuse_seeds && ov.reuse_timing) && isempty(gcp("nocreate"))
    parpool("Processes");
end
use_par = ~isempty(gcp("nocreate"));
m_cdf = "VCP 212 SAILSx3 AVSET 6.4";
mc = vcp_mode_cases_rev1();
mc = mc(2:end);                                     % c{1} is the working-paper baseline, not a WSR-88D mode
labs = cellfun(@(x) string(x.label), mc);

% ---- A. Per-hour statistics for every mode (S1, S4)
T_all = readtable(fullfile(ov.in_dir, "ts_all_modes_rev8.csv"), "TextType", "string");
s0 = nexrad_asr11_mag_rev1(tiny);
t_rev_atc = 360/(6*s0.sys_asr11.scan.rate_rpm);    % ATC revolution, s
[scans, names] = nexrad_mode_scans_rev1();
n05 = zeros(height(T_all), 1);
for i = 1:height(T_all)
    s = scans{names == T_all.mode(i)};
    n05(i) = sum(s.n_rev(abs(s.tilt_deg - 0.5) < 1e-6));
end
vph = 3600./T_all.t_volume_s;                       % volumes per hour
P = table(T_all.mode, T_all.el_deg, T_all.t_volume_s, n05, ...
    T_all.vol_pct_60/100.*vph, T_all.cuts_hit_60_pct/100.*n05.*vph, T_all.rev_pct_60/100*3600/t_rev_atc, ...
    T_all.vol_pct_60, T_all.cuts_hit_60_pct, T_all.rev_pct_60, ...
    'VariableNames', ["mode", "el_deg", "t_volume_s", "cuts05_per_volume", "vols_per_h_60", ...
    "cuts_per_h_60", "revs_per_h_60", "vol_pct_60", "cuts_hit_60_pct", "rev_pct_60"]);
mets = ["vols_per_h_60", "cuts_per_h_60", "revs_per_h_60", "vol_pct_60", "cuts_hit_60_pct", "rev_pct_60"];
for q = mets
    P.("rank_" + q) = zeros(height(P), 1);
    for a = [0 0.5]
        k = P.el_deg == a;
        [~, ~, r] = unique(-P.(q)(k));              % dense rank, 1 = largest
        P.("rank_" + q)(k) = r;
    end
end
writetable(P, fullfile(out_dir, "per_hour_modes_rev9.csv"));
EV = {};
for a = [0 0.5]
    Pa = P(P.el_deg == a, :);
    for q = mets
        [mx, im] = max(Pa.(q));
        c = Pa(Pa.mode == m_cdf, :);
        EV(end + 1, :) = {a, q, mx, Pa.mode(im), c.(q), c.("rank_" + q), c.(q)/mx}; %#ok<SAGROW>
    end
end
writetable(cell2table(EV, 'VariableNames', ["el_deg", "metric", "max_over_modes", "worst_mode", ...
    "at_cdf_worst_mode", "rank_of_cdf_worst", "ratio_to_max"]), fullfile(out_dir, "per_hour_envelope_rev9.csv"));

% ---- B. Run statistics over the 216 modes (S3), one-hour trials
RM = {};
for a = [0 0.5]
    Ta = T_all(T_all.el_deg == a, :);
    for q = ["p95_run_60_vols", "share_runs_ge2_60_pct", "max_run_60_vols"]
        v = Ta.(q);
        RM(end + 1, :) = {a, q, min(v), median(v), max(v), sum(v == min(v)), numel(v)}; %#ok<SAGROW>
    end
end
writetable(cell2table(RM, 'VariableNames', ["el_deg", "statistic", "min_over_modes", "median_over_modes", ...
    "max_over_modes", "n_modes_at_min", "n_modes"]), fullfile(out_dir, "run_modes_rev9.csv"));

% ---- C. Seeds with censored runs flagged (S2)
modesC = [T_all.mode(T_all.el_deg == 0.5 & T_all.vol_pct_60 == max(T_all.vol_pct_60(T_all.el_deg == 0.5))); m_cdf];
modesC = unique(modesC, "stable");
models = ["per_trial", "random_walk", "random_walk"];
steps = [NaN 0.01 0.05];
[MM, DD, SS] = ndgrid(1:numel(modesC), 1:numel(models), 1:ov.n_seeds);
nC = numel(MM);
rowsC = cell(nC, 1);
rowsL = cell(nC, 1);
dt_ = ov.drift_trials;
dh = ov.drift_hours;
if ov.reuse_seeds
    rowsC = {readtable(fullfile(out_dir, "seeds_runs_rev9.csv"), "TextType", "string")};
    rowsL = {readtable(fullfile(out_dir, "lock_trials_rev9.csv"), "TextType", "string")};
elseif use_par
    parfor i = 1:nC
        [rowsC{i}, rowsL{i}] = seed_runs(mc{labs == modesC(MM(i))}, tiny, models(DD(i)), steps(DD(i)), ...
            20261200 + SS(i), dt_, dh); %#ok<PFBNS>
    end
else
    for i = 1:nC
        [rowsC{i}, rowsL{i}] = seed_runs(mc{labs == modesC(MM(i))}, tiny, models(DD(i)), steps(DD(i)), ...
            20261200 + SS(i), dt_, dh);
    end
end
T_seed = vertcat(rowsC{:});
writetable(T_seed, fullfile(out_dir, "seeds_runs_rev9.csv"));
T_lock = vertcat(rowsL{:});
writetable(T_lock, fullfile(out_dir, "lock_trials_rev9.csv"));
step_key = T_seed.step_deg_s_per_rev;
step_key(isnan(step_key)) = 0;                      % findgroups drops NaN keys (the fixed-offset model)
G = findgroups(T_seed.mode, T_seed.drift_model, step_key, T_seed.el_deg, T_seed.threshold_dBi);
S = splitapply(@(md, dm, st, el, th, vp, med, p95, mxr, mxu, nce, mxc, s2) {md(1), dm(1), st(1), el(1), th(1), ...
    numel(vp), min(vp), max(vp), prctile(med, 50), min(p95), max(p95), max(mxr), max(mxu), sum(nce), max(mxc), ...
    min(s2), max(s2)}, T_seed.mode, T_seed.drift_model, T_seed.step_deg_s_per_rev, T_seed.el_deg, ...
    T_seed.threshold_dBi, T_seed.vol_pct, T_seed.median_run_vols, T_seed.p95_run_vols, T_seed.max_run_vols, ...
    T_seed.max_uncensored_run_vols, T_seed.n_censored_runs, T_seed.max_censored_run_vols, ...
    T_seed.share_runs_ge2_pct, G);
T_ss = cell2table(S, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "el_deg", "threshold_dBi", ...
    "n_seeds", "vol_pct_min", "vol_pct_max", "median_run_median", "p95_run_min", "p95_run_max", ...
    "max_run", "max_uncensored_run", "n_censored_runs", "max_censored_run", "share_ge2_min", "share_ge2_max"]);
writetable(T_ss, fullfile(out_dir, "seeds_summary_rev9.csv"));

% Longest run per fixed-offset trial (S2): the share of 24 h trials (one
% rate draw each) whose longest run reaches N volumes, how many of those
% the trial end cut, and the smallest |phase slip| among them.
FL = {};
for m = modesC'
    t_v = T_all.t_volume_s(find(T_all.mode == m, 1));
    for a = [0 0.5]
        for th = [60 70]
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            for N = [6 10 20 40]
                k = L.longest_run_vols >= N;
                FL(end + 1, :) = {m, a, th, N, N*t_v/3600, 100*mean(k), nnz(k), nnz(k & L.longest_censored), ...
                    min([Inf; abs(L.slip_deg_per_vol(k))]), height(L)}; %#ok<SAGROW>
            end
        end
    end
end
writetable(cell2table(FL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "n_vols", "hours", ...
    "trials_with_run_ge_N_pct", "n_trials_ge_N", "n_censored_ge_N", "min_abs_slip_deg", "n_trials"]), ...
    fullfile(out_dir, "lock_fit_rev9.csv"));

% ---- D. Time-series cost at several modes (S10)
tm = ["VCP 212 SAILSx3 AVSET 6.4", "VCP 34 SAILSx1", "VCP 35 SAILSx1", "VCP 32", "VCP 31", "VCP 12"];
tm = tm(ismember(tm, labs));
tm = tm(1:min(numel(tm), ov.n_modes_time));
rst = RandStream("mrg32k3a", "Seed", 3);
eb = 2.5*rand(rst, 32, 1) - 0.5;
TT = {};
if ov.reuse_timing
    tm = strings(1, 0);                              % keep the timing tables already written
end
for m = tm
    w = nexrad_asr11_mag_rev1(merge(tiny, rmfield(mc{labs == m}, "label")));
    lopt = struct("seed", 20261105, "n_trials", ov.time_hours, "duration_s", 3600);
    t_lib = timeit(@() scan_event_library_rev1(w.scan_max.clock_spec, lopt));
    lib = scan_event_library_rev1(w.scan_max.clock_spec, lopt);
    t_one = timeit(@() mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0));
    t_b = timeit(@() mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, eb, eb))/32;
    TT(end + 1, :) = {m, w.scan_max.t_volume_s, t_lib, t_one, t_b}; %#ok<SAGROW>
end
if ~ov.reuse_timing
    writetable(cell2table(TT, 'VariableNames', ["mode", "t_volume_s", "library_s", "one_pair_s", ...
        "per_pair_batch32_s"]), fullfile(out_dir, "timing_modes_rev9.csv"));
    fid = fopen(fullfile(out_dir, "timing_rev9_machine.txt"), "w");
    fprintf(fid, "%s; %s; computational threads %d; %s; parallel pool for seeds: %d; %d h per timed time series\n", ...
        computer, version, maxNumCompThreads, string(datetime("now")), use_par, ov.time_hours);
    fclose(fid);
end

% ---- E. Figure 13 and the lock figure, drawn from the tables just written
report_figures_rev9(out_dir);

% =========================================================================
function o = merge(a, b)
o = a;
f = fieldnames(b);
for i = 1:numel(f)
    o.(f{i}) = b.(f{i});
end
end

function [T, L] = seed_runs(mcase, tiny, model, step, seed, n_tr, hours)
% Run statistics of one seed of n_tr trials of `hours` each, and, for the
% fixed-offset model, the rate, phase slip and longest run of each trial.
lab = string(mcase.label);
w = nexrad_asr11_mag_rev1(merge(tiny, rmfield(mcase, "label")));
spec = w.scan_max.clock_spec;
spec.jitter_mode = model;
if model == "random_walk"
    spec.drift_step_a = step;
end
dur = 3600*hours;
lib = scan_event_library_rev1(spec, struct("seed", seed, "n_trials", n_tr, "duration_s", dur));
R = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, [0; 0.5], [0; 0.5], ...
    struct("persist_thr_dBi", [60 70], "radial_thr_dBi", 60, "keep_scans", true));
t_vol = w.scan_max.t_volume_s;
slip = NaN(n_tr, 1);
rate = NaN(n_tr, 1);
if model == "per_trial"
    for j = 1:n_tr
        C = scan_clock_rev1(spec, seed, lib.trials(j), dur);
        rate(j) = C.r;
        slip(j) = mod(C.r*t_vol + 180, 360) - 180;  % ATC azimuth shift per WSR-88D volume, deg
    end
end
rows = {};
rowsL = {};
for a = 1:2
    for th = [60 70]
        [rs, tr] = run_stats(R(a).nexrad_volume, th);
        V = R(a).nexrad_volume(R(a).nexrad_volume.complete, :);
        rows(end + 1, :) = {lab, model, step, seed, R(a).el_n_deg, th, 100*mean(V.max_g21_dBi >= th), ...
            rs(1), rs(2), rs(3), rs(4), rs(5), rs(6), rs(7), rs(8), rs(9)}; %#ok<AGROW>
        if model == "per_trial"
            for j = 1:n_tr
                rowsL(end + 1, :) = {lab, seed, lib.trials(j), R(a).el_n_deg, th, rate(j), slip(j), ...
                    tr.longest(j), tr.censored(j)}; %#ok<AGROW>
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "seed", "el_deg", ...
    "threshold_dBi", "vol_pct", "median_run_vols", "p95_run_vols", "n_runs", "max_run_vols", ...
    "share_runs_ge2_pct", "share_runs_ge3_pct", "max_uncensored_run_vols", "n_censored_runs", ...
    "max_censored_run_vols"]);
L = cell2table(cell(0, 9), 'VariableNames', ["mode", "seed", "trial", "el_deg", "threshold_dBi", ...
    "rate_deg_s", "slip_deg_per_vol", "longest_run_vols", "longest_censored"]);
if ~isempty(rowsL)
    L = cell2table(rowsL, 'VariableNames', L.Properties.VariableNames);
end
end

function [s, tr] = run_stats(Vt, thr)
% Runs of consecutive complete volumes reaching thr, per trial. A run that
% touches the first or last complete volume of its trial is censored: the
% trial cut it, so its length is a lower bound.
%   s = [median run, 95th-percentile run, n runs, longest run, share of
%        high volumes in runs of 2 or more (%), in runs of 3 or more (%),
%        longest uncensored run, n censored runs, longest censored run]
%   tr.longest, tr.censored: longest run of each trial and whether it is
%        censored
V = Vt(Vt.complete, :);
tu = unique(V.trial)';
len = [];
cen = false(0, 1);
tr = struct("longest", zeros(numel(tu), 1), "censored", false(numel(tu), 1));
for jj = 1:numel(tu)
    hit = V.max_g21_dBi(V.trial == tu(jj)) >= thr;
    d = diff([0; hit(:); 0]);
    b = find(d == 1);
    e = find(d == -1) - 1;
    l = e - b + 1;
    c = b == 1 | e == numel(hit);
    len = [len; l]; %#ok<AGROW>
    cen = [cen; c]; %#ok<AGROW>
    if ~isempty(l)
        [tr.longest(jj), im] = max(l);
        tr.censored(jj) = c(im);
    end
end
if isempty(len)
    s = [NaN NaN 0 0 0 0 0 0 0];
    return
end
mu = max([0; len(~cen)]);
mcn = max([0; len(cen)]);
s = [median(len), prctile(len, 95), numel(len), max(len), 100*sum(len(len >= 2))/sum(len), ...
    100*sum(len(len >= 3))/sum(len), mu, nnz(cen), mcn];
end
