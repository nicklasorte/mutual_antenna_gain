%REPORT_NUMBERS_REV12 Tables and figure for methodology report rev12.
%   Answers the independent recheck of revision 11
%   (docs/rev11_independent_review_rev1.md, findings V1 to V10). Reruns the
%   20-seed study of revisions 8 to 11 (same seeds; periodic runs now break
%   ties as documented, finding V7) and adds:
%     - the lock behind each longest periodic run: the drift over a
%       multiple of its period up to 16 volumes (lock_drift_multiple_rev1,
%       finding V2);
%     - the high-volume rate inside each trial's longest periodic run (V3);
%     - the chance comparison with 20,000 null trials, each seed's pooled
%       share of high volumes, and 95% Clopper-Pearson intervals (V4, V5);
%     - the chance of the random-walk periodic maxima (V6).
%   Reads ts_all_modes_rev8.csv for volume times.
%
%     seeds_runs_rev12.csv          per seed: runs, episodes (3, 5, 8 dB),
%                                   periodic runs, with censoring
%     lock_trials_rev12.csv         per fixed-offset trial, as in revision
%                                   11, plus the lock behind the periodic
%                                   run (drift, k), the run in hours, and its
%                                   high volumes per hour
%     lock_k16_rev12.csv            periodic runs of 4, 6.5, and 8 hours or
%                                   more: their lock drift, and how many
%                                   trials have a lock that tight
%     lock_fit_rev12.csv            share of rate draws whose longest
%                                   periodic run reaches 1, 2, 4, or 8 hours,
%                                   with 95% intervals, against chance
%     lock_walk_null_rev12.csv      chance of the random-walk periodic maxima
%     fig_lock_k16_rev12.png        longest periodic run against its
%                                   lock drift
%
%   Set rev12_overrides (struct: out_dir, n_seeds, drift_trials,
%   drift_hours, n_null, parallel_flag, reuse_seeds) for a small run;
%   reuse_seeds = true reads seeds_runs_rev12.csv and lock_trials_rev12.csv.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also report_numbers_rev11, periodic_runs_rev1, lock_drift_multiple_rev1.

if ~exist("rev12_overrides", "var"), rev12_overrides = struct(); end
here = fileparts(which("report_numbers_rev12"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "in_dir", fullfile(here, "docs", "figures"), ...
    "n_seeds", 20, "drift_trials", 10, "drift_hours", 24, "parallel_flag", 1, "reuse_seeds", false, "n_null", 20000);
f = fieldnames(rev12_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev12_overrides.(f{i});
end
out_dir = ov.out_dir;
if ~isfolder(out_dir), mkdir(out_dir); end
tiny = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1);
if ov.parallel_flag == 1 && ~ov.reuse_seeds && isempty(gcp("nocreate"))
    parpool("Processes");
end
use_par = ~isempty(gcp("nocreate"));
m_cdf = "VCP 212 SAILSx3 AVSET 6.4";
mc = vcp_mode_cases_rev1();
mc = mc(2:end);                                     % c{1} is the working-paper baseline, not a WSR-88D mode
labs = cellfun(@(x) string(x.label), mc);
T_all = readtable(fullfile(ov.in_dir, "ts_all_modes_rev8.csv"), "TextType", "string");
modesC = unique([T_all.mode(T_all.el_deg == 0.5 & T_all.vol_pct_60 == max(T_all.vol_pct_60(T_all.el_deg == 0.5))); ...
    m_cdf], "stable");
models = ["per_trial", "random_walk", "random_walk"];
steps = [NaN 0.01 0.05];
[MM, DD, SS] = ndgrid(1:numel(modesC), 1:numel(models), 1:ov.n_seeds);
nC = numel(MM);
rowsC = cell(nC, 1);
rowsL = cell(nC, 1);
dt_ = ov.drift_trials;
dh = ov.drift_hours;
if ov.reuse_seeds
    rowsC = {readtable(fullfile(out_dir, "seeds_runs_rev12.csv"), "TextType", "string")};
    rowsL = {readtable(fullfile(out_dir, "lock_trials_rev12.csv"), "TextType", "string")};
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
writetable(T_seed, fullfile(out_dir, "seeds_runs_rev12.csv"));
T_lock = vertcat(rowsL{:});
writetable(T_lock, fullfile(out_dir, "lock_trials_rev12.csv"));

% Per-trial additions: the lock behind each trial's longest periodic run,
% the smallest drift from whole ATC revolutions over a multiple of the
% run's period up to 16 volumes (V2), and the high-volume rate inside the
% run (V3).
[dk, kk] = lock_drift_multiple_rev1(T_lock.slip_deg_per_vol, T_lock.per_k, 16);
T_lock.lock16_drift_deg = dk;
T_lock.lock16_k = kk;
t_v_row = zeros(height(T_lock), 1);
for m = modesC'
    t_v_row(T_lock.mode == m) = T_all.t_volume_s(find(T_all.mode == m, 1));
end
T_lock.per_hours = T_lock.per_vols.*t_v_row/3600;
T_lock.per_high_per_h = T_lock.per_hits./max(T_lock.per_vols, 1)*3600./t_v_row;   % high volumes per hour inside the run
writetable(T_lock, fullfile(out_dir, "lock_trials_rev12.csv"));

% Long periodic runs and their lock (V2): for each duration, the largest
% lock drift among the trials that reach it, how many trials have a lock
% that tight, and how many of those reach it.
K2 = {};
for m = modesC'
    for a = [0 0.5]
        for th = [60 70]
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            for H = [4 6.5 8]
                lg = L.per_hours >= H;
                if ~any(lg)
                    continue
                end
                dmax = max(L.lock16_drift_deg(lg));
                tight = L.lock16_drift_deg <= dmax;
                K2(end + 1, :) = {m, a, th, H, nnz(lg), dmax, median(L.lock16_drift_deg(lg)), nnz(tight), ...
                    100*nnz(lg)/nnz(tight), height(L), nnz(lg & L.m4_deg < 8), ...
                    median(L.per_high_per_h(lg))}; %#ok<SAGROW>
            end
        end
    end
end
if isempty(K2), K2 = cell(0, 12); end
writetable(cell2table(K2, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_long", ...
    "max_lock_drift_long_deg", "median_lock_drift_long_deg", "n_trials_that_tight", "pct_of_those_long", ...
    "n_trials", "n_long_m4_below_8", "median_high_per_h_long"]), fullfile(out_dir, "lock_k16_rev12.csv"));

% Chance: independent volumes, each null trial drawing the pooled share of
% high volumes of one seed in turn (the seeds' shares, not each trial's;
% recheck V4), same number of complete volumes, same statistics. Each null
% trial has its own RandStream substream, so serial and parallel agree.
n_null = ov.n_null;
FL = {};
WL = {};
cell_id = 0;
for m = modesC'
    t_v = T_all.t_volume_s(find(T_all.mode == m, 1));
    n_v = floor(3600*dh/t_v) - 1;                   % complete volumes per trial
    for a = [0 0.5]
        for th = [60 70]
            cell_id = cell_id + 1;
            sel = T_seed.mode == m & T_seed.drift_model == "per_trial" & T_seed.el_deg == a & T_seed.threshold_dBi == th;
            pv = T_seed.vol_pct(sel)/100;
            pmax = zeros(n_null, 1);
            base = 1e6*cell_id;
            if use_par
                parfor i = 1:n_null
                    pmax(i) = null_trial(pv, i, n_v, base);
                end
            else
                for i = 1:n_null
                    pmax(i) = null_trial(pv, i, n_v, base);
                end
            end
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            n_m = height(L);
            for H = [1 2 4 8]
                N = ceil(H*3600/t_v);
                xm = nnz(L.per_vols >= N);
                xn = nnz(pmax >= N);
                [pm, cm] = binofit(xm, n_m);
                [pn, cn] = binofit(xn, n_null);
                FL(end + 1, :) = {m, a, th, H, N, 100*mean(L.longest_run_vols >= N), 100*mean(L.ep5_vols >= N), ...
                    100*pm, 100*cm(1), 100*cm(2), nnz(L.per_vols >= N & L.per_censored == 1), n_m, ...
                    100*pn, 100*cn(1), 100*cn(2), n_null, pm/pn, cm(1)/cn(2), cm(2)/cn(1)}; %#ok<SAGROW>
            end
            % Random-walk maxima against chance (V6): the chance that at
            % least one of the 200 trials of a walk model reaches the
            % largest periodic run it gave.
            for st = [0.01 0.05]
                w = T_seed.mode == m & T_seed.drift_model == "random_walk" & T_seed.step_deg_s_per_rev == st & ...
                    T_seed.el_deg == a & T_seed.threshold_dBi == th;
                obs = max(T_seed.max_periodic_vols(w));
                p1 = mean(pmax >= obs);
                n_tr = 10*nnz(w);
                WL(end + 1, :) = {m, a, th, st, obs, obs*t_v/3600, p1, 1 - (1 - p1)^n_tr, n_tr}; %#ok<SAGROW>
            end
        end
    end
end
writetable(cell2table(FL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_vols", ...
    "trials_run_ge_pct", "trials_episode5_ge_pct", "trials_periodic_ge_pct", "periodic_lo_pct", ...
    "periodic_hi_pct", "n_censored_periodic", "n_trials", "null_periodic_ge_pct", "null_lo_pct", ...
    "null_hi_pct", "n_null_trials", "ratio", "ratio_lo", "ratio_hi"]), fullfile(out_dir, "lock_fit_rev12.csv"));
writetable(cell2table(WL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "step_deg_s_per_rev", ...
    "max_periodic_vols", "max_periodic_hours", "null_trial_ge_pct", "chance_any_of_trials", "n_trials"]), ...
    fullfile(out_dir, "lock_walk_null_rev12.csv"));

report_figures_rev12(out_dir);

% -------------------------------------------------------------------------
function p = null_trial(pv, i, n_v, base)
% Longest periodic run (k = 1 to 4) of one null trial of n_v independent
% volumes, high with the share of seed 1 + mod(i - 1, numel(pv)).
rs = RandStream("mrg32k3a", "Seed", 20261112);
rs.Substream = base + i;
g = 50 + 11*(rand(rs, n_v, 1) < pv(1 + mod(i - 1, numel(pv))));   % 61 dBi when high, 50 dBi otherwise
V = table(ones(n_v, 1), true(n_v, 1), g, 'VariableNames', ["trial" "complete" "max_g21_dBi"]);
pr = periodic_runs_rev1(V, 60, 4);
p = pr(1);
end

% =========================================================================
function o = merge(a, b)
o = a;
f = fieldnames(b);
for i = 1:numel(f)
    o.(f{i}) = b.(f{i});
end
end

function [T, L] = seed_runs(mcase, tiny, model, step, seed, n_tr, hours)
% Runs, episodes at three margins, and periodic runs of one seed of n_tr
% trials of `hours` each; per-trial values for the fixed-offset model.
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
[slip, rate] = deal(NaN(n_tr, 1));
if model == "per_trial"
    for j = 1:n_tr
        C = scan_clock_rev1(spec, seed, lib.trials(j), dur);
        rate(j) = C.r;
        slip(j) = mod(C.r*t_vol + 180, 360) - 180;  % ATC azimuth shift per WSR-88D volume, deg
    end
end
[m4, kb] = lock_drift_m4_rev1(slip);
rows = {};
rowsL = {};
for a = 1:2
    Vt = R(a).nexrad_volume;
    for th = [60 70]
        [r5, t5] = run_stats_rev1(Vt, th, 5);
        [r3, t3] = run_stats_rev1(Vt, th, 3);
        [r8, t8] = run_stats_rev1(Vt, th, 8);
        [pr, tp] = periodic_runs_rev1(Vt, th, 4);
        V = Vt(Vt.complete, :);
        rows(end + 1, :) = {lab, model, step, seed, R(a).el_n_deg, th, 100*mean(V.max_g21_dBi >= th), ...
            r5(1), r5(2), r5(4), r5(7), r5(9), r3(10), r5(10), r5(11), r8(10), pr(1), pr(2), pr(4), pr(5)}; %#ok<AGROW>
        if model == "per_trial"
            for j = 1:n_tr
                rowsL(end + 1, :) = {lab, seed, lib.trials(j), R(a).el_n_deg, th, rate(j), slip(j), m4(j), kb(j), ...
                    t5.longest(j), t5.censored(j), t3.ep_span(j), t5.ep_span(j), t5.ep_censored(j), t8.ep_span(j), ...
                    tp.dur(j), tp.k(j), tp.hits(j), tp.censored(j)}; %#ok<AGROW>
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "seed", "el_deg", ...
    "threshold_dBi", "vol_pct", "median_run_vols", "p95_run_vols", "max_run_vols", "max_uncensored_run_vols", ...
    "max_censored_run_vols", "max_episode3_vols", "max_episode5_vols", "max_uncensored_episode5_vols", ...
    "max_episode8_vols", "max_periodic_vols", "max_periodic_k", "max_uncensored_periodic_vols", ...
    "max_censored_periodic_vols"]);
vn = ["mode", "seed", "trial", "el_deg", "threshold_dBi", "rate_deg_s", "slip_deg_per_vol", "m4_deg", "m4_k", ...
    "longest_run_vols", "longest_censored", "ep3_vols", "ep5_vols", "ep5_censored", "ep8_vols", "per_vols", ...
    "per_k", "per_hits", "per_censored"];
L = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
if ~isempty(rowsL)
    L = cell2table(rowsL, 'VariableNames', vn);
end
end
