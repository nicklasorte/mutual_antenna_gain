%REPORT_NUMBERS_REV11 Tables and figure for methodology report rev11.
%   Answers the independent recheck of revision 10
%   (docs/rev10_independent_review_rev1.md, findings U1 to U11). Reruns the
%   20-seed study of revisions 8 to 10 (same seeds) and adds periodic runs
%   (periodic_runs_rev1: high volumes recurring every k = 1 to 4 volumes,
%   finding U1), lock episodes at margins of 3, 5, and 8 dB (U2), and the
%   longest run, episode, and periodic run from separate trials (U3).
%   Reads ts_all_modes_rev8.csv for volume times.
%
%     seeds_runs_rev11.csv       per seed: runs, episodes (3, 5, 8 dB),
%                                periodic runs, with censoring
%     lock_trials_rev11.csv      per fixed-offset trial: rate, slip, m4 and
%                                its k, longest run, episodes at each
%                                margin, longest periodic run and its k
%     lock_fit_rev11.csv         share of rate draws whose longest run,
%                                episode (5 dB), or periodic run reaches 1,
%                                2, 4, or 8 hours
%     lock_top_trials_rev11.csv  the five trials with the longest periodic
%                                run, episode, and run for each mode, path,
%                                and level (one row per trial)
%     lock_m4_rev11.csv          share of trials with a periodic run of 14
%                                volumes or more by m4 bin, and the largest
%                                m4 among them (U4)
%     lock_null_rev11.csv        the same shares for independent volumes
%                                (each volume high with the trial's own
%                                probability, order shuffled away): what
%                                chance alone gives, for comparison
%     fig_lock_periodic_rev11.png  longest periodic run against m4
%
%   Set rev11_overrides (struct: out_dir, n_seeds, drift_trials,
%   drift_hours, parallel_flag, reuse_seeds) for a small run;
%   reuse_seeds = true reads seeds_runs_rev11.csv and lock_trials_rev11.csv
%   instead of rerunning the seeds.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also report_numbers_rev10, periodic_runs_rev1, run_stats_rev1,
%   lock_drift_m4_rev1.

if ~exist("rev11_overrides", "var"), rev11_overrides = struct(); end
here = fileparts(which("report_numbers_rev11"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "in_dir", fullfile(here, "docs", "figures"), ...
    "n_seeds", 20, "drift_trials", 10, "drift_hours", 24, "parallel_flag", 1, "reuse_seeds", false);
f = fieldnames(rev11_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev11_overrides.(f{i});
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
    rowsC = {readtable(fullfile(out_dir, "seeds_runs_rev11.csv"), "TextType", "string")};
    rowsL = {readtable(fullfile(out_dir, "lock_trials_rev11.csv"), "TextType", "string")};
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
writetable(T_seed, fullfile(out_dir, "seeds_runs_rev11.csv"));
T_lock = vertcat(rowsL{:});
writetable(T_lock, fullfile(out_dir, "lock_trials_rev11.csv"));

% Share of rate draws (fixed-offset trials) reaching 1, 2, or 4 hours.
FL = {};
TP = {};
MB = {};
for m = modesC'
    t_v = T_all.t_volume_s(find(T_all.mode == m, 1));
    for a = [0 0.5]
        for th = [60 70]
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            for H = [1 2 4 8]
                N = ceil(H*3600/t_v);
                FL(end + 1, :) = {m, a, th, H, N, 100*mean(L.longest_run_vols >= N), ...
                    100*mean(L.ep5_vols >= N), 100*mean(L.per_vols >= N), ...
                    nnz(L.per_vols >= N & L.per_censored == 1), height(L)}; %#ok<SAGROW>
            end
            % Five longest trials by each statistic (one row per trial).
            for q = ["per_vols", "ep5_vols", "longest_run_vols"]
                [~, o] = sort(L.(q), "descend");
                for r = 1:min(5, height(L))
                    x = L(o(r), :);
                    TP(end + 1, :) = {m, a, th, q, r, x.seed, x.trial, x.m4_deg, x.m4_k, x.(q), x.(q)*t_v/3600, ...
                        x.per_k, x.longest_censored, x.ep5_censored, x.per_censored}; %#ok<SAGROW>
                end
            end
            % Periodic runs of 14 volumes or more by m4 bin (U4).
            edges = [0 2 4 8 12 16 45 180];
            for b = 1:numel(edges) - 1
                k = L.m4_deg >= edges(b) & L.m4_deg < edges(b + 1);
                MB(end + 1, :) = {m, a, th, edges(b), edges(b + 1), nnz(k), nnz(k & L.per_vols >= 14), ...
                    nnz(k & L.longest_run_vols >= 14), max([0; L.per_vols(k)])}; %#ok<SAGROW>
            end
        end
    end
end
writetable(cell2table(FL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_vols", ...
    "trials_run_ge_pct", "trials_episode5_ge_pct", "trials_periodic_ge_pct", "n_censored_periodic", ...
    "n_trials"]), fullfile(out_dir, "lock_fit_rev11.csv"));
writetable(cell2table(TP, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "statistic", "rank", "seed", ...
    "trial", "m4_deg", "m4_k", "vols", "hours", "per_k", "run_censored", "ep5_censored", "per_censored"]), ...
    fullfile(out_dir, "lock_top_trials_rev11.csv"));
writetable(cell2table(MB, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "m4_lo_deg", "m4_hi_deg", ...
    "n_trials", "n_periodic_ge14", "n_run_ge14", "max_periodic_vols"]), fullfile(out_dir, "lock_m4_rev11.csv"));

% Chance alone: independent volumes with each trial's share of high volumes
% (Bernoulli), same trial length, same statistics. 200 null trials per
% fixed-offset trial set; RandStream substreams keep it reproducible.
rsn = RandStream("mrg32k3a", "Seed", 20261111);
NL = {};
for m = modesC'
    t_v = T_all.t_volume_s(find(T_all.mode == m, 1));
    n_v = floor(3600*dh/t_v) - 1;                   % complete volumes per trial
    for a = [0 0.5]
        for th = [60 70]
            sel = T_seed.mode == m & T_seed.drift_model == "per_trial" & T_seed.el_deg == a & T_seed.threshold_dBi == th;
            pv = T_seed.vol_pct(sel)/100;           % share of high volumes per seed
            n_null = 200;
            [rmax, pmax] = deal(zeros(n_null, 1));
            for i = 1:n_null
                p_i = pv(1 + mod(i - 1, numel(pv)));
                g = 50 + 11*(rand(rsn, n_v, 1) < p_i);   % 61 dBi when high, 50 dBi otherwise
                V = table(ones(n_v, 1), true(n_v, 1), g, 'VariableNames', ["trial" "complete" "max_g21_dBi"]);
                r = run_stats_rev1(V, 60, 5);
                pr = periodic_runs_rev1(V, 60, 4);
                rmax(i) = r(4);
                pmax(i) = pr(1);
            end
            for H = [1 2 4 8]
                N = ceil(H*3600/t_v);
                NL(end + 1, :) = {m, a, th, H, N, 100*mean(rmax >= N), 100*mean(pmax >= N), ...
                    median(pmax), max(pmax), n_null}; %#ok<SAGROW>
            end
        end
    end
end
writetable(cell2table(NL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_vols", ...
    "null_run_ge_pct", "null_periodic_ge_pct", "null_periodic_median_vols", "null_periodic_max_vols", ...
    "n_null_trials"]), fullfile(out_dir, "lock_null_rev11.csv"));

report_figures_rev11(out_dir);

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
