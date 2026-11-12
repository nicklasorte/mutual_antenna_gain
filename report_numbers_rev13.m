%REPORT_NUMBERS_REV13 Tables and figure for methodology report rev13.
%   Answers the independent recheck of revision 12
%   (docs/rev12_independent_review_rev1.md, findings W1 to W14). Reruns the
%   20-seed study of revisions 8 to 12 (same seeds) and saves per trial the
%   share of high volumes, the number of complete volumes, every high
%   volume inside the longest periodic run, and periodic runs with periods
%   up to 16 volumes. Adds:
%     - high volumes per hour inside each long periodic run, all of them,
%       not only the repeating ones (W1);
%     - a chance comparison that keeps each trial's own share of high
%       volumes and length (W2);
%     - the lock from the slip alone, any k from 1 to 16, beside the lock at
%       the run's own period (W3);
%     - periodic runs with periods up to 16 volumes, and the trials with a
%       tight lock over 5 to 16 volumes (W4).
%   Reads ts_all_modes_rev8.csv for volume times.
%
%     seeds_runs_rev13.csv        per seed: runs, episodes, periodic runs
%                                 (periods up to 4 and up to 16)
%     lock_trials_rev13.csv       per fixed-offset trial: share and number
%                                 of complete volumes, slip, m4, both locks,
%                                 periodic runs (4 and 16), rates
%     lock_k16_rev13.csv          long periodic runs: both locks, how many
%                                 trials have locks that tight, and rates
%     lock_long_k_rev13.csv       trials with a lock over 5 to 16 volumes
%                                 within 1 deg but a periodic run under 20
%                                 volumes at periods up to 4
%     lock_fit_rev13.csv          share of rate draws whose longest periodic
%                                 run reaches 1, 2, 4, or 8 hours (periods up
%                                 to 4 and to 16) against chance, with 95%
%                                 Clopper-Pearson intervals
%     lock_walk_null_rev13.csv    chance of the random-walk periodic maxima
%     fig_lock_any16_rev13.png    longest periodic run against the lock
%                                 from the slip alone
%
%   Set rev13_overrides (struct: out_dir, n_seeds, drift_trials,
%   drift_hours, n_null, parallel_flag, reuse_seeds) for a small run;
%   reuse_seeds = true reads seeds_runs_rev13.csv and lock_trials_rev13.csv.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also report_numbers_rev12, periodic_runs_rev1, lock_drift_multiple_rev1.

if ~exist("rev13_overrides", "var"), rev13_overrides = struct(); end
here = fileparts(which("report_numbers_rev13"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "in_dir", fullfile(here, "docs", "figures"), ...
    "n_seeds", 20, "drift_trials", 10, "drift_hours", 24, "parallel_flag", 1, "reuse_seeds", false, "n_null", 20000);
f = fieldnames(rev13_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev13_overrides.(f{i});
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
    rowsC = {readtable(fullfile(out_dir, "seeds_runs_rev13.csv"), "TextType", "string")};
    rowsL = {readtable(fullfile(out_dir, "lock_trials_rev13.csv"), "TextType", "string")};
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
writetable(T_seed, fullfile(out_dir, "seeds_runs_rev13.csv"));
T_lock = vertcat(rowsL{:});
writetable(T_lock, fullfile(out_dir, "lock_trials_rev13.csv"));

% Per-trial additions (W1, W3, W4): the lock behind the period-4 run (drift
% over multiples of its period up to 16 volumes, as in revision 12), the
% lock from the slip alone (any k from 1 to 16; what a rate tolerance
% gives), the run in hours, and its high volumes per hour two ways: the
% repeating volumes only, and every high volume inside the run's span.
[dk, kk] = lock_drift_multiple_rev1(T_lock.slip_deg_per_vol, T_lock.per_k, 16);
T_lock.lock_run_drift_deg = dk;
T_lock.lock_run_k = kk;
[da, ka] = lock_drift_multiple_rev1(T_lock.slip_deg_per_vol, ones(height(T_lock), 1), 16);
T_lock.lock_any16_drift_deg = da;
T_lock.lock_any16_k = ka;
t_v_row = zeros(height(T_lock), 1);
for m = modesC'
    t_v_row(T_lock.mode == m) = T_all.t_volume_s(find(T_all.mode == m, 1));
end
T_lock.per_hours = T_lock.per_vols.*t_v_row/3600;
T_lock.per16_hours = T_lock.per16_vols.*t_v_row/3600;
T_lock.per_pattern_per_h = T_lock.per_hits./max(T_lock.per_vols, 1)*3600./t_v_row;
T_lock.per_span_per_h = T_lock.per_span_high./max(T_lock.per_vols, 1)*3600./t_v_row;
T_lock.trial_per_h = T_lock.trial_vol_pct/100*3600./t_v_row;
writetable(T_lock, fullfile(out_dir, "lock_trials_rev13.csv"));

% Long periodic runs: their locks, measured two ways, and their rates (W1,
% W3). "Tight" counts trials whose lock is at or within the loosest lock
% among the long runs.
K2 = {};
for m = modesC'
    for a = [0 0.5]
        for th = [60 70]
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            avg_h = mean(L.trial_per_h);
            for H = [4 6.5 8]
                lg = L.per_hours >= H;
                if ~any(lg)
                    continue
                end
                dr = max(L.lock_run_drift_deg(lg));
                dA = max(L.lock_any16_drift_deg(lg));
                K2(end + 1, :) = {m, a, th, H, nnz(lg), dr, nnz(L.lock_run_drift_deg <= dr), dA, ...
                    nnz(L.lock_any16_drift_deg <= dA), height(L), nnz(lg & L.m4_deg < 8), avg_h, ...
                    median(L.per_pattern_per_h(lg)), median(L.per_span_per_h(lg)), ...
                    nnz(lg & L.per_span_per_h > avg_h)}; %#ok<SAGROW>
            end
        end
    end
end
if isempty(K2), K2 = cell(0, 15); end
writetable(cell2table(K2, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_long", ...
    "max_lock_run_drift_deg", "n_tight_run", "max_lock_any16_drift_deg", "n_tight_any16", "n_trials", ...
    "n_long_m4_below_8", "mode_avg_high_per_h", "median_pattern_per_h_long", "median_span_per_h_long", ...
    "n_long_span_above_avg"]), fullfile(out_dir, "lock_k16_rev13.csv"));

% Trials with a tight lock over 5 to 16 volumes but a short period-4 run
% (W4), and what periods up to 16 find in them.
W4 = {};
for m = modesC'
    L = T_lock(T_lock.mode == m & T_lock.el_deg == 0.5 & T_lock.threshold_dBi == 60, :);
    k = L.lock_any16_drift_deg <= 1 & L.lock_any16_k >= 5 & L.per_vols < 20;
    for r = find(k)'
        x = L(r, :);
        W4(end + 1, :) = {m, x.seed, x.trial, x.lock_any16_k, x.lock_any16_drift_deg, x.per_vols, x.per_k, ...
            x.per16_vols, x.per16_k, x.per16_hits, x.per16_censored, x.per16_hours}; %#ok<SAGROW>
    end
end
if isempty(W4), W4 = cell(0, 12); end
writetable(cell2table(W4, 'VariableNames', ["mode", "seed", "trial", "lock_k", "lock_drift_deg", "per_vols", ...
    "per_k", "per16_vols", "per16_k", "per16_hits", "per16_censored", "per16_hours"]), ...
    fullfile(out_dir, "lock_long_k_rev13.csv"));

% Chance (W2): independent volumes; each null trial draws the share of
% high volumes of one real fixed-offset trial in turn (each trial's own
% share, so the spread between trials is kept), with that trial's number
% of complete volumes. Periodic runs with periods up to 4 and up to 16.
n_null = ov.n_null;
FL = {};
WL = {};
cell_id = 0;
for m = modesC'
    t_v = T_all.t_volume_s(find(T_all.mode == m, 1));
    for a = [0 0.5]
        for th = [60 70]
            cell_id = cell_id + 1;
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            pv = L.trial_vol_pct/100;
            nv = L.trial_n_vols;
            [p4, p16] = deal(zeros(n_null, 1));
            base = 1e6*cell_id;
            if use_par
                parfor i = 1:n_null
                    [p4(i), p16(i)] = null_trial(pv, nv, i, base);
                end
            else
                for i = 1:n_null
                    [p4(i), p16(i)] = null_trial(pv, nv, i, base);
                end
            end
            n_m = height(L);
            for H = [1 2 4 8]
                N = ceil(H*3600/t_v);
                for kmax = [4 16]
                    if kmax == 4, mv = L.per_vols; nl = p4; else, mv = L.per16_vols; nl = p16; end
                    xm = nnz(mv >= N);
                    xn = nnz(nl >= N);
                    [pm, cm] = binofit(xm, n_m);
                    [pn, cn] = binofit(xn, n_null);
                    FL(end + 1, :) = {m, a, th, H, N, kmax, 100*pm, 100*cm(1), 100*cm(2), n_m, ...
                        100*pn, 100*cn(1), 100*cn(2), n_null, pm/pn, cm(1)/cn(2), cm(2)/cn(1)}; %#ok<SAGROW>
                end
            end
            % Random-walk maxima against chance (periods up to 4).
            for st = [0.01 0.05]
                w = T_seed.mode == m & T_seed.drift_model == "random_walk" & T_seed.step_deg_s_per_rev == st & ...
                    T_seed.el_deg == a & T_seed.threshold_dBi == th;
                obs = max(T_seed.max_periodic_vols(w));
                p1 = mean(p4 >= obs);
                n_tr = 10*nnz(w);
                WL(end + 1, :) = {m, a, th, st, obs, obs*t_v/3600, p1, 1 - (1 - p1)^n_tr, n_tr}; %#ok<SAGROW>
            end
        end
    end
end
writetable(cell2table(FL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_vols", ...
    "max_period", "trials_periodic_ge_pct", "periodic_lo_pct", "periodic_hi_pct", "n_trials", ...
    "null_periodic_ge_pct", "null_lo_pct", "null_hi_pct", "n_null_trials", "ratio", "ratio_lo", "ratio_hi"]), ...
    fullfile(out_dir, "lock_fit_rev13.csv"));
writetable(cell2table(WL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "step_deg_s_per_rev", ...
    "max_periodic_vols", "max_periodic_hours", "null_trial_ge_pct", "chance_any_of_trials", "n_trials"]), ...
    fullfile(out_dir, "lock_walk_null_rev13.csv"));

report_figures_rev13(out_dir);

% -------------------------------------------------------------------------
function [p4, p16] = null_trial(pv, nv, i, base)
% Longest periodic run (periods up to 4 and up to 16) of one null trial of
% independent volumes with the share and length of real trial
% 1 + mod(i - 1, numel(pv)).
j = 1 + mod(i - 1, numel(pv));
rs = RandStream("mrg32k3a", "Seed", 20261113);
rs.Substream = base + i;
n_v = nv(j);
g = 50 + 11*(rand(rs, n_v, 1) < pv(j));          % 61 dBi when high, 50 dBi otherwise
V = table(ones(n_v, 1), true(n_v, 1), g, 'VariableNames', ["trial" "complete" "max_g21_dBi"]);
a = periodic_runs_rev1(V, 60, 4);
b = periodic_runs_rev1(V, 60, 16);
p4 = a(1);
p16 = b(1);
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
        [pr16, tp16] = periodic_runs_rev1(Vt, th, 16);
        V = Vt(Vt.complete, :);
        rows(end + 1, :) = {lab, model, step, seed, R(a).el_n_deg, th, 100*mean(V.max_g21_dBi >= th), ...
            r5(1), r5(2), r5(4), r5(7), r5(9), r3(10), r5(10), r5(11), r8(10), pr(1), pr(2), pr(4), pr(5), ...
            pr16(1), pr16(2)}; %#ok<AGROW>
        if model == "per_trial"
            tu = unique(V.trial)';
            for j = 1:n_tr
                hit = V.max_g21_dBi(V.trial == tu(j)) >= th;
                sh = 0;
                if tp.first(j) > 0
                    sh = sum(hit(tp.first(j):tp.last(j)));   % every high volume inside the run's span
                end
                rowsL(end + 1, :) = {lab, seed, lib.trials(j), R(a).el_n_deg, th, 100*mean(hit), numel(hit), ...
                    rate(j), slip(j), m4(j), kb(j), t5.longest(j), t5.censored(j), t3.ep_span(j), t5.ep_span(j), ...
                    t5.ep_censored(j), t8.ep_span(j), tp.dur(j), tp.k(j), tp.hits(j), tp.censored(j), sh, ...
                    tp16.dur(j), tp16.k(j), tp16.hits(j), tp16.censored(j)}; %#ok<AGROW>
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "seed", "el_deg", ...
    "threshold_dBi", "vol_pct", "median_run_vols", "p95_run_vols", "max_run_vols", "max_uncensored_run_vols", ...
    "max_censored_run_vols", "max_episode3_vols", "max_episode5_vols", "max_uncensored_episode5_vols", ...
    "max_episode8_vols", "max_periodic_vols", "max_periodic_k", "max_uncensored_periodic_vols", ...
    "max_censored_periodic_vols", "max_periodic16_vols", "max_periodic16_k"]);
vn = ["mode", "seed", "trial", "el_deg", "threshold_dBi", "trial_vol_pct", "trial_n_vols", "rate_deg_s", ...
    "slip_deg_per_vol", "m4_deg", "m4_k", "longest_run_vols", "longest_censored", "ep3_vols", "ep5_vols", ...
    "ep5_censored", "ep8_vols", "per_vols", "per_k", "per_hits", "per_censored", "per_span_high", ...
    "per16_vols", "per16_k", "per16_hits", "per16_censored"];
L = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
if ~isempty(rowsL)
    L = cell2table(rowsL, 'VariableNames', vn);
end
end
