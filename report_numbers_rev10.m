%REPORT_NUMBERS_REV10 Tables and figures for methodology report rev10.
%   Answers the independent recheck of revision 9
%   (docs/rev9_independent_review_rev1.md, findings T1 to T14) and writes
%   every evidence file that revision 10 adds to docs/figures. Runs under
%   the revision 6 defaults (ROC specification, 216 modes, ASR-11, exact
%   CDF). Reads ts_all_modes_rev8.csv (report_numbers_rev8).
%
%     per_hour_modes_rev10.csv    per-hour counts at 60 dBi for each mode
%                                 and path, with competition ranks (1 =
%                                 largest, ties share the best rank) and
%                                 0.5 deg revisits per hour (T4, T6)
%     per_hour_envelope_rev10.csv largest of each statistic, its mode, and
%                                 the rank and ratio of the CDF-worst mode
%     per_hour_corr_rev10.csv     Spearman correlation of each per-hour
%                                 count with the 0.5 deg revisits per hour
%     seeds_runs_rev10.csv        20 seeds x ten 24 h trials, three drift
%                                 models, two modes: runs and lock episodes
%                                 (stretches of volumes all within 5 dB of
%                                 the level that hold a volume at it),
%                                 censored flags (T1, T14)
%     seeds_summary_rev10.csv     spread over seeds
%     lock_trials_rev10.csv       fixed-offset trials: rate, phase slip per
%                                 volume, the k-volume slip m4 = min over
%                                 k = 1..4 of |k x slip| (modulo 360), the k
%                                 that sets it, longest run and episode (T2)
%     lock_m4_rev10.csv           longest runs and episodes binned by m4
%     lock_fit_rev10.csv          share of rate draws whose longest run or
%                                 episode reaches N volumes or H hours (T3)
%     timing_modes_rev10.csv      time series cost at eight modes, three
%                                 repeats each, no parallel pool (T8)
%     timing_rev10_machine.txt    machine and settings of the timings
%     fig_lock_vs_m4_rev10.png    longest run against m4
%
%   Parts: rev10_overrides.part = "all" (default), "seeds" (tables A and C,
%   with a pool), or "timing" (table D only; run it with no pool open so
%   the timings don't share the machine with workers). Set rev10_overrides
%   (struct: out_dir, n_seeds, drift_trials, drift_hours, time_hours,
%   n_repeats, modes_time, margin_dB, parallel_flag, part) for a small run.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also report_numbers_rev9, run_stats_rev1, report_figures_rev10.

if ~exist("rev10_overrides", "var"), rev10_overrides = struct(); end
here = fileparts(which("report_numbers_rev10"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "in_dir", fullfile(here, "docs", "figures"), ...
    "n_seeds", 20, "drift_trials", 10, "drift_hours", 24, "time_hours", 120, "n_repeats", 3, ...
    "modes_time", ["VCP 212 SAILSx3 AVSET 6.4", "VCP 12 SAILSx3 AVSET 6.4", ...
    "VCP 12 SAILSx3 AVSET 10.0, CD 12.524", "VCP 112 SAILSx1 AVSET 10.0", "VCP 34 SAILSx1", ...
    "VCP 35 SAILSx1", "VCP 32", "VCP 31"], "margin_dB", 5, "parallel_flag", 1, "part", "all");
f = fieldnames(rev10_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev10_overrides.(f{i});
end
out_dir = ov.out_dir;
if ~isfolder(out_dir), mkdir(out_dir); end
tiny = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1);
m_cdf = "VCP 212 SAILSx3 AVSET 6.4";
mc = vcp_mode_cases_rev1();
mc = mc(2:end);                                     % c{1} is the working-paper baseline, not a WSR-88D mode
labs = cellfun(@(x) string(x.label), mc);
do_seeds = any(ov.part == ["all", "seeds"]);
do_time = any(ov.part == ["all", "timing"]);

% ---- D first. Timings, no pool, before any pool opens (T8)
if do_time
    if ~isempty(gcp("nocreate"))
        error("report_numbers_rev10:poolOpen", "Close the parallel pool before the timing part.");
    end
    rst = RandStream("mrg32k3a", "Seed", 3);
    eb = 2.5*rand(rst, 32, 1) - 0.5;
    tm = ov.modes_time(ismember(ov.modes_time, labs));
    TT = {};
    for m = tm
        w = nexrad_asr11_mag_rev1(merge(tiny, rmfield(mc{labs == m}, "label")));
        lopt = struct("seed", 20261105, "n_trials", ov.time_hours, "duration_s", 3600);
        t = zeros(ov.n_repeats, 3);
        for r = 1:ov.n_repeats
            tic; lib = scan_event_library_rev1(w.scan_max.clock_spec, lopt); t(r, 1) = toc;
            tic; mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0); t(r, 2) = toc;
            tic; mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, eb, eb); t(r, 3) = toc/32;
        end
        TT(end + 1, :) = {m, w.scan_max.t_volume_s, median(t(:, 1)), min(t(:, 1)), max(t(:, 1)), ...
            median(t(:, 2)), min(t(:, 2)), max(t(:, 2)), median(t(:, 3)), min(t(:, 3)), max(t(:, 3))}; %#ok<SAGROW>
    end
    writetable(cell2table(TT, 'VariableNames', ["mode", "t_volume_s", "library_s", "library_min_s", ...
        "library_max_s", "one_pair_s", "one_pair_min_s", "one_pair_max_s", "per_pair_batch32_s", ...
        "per_pair_batch32_min_s", "per_pair_batch32_max_s"]), fullfile(out_dir, "timing_modes_rev10.csv"));
    fid = fopen(fullfile(out_dir, "timing_rev10_machine.txt"), "w");
    fprintf(fid, "%s; %s; computational threads %d; %s; parallel pool open: 0; %d repeats, median; " + ...
        "%d h per timed time series; single wall-clock runs (tic/toc)\n", computer, version, ...
        maxNumCompThreads, string(datetime("now")), ov.n_repeats, ov.time_hours);
    fclose(fid);
end
if ~do_seeds
    return
end
if ov.parallel_flag == 1 && isempty(gcp("nocreate"))
    parpool("Processes");
end
use_par = ~isempty(gcp("nocreate"));

% ---- A. Per-hour statistics with competition ranks (T4, T6)
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
P = table(T_all.mode, T_all.el_deg, T_all.t_volume_s, n05, n05.*vph, ...
    T_all.vol_pct_60/100.*vph, T_all.cuts_hit_60_pct/100.*n05.*vph, T_all.rev_pct_60/100*3600/t_rev_atc, ...
    T_all.vol_pct_60, T_all.cuts_hit_60_pct, T_all.rev_pct_60, ...
    'VariableNames', ["mode", "el_deg", "t_volume_s", "cuts05_per_volume", "revisits05_per_h", ...
    "vols_per_h_60", "cuts_per_h_60", "revs_per_h_60", "vol_pct_60", "cuts_hit_60_pct", "rev_pct_60"]);
mets = ["revisits05_per_h", "vols_per_h_60", "cuts_per_h_60", "revs_per_h_60", "vol_pct_60", ...
    "cuts_hit_60_pct", "rev_pct_60"];
for q = mets
    P.("rank_" + q) = zeros(height(P), 1);
    for a = [0 0.5]
        k = find(P.el_deg == a);
        v = P.(q)(k);
        P.("rank_" + q)(k) = 1 + sum(v' > v + 1e-12*abs(v), 2);   % competition rank: 1 + modes strictly larger
    end
end
writetable(P, fullfile(out_dir, "per_hour_modes_rev10.csv"));
EV = {};
CR = {};
for a = [0 0.5]
    Pa = P(P.el_deg == a, :);
    for q = mets
        [mx, im] = max(Pa.(q));
        c = Pa(Pa.mode == m_cdf, :);
        EV(end + 1, :) = {a, q, mx, Pa.mode(im), c.(q), c.("rank_" + q), c.(q)/mx}; %#ok<SAGROW>
    end
    for q = ["vols_per_h_60", "cuts_per_h_60", "revs_per_h_60"]
        CR(end + 1, :) = {a, q, corr(Pa.revisits05_per_h, Pa.(q), "Type", "Spearman"), height(Pa)}; %#ok<SAGROW>
    end
end
writetable(cell2table(EV, 'VariableNames', ["el_deg", "metric", "max_over_modes", "worst_mode", ...
    "at_cdf_worst_mode", "rank_of_cdf_worst", "ratio_to_max"]), fullfile(out_dir, "per_hour_envelope_rev10.csv"));
writetable(cell2table(CR, 'VariableNames', ["el_deg", "count", "spearman_vs_revisits", "n_modes"]), ...
    fullfile(out_dir, "per_hour_corr_rev10.csv"));

% ---- C. Seeds: runs, lock episodes, k-volume slip (T1, T2, T3, T14)
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
gv = ov.margin_dB;
if use_par
    parfor i = 1:nC
        [rowsC{i}, rowsL{i}] = seed_runs(mc{labs == modesC(MM(i))}, tiny, models(DD(i)), steps(DD(i)), ...
            20261200 + SS(i), dt_, dh, gv); %#ok<PFBNS>
    end
else
    for i = 1:nC
        [rowsC{i}, rowsL{i}] = seed_runs(mc{labs == modesC(MM(i))}, tiny, models(DD(i)), steps(DD(i)), ...
            20261200 + SS(i), dt_, dh, gv);
    end
end
T_seed = vertcat(rowsC{:});
writetable(T_seed, fullfile(out_dir, "seeds_runs_rev10.csv"));
T_lock = vertcat(rowsL{:});
writetable(T_lock, fullfile(out_dir, "lock_trials_rev10.csv"));
step_key = T_seed.step_deg_s_per_rev;
step_key(isnan(step_key)) = 0;                      % findgroups drops NaN keys (the fixed-offset model)
G = findgroups(T_seed.mode, T_seed.drift_model, step_key, T_seed.el_deg, T_seed.threshold_dBi);
S = splitapply(@(md, dm, st, el, th, vp, p95, mxr, mxu, mxc, eps_, epu, epc) {md(1), dm(1), st(1), el(1), th(1), ...
    numel(vp), min(p95), max(p95), max(mxr), max(mxu), max(mxc), max(eps_), max(epu), max(epc)}, ...
    T_seed.mode, T_seed.drift_model, T_seed.step_deg_s_per_rev, T_seed.el_deg, T_seed.threshold_dBi, ...
    T_seed.vol_pct, T_seed.p95_run_vols, T_seed.max_run_vols, T_seed.max_uncensored_run_vols, ...
    T_seed.max_censored_run_vols, T_seed.max_episode_vols, T_seed.max_uncensored_episode_vols, ...
    T_seed.max_censored_episode_vols, G);
T_ss = cell2table(S, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "el_deg", "threshold_dBi", ...
    "n_seeds", "p95_run_min", "p95_run_max", "max_run", "max_uncensored_run", "max_censored_run", ...
    "max_episode", "max_uncensored_episode", "max_censored_episode"]);
writetable(T_ss, fullfile(out_dir, "seeds_summary_rev10.csv"));

% Longest run and episode against the k-volume slip m4, and the share of
% rate draws reaching N volumes or H hours.
FB = {};
FL = {};
edges = [0 2 4 8 16 45 180];
for m = modesC'
    t_v = T_all.t_volume_s(find(T_all.mode == m, 1));
    for a = [0 0.5]
        for th = [60 70]
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            for b = 1:numel(edges) - 1
                k = L.m4_deg >= edges(b) & L.m4_deg < edges(b + 1);
                FB(end + 1, :) = {m, a, th, edges(b), edges(b + 1), nnz(k), median(L.longest_run_vols(k)), ...
                    max([0; L.longest_run_vols(k)]), max([0; L.ep_span_vols(k)]), ...
                    100*mean(L.longest_run_vols(k) >= 10)}; %#ok<SAGROW>
            end
            for N = [6 10 14 20 40]
                kr = L.longest_run_vols >= N;
                ke = L.ep_span_vols >= N;
                FL(end + 1, :) = {m, a, th, "volumes", N, N*t_v/3600, 100*mean(kr), nnz(kr & L.longest_censored == 1), ...
                    100*mean(ke), nnz(ke & L.ep_censored == 1), height(L)}; %#ok<SAGROW>
            end
            for H = [1 2 4]
                N = ceil(H*3600/t_v);
                kr = L.longest_run_vols >= N;
                ke = L.ep_span_vols >= N;
                FL(end + 1, :) = {m, a, th, "hours", N, H, 100*mean(kr), nnz(kr & L.longest_censored == 1), ...
                    100*mean(ke), nnz(ke & L.ep_censored == 1), height(L)}; %#ok<SAGROW>
            end
        end
    end
end
writetable(cell2table(FB, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "m4_lo_deg", "m4_hi_deg", ...
    "n_trials", "median_longest_run", "max_longest_run", "max_episode", "trials_run_ge10_pct"]), ...
    fullfile(out_dir, "lock_m4_rev10.csv"));
writetable(cell2table(FL, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "by", "n_vols", "hours", ...
    "trials_run_ge_pct", "n_censored_runs", "trials_episode_ge_pct", "n_censored_episodes", "n_trials"]), ...
    fullfile(out_dir, "lock_fit_rev10.csv"));

report_figures_rev10(out_dir);

% =========================================================================
function o = merge(a, b)
o = a;
f = fieldnames(b);
for i = 1:numel(f)
    o.(f{i}) = b.(f{i});
end
end

function [T, L] = seed_runs(mcase, tiny, model, step, seed, n_tr, hours, margin_dB)
% Runs and episodes of one seed of n_tr trials of `hours` each, and, for
% the fixed-offset model, the rate, slips and longest run of each trial.
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
[slip, rate, m4, kb] = deal(NaN(n_tr, 1));
if model == "per_trial"
    for j = 1:n_tr
        C = scan_clock_rev1(spec, seed, lib.trials(j), dur);
        rate(j) = C.r;
        slip(j) = mod(C.r*t_vol + 180, 360) - 180;  % ATC azimuth shift per WSR-88D volume, deg
        [m4(j), kb(j)] = min(abs(mod((1:4)*slip(j) + 180, 360) - 180));   % k volumes near whole revolutions
    end
end
rows = {};
rowsL = {};
for a = 1:2
    for th = [60 70]
        [rs, tr] = run_stats_rev1(R(a).nexrad_volume, th, margin_dB);
        V = R(a).nexrad_volume(R(a).nexrad_volume.complete, :);
        rows(end + 1, :) = [{lab, model, step, seed, R(a).el_n_deg, th, 100*mean(V.max_g21_dBi >= th)}, ...
            num2cell(rs)]; %#ok<AGROW>
        if model == "per_trial"
            for j = 1:n_tr
                rowsL(end + 1, :) = {lab, seed, lib.trials(j), R(a).el_n_deg, th, rate(j), slip(j), m4(j), kb(j), ...
                    tr.longest(j), tr.censored(j), tr.ep_span(j), tr.ep_high(j), tr.ep_censored(j)}; %#ok<AGROW>
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "seed", "el_deg", ...
    "threshold_dBi", "vol_pct", "median_run_vols", "p95_run_vols", "n_runs", "max_run_vols", ...
    "share_runs_ge2_pct", "share_runs_ge3_pct", "max_uncensored_run_vols", "n_censored_runs", ...
    "max_censored_run_vols", "max_episode_vols", "max_uncensored_episode_vols", "max_censored_episode_vols"]);
vn = ["mode", "seed", "trial", "el_deg", "threshold_dBi", "rate_deg_s", "slip_deg_per_vol", "m4_deg", ...
    "m4_k", "longest_run_vols", "longest_censored", "ep_span_vols", "ep_high_vols", "ep_censored"];
L = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
if ~isempty(rowsL)
    L = cell2table(rowsL, 'VariableNames', vn);
end
end
