%REPORT_NUMBERS_REV8 Tables and figures for methodology report rev8.
%   Answers the independent review of revision 7
%   (docs/rev7_independent_review_rev1.md) and writes every evidence file
%   that revision 8 cites to docs/figures (finding Q7: no cited number
%   without a script). Runs under the revision 6 defaults (ROC
%   specification, 216 modes, ASR-11, exact CDF) unless a table says
%   otherwise.
%
%     takeoff_worst_rev8.csv     largest envelope percentile over takeoff
%                                angles -1 to 3 deg on a 0.05 deg grid, the
%                                range of angles within 0.01 dB of it (Q11)
%     worst_modes_rev8.csv       mode that sets each envelope percentile
%                                and how many modes tie with it (Q9)
%     ts_all_modes_rev8.csv      120 h time series for every one of the 216
%                                modes at 0/0 and 0.5/0.5 deg: persistence,
%                                run statistics and harm metric (Q1)
%     ts_envelope_rev8.csv       largest value of each time-series metric
%                                over the 216 modes and the mode setting it
%     vol_exceed_modes_rev8.csv  per-volume exceedance curve, envelope over
%                                the modes and the CDF-worst mode
%     seeds_runs_rev8.csv        run statistics for each of 20 seeds of ten
%                                24 h trials, three drift models, the
%                                persistence-worst and CDF-worst modes (Q2)
%     seeds_summary_rev8.csv     spread of those statistics over seeds (Q3)
%     noaa_case_rev8.csv         NOAA case: envelope over 216 modes, value
%                                at 0/0 and largest over the takeoff surface,
%                                with and without a supplemental cut (Q4, Q5)
%     worked_example_rev8.csv    worst WSR-88D takeoff angle example (0.5/0
%                                deg): CDF, time series at the CDF-worst and
%                                persistence-worst modes
%     ref_levels_rev8.csv        reference levels at 0.5/0 deg with a 120 h
%                                revolution share (Q6, Q7)
%     bootstrap_rev8.csv         95% interval of the share of volumes
%                                reaching each level, headline run
%     trial_length_rev8.csv      1 h against 12 h trials, same model
%     timing_rev8.csv            timings and the machine they ran on (Q12)
%     fig12_volume_exceedance.png, fig13_runs_vs_drift.png, fig14_harm.png
%
%   The time-series sweeps run in parfor when a pool is open (or
%   rev8_overrides.parallel_flag = 1 opens one). Set rev8_overrides
%   (struct: out_dir, n_trials, duration_s, n_seeds, drift_trials,
%   drift_hours, grid_step, n_modes, parallel_flag) before running for a
%   small run.
%
%   See also report_numbers_rev7, mag_pair_timeseries_rev1, mag_pair_cdf_rev1.

if ~exist("rev8_overrides", "var"), rev8_overrides = struct(); end
here = fileparts(which("report_numbers_rev8"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "n_trials", 120, "duration_s", 3600, ...
    "n_seeds", 20, "drift_trials", 10, "drift_hours", 24, "grid_step", 0.05, "n_modes", Inf, ...
    "parallel_flag", 1);
f = fieldnames(rev8_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev8_overrides.(f{i});
end
out_dir = ov.out_dir;
if ~isfolder(out_dir), mkdir(out_dir); end
pct = [50 90 99 99.9 99.99 99.999];
tiny = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1);
if ov.parallel_flag == 1 && isempty(gcp("nocreate"))
    parpool("Processes");
end
use_par = ~isempty(gcp("nocreate"));

% ---- Tables for the CDF parts
cases = ["roc", "roc_m3", "roc_m6", "nssl", "roc_asr9", "noaa"];
copt = {struct(), struct("nexrad_sl_margin_dB", 3), struct("nexrad_sl_margin_dB", 6), ...
    struct("nexrad_pattern", "nssl_measured"), struct("atc_radar", "asr9"), struct("asr_remote_sl_dBi", 6.5)};
[scans, names] = nexrad_mode_scans_rev1();
g_el = -1:0.05:3;
Tn = cell(numel(cases), 1);
Ta = cell(numel(cases), 1);
for k = 1:numel(cases)
    o = merge(tiny, copt{k});
    s = nexrad_asr11_mag_rev1(o);
    Tn{k} = gain_dist_table_rev1(s.sys_nexrad.tx, scans, struct("mode_names", names, "el_grid_deg", g_el));
    Ta{k} = gain_dist_table_rev1(s.sys_asr11.tx, s.sys_asr11.scan, struct("el_grid_deg", g_el));
end

% ---- A. Largest value over the takeoff surface (Q11) and worst modes (Q9)
e = (-1:ov.grid_step:3)';
[EN, EA] = ndgrid(e, e);
W = cell(numel(cases), 1);
for k = 1:numel(cases)
    R = mag_pair_cdf_rev1(Ta{k}, Tn{k}, EA(:), EN(:), struct("pct", pct));
    P = R.pct{:, 3:2 + numel(pct)};
    z = P(EN(:) == 0 & EA(:) == 0, :);
    rows = cell(numel(pct), 1);
    for j = 1:numel(pct)
        mx = max(P(:, j));
        at = abs(P(:, j) - mx) < 0.01;
        rows{j} = {cases(k), pct(j), z(j), mx, min(EN(at)), max(EN(at)), min(EA(at)), max(EA(at)), nnz(at)};
    end
    W{k} = cell2table(vertcat(rows{:}), 'VariableNames', ["case", "pct", "at_0_0_dBi", "worst_dBi", ...
        "el_nexrad_min_deg", "el_nexrad_max_deg", "el_atc_min_deg", "el_atc_max_deg", "n_grid_nodes"]);
end
writetable(vertcat(W{:}), fullfile(out_dir, "takeoff_worst_rev8.csv"));

ang = [0 0; 0.5 0.5; 0.5 0];
WM = {};
for k = 1:numel(cases)
    R = mag_pair_cdf_rev1(Ta{k}, Tn{k}, ang(:, 2), ang(:, 1), struct("pct", pct, "per_mode", true));
    for a = 1:size(ang, 1)
        for j = 1:numel(pct)
            qm = squeeze(R.mode_pct(a, :, j));
            env = R.pct{a, 2 + j};
            WM(end + 1, :) = {cases(k), ang(a, 1), ang(a, 2), pct(j), env, ...
                extractAfter(R.worst_mode(a, j), " | "), nnz(abs(qm - env) < 1e-9)}; %#ok<SAGROW>
        end
    end
end
T_wm = cell2table(WM, 'VariableNames', ["case", "el_nexrad_deg", "el_atc_deg", "pct", "env_dBi", ...
    "worst_mode", "n_modes_tied"]);
writetable(T_wm, fullfile(out_dir, "worst_modes_rev8.csv"));

% ---- B. Time series for every mode (Q1), ROC, ASR-11, 0/0 and 0.5/0.5
mc = vcp_mode_cases_rev1();
mc = mc(2:end);                                     % c{1} is the working-paper baseline, not a WSR-88D mode
labs = cellfun(@(x) string(x.label), mc);
m_cdf = "VCP 212 SAILSx3 AVSET 6.4";
idxB = 1:numel(mc);
if isfinite(ov.n_modes)                             % small runs: a few modes plus the CDF-worst mode
    idxB = unique([1:ov.n_modes - 1, find(labs == m_cdf, 1)]);
end
n_m = numel(idxB);
thr = [40 45 50 55 60 65 70];
gx = (40:0.1:80)';
rowsB = cell(n_m, 1);
exB = cell(n_m, 1);
tsopt = struct("persist_thr_dBi", thr, "radial_thr_dBi", [50 60 70], "keep_scans", true);
nt = ov.n_trials;
ds = ov.duration_s;
if use_par
    parfor m = 1:n_m
        [rowsB{m}, exB{m}] = mode_ts(mc{idxB(m)}, tiny, nt, ds, tsopt, gx); %#ok<PFBNS>
    end
else
    for m = 1:n_m
        [rowsB{m}, exB{m}] = mode_ts(mc{idxB(m)}, tiny, nt, ds, tsopt, gx);
    end
end
T_all = vertcat(rowsB{:});
writetable(T_all, fullfile(out_dir, "ts_all_modes_rev8.csv"));
metrics = ["vol_pct_50", "vol_pct_55", "vol_pct_60", "vol_pct_70", "rev_pct_60", "max_run_60_vols", ...
    "share_runs_ge2_60_pct", "cuts_hit_50_pct", "cuts_hit_60_pct", "cuts_hit_70_pct", "max_radials_60"];
EV = {};
for a = [0 0.5]
    Ta_ = T_all(T_all.el_deg == a, :);
    for q = metrics
        [mx, im] = max(Ta_.(q));
        EV(end + 1, :) = {a, q, mx, Ta_.mode(im), Ta_.(q)(Ta_.mode == m_cdf)}; %#ok<SAGROW>
    end
end
T_env = cell2table(EV, 'VariableNames', ["el_deg", "metric", "max_over_modes", "worst_mode", "at_cdf_worst_mode"]);
writetable(T_env, fullfile(out_dir, "ts_envelope_rev8.csv"));
X = cat(3, exB{:});                                  % gains x [0, 0.5] x modes
i_cdf = find(labs(idxB) == m_cdf, 1);
T_ex = table(gx, max(X(:, 1, :), [], 3), squeeze(X(:, 1, i_cdf)), max(X(:, 2, :), [], 3), squeeze(X(:, 2, i_cdf)), ...
    'VariableNames', ["gain_dBi", "env_0_pct", "cdf_worst_0_pct", "env_05_pct", "cdf_worst_05_pct"]);
writetable(T_ex, fullfile(out_dir, "vol_exceed_modes_rev8.csv"));

% ---- C. Run lengths over seeds (Q2, Q3)
m_vol = T_env.worst_mode(T_env.el_deg == 0.5 & T_env.metric == "vol_pct_60");
modesC = unique([m_vol, m_cdf], "stable");
models = ["per_trial", "random_walk", "random_walk"];
steps = [NaN 0.01 0.05];
[MM, DD, SS] = ndgrid(1:numel(modesC), 1:numel(models), 1:ov.n_seeds);
nC = numel(MM);
rowsC = cell(nC, 1);
dt_ = ov.drift_trials;
dh = ov.drift_hours;
if use_par
    parfor i = 1:nC
        rowsC{i} = seed_runs(mc{labs == modesC(MM(i))}, tiny, models(DD(i)), steps(DD(i)), ...
            20261200 + SS(i), dt_, dh); %#ok<PFBNS>
    end
else
    for i = 1:nC
        rowsC{i} = seed_runs(mc{labs == modesC(MM(i))}, tiny, models(DD(i)), steps(DD(i)), ...
            20261200 + SS(i), dt_, dh);
    end
end
T_seed = vertcat(rowsC{:});
writetable(T_seed, fullfile(out_dir, "seeds_runs_rev8.csv"));
step_key = T_seed.step_deg_s_per_rev;
step_key(isnan(step_key)) = 0;                      % findgroups drops NaN keys (the fixed-offset model)
G = findgroups(T_seed.mode, T_seed.drift_model, step_key, T_seed.el_deg, T_seed.threshold_dBi);
S = splitapply(@(md, dm, st, el, th, vp, med, p95, mxr, s2, s3) {md(1), dm(1), st(1), el(1), th(1), numel(vp), ...
    min(vp), max(vp), prctile(med, 50), min(p95), max(p95), min(mxr), prctile(mxr, 50), max(mxr), ...
    min(s2), max(s2), min(s3), max(s3)}, T_seed.mode, T_seed.drift_model, T_seed.step_deg_s_per_rev, ...
    T_seed.el_deg, T_seed.threshold_dBi, T_seed.vol_pct, T_seed.median_run_vols, T_seed.p95_run_vols, ...
    T_seed.max_run_vols, T_seed.share_runs_ge2_pct, T_seed.share_runs_ge3_pct, G);
T_ss = cell2table(S, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "el_deg", "threshold_dBi", ...
    "n_seeds", "vol_pct_min", "vol_pct_max", "median_run_median", "p95_run_min", "p95_run_max", ...
    "max_run_min", "max_run_median", "max_run_max", "share_ge2_min", "share_ge2_max", "share_ge3_min", "share_ge3_max"]);
writetable(T_ss, fullfile(out_dir, "seeds_summary_rev8.csv"));

% ---- D. NOAA case over the takeoff surface (Q4, Q5)
noaa = {"baseline", 1, []; "noaa", 6, []; "noaa_cut_0.2", 6, 0.2; "noaa_cut_-0.2", 6, -0.2};
ND = {};
for r = 1:size(noaa, 1)
    k = noaa{r, 2};
    Tn_r = Tn{k};
    if ~isempty(noaa{r, 3})
        [sc, nm] = nexrad_mode_scans_rev1("all", low_el_deg=noaa{r, 3}, low_el_mode="add");
        s = nexrad_asr11_mag_rev1(merge(tiny, copt{k}));
        Tn_r = gain_dist_table_rev1(s.sys_nexrad.tx, sc, struct("mode_names", nm, "el_grid_deg", g_el));
    end
    R = mag_pair_cdf_rev1(Ta{k}, Tn_r, EA(:), EN(:), struct("pct", pct));
    P = [R.pct{:, 3:2 + numel(pct)}, R.pct.max_on_path_dBi];
    z = P(EN(:) == 0 & EA(:) == 0, :);
    h = P(EN(:) == 0.5 & EA(:) == 0.5, :);
    [mx, im] = max(P, [], 1);
    ND(end + 1, :) = [{noaa{r, 1}, "at_0_0"}, num2cell(z), {0, 0}]; %#ok<SAGROW>
    ND(end + 1, :) = [{noaa{r, 1}, "at_0.5_0.5"}, num2cell(h), {0.5, 0.5}]; %#ok<SAGROW>
    ND(end + 1, :) = [{noaa{r, 1}, "largest"}, num2cell(mx), {NaN, NaN}]; %#ok<SAGROW>
    ND(end + 1, :) = [{noaa{r, 1}, "largest_at_el_nexrad"}, num2cell(reshape(EN(im), 1, [])), {NaN, NaN}]; %#ok<SAGROW>
    ND(end + 1, :) = [{noaa{r, 1}, "largest_at_el_atc"}, num2cell(reshape(EA(im), 1, [])), {NaN, NaN}]; %#ok<SAGROW>
end
T_noaa = cell2table(ND, 'VariableNames', ["case", "row", "p50", "p90", "p99", "p99_9", "p99_99", "p99_999", ...
    "max_on_path_dBi", "el_nexrad_deg", "el_atc_deg"]);
writetable(T_noaa, fullfile(out_dir, "noaa_case_rev8.csv"));

% ---- E. Worked example at the worst WSR-88D takeoff angle (0.5/0 deg)
EXR = {};
for a = [0 0.5]
    for m = unique([m_cdf, m_vol], "stable")
        q = merge(merge(tiny, rmfield(mc{labs == m}, "label")), struct("el_nexrad_deg", a, "el_atc_deg", 0));
        w = nexrad_asr11_mag_rev1(q);
        lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 20261105, ...
            "n_trials", ov.n_trials, "duration_s", ov.duration_s));
        Rt = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, a, 0, tsopt);
        P = Rt.persistence_volume;
        H = Rt.radial_summary;
        rs = run_stats(Rt.nexrad_volume, 60);
        g = w.cdf.pct21.gain_dBi';
        EXR(end + 1, :) = {a, 0, m, g(3), g(4), g(5), g(6), w.ref_levels.mag_dBi(1), ...
            P.vol_pct(P.threshold_dBi == 60), P.vol_pct(P.threshold_dBi == 70), rs(4), rs(5), ...
            H.cuts_with_hit_pct(2), H.max_radials_per_cut(2)}; %#ok<SAGROW>
    end
end
T_wex = cell2table(EXR, 'VariableNames', ["el_nexrad_deg", "el_atc_deg", "mode", "p99", "p99_9", "p99_99", ...
    "p99_999", "mbmb_on_path_dBi", "vol_pct_60", "vol_pct_70", "max_run_60_vols", "share_runs_ge2_60_pct", ...
    "cuts_hit_60_pct", "max_radials_60"]);
writetable(T_wex, fullfile(out_dir, "worked_example_rev8.csv"));
q = merge(merge(tiny, rmfield(mc{labs == m_cdf}, "label")), struct("el_nexrad_deg", 0.5, "el_atc_deg", 0, ...
    "ts_method", "events", "n_trials", ov.n_trials, "duration_s", ov.duration_s, "seed", 20261105));
w = nexrad_asr11_mag_rev1(q);
writetable(w.ref_levels, fullfile(out_dir, "ref_levels_rev8.csv"));

% ---- F. Bootstrap and trial length (headline run, ROC, ASR-11, 0/0)
q = merge(tiny, rmfield(mc{labs == m_cdf}, "label"));
w = nexrad_asr11_mag_rev1(q);
lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 20261105, "n_trials", ov.n_trials, ...
    "duration_s", ov.duration_s));
Rh = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, [0; 0.5], [0; 0.5], tsopt);
thr_b = [50 55 60 70];
V = Rh(1).nexrad_volume(Rh(1).nexrad_volume.complete, :);
tr = unique(V.trial);
hit = zeros(numel(tr), numel(thr_b));
nv = zeros(numel(tr), 1);
for j = 1:numel(tr)
    v = V.max_g21_dBi(V.trial == tr(j));
    nv(j) = numel(v);
    hit(j, :) = sum(v >= thr_b, 1);
end
rsb = RandStream("mrg32k3a", "Seed", 21);
B = randi(rsb, numel(tr), numel(tr), 2000);
st = zeros(2000, numel(thr_b));
for b = 1:2000
    st(b, :) = sum(hit(B(:, b), :), 1)/sum(nv(B(:, b)));
end
writetable(table(thr_b', 100*sum(hit, 1)'/sum(nv), 100*prctile(st, 2.5)', 100*prctile(st, 97.5)', ...
    'VariableNames', ["threshold_dBi", "vol_pct", "ci95_lo_pct", "ci95_hi_pct"]), fullfile(out_dir, "bootstrap_rev8.csv"));
lib12 = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 20261105, "n_trials", 10, ...
    "duration_s", 12*3600));
R12 = mag_pair_timeseries_rev1(lib12, w.sys_nexrad.tx, w.sys_asr11.tx, [0; 0.5], [0; 0.5], tsopt);
TL = {};
for a = 1:2
    for th = [55 60 70]
        for L = 1:2
            if L == 1, Rr = Rh(a); h_ = 1; else, Rr = R12(a); h_ = 12; end
            P = Rr.persistence_volume;
            TL(end + 1, :) = {h_, Rr.el_n_deg, th, P.vol_pct(P.threshold_dBi == th), ...
                P.max_run_vols(P.threshold_dBi == th), P.max_run_s(P.threshold_dBi == th)/60, ...
                P.censored_runs(P.threshold_dBi == th)}; %#ok<SAGROW>
        end
    end
end
writetable(cell2table(TL, 'VariableNames', ["trial_h", "el_deg", "threshold_dBi", "vol_pct", "max_run_vols", ...
    "max_run_min", "censored_runs"]), fullfile(out_dir, "trial_length_rev8.csv"));

% ---- G. Timings (Q12)
rst = RandStream("mrg32k3a", "Seed", 3);
n_p = 500;
ea = -5 + 15*rand(rst, n_p, 1);
en = -5 + 15*rand(rst, n_p, 1);
g5 = -5:0.05:10;
s = nexrad_asr11_mag_rev1(tiny);
tic;
Tn5 = gain_dist_table_rev1(s.sys_nexrad.tx, scans, struct("mode_names", names, "el_grid_deg", g5));
t_tab = toc;
Ta5 = gain_dist_table_rev1(s.sys_asr11.tx, s.sys_asr11.scan, struct("el_grid_deg", g5));
t_bis = timeit(@() mag_pair_cdf_rev1(Ta5, Tn5, ea, en))/n_p;
t_full = timeit(@() mag_pair_cdf_rev1(Ta5, Tn5, ea, en, struct("pct_search", "full")))/n_p;
t_one = timeit(@() mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0));
eb = 2.5*rand(rst, 32, 1) - 0.5;
t_b = timeit(@() mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, eb, eb))/32;
t_lib = timeit(@() scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 1, "n_trials", ov.n_trials, ...
    "duration_s", ov.duration_s)));
writetable(table(["WSR-88D table, 216 modes, 301 takeoff angles"; "Pair CDF, bisection (per pair, 500 pairs)"; ...
    "Pair CDF, full search (per pair, 500 pairs)"; "Time series, one pair"; "Time series, per pair in a batch of 32"; ...
    "Event library, one pair of modes"], [t_tab; t_bis; t_full; t_one; t_b; t_lib], ...
    'VariableNames', ["item", "seconds"]), fullfile(out_dir, "timing_rev8.csv"));
fid = fopen(fullfile(out_dir, "timing_rev8_machine.txt"), "w");
fprintf(fid, "%s; %s; computational threads %d; %s; parallel pool for sweeps: %d; %s\n", computer, version, ...
    maxNumCompThreads, string(datetime("now")), use_par, sprintf("%d h per time series", ov.n_trials*ov.duration_s/3600));
fclose(fid);

% ---- H. Figures 12 to 14, drawn from the tables just written
report_figures_rev8(out_dir);

% =========================================================================
function o = merge(a, b)
o = a;
f = fieldnames(b);
for i = 1:numel(f)
    o.(f{i}) = b.(f{i});
end
end

function [row, ex] = mode_ts(mcase, tiny, nt, ds, tsopt, gx)
% 120 h time series of one scan mode at 0/0 and 0.5/0.5 deg.
lab = string(mcase.label);
w = nexrad_asr11_mag_rev1(merge(tiny, rmfield(mcase, "label")));
lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 20261105, "n_trials", nt, "duration_s", ds));
R = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, [0; 0.5], [0; 0.5], tsopt);
rows = cell(2, 1);
ex = zeros(numel(gx), 2);
for a = 1:2
    P = R(a).persistence_volume;
    PR = R(a).persistence;
    H = R(a).radial_summary;
    rs = run_stats(R(a).nexrad_volume, 60);
    v = R(a).nexrad_volume.max_g21_dBi(R(a).nexrad_volume.complete);
    ex(:, a) = 100*mean(v(:)' >= gx, 2);
    pv = @(t) P.vol_pct(P.threshold_dBi == t);
    rows{a} = table(lab, R(a).el_n_deg, w.scan_max.t_volume_s, pv(50), pv(55), pv(60), pv(70), ...
        PR.rev_pct(PR.threshold_dBi == 60), P.max_run_vols(P.threshold_dBi == 60), ...
        P.censored_runs(P.threshold_dBi == 60), rs(2), rs(5), ...
        H.cuts_with_hit_pct(1), H.cuts_with_hit_pct(2), H.cuts_with_hit_pct(3), H.max_radials_per_cut(2), ...
        'VariableNames', ["mode", "el_deg", "t_volume_s", "vol_pct_50", "vol_pct_55", "vol_pct_60", "vol_pct_70", ...
        "rev_pct_60", "max_run_60_vols", "censored_runs_60", "p95_run_60_vols", "share_runs_ge2_60_pct", ...
        "cuts_hit_50_pct", "cuts_hit_60_pct", "cuts_hit_70_pct", "max_radials_60"]);
end
row = vertcat(rows{:});
end

function T = seed_runs(mcase, tiny, model, step, seed, n_tr, hours)
% Run statistics of one seed of n_tr trials of `hours` each.
lab = string(mcase.label);
w = nexrad_asr11_mag_rev1(merge(tiny, rmfield(mcase, "label")));
spec = w.scan_max.clock_spec;
spec.jitter_mode = model;
if model == "random_walk"
    spec.drift_step_a = step;
end
lib = scan_event_library_rev1(spec, struct("seed", seed, "n_trials", n_tr, "duration_s", 3600*hours));
R = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, [0; 0.5], [0; 0.5], ...
    struct("persist_thr_dBi", [60 70], "radial_thr_dBi", 60, "keep_scans", true));
rows = {};
for a = 1:2
    for th = [60 70]
        rs = run_stats(R(a).nexrad_volume, th);
        V = R(a).nexrad_volume(R(a).nexrad_volume.complete, :);
        rows(end + 1, :) = {lab, model, step, seed, R(a).el_n_deg, th, 100*mean(V.max_g21_dBi >= th), ...
            rs(1), rs(2), rs(3), rs(4), rs(5), rs(6)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', ["mode", "drift_model", "step_deg_s_per_rev", "seed", "el_deg", ...
    "threshold_dBi", "vol_pct", "median_run_vols", "p95_run_vols", "n_runs", "max_run_vols", ...
    "share_runs_ge2_pct", "share_runs_ge3_pct"]);
end

function s = run_stats(Vt, thr)
% [median run, 95th-percentile run, n runs, longest run, share of high
% volumes in runs of 2 or more (%), share in runs of 3 or more (%)] of
% consecutive complete volumes reaching thr, per trial.
V = Vt(Vt.complete, :);
len = [];
for j = unique(V.trial)'
    hit = V.max_g21_dBi(V.trial == j) >= thr;
    d = diff([0; hit(:); 0]);
    len = [len; find(d == -1) - find(d == 1)]; %#ok<AGROW>
end
if isempty(len)
    s = [NaN NaN 0 0 0 0];
    return
end
s = [median(len), prctile(len, 95), numel(len), max(len), 100*sum(len(len >= 2))/sum(len), ...
    100*sum(len(len >= 3))/sum(len)];
end
