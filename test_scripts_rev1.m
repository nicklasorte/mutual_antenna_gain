function tests = test_scripts_rev1
%TEST_SCRIPTS_REV1 Tests for the example scripts, the C# crosscheck, the
%   report figure script, merge_overrides_rev1 and peak_label_rev1. Each script runs end to
%   end with small runs set through a workspace struct example_overrides,
%   so every line of each script executes. report_figures_rev5 writes to a
%   temporary folder, not docs/figures.
%   Run with runtests("test_scripts_rev1").
tests = functiontests(localfunctions);
end

function setup(tc)
tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
    ["radar_coupling:commensurateRates", "radar_coupling:equalRates"]));
end

function teardown(~)
close all force
end

% =========================================================================
function testMergeOverrides(tc)
o = struct("seed", 1, "n_samples", 5e6, "duration_s", 3600, "vcp", 212);
ov = struct("n_samples", 10, "duration_s", 5, "n_trials", 2, "dt_s", 0.001, "out_dir", "x", "vcp", 35);
m = merge_overrides_rev1(o, ov);
verifyEqual(tc, [m.n_samples m.duration_s m.n_trials m.dt_s], [10 5 2 0.001]);
verifyEqual(tc, m.vcp, 212, "only run-size fields are copied");
verifyFalse(tc, isfield(m, "out_dir"));
verifyEqual(tc, merge_overrides_rev1(o, struct()), o);
end

function testExampleRadarMutualGain(tc)
W = run_script("example_radar_mutual_gain_rev1", small());
verifyEqual(tc, numel(W.res), 3);
for i = 1:3
    verifyEqual(tc, W.res{i}.opts.n_samples, small().n_samples);
    verifyTrue(tc, all(isfinite(W.res{i}.pct21.gain_dBi)));
end
verifyEqual(tc, W.ts1.opts.mode, "time_series");
verifyTrue(tc, all(isfinite(W.sy.sync.max_g21_dBi)));
% Synchronized with offset 0: both beams point the same way, so one
% always looks away and the pair never reaches main beam to main beam.
verifyLessThan(tc, W.sy.sync.max_g21_dBi, 34 + 34 - 10);
end

function testCrosscheckCsharpPasses(tc)
W = run_script("crosscheck_csharp_rev1", struct());
T = W.T;
chk = T.pass ~= "info";
verifyGreaterThan(tc, sum(chk), 10);
verifyTrue(tc, all(T.pass(chk) == "true"), join(T.check(chk & T.pass ~= "true"), "; "));
end

function testExampleNexradAsr11Sweep(tc)
W = run_script("example_nexrad_asr11_rev1", small());
lab = cellfun(@(c) string(c.label), W.c);
verifyEqual(tc, numel(lab), 75, "74 VCP modes and the paper baseline");
verifyEqual(tc, numel(unique(lab)), numel(lab), "labels are unique");
verifyEqual(tc, height(W.T), 75);
verifyEqual(tc, height(W.V), 74);
% Every allowed combination is present; SAILS never goes with MRLE.
verifyTrue(tc, any(lab == "VCP 212 SAILSx3 AVSET 6.4"));
verifyTrue(tc, any(lab == "VCP 215 MRLEx2 AVSET 6.4, CD 12.524"));
verifyTrue(tc, any(lab == "VCP 31") && ~any(startsWith(lab, "VCP 31 ")));
verifyFalse(tc, any(contains(lab, "SAILS") & contains(lab, "MRLE")));
verifyFalse(tc, any(contains(lab, "VCP 212") & contains(lab, "CD 12.524")), ...
    "VCP 212 already runs its Doppler cut at 12.524 deg/s");
% The envelope bounds every mode at every percentile.
E = W.env.pct;
for k = 1:height(E)
    row = E{k, 4:end};
    verifyGreaterThanOrEqual(tc, E.envelope_dBi(k), max(row) - 1e-9);
end
verifyEqual(tc, W.V.case(W.i999), W.env.pct.worst_mode(W.env.pct.pct == 99.9));
% Results are stored lean: samples only for the paper baseline, VCP 215
% and the worst mode (rerun with the same seed), in single; g12 dropped.
i215 = find(W.T.case == "VCP 215");
iw = 1 + W.i999;
for i = 1:numel(W.res)
    verifyEmpty(tc, W.res{i}.cdf.g12_dBi);
    if ismember(i, [1 i215 iw])
        verifyClass(tc, W.res{i}.cdf.g21_dBi, "single");
    else
        verifyEmpty(tc, W.res{i}.cdf.g21_dBi);
    end
end
% The rerun reproduces the worst mode's sweep samples exactly.
verifyEqual(tc, prctile(double(W.res{iw}.cdf.g21_dBi), W.summ{iw}.pct(:)), W.summ{iw}.q_dBi, "AbsTol", 1e-5);
end

function testExampleSensitivity(tc)
W = run_script("example_sensitivity_rev1", small());
T = W.Tsens;
verifyEqual(tc, height(T), 19);
verifyEqual(tc, T.group(1), "Baseline");
verifyTrue(tc, all(T.max_run_45_vols >= 0));
verifyEqual(tc, T.d99_9_dB(1), 0);
verifyTrue(tc, all(ismember(["NEXRAD pattern", "NEXRAD back lobe", "ASR-11 pattern", ...
    "ASR-11 remote sidelobe", "ASR-11 tilt", "Path elevation", "ASR-11 rate tolerance", ...
    "NEXRAD sector"], T.group)));
% Same seed for every case: a higher ASR-11 floor raises every sample.
up = T.label == "+6.5 dBi";
dn = T.label == "+0.5 dBi";
verifyGreaterThanOrEqual(tc, T.p99_9_dBi(up), T.p99_9_dBi(1));
verifyLessThanOrEqual(tc, T.p99_9_dBi(dn), T.p99_9_dBi(1));
% Rate tolerance changes only the time series, not the random-mode CDF.
rt = T.group == "ASR-11 rate tolerance";
verifyEqual(tc, T.p99_9_dBi(rt), T.p99_9_dBi(1) + [0; 0], "AbsTol", 1e-9);
end

function testExampleDrift(tc)
W = run_script("example_asr11_drift_rev1", struct());
T = W.Tdrift;
verifyEqual(tc, height(T), 7);
verifyEqual(tc, T.ratio, ["1:1"; "2:1"; "6:1"; "74:1"; "1:1"; "1:1"; "74:1"]);
verifyEqual(tc, T.nominal_episode_s(1:2), [Inf; Inf]);
verifyEqual(tc, T.nominal_episode_s(3), 360/12.524, "RelTol", 1e-12, "capped at one NEXRAD cut");
verifyEqual(tc, T{3, ["p_gt_1min", "p_gt_10min", "p_gt_1h", "p_gt_1d", "p_gt_30d"]}, zeros(1, 5));
verifyTrue(tc, issorted(flip(T{1, 7:end})), "longer episodes are less likely");
% Volume lock: 74 ASR-11 revolutions ~ one 353.7 s volume. The ASR-11
% offset is uniform on +/-0.5 deg/s and delta = delta0 + d1, so for a
% window inside the range P(T_ep > t) = 2*(W/t)/1.
W = T.window_deg(4);
verifyEqual(tc, W, 1.45 + 74*0.925*(360/W_vol())/12.524, "RelTol", 1e-9);
verifyEqual(tc, T.p_gt_1h(4), 2*W/3600, "RelTol", 1e-9);
% ASR-9 rows (M.1464-2 Radar B, 1.3 deg): equal rates give a triangular
% drift on +/-1 deg/s, so P(T_ep > t) = 1 - (1 - W/t)^2.
verifyEqual(tc, T.window_deg(5:6), [2.6; 2.75], "RelTol", 1e-12);
verifyEqual(tc, T.p_gt_1h(5:6), 1 - (1 - [2.6; 2.75]/3600).^2, "RelTol", 1e-9);
verifyEqual(tc, T.window_deg(7), 1.3 + 74*0.925*(360/W_vol())/12.524, "RelTol", 1e-9);
verifyLessThan(tc, T.p_gt_1h(7), T.p_gt_1h(4), "narrower ASR-9 beam, shorter lock");
end

function testExampleAsr9(tc)
W = run_script("example_asr9_rev1", small());
verifyEqual(tc, height(W.T9), 74, "every VCP mode, no paper baseline");
verifyEqual(tc, W.asr9_base.atc_radar, "asr9");
verifyEqual(tc, W.w9.sys_asr11.name, "ASR-9");
verifyEqual(tc, W.w11.sys_asr11.name, "ASR-11");
verifyEqual(tc, W.cmp9_label, W.T9.case(W.j9));
verifyEqual(tc, W.T9.p99_9(W.j9), max(W.T9.p99_9));
verifyEqual(tc, W.env9.pct.envelope_dBi', max(W.T9{:, ["p99", "p99_9", "p99_99", "p99_999"]}, [], 1), ...
    "AbsTol", 0.05, "the envelope is the largest mode at each percentile");
C = W.Tcmp9;
verifyEqual(tc, height(C), 15);
verifyEqual(tc, C.asr9_minus_asr11(1), -0.5, "AbsTol", 1e-9, "peak gains 33.5 and 34 dBi");
verifyEqual(tc, C.asr9_minus_asr11(2), -0.5, "AbsTol", 1e-6, "main beam on the path");
verifyEqual(tc, C.asr9_minus_asr11(4), 0, "AbsTol", 1e-6, "same +3.5 dBi remote sidelobe");
verifyEqual(tc, C.asr9_minus_asr11(6), -0.5, "AbsTol", 1e-6);
verifyLessThan(tc, C.asr9(end), C.asr11(end), "narrower main-beam region");
verifyEqual(tc, W.w9.opts.n_trials, 1, "example_overrides sets the comparison trials");
end

function testPeakLabelBranches(tc)
% Every branch of the Figure 5 label: main beam = within 3 dB of the
% antenna's peak on the path (NEXRAD 41.99 dBi at 0.5 deg; ASR-11 34 dBi).
P = table([0.5; 0.9], [41.99; 34.14], 'VariableNames', ["el_deg", "nexrad_gain_on_path_dBi"]);
mk = @(el, gn, ga) table(gn + ga, el, gn, ga, 0.7, 0.3, 'VariableNames', ["max_g21_dBi", ...
    "nexrad_el_at_max_deg", "nexrad_gain_at_max_dBi", "asr11_gain_at_max_dBi", ...
    "nexrad_az_off_at_max_deg", "asr11_az_off_at_max_deg"]);
verifyTrue(tc, startsWith(peak_label_rev1(mk(0.5, 41, 33.8), P, 34), "74.8 dBi: both main beams"));
verifyTrue(tc, contains(peak_label_rev1(mk(0.9, 26.5, 33.6), P, 34), "ASR-11 main beam + NEXRAD 0.9 deg cut, beam 0.7 deg"));
verifyTrue(tc, contains(peak_label_rev1(mk(0.5, 41, 20), P, 34), "NEXRAD main beam (0.5 deg cut) + ASR-11 beam 0.3 deg"));
verifyTrue(tc, endsWith(peak_label_rev1(mk(0.9, 20, 20), P, 34), "both off axis"));
% A cut missing from the plateau table is never NEXRAD main beam.
verifyTrue(tc, endsWith(peak_label_rev1(mk(1.8, 10, 20), P, 34), "both off axis"));
end

function testReportFiguresWritesEveryOutput(tc)
d = fullfile(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder, "figs");
ov = small();
ov.out_dir = d;                                 % created by the script
ov.n_h = 3;
ov.long_duration_s = 1200;                      % 1 h in all: peaks above the plateau
ov.n_list = [1e4 2e4];
run_script("report_figures_rev5", ov);
png = "fig" + ["1_antenna_patterns", "2_vcp_schedule", "3_low_tilt_driver", "4_exceedance", ...
    "5_max_per_rev", "6_peak_anatomy", "7_mbmb_vs_statistical", "8_sensitivity", "9_asr9_comparison", ...
    "_a1_call_tree"] + ".png";
csv = ["vcp_sweep", "ref_levels", "by_elevation", "asr_rev_pct", "sensitivity", "plateaus", ...
    "asr_rev_pct_120h", "lobe_pairs", "persistence_120h", "persistence_volume_120h", "peaks_above_plateau_120h", ...
    "rev_exceedance_120h", "convergence", "mode_envelope", "sensitivity_sweep", ...
    "alignment_dwell", "cdf_array", "mbmb_vs_statistical", "case1_percentiles", "vcp_table", ...
    "worst_mode_weights", "nexrad_gain_by_cut", "asr11_gain_cut", "call_tree", ...
    "asr9_" + ["vcp_sweep", "mode_envelope", "comparison", "ref_levels", "lobe_pairs", ...
    "persistence", "persistence_volume", "sensitivity", "assumptions"]] + "_rev5.csv";
for f = [png, csv]
    verifyTrue(tc, isfile(fullfile(d, f)), f);
end
S = readtable(fullfile(d, "vcp_sweep_rev5.csv"), "TextType", "string", "VariableNamingRule", "preserve");
verifyEqual(tc, height(S), 75);
C = readtable(fullfile(d, "convergence_rev5.csv"));
verifyEqual(tc, unique(C.n_samples)', ov.n_list);
L = readtable(fullfile(d, "lobe_pairs_rev5.csv"));
verifyEqual(tc, sum(L.time_pct), 100, "AbsTol", 1e-6);
% Typed report tables match their sources: NEXRAD 41.99 dBi at 0.5 deg
% (Step 2), ASR-11 34 dBi on the path (Step 3), worst-mode shares sum to 100.
Nb = readtable(fullfile(d, "nexrad_gain_by_cut_rev5.csv"));
verifyEqual(tc, Nb.nexrad_gain_on_path_dBi(Nb.el_deg == 0.5), 45.5 - 12*(0.5/0.925)^2, "AbsTol", 1e-3);
Ab = readtable(fullfile(d, "asr11_gain_cut_rev5.csv"));
verifyEqual(tc, Ab.asr11_gain_dBi(1), 34, "AbsTol", 1e-9);
Wb = readtable(fullfile(d, "worst_mode_weights_rev5.csv"));
verifyEqual(tc, sum(Wb.share_pct), 100, "AbsTol", 1e-9);
Vt = readtable(fullfile(d, "vcp_table_rev5.csv"));
verifyEqual(tc, height(Vt), 8);
Ct = readtable(fullfile(d, "call_tree_rev5.csv"), "TextType", "string");
verifyTrue(tc, any(Ct.caller == "report_figures_rev5" & Ct.callee == "example_nexrad_asr11_rev1"));
verifyTrue(tc, any(Ct.caller == "nexrad_asr11_mag_rev1" & Ct.callee == "radar_mutual_gain_cdf_rev1"));
A9 = readtable(fullfile(d, "asr9_assumptions_rev5.csv"), "TextType", "string");
verifyEqual(tc, A9.value(A9.parameter == "ASR-9 peak gain"), 33.5);
P = readtable(fullfile(d, "peaks_above_plateau_120h_rev5.csv"));
verifyGreaterThan(tc, height(P), 0, "the labeled-peak code runs");
end

function testReportNumbersRev6(tc)
% Small run of the Rev 6 tables: every file written, five cases, and the
% Rev 6 defaults at a horizontal path (ROC specification far level 5.5 dBi
% plus ASR-11 remote 3.5 dBi = 9 dBi at 50%; NSSL -9.5 + 3.5 = -6 dBi).
d = fullfile(tempname, "r6");
rev6_overrides = struct("out_dir", d, "n_trials", 1, "duration_s", 120); %#ok<NASGU>
evalc("report_numbers_rev6");
f = ["headline", "noaa_tilt", "timeseries", "harm", "lookup_pair", "ref_levels"] + "_rev6.csv";
for k = f
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
H = readtable(fullfile(d, "headline_rev6.csv"), "TextType", "string", "VariableNamingRule", "preserve");
verifyEqual(tc, unique(H.case)', sort(["roc_asr11" "roc_asr9" "nssl_asr11" "nssl_asr9" "noaa_asr11"]));
verifyEqual(tc, H.env216_dBi(H.case == "roc_asr11" & H.pct == 50), 9, "AbsTol", 0.35);
verifyEqual(tc, H.env216_dBi(H.case == "nssl_asr11" & H.pct == 50), -6, "AbsTol", 0.35);
verifyTrue(tc, all(H.env216_dBi >= H.env74_dBi), "216 modes include the 74");
M = readtable(fullfile(d, "harm_rev6.csv"), "TextType", "string", "VariableNamingRule", "preserve");
verifyEqual(tc, unique(M.radial_deg)', [0.5 1]);
rmdir(fileparts(d), "s");
end

function testReportNumbersRev8(tc)
% Small run of the Rev 8 tables and figures: every cited file is written
% by the script (finding Q7), the time-series envelope is at or above the
% CDF-worst mode (Q1), the reference levels see the cut on the path (Q6),
% and the NOAA 50% value is the far level 5.5 dBi plus the +6.5 dBi
% remote sidelobe.
d = fullfile(tempname, "r8");
rev8_overrides = struct("out_dir", d, "n_trials", 1, "duration_s", 900, "n_seeds", 1, ...
    "drift_trials", 1, "drift_hours", 0.25, "grid_step", 0.5, "n_modes", 3, "parallel_flag", 0); %#ok<NASGU>
evalc("report_numbers_rev8");
f = ["takeoff_worst", "worst_modes", "ts_all_modes", "ts_envelope", "vol_exceed_modes", "seeds_runs", ...
    "seeds_summary", "noaa_case", "worked_example", "ref_levels", "bootstrap", "trial_length", "timing"] + "_rev8.csv";
for k = [f, "timing_rev8_machine.txt", "fig12_volume_exceedance.png", "fig13_runs_vs_drift.png", "fig14_harm.png"]
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
E = readtable(fullfile(d, "ts_envelope_rev8.csv"), "TextType", "string");
verifyTrue(tc, all(E.max_over_modes >= E.at_cdf_worst_mode), "envelope below the CDF-worst mode");
R = readtable(fullfile(d, "ref_levels_rev8.csv"), "TextType", "string");
verifyEqual(tc, R.mag_dBi(1), 79.5, "AbsTol", 1e-6);
N = readtable(fullfile(d, "noaa_case_rev8.csv"), "TextType", "string", "VariableNamingRule", "preserve");
verifyEqual(tc, N.p50(N.case == "noaa" & N.row == "at_0_0"), 12, "AbsTol", 0.35);
W = readtable(fullfile(d, "worst_modes_rev8.csv"), "TextType", "string", "VariableNamingRule", "preserve");
verifyTrue(tc, all(W.n_modes_tied >= 1));
rmdir(fileparts(d), "s");
end

function testReportNumbersRev9(tc)
% Small run of the Rev 9 tables and figures (independent recheck S1, S2,
% S3, S10): every file is written, per-hour volumes equal the share times
% volumes per hour, and a run touching a trial end is flagged censored.
d = fullfile(tempname, "r9");
rev9_overrides = struct("out_dir", d, "n_seeds", 1, "drift_trials", 2, "drift_hours", 0.5, ...
    "n_modes_time", 1, "time_hours", 1, "parallel_flag", 0); %#ok<NASGU>
evalc("report_numbers_rev9");
f = ["per_hour_modes", "per_hour_envelope", "run_modes", "seeds_runs", "seeds_summary", "lock_trials", ...
    "lock_fit", "timing_modes"] + "_rev9.csv";
for k = [f, "timing_rev9_machine.txt", "fig13_runs_vs_drift_rev9.png", "fig_lock_vs_slip_rev9.png"]
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
P = readtable(fullfile(d, "per_hour_modes_rev9.csv"), "TextType", "string");
verifyEqual(tc, height(P), 432, "216 modes at two paths");
verifyEqual(tc, P.vols_per_h_60, P.vol_pct_60/100*3600./P.t_volume_s, "RelTol", 1e-12);
k = P.mode == "VCP 212 SAILSx3 AVSET 6.4";
verifyEqual(tc, P.cuts05_per_volume(k), [8; 8]);
S = readtable(fullfile(d, "seeds_runs_rev9.csv"), "TextType", "string");
verifyTrue(tc, all(max(S.max_uncensored_run_vols, S.max_censored_run_vols) == S.max_run_vols));
L = readtable(fullfile(d, "lock_trials_rev9.csv"), "TextType", "string");
verifyTrue(tc, all(abs(L.slip_deg_per_vol) <= 180));
rmdir(fileparts(d), "s");
end

function testReportNumbersRev10(tc)
% Small run of the Rev 10 tables (independent recheck T1 to T8): every
% file is written, ranks are competition ranks, m4 is the smallest of
% |k x slip| over k = 1 to 4, and an episode spans at least its longest run.
d = fullfile(tempname, "r10");
rev10_overrides = struct("out_dir", d, "n_seeds", 1, "drift_trials", 2, "drift_hours", 0.5, ...
    "time_hours", 1, "n_repeats", 1, "modes_time", "VCP 212 SAILSx3 AVSET 6.4", "parallel_flag", 0); %#ok<NASGU>
if ~isempty(gcp("nocreate"))
    % The timing part refuses an open pool (its guard, T8); a batch job
    % keeps its pool open, so check the guard and run the other parts
    % (recheck of revision 13, Y8).
    verifyError(tc, @() run_rev10_in(rev10_overrides), "report_numbers_rev10:poolOpen");
    rev10_overrides.part = "seeds";
    run_rev10_in(rev10_overrides);
    f = ["per_hour_modes", "per_hour_envelope", "per_hour_corr", "seeds_runs", "seeds_summary", "lock_trials", ...
        "lock_m4", "lock_fit"] + "_rev10.csv";
    for k = [f, "fig_lock_vs_m4_rev10.png"]
        verifyTrue(tc, isfile(fullfile(d, k)), k);
    end
    L = readtable(fullfile(d, "lock_trials_rev10.csv"), "TextType", "string");
    m4 = min(abs(mod((1:4).*L.slip_deg_per_vol + 180, 360) - 180), [], 2);
    verifyEqual(tc, L.m4_deg, m4, "AbsTol", 1e-9);
    rmdir(fileparts(d), "s");
    return
end
evalc("report_numbers_rev10");
f = ["per_hour_modes", "per_hour_envelope", "per_hour_corr", "seeds_runs", "seeds_summary", "lock_trials", ...
    "lock_m4", "lock_fit", "timing_modes"] + "_rev10.csv";
for k = [f, "timing_rev10_machine.txt", "fig_lock_vs_m4_rev10.png"]
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
P = readtable(fullfile(d, "per_hour_modes_rev10.csv"), "TextType", "string");
k = P.el_deg == 0.5;
v = P.cuts_hit_60_pct(k);
verifyEqual(tc, P.rank_cuts_hit_60_pct(k), 1 + sum(v' > v, 2), "competition rank");
L = readtable(fullfile(d, "lock_trials_rev10.csv"), "TextType", "string");
m4 = min(abs(mod((1:4).*L.slip_deg_per_vol + 180, 360) - 180), [], 2);
verifyEqual(tc, L.m4_deg, m4, "AbsTol", 1e-9);
verifyTrue(tc, all(L.ep_span_vols >= L.longest_run_vols));
rmdir(fileparts(d), "s");
end

function testReportNumbersRev11(tc)
% Small run of the Rev 11 tables (independent recheck U1 to U3): every file
% is written, the periodic run is at least the run (k = 1 is a run), and
% the episodes grow with the margin.
d = fullfile(tempname, "r11");
rev11_overrides = struct("out_dir", d, "n_seeds", 1, "drift_trials", 2, "drift_hours", 0.5, ...
    "parallel_flag", 0); %#ok<NASGU>
evalc("report_numbers_rev11");
f = ["seeds_runs", "lock_trials", "lock_fit", "lock_top_trials", "lock_m4"] + "_rev11.csv";
for k = [f, "fig_lock_periodic_rev11.png"]
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
L = readtable(fullfile(d, "lock_trials_rev11.csv"), "TextType", "string");
verifyTrue(tc, all(L.per_vols >= L.longest_run_vols));
verifyTrue(tc, all(L.ep3_vols <= L.ep5_vols & L.ep5_vols <= L.ep8_vols));
rmdir(fileparts(d), "s");
end

function testReportNumbersRev12(tc)
% Small run of the Rev 12 tables (independent recheck V2 to V6, V9): every
% file is written, including the chance tables, and the intervals bracket
% their estimates.
d = fullfile(tempname, "r12");
rev12_overrides = struct("out_dir", d, "n_seeds", 1, "drift_trials", 2, "drift_hours", 0.5, ...
    "n_null", 200, "parallel_flag", 0); %#ok<NASGU>
evalc("report_numbers_rev12");
f = ["seeds_runs", "lock_trials", "lock_k16", "lock_fit", "lock_walk_null"] + "_rev12.csv";
for k = [f, "fig_lock_k16_rev12.png"]
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
F = readtable(fullfile(d, "lock_fit_rev12.csv"), "TextType", "string");
verifyEqual(tc, unique(F.n_null_trials), 200);
% Recheck W14: the interval method against a hand value. Clopper-Pearson
% for 5 of 200: beta quantiles 0.025 of Beta(5, 196) and 0.975 of
% Beta(6, 195), worked outside MATLAB (scipy.stats.beta).
[~, ci] = binofit(5, 200);
verifyEqual(tc, ci, [0.0081662 0.0573744], "AbsTol", 1e-6);
rmdir(fileparts(d), "s");
end

function testReportNumbersRev13(tc)
% Small run of the Rev 13 tables (independent recheck W1 to W4): every file
% is written, the run span holds at least the repeating high volumes, and
% periods up to 16 find runs at least as long as periods up to 4.
d = fullfile(tempname, "r13");
rev13_overrides = struct("out_dir", d, "n_seeds", 1, "drift_trials", 2, "drift_hours", 0.5, ...
    "n_null", 100, "parallel_flag", 0); %#ok<NASGU>
evalc("report_numbers_rev13");
f = ["seeds_runs", "lock_trials", "lock_k16", "lock_long_k", "lock_fit", "lock_walk_null"] + "_rev13.csv";
for k = [f, "fig_lock_any16_rev13.png"]
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
L = readtable(fullfile(d, "lock_trials_rev13.csv"), "TextType", "string");
verifyTrue(tc, all(L.per_span_high >= L.per_hits));
verifyTrue(tc, all(L.per16_vols >= L.per_vols));
verifyTrue(tc, all(L.lock_any16_drift_deg <= L.lock_run_drift_deg + 1e-9 | isnan(L.lock_run_drift_deg)));
rmdir(fileparts(d), "s");
end

function testReportNumbersRev14(tc)
% Small run of the Rev 14 tables (independent recheck Y1 to Y10) on the
% committed Rev 13 per-trial file: the between share matches a hand count,
% the threshold table is cumulative, the bootstrap brackets the grid
% chance, and the walk table is in percent.
src = fullfile(fileparts(which("report_numbers_rev14")), "docs", "figures");
d = fullfile(tempname, "r14");
run_rev14_in(struct("out_dir", d, "in_dir", src, "n_grid_sims", 200, "p_step", 0.05, ...
    "n_boot", 200, "parallel_flag", 0, "cells", [0.5 60]));
f = ["lock_offlattice", "lock_share_drift", "lock_null_grid", "lock_null_boot", "lock_long_k", "lock_walk_null"] + "_rev14.csv";
for k = [f, "fig_lock_any16_rev14.png"]
    verifyTrue(tc, isfile(fullfile(d, k)), k);
end
L = readtable(fullfile(src, "lock_trials_rev13.csv"), "TextType", "string");
L = L(L.mode == "VCP 212 SAILSx3 AVSET 6.4" & L.el_deg == 0.5 & L.threshold_dBi == 60, :);
A = readtable(fullfile(src, "ts_all_modes_rev8.csv"), "TextType", "string");
tv = A.t_volume_s(find(A.mode == "VCP 212 SAILSx3 AVSET 6.4", 1));
lg = L.per_vols*tv/3600 >= 4;
bs = (L.per_span_high(lg) - L.per_hits(lg))./(L.per_vols(lg) - L.per_hits(lg));
O = readtable(fullfile(d, "lock_offlattice_rev14.csv"), "TextType", "string");
o = O(O.mode == "VCP 212 SAILSx3 AVSET 6.4" & O.el_deg == 0.5 & O.threshold_dBi == 60 & O.hours == 4, :);
verifyEqual(tc, o.n_long, nnz(lg));
verifyEqual(tc, o.median_between_share_long, median(bs, "omitnan"), "AbsTol", 1e-12);
S = readtable(fullfile(d, "lock_share_drift_rev14.csv"), "TextType", "string");
s = S(S.mode == "VCP 212 SAILSx3 AVSET 6.4" & S.hours == 4, :);
verifyTrue(tc, all(diff(s.n_tight) >= 0) && all(diff(s.n_tight_long) >= 0), "cumulative in the threshold");
verifyEqual(tc, s.n_tight(end), height(L), "every slip has a lock within 360/17 deg");
B = readtable(fullfile(d, "lock_null_boot_rev14.csv"), "TextType", "string");
verifyTrue(tc, all(B.null_boot_lo_pct <= B.null_grid_pct + 1e-12 & B.null_grid_pct <= B.null_boot_hi_pct + 1e-12));
verifyEqual(tc, B.ratio_lo, B.model_lo_pct./B.null_boot_hi_pct, "RelTol", 1e-12);
W0 = readtable(fullfile(src, "lock_walk_null_rev13.csv"), "TextType", "string");
W = readtable(fullfile(d, "lock_walk_null_rev14.csv"), "TextType", "string");
verifyEqual(tc, W.null_trial_ge_pct, 100*W0.null_trial_ge_pct, "AbsTol", 1e-12);
verifyEqual(tc, W.chance_any_of_trials_pct, 100*(1 - (1 - W0.null_trial_ge_pct).^W0.n_trials), "AbsTol", 1e-9);
rmdir(fileparts(d), "s");
end

function run_rev14_in(rev14_overrides) %#ok<INUSD>
% Run report_numbers_rev14 in its own workspace, so its variables don't
% overwrite the test's.
evalc("report_numbers_rev14");
end

function run_rev10_in(rev10_overrides) %#ok<INUSD>
% Run report_numbers_rev10 in its own workspace with the given overrides.
evalc("report_numbers_rev10");
end

% =========================================================================
function v = W_vol()
v = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4).volume_s;
end

function o = small()
o = struct("n_samples", 2e4, "duration_s", 60, "n_trials", 1);
end

function W = run_script(name, ov)
% Run a script in this function's workspace and return its variables.
% Scripts clear their workspace except example_overrides.
example_overrides = ov; %#ok<NASGU>
evalc(name);
v = who;
W = struct();
for i = 1:numel(v)
    W.(v{i}) = eval(v{i});
end
end
