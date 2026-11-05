function tests = test_scripts_rev1
%TEST_SCRIPTS_REV1 Tests for the example scripts, the C# crosscheck, the
%   report figure script, merge_overrides_rev1 and peak_label_rev1. Each script runs end to
%   end with small runs set through a workspace struct example_overrides,
%   so every line of each script executes. report_figures_rev4 writes to a
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
run_script("report_figures_rev4", ov);
png = "fig" + ["1_antenna_patterns", "2_vcp_schedule", "3_low_tilt_driver", "4_exceedance", ...
    "5_max_per_rev", "6_peak_anatomy", "7_mbmb_vs_statistical", "8_sensitivity", "9_asr9_comparison", ...
    "_a1_call_tree"] + ".png";
csv = ["vcp_sweep", "ref_levels", "by_elevation", "asr_rev_pct", "sensitivity", "plateaus", ...
    "asr_rev_pct_120h", "lobe_pairs", "persistence_120h", "persistence_volume_120h", "peaks_above_plateau_120h", ...
    "rev_exceedance_120h", "convergence", "mode_envelope", "sensitivity_sweep", ...
    "alignment_dwell", "cdf_array", "mbmb_vs_statistical", "case1_percentiles", "vcp_table", ...
    "worst_mode_weights", "nexrad_gain_by_cut", "asr11_gain_cut", "call_tree", ...
    "asr9_" + ["vcp_sweep", "mode_envelope", "comparison", "ref_levels", "lobe_pairs", ...
    "persistence", "persistence_volume", "sensitivity", "assumptions"]] + "_rev4.csv";
for f = [png, csv]
    verifyTrue(tc, isfile(fullfile(d, f)), f);
end
S = readtable(fullfile(d, "vcp_sweep_rev4.csv"), "TextType", "string", "VariableNamingRule", "preserve");
verifyEqual(tc, height(S), 75);
C = readtable(fullfile(d, "convergence_rev4.csv"));
verifyEqual(tc, unique(C.n_samples)', ov.n_list);
L = readtable(fullfile(d, "lobe_pairs_rev4.csv"));
verifyEqual(tc, sum(L.time_pct), 100, "AbsTol", 1e-6);
% Typed report tables match their sources: NEXRAD 41.99 dBi at 0.5 deg
% (Step 2), ASR-11 34 dBi on the path (Step 3), worst-mode shares sum to 100.
Nb = readtable(fullfile(d, "nexrad_gain_by_cut_rev4.csv"));
verifyEqual(tc, Nb.nexrad_gain_on_path_dBi(Nb.el_deg == 0.5), 45.5 - 12*(0.5/0.925)^2, "AbsTol", 1e-3);
Ab = readtable(fullfile(d, "asr11_gain_cut_rev4.csv"));
verifyEqual(tc, Ab.asr11_gain_dBi(1), 34, "AbsTol", 1e-9);
Wb = readtable(fullfile(d, "worst_mode_weights_rev4.csv"));
verifyEqual(tc, sum(Wb.share_pct), 100, "AbsTol", 1e-9);
Vt = readtable(fullfile(d, "vcp_table_rev4.csv"));
verifyEqual(tc, height(Vt), 8);
Ct = readtable(fullfile(d, "call_tree_rev4.csv"), "TextType", "string");
verifyTrue(tc, any(Ct.caller == "report_figures_rev4" & Ct.callee == "example_nexrad_asr11_rev1"));
verifyTrue(tc, any(Ct.caller == "nexrad_asr11_mag_rev1" & Ct.callee == "radar_mutual_gain_cdf_rev1"));
A9 = readtable(fullfile(d, "asr9_assumptions_rev4.csv"), "TextType", "string");
verifyEqual(tc, A9.value(A9.parameter == "ASR-9 peak gain"), 33.5);
P = readtable(fullfile(d, "peaks_above_plateau_120h_rev4.csv"));
verifyGreaterThan(tc, height(P), 0, "the labeled-peak code runs");
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
