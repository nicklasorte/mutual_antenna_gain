function tests = test_scripts_rev1
%TEST_SCRIPTS_REV1 Tests for the example scripts, the C# crosscheck, the
%   report figure script and merge_overrides_rev1. Each script runs end to
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
% Results are stored lean: g21 in single, g12 dropped.
verifyClass(tc, W.res{2}.cdf.g21_dBi, "single");
verifyEmpty(tc, W.res{2}.cdf.g12_dBi);
end

function testExampleSensitivity(tc)
W = run_script("example_sensitivity_rev1", small());
T = W.Tsens;
verifyEqual(tc, height(T), 19);
verifyEqual(tc, T.group(1), "Baseline");
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
verifyEqual(tc, height(T), 3);
verifyEqual(tc, T.ratio, ["1:1"; "2:1"; "6:1"]);
verifyEqual(tc, T.nominal_episode_s(1:2), [Inf; Inf]);
verifyEqual(tc, T.nominal_episode_s(3), 360/12.524, "RelTol", 1e-12, "capped at one NEXRAD cut");
verifyEqual(tc, T{3, ["p_gt_1min", "p_gt_10min", "p_gt_1h", "p_gt_1d", "p_gt_30d"]}, zeros(1, 5));
verifyTrue(tc, issorted(flip(T{1, 7:end})), "longer episodes are less likely");
end

function testReportFiguresWritesEveryOutput(tc)
d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
ov = small();
ov.out_dir = d;
ov.n_h = 2;
ov.n_list = [1e4 2e4];
run_script("report_figures_rev4", ov);
png = "fig" + ["1_antenna_patterns", "2_vcp_schedule", "3_low_tilt_driver", "4_exceedance", ...
    "5_max_per_rev", "6_peak_anatomy", "7_mbmb_vs_statistical", "8_sensitivity"] + ".png";
csv = ["vcp_sweep", "ref_levels", "by_elevation", "asr_rev_pct", "sensitivity", "plateaus", ...
    "asr_rev_pct_120h", "lobe_pairs", "persistence_120h", "peaks_above_plateau_120h", ...
    "rev_exceedance_120h", "convergence", "mode_envelope", "sensitivity_sweep", ...
    "alignment_dwell", "cdf_array", "mbmb_vs_statistical"] + "_rev4.csv";
for f = [png, csv]
    verifyTrue(tc, isfile(fullfile(d, f)), f);
end
S = readtable(fullfile(d, "vcp_sweep_rev4.csv"), "TextType", "string");
verifyEqual(tc, height(S), 75);
C = readtable(fullfile(d, "convergence_rev4.csv"));
verifyEqual(tc, unique(C.n_samples)', ov.n_list);
L = readtable(fullfile(d, "lobe_pairs_rev4.csv"));
verifyEqual(tc, sum(L.time_pct), 100, "AbsTol", 1e-6);
end

% =========================================================================
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
