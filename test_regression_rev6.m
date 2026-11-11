function tests = test_regression_rev6
%TEST_REGRESSION_REV6 Pinned revision 6 outputs (guard for the speed work).
%   The pair CDFs (envelope percentiles, worst modes, largest MAG on the
%   path) and the time series (every per-scan maximum and harm count of
%   four 1 h trials) must equal test_data/rev6_*_reference_rev1.csv,
%   written by make_rev6_reference_rev1 at commit a33e806 before any speed
%   change (docs/rev6_redteam_plan_rev1.md, V1 to V8).
%
%   Run: results = runtests("test_regression_rev6");
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(which("test_regression_rev6"));
tc.TestData.pair = readtable(fullfile(here, "test_data", "rev6_pair_reference_rev1.csv"), "TextType", "string");
tc.TestData.ts = readtable(fullfile(here, "test_data", "rev6_ts_reference_rev1.csv"), "TextType", "string");
[tc.TestData.new_pair, tc.TestData.new_ts] = rev6_reference_rev1();
end

function testPairCdfUnchanged(tc)
R = tc.TestData.pair;
N = tc.TestData.new_pair;
verifyEqual(tc, N.mag_dBi, R.mag_dBi, "AbsTol", 1e-9);
verifyEqual(tc, N.worst_mode, R.worst_mode);
verifyEqual(tc, N.max_on_path_dBi, R.max_on_path_dBi, "AbsTol", 1e-9);
end

function testTimeSeriesUnchanged(tc)
R = tc.TestData.ts;
N = tc.TestData.new_ts;
verifyEqual(tc, height(N), height(R));
verifyEqual(tc, N.kind, R.kind);
verifyEqual(tc, [N.trial N.scan], [R.trial R.scan]);
verifyEqual(tc, N.complete, logical(R.complete));
verifyEqual(tc, N.max_g21_dBi, R.max_g21_dBi, "AbsTol", 1e-9);
verifyEqual(tc, N.max_g12_dBi, R.max_g12_dBi, "AbsTol", 1e-9);
end
