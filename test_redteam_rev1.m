function tests = test_redteam_rev1
%TEST_REDTEAM_REV1 Code-path tests for redteam_runs_rev1 (quick runs only).
%   The quick runs check that every attack runs, writes its evidence and
%   its register row, and that the register assembles. They are not
%   evidence: the full runs are in docs/figures/redteam.
%
%   Run: results = runtests("test_redteam_rev1");
tests = functiontests(localfunctions);
end

function setupOnce(tc)
d = fullfile(tempname, "repo");
mkdir(fullfile(d, "docs"));
tc.TestData.repo = d;
tc.TestData.opts = struct("repo_dir", d, "quick", true);
end

function teardownOnce(tc)
rmdir(fileparts(tc.TestData.repo), "s");
end

function testEveryAttackRunsAndWritesRow(tc)
o = tc.TestData.opts;
R = redteam_runs_rev1("all", o);
ids = ["A1" "A2" "A3" "A4" "B1" "B2" "B4" "B5" "C1" "C2" "C3" "C4" "C5" "D2" "F1" "NOAA"];
verifyEqual(tc, R.id', ids);
verifyTrue(tc, all(R.quick));
out = fullfile(o.repo_dir, "docs", "figures", "redteam");
for id = [ids "B3"]
    verifyTrue(tc, isfile(fullfile(out, id + "_row.csv")), id + " row");
end
for id = ids
    verifyTrue(tc, isfile(fullfile(out, id + "_rev1.csv")), id + " evidence");
end
verifyTrue(tc, all(strlength(R.key_result) > 0));
end

function testRegisterAssembles(tc)
o = tc.TestData.opts;
redteam_runs_rev1(["C3" "C5"], o);
T = verifyWarning(tc, @() redteam_runs_rev1("register", o), "redteam:quickRows");
verifyTrue(tc, all(ismember(["C3" "C5" "D1" "E1" "F2" "F3" "L2" "G1"], T.id)));
verifyTrue(tc, isfile(fullfile(o.repo_dir, "docs", "redteam_register_rev1.csv")));
verifyFalse(tc, any(ismissing(T.result)));
end

function testUnknownAttack(tc)
verifyError(tc, @() redteam_runs_rev1("Z9", tc.TestData.opts), "redteam:unknownAttack");
verifyError(tc, @() redteam_runs_rev1("C5", struct("bogus", 1)), "redteam:unknownOption");
end

function testRev6InputsQuickRun(tc)
% The revision 6 inputs (ROC specification, 216 modes) run the attacks into
% their own folder and refuse to touch the revision 5 register.
d = fullfile(tc.TestData.repo, "rev6");
o = struct("repo_dir", tc.TestData.repo, "out_dir", d, "quick", true, "inputs", "rev6");
R = redteam_runs_rev1(["A3" "C4" "NOAA"], o);
verifyEqual(tc, R.id', ["A3" "C4" "NOAA"]);
verifyTrue(tc, isfile(fullfile(d, "A3_rev1.csv")));
verifyError(tc, @() redteam_runs_rev1("register", o), "redteam:registerRev5");
verifyError(tc, @() redteam_runs_rev1("A3", struct("inputs", "rev9", "quick", true)), "redteam:badInputs");
end
