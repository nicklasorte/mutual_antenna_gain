function R = redteam_fault_seed_rev1(repo_dir, opts)
%REDTEAM_FAULT_SEED_REV1 Check that the tests catch deliberate bugs (red team F2).
%   R = redteam_fault_seed_rev1(repo_dir) copies the repository to a
%   temporary folder once per fault, plants one known bug by an exact text
%   replacement, runs the test files, and records whether any test fails.
%   A clean copy runs first as a control and must pass, or the function
%   stops. A fault that no test catches is a gap in the tests, not in the
%   model: add a test that would have caught it.
%
%   Plan: docs/mag_optimization_redteam_plan_rev1.md, attack F2.
%
%   Tests whose names contain "parallel" are left out (a batch worker
%   cannot start a pool); run them in the client session.
%
%   Run it in a fresh MATLAB process (batch) whose path does not hold
%   repo_dir, so the temporary copy is the one that runs:
%     j = batch(parcluster("Processes"), @redteam_fault_seed_rev1, 1, {pwd}, ...
%         "CurrentFolder", pwd, "AutoAddClientPath", false);
%
%   repo_dir  repository folder
%   opts      struct, all fields optional:
%     .tests    test files to run (default: both models' unit tests and
%               test_regression_rev5; test_scripts_rev1 and
%               test_validation_rev1 are left out for run time)
%     .out_csv  results file (default docs/figures/fault_seed_rev1.csv)
%
%   R table: id, file, fault, caught, failed_tests.
%
%   See also test_gain_slice_rev1, test_mag_lookup_rev1, test_regression_rev5.

arguments
    repo_dir (1, 1) string
    opts (1, 1) struct = struct()
end
if ~isfield(opts, "tests")
    opts.tests = ["test_gain_slice_rev1", "test_mag_lookup_rev1", "test_regression_rev5", ...
        "test_radar_mutual_gain_cdf_rev1", "test_nexrad_asr11_rev1"];
end
if ~isfield(opts, "out_csv")
    opts.out_csv = fullfile(repo_dir, "docs", "figures", "fault_seed_rev1.csv");
end
F = faults();
here = pwd;
cleanup = onCleanup(@() cd(here));

% ---- Control: the unmodified copy must pass.
[ok, failed] = run_copy(repo_dir, "", "", "", opts.tests);
if ~ok
    error("fault_seed:controlFailed", "The unmodified copy fails: %s. Fix the tests first.", strjoin(failed, ", "));
end

n = height(F);
caught = false(n, 1);
failed_tests = strings(n, 1);
for i = 1:n
    [ok, failed] = run_copy(repo_dir, F.file(i), F.old(i), F.new(i), opts.tests);
    caught(i) = ~ok;
    failed_tests(i) = strjoin(failed, "; ");
    fprintf("%-3d %-6s %s\n", F.id(i), string(caught(i)), F.fault(i));
end
R = table(F.id, F.file, F.fault, caught, failed_tests, ...
    'VariableNames', ["id", "file", "fault", "caught", "failed_tests"]);
writetable(R, opts.out_csv);
fprintf("%d of %d faults caught. Wrote %s\n", sum(caught), n, opts.out_csv);
end

% =========================================================================
function [ok, failed] = run_copy(repo_dir, file, old, new, tests)
tmp = tempname;
mkdir(tmp);
del = onCleanup(@() rmdir(tmp, "s"));
copyfile(fullfile(repo_dir, "*.m"), tmp);
mkdir(fullfile(tmp, "test_data"));
copyfile(fullfile(repo_dir, "test_data", "*"), fullfile(tmp, "test_data"));
mkdir(fullfile(tmp, "docs", "figures"));
copyfile(fullfile(repo_dir, "docs", "figures", "vcp_sweep_rev5.csv"), fullfile(tmp, "docs", "figures"));
if file ~= ""
    f = fullfile(tmp, file);
    txt = string(fileread(f));
    k = count(txt, old);
    if k ~= 1
        error("fault_seed:badFault", "Fault text must occur once in %s (found %d): %s", file, k, old);
    end
    writelines(replace(txt, old, new), f);
end
cd(tmp);
clear functions; %#ok<CLFUNC> new copy of every function
rehash;
% Parallel tests need a pool, which a batch worker cannot start; they are
% checked in the client session instead. Warnings stay on: some tests
% expect them.
suite = testsuite(cellstr(tests));
keep = ~contains(lower(string({suite.Name})), "parallel");
r = run(suite(keep));
cd(fileparts(tmp));
bad = [r.Failed];                                % filtered tests are not failures
ok = ~any(bad);
failed = string({r(bad).Name});
end

% -------------------------------------------------------------------------
function F = faults()
% Each fault is one exact text replacement that a real bug could make.
c = {
1, "gain_slice_rev1.m", "Closed form: wrong sign of the elevation term in the azimuth share", ...
    "phi = acosd(min(max((cosd(gam) - A)/B, -1), 1));", "phi = acosd(min(max((cosd(gam) + A)/B, -1), 1));"
2, "gain_slice_rev1.m", "Closed form: gains rounded down instead of up", ...
    "k = ceil(gu(keep)/dg - 1e-9);", "k = floor(gu(keep)/dg - 1e-9);"
3, "gain_slice_rev1.m", "Separable: elevation term with the wrong sign", ...
    "s = g0 - S.g_max;", "s = S.g_max - g0;"
4, "gain_slice_rev1.m", "Sector: nearest offset to the path always 0 (sector always sees the main beam)", ...
    "d = min(abs([lo, hi, lo - 360, hi - 360]));", "d = 0;"
5, "mag_pair_cdf_rev1.m", "Pair: receiver evaluated at the transmitter's takeoff angle", ...
    "[pr, k0r, gr] = slices_at(sd_r, er, opts.interp);", "[pr, k0r, gr] = slices_at(sd_r, et, opts.interp);"
6, "mag_pair_cdf_rev1.m", "Pair: interpolation weights swapped (a and 1 - a)", ...
    "q = (1 - a)*(count_less(F0, u) + 1) + a*(count_less(F1, u) + 1);", ...
    "q = a*(count_less(F0, u) + 1) + (1 - a)*(count_less(F1, u) + 1);"
7, "mag_pair_cdf_rev1.m", "Pair: envelope takes the smallest mode instead of the largest", ...
    "[q_env, worst] = max(q_mode, [], 2);", "[q_env, worst] = min(q_mode, [], 2);"
8, "gain_dist_table_rev1.m", "Table: mode time weights not normalized", ...
    "W(m, :) = W(m, :)/sum(W(m, :));", "W(m, :) = W(m, :)/max(W(m, :));"
9, "antenna_gain_rev1.m", "Pattern: wrong sign in the off-axis angle (spherical law of cosines)", ...
    "cos_g = cosd(el).*cosd(tilt).*cosd(az_off) + sind(el).*sind(tilt);", ...
    "cos_g = cosd(el).*cosd(tilt).*cosd(az_off) - sind(el).*sind(tilt);"
10, "radar_mutual_gain_cdf_rev1.m", "Original model: time series starts every trial at the same point of the volume", ...
    "t0 = t_cycle*rand(rs);", "t0 = 0*rand(rs);"
};
F = cell2table(c, 'VariableNames', ["id", "file", "fault", "old", "new"]);
F.file = string(F.file);
F.fault = string(F.fault);
F.old = string(F.old);
F.new = string(F.new);
end
