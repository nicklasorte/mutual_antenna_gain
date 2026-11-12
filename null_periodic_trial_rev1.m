function [p4, p16, j] = null_periodic_trial_rev1(pv, nv, i, base, seed)
%NULL_PERIODIC_TRIAL_REV1 Longest periodic run of one null trial of independent volumes.
%   [p4, p16, j] = null_periodic_trial_rev1(pv, nv, i, base, seed) draws
%   null trial i: nv(j) independent WSR-88D volumes, each high with
%   probability pv(j), where j = 1 + mod(i - 1, numel(pv)) cycles through
%   the model trials. It returns the longest periodic run (periodic_runs_rev1)
%   with periods up to 4 (p4) and up to 16 (p16), in volumes, and j.
%   With 20,000 null trials and 200 model trials, each model trial's share
%   and length is used 100 times (independent recheck of revision 13,
%   findings W2 and Y11).
%
%   pv    shares of high volumes, one per model trial (0 to 1)
%   nv    numbers of complete volumes, one per model trial
%   i     null trial index (positive integer)
%   base  substream offset, so each cell of a study uses its own streams
%   seed  RandStream seed (default 20261113, as in report_numbers_rev13)
%
%   Draws use RandStream("mrg32k3a") with substream base + i, so a parfor
%   over i gives the same results as a for loop.
%
%   Example:
%     [p4, ~, j] = null_periodic_trial_rev1([0 1], [10 12], 2, 0);  % 12, 2: all high
%
%   See also periodic_runs_rev1, report_numbers_rev14.

arguments
    pv (:,1) double {mustBeInRange(pv, 0, 1)}
    nv (:,1) double {mustBeInteger, mustBePositive}
    i (1,1) double {mustBeInteger, mustBePositive}
    base (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    seed (1,1) double {mustBeInteger, mustBeNonnegative} = 20261113
end
if numel(pv) ~= numel(nv)
    error("null_periodic_trial_rev1:sizeMismatch", "pv and nv must have the same number of elements.");
end
j = 1 + mod(i - 1, numel(pv));
rs = RandStream("mrg32k3a", "Seed", seed);
rs.Substream = base + i;
n_v = nv(j);
g = 50 + 11*(rand(rs, n_v, 1) < pv(j));          % 61 dBi when high, 50 dBi otherwise
V = table(ones(n_v, 1), true(n_v, 1), g, 'VariableNames', ["trial" "complete" "max_g21_dBi"]);
a = periodic_runs_rev1(V, 60, 4);
p4 = a(1);
p16 = NaN;
if nargout > 1
    b = periodic_runs_rev1(V, 60, 16);
    p16 = b(1);
end
end
