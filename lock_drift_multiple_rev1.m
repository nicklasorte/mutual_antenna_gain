function [drift, k] = lock_drift_multiple_rev1(slip_deg, period, k_max)
%LOCK_DRIFT_MULTIPLE_REV1 Lock behind a periodic run: drift over multiples of its period.
%   [drift, k] = lock_drift_multiple_rev1(slip_deg, period, k_max) returns,
%   for each ATC phase slip per WSR-88D volume (deg) and the period of a
%   periodic run (volumes), the smallest drift of the ATC azimuth from a
%   whole number of revolutions over k = period, 2 period, ... up to k_max
%   volumes (default 16), and the k that gives it. A periodic run every p
%   volumes can come from a lock over a multiple of p: over k volumes the
%   ATC comes back to nearly the same azimuth, and the volumes between
%   cycle through the other alignments that also reach the level
%   (independent recheck of revision 11, finding V2). A period of 0 or NaN
%   gives NaN. On a tie, the smallest k.
%
%   Example:
%     [d, k] = lock_drift_multiple_rev1(-159.997, 3)   % 9 x -159.997 = -1439.973: d = 0.027, k = 9
%
%   See also lock_drift_m4_rev1, periodic_runs_rev1.

arguments
    slip_deg double
    period double
    k_max (1,1) double {mustBeInteger, mustBePositive} = 16
end
drift = NaN(size(slip_deg));
k = NaN(size(slip_deg));
p = period .* ones(size(slip_deg));
for i = 1:numel(slip_deg)
    if ~(p(i) >= 1) || isnan(slip_deg(i))
        continue
    end
    kk = p(i):p(i):k_max;
    if isempty(kk)
        continue
    end
    d = abs(mod(slip_deg(i)*kk + 180, 360) - 180);
    [drift(i), j] = min(d);
    k(i) = kk(j);
end
end
