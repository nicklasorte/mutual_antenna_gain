function [m4, k] = lock_drift_m4_rev1(slip_deg, k_max)
%LOCK_DRIFT_M4_REV1 Drift from whole ATC revolutions over 1 to k_max volumes.
%   [m4, k] = lock_drift_m4_rev1(slip_deg, k_max) returns, for each ATC
%   phase slip per WSR-88D volume (deg, any value), the smallest drift of
%   the ATC azimuth from a whole number of revolutions over k = 1 to k_max
%   volumes (default 4), and the k that gives it. A small m4 means the ATC
%   comes back to nearly the same azimuth, relative to the WSR-88D scan,
%   every k volumes. On a tie, the smallest k.
%
%   Example:
%     [m, k] = lock_drift_m4_rev1([-119 91 3 179.5])   % m = [3 4 3 1], k = [3 4 1 2]
%
%   See also periodic_runs_rev1, report_numbers_rev11.

arguments
    slip_deg double
    k_max (1,1) double {mustBeInteger, mustBePositive} = 4
end
d = abs(mod(slip_deg(:).*(1:k_max) + 180, 360) - 180);   % |drift| over k volumes, deg
[m4, k] = min(d, [], 2);
m4 = reshape(m4, size(slip_deg));
k = reshape(k, size(slip_deg));
end
