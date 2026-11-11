function [s, tr] = periodic_runs_rev1(Vt, thr, k_max)
%PERIODIC_RUNS_REV1 Periodic runs of high WSR-88D volumes (locks over k volumes).
%   [s, tr] = periodic_runs_rev1(Vt, thr, k_max) finds, in each trial, the
%   longest periodic run: volumes v, v + k, v + 2k, ... that all reach thr,
%   for a period k of 1 to k_max volumes (default 4). A lock over k volumes
%   brings the ATC radar back to nearly the same azimuth every k volumes, so
%   it can repeat high coupling every k volumes for hours while the volumes
%   between stay low; runs and episodes (run_stats_rev1), which need
%   consecutive volumes, count that as length 1 (independent recheck of
%   revision 10, finding U1). k = 1 gives the ordinary run.
%
%   The duration of a periodic run of h high volumes at period k is
%   (h - 1) k + 1 volumes, first high volume to last. A periodic run is
%   censored when the trial cut it: its first high volume is within k of
%   the first complete volume, or its last within k of the last.
%
%   Vt   table with columns trial, complete, max_g21_dBi (one row per
%        volume, in time order within each trial), as in
%        mag_pair_timeseries_rev1 output .nexrad_volume
%
%   s    row: [longest duration (volumes), its period k, high volumes in
%        it, longest uncensored duration, longest censored duration]
%   tr   struct, one entry per trial: dur, k, hits, censored (of the
%        longest periodic run; on a tie, the shortest period, and censored
%        when any tied run is censored)
%
%   Example:
%     g = 50*ones(1, 15); g([5 8 11]) = 61;        % every third volume
%     V = table(ones(15, 1), true(15, 1), g(:), ...
%         'VariableNames', ["trial" "complete" "max_g21_dBi"]);
%     s = periodic_runs_rev1(V, 60);   % [7 3 3 7 0]: 3 hits over 7 volumes, k = 3
%
%   See also run_stats_rev1, report_numbers_rev11.

arguments
    Vt table
    thr (1,1) double
    k_max (1,1) double {mustBeInteger, mustBePositive} = 4
end
V = Vt(Vt.complete, :);
tu = unique(V.trial)';
n_t = numel(tu);
z = zeros(n_t, 1);
tr = struct("dur", z, "k", z, "hits", z, "censored", false(n_t, 1));
all_dur = [];
all_cen = false(0, 1);
for jj = 1:n_t
    hit = V.max_g21_dBi(V.trial == tu(jj)) >= thr;
    n = numel(hit);
    best = [0 0 0 0];                               % [duration, k, hits, censored]
    for k = 1:k_max
        for r = 1:k
            idx = r:k:n;
            h = hit(idx);
            d = diff([0; h(:); 0]);
            b = find(d == 1);
            e = find(d == -1) - 1;
            if isempty(b)
                continue
            end
            nh = e - b + 1;
            dur = (nh - 1)*k + 1;
            first = idx(b)';
            last = idx(e)';
            cen = first <= k | last > n - k;        % the lattice point before or after lies outside the trial
            all_dur = [all_dur; dur(:)]; %#ok<AGROW>
            all_cen = [all_cen; cen(:)]; %#ok<AGROW>
            [dm, im] = max(dur);
            cm = any(cen(dur == dm));
            if dm > best(1)
                best = [dm k nh(im) cm];
            elseif dm == best(1) && dm > 0
                best(4) = best(4) || cm;            % a tie keeps the shortest period, censored if any tied run is (recheck V7)
            end
        end
    end
    tr.dur(jj) = best(1);
    tr.k(jj) = best(2);
    tr.hits(jj) = best(3);
    tr.censored(jj) = best(4) == 1;
end
if isempty(all_dur)
    s = [0 0 0 0 0];
    return
end
[~, jm] = max(tr.dur);
mx = @(v) max([0; v(:)]);
s = [tr.dur(jm), tr.k(jm), tr.hits(jm), mx(all_dur(~all_cen)), mx(all_dur(all_cen))];
end
