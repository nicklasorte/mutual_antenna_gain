function [s, tr] = run_stats_rev1(Vt, thr, margin_dB)
%RUN_STATS_REV1 Runs and lock episodes of consecutive high WSR-88D volumes.
%   [s, tr] = run_stats_rev1(Vt, thr, margin_dB) counts, in each trial,
%   runs of consecutive complete volumes whose maximum MAG reaches thr, and
%   lock episodes: stretches of consecutive volumes that all reach
%   thr - margin_dB (default 5 dB) and hold at least one volume at thr. An
%   episode joins runs at thr split by short dips that stay near the level,
%   as in a near-lock, without joining unrelated runs across volumes far
%   below it (which a gap-length rule does when 30% to 50% of volumes are
%   high). A run or episode that touches the first or last complete volume
%   of its trial is censored: the trial cut it, so its length is a lower
%   bound.
%
%   Vt   table with columns trial, complete, max_g21_dBi (one row per
%        volume, in time order within each trial), as returned in
%        mag_pair_timeseries_rev1 output .nexrad_volume
%
%   s    row: [median run, 95th-percentile run, n runs, longest run,
%        share of high volumes in runs of 2 or more (%), in runs of 3 or
%        more (%), longest uncensored run, n censored runs, longest
%        censored run, longest episode (volumes), longest uncensored
%        episode, longest censored episode]
%   tr   struct, one entry per trial: longest (longest run), censored
%        (true when any run of that length is censored), ep_span (longest
%        episode span), ep_high (high volumes in that episode), ep_censored
%
%   Independent recheck of revision 9: T1 (one near-lock split into
%   several runs by short dips) and T14 (censored flag on ties).
%
%   Example:
%     V = table([1;1;1;1;1;1], true(6,1), [50;61;57;61;61;61], ...
%         'VariableNames', ["trial" "complete" "max_g21_dBi"]);
%     s = run_stats_rev1(V, 60, 5);   % runs 1 and 3; one episode of 5, censored
%
%   See also mag_pair_timeseries_rev1, report_numbers_rev10.

arguments
    Vt table
    thr (1,1) double
    margin_dB (1,1) double {mustBeNonnegative} = 5
end
V = Vt(Vt.complete, :);
tu = unique(V.trial)';
n_t = numel(tu);
len = [];
cen = false(0, 1);
z = zeros(n_t, 1);
tr = struct("longest", z, "censored", false(n_t, 1), "ep_span", z, "ep_high", z, "ep_censored", false(n_t, 1));
ep_len = [];
ep_cen = false(0, 1);
for jj = 1:n_t
    g_v = V.max_g21_dBi(V.trial == tu(jj));
    hit = g_v >= thr;
    n = numel(hit);
    d = diff([0; hit(:); 0]);
    b = find(d == 1);
    e = find(d == -1) - 1;
    if isempty(b)
        continue
    end
    l = e - b + 1;
    c = b == 1 | e == n;
    len = [len; l]; %#ok<AGROW>
    cen = [cen; c]; %#ok<AGROW>
    tr.longest(jj) = max(l);
    tr.censored(jj) = any(c(l == max(l)));
    % Episodes: stretches at thr - margin_dB that hold a volume at thr.
    lo = g_v >= thr - margin_dB;
    dl = diff([0; lo(:); 0]);
    eb = find(dl == 1);
    ee = find(dl == -1) - 1;
    cs = [0; cumsum(hit(:))];
    eh = cs(ee + 1) - cs(eb);                       % high volumes in each stretch
    keep = eh > 0;
    [eb, ee, eh] = deal(eb(keep), ee(keep), eh(keep));
    sp = ee - eb + 1;
    ec = eb == 1 | ee == n;
    ep_len = [ep_len; sp]; %#ok<AGROW>
    ep_cen = [ep_cen; ec]; %#ok<AGROW>
    [tr.ep_span(jj), im] = max(sp);
    tr.ep_high(jj) = eh(im);
    tr.ep_censored(jj) = any(ec(sp == max(sp)));
end
if isempty(len)
    s = [NaN NaN 0 0 0 0 0 0 0 0 0 0];
    return
end
mx = @(v) max([0; v(:)]);
s = [median(len), prctile(len, 95), numel(len), max(len), 100*sum(len(len >= 2))/sum(len), ...
    100*sum(len(len >= 3))/sum(len), mx(len(~cen)), nnz(cen), mx(len(cen)), ...
    max(ep_len), mx(ep_len(~ep_cen)), mx(ep_len(ep_cen))];
end
