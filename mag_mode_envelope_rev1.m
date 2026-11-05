function env = mag_mode_envelope_rev1(g_list, names, pct)
%MAG_MODE_ENVELOPE_REV1 MAG bound that holds for any mix of scan modes.
%   env = mag_mode_envelope_rev1(g_list, names, pct) takes MAG samples for
%   several scan modes (for example, every WSR-88D VCP and option) and
%   returns the worst-mode envelope. The envelope bounds the MAG of any
%   mixture of the modes, whatever fraction of time a site spends in each.
%
%   Why it is a bound: for mode weights w_i >= 0 that sum to 1, the
%   mixture exceedance is P_mix(G > x) = sum_i w_i P_i(G > x), which is
%   at most max_i P_i(G > x). At x = max_i q_i(p), where q_i(p) is the p
%   percentile of mode i, every P_i(G > x) <= 1 - p/100, so the mixture's
%   p percentile is at most max_i q_i(p). The bound needs no mode
%   statistics, which NOAA reported are site-dependent, operator-selected
%   and not archived in a usable form (2.7 GHz Retuning LOE meeting,
%   November 4, 2026).
%
%   g_list  cell array of MAG sample vectors (dBi), one per mode
%   names   string array of mode names, one per cell
%   pct     percentiles to report (default [99 99.9 99.99 99.999])
%
%   env struct:
%     .pct     table: pct, envelope_dBi, worst_mode, and one column per mode
%     .exceed  table on a 0.1 dB grid: gain_dBi, envelope_exceed_pct,
%              worst_mode (the mode with the largest exceedance at that gain)
%
%   See also nexrad_asr11_mag_rev1, example_nexrad_asr11_rev1.

arguments
    g_list (1,:) cell
    names (1,:) string
    pct (1,:) double = [99 99.9 99.99 99.999]
end
if numel(g_list) ~= numel(names) || isempty(g_list)
    error("radar_coupling:envelopeInput", "g_list and names must be nonempty and the same length.");
end
if any(pct <= 0 | pct >= 100)
    error("radar_coupling:envelopeInput", "pct must lie strictly between 0 and 100.");
end

% ---- percentile envelope
n = numel(g_list);
q = zeros(numel(pct), n);
for i = 1:n
    q(:, i) = prctile(g_list{i}(:), pct(:));
end
[q_env, i_worst] = max(q, [], 2);
P = table(pct(:), q_env, names(i_worst)', 'VariableNames', ["pct", "envelope_dBi", "worst_mode"]);
P = [P, array2table(q, "VariableNames", matlab.lang.makeValidName(names))];

% ---- exceedance envelope on a common 0.1 dB grid
lo = floor(min(cellfun(@min, g_list))*10)/10 - 0.05;
hi = ceil(max(cellfun(@max, g_list))*10)/10 + 0.05;
edges = lo:0.1:(hi + 0.1);
ex = zeros(numel(edges) - 1, n);
for i = 1:n
    c = histcounts(g_list{i}, edges);
    ex(:, i) = 100*(1 - cumsum(c(:))/numel(g_list{i}));
end
[ex_env, j_worst] = max(ex, [], 2);
E = table(edges(2:end)', ex_env, names(j_worst)', ...
    'VariableNames', ["gain_dBi", "envelope_exceed_pct", "worst_mode"]);

env = struct("pct", P, "exceed", E);
end
