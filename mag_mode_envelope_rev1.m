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
%   g_list  cell array, one per mode, of MAG sample vectors (dBi) or of
%           mag_mode_summary_rev1 structs (all with the same pct), so a
%           long sweep need not keep every sample
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

is_sum = cellfun(@isstruct, g_list);
if any(is_sum) && ~all(is_sum)
    error("radar_coupling:envelopeInput", "g_list must hold all samples or all summaries.");
end
if all(is_sum)
    env = from_summaries(g_list, names, pct);
    return
end

% ---- percentile envelope
n = numel(g_list);
q = zeros(numel(pct), n);
for i = 1:n
    q(:, i) = prctile(double(g_list{i}(:)), pct(:));
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

% -------------------------------------------------------------------------
function env = from_summaries(S, names, pct)
% Same envelope from mag_mode_summary_rev1 structs.
n = numel(S);
q = zeros(numel(pct), n);
ex = zeros(numel(S{1}.exceed), n);
for i = 1:n
    [ok, k] = ismember(pct, S{i}.pct);
    if ~all(ok) || numel(S{i}.exceed) ~= size(ex, 1)
        error("radar_coupling:envelopeInput", "Every summary needs the requested pct and the same grid.");
    end
    q(:, i) = S{i}.q_dBi(k);
    ex(:, i) = S{i}.exceed;
end
[q_env, i_worst] = max(q, [], 2);
P = table(pct(:), q_env, names(i_worst)', 'VariableNames', ["pct", "envelope_dBi", "worst_mode"]);
P = [P, array2table(q, "VariableNames", matlab.lang.makeValidName(names))];
% Trim the grid to where some mode has samples above the bin.
edges = S{1}.edges;
keep = find(min(ex, [], 2) < 100, 1) - 1;      % first bin with samples in any mode
last = find(max(ex, [], 2) > 0, 1, "last") + 1;
r = max(keep, 1):min(last, size(ex, 1));
[ex_env, j_worst] = max(ex(r, :), [], 2);
E = table(edges(r + 1), ex_env, names(j_worst)', ...
    'VariableNames', ["gain_dBi", "envelope_exceed_pct", "worst_mode"]);
env = struct("pct", P, "exceed", E);
end
