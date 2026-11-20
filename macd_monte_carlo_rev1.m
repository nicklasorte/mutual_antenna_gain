function out = macd_monte_carlo_rev1(ant1, ant2, opts)
%MACD_MONTE_CARLO_REV1 Monte Carlo check of the MACD convolution.
%   out = macd_monte_carlo_rev1(ant1, ant2) draws random, independent
%   pointings of the two antennas and counts the MAG of each draw, as a
%   check on the convolution of macd_rev1. Each draw takes:
%
%     - an azimuth for each antenna, uniform over its scan (no step: the
%       pointing is continuous);
%     - a mode and cut for each antenna, with probability modes.prob;
%     - the gain of each antenna toward the other from its gain_fn, and
%       their sum, the MAG.
%
%   The MAGs go into counts on a fine grid (bin_dB), so billions of draws
%   need no more memory than one chunk. Percentiles and their 95%
%   intervals come from the counts: the interval is the range of MAG
%   between the ranks n q -/+ 1.96 sqrt(n q (1 - q)), the normal
%   approximation to the binomial distribution of the number of draws
%   below the true percentile.
%
%   ant1, ant2   antenna structs as in macd_rev1 (g_max_dBi, gain_fn,
%                scan_width_deg, mainbeam_az_deg, bearing_deg, modes)
%
%   opts (name-value):
%     n_draws        total draws (1e8); rounded up to whole chunks
%     chunk          draws per chunk (1e7)
%     checkpoints    draw counts at which to record the percentiles
%                    ([] = n_draws only); each is rounded to whole chunks
%     bin_dB         count grid (0.001)
%     pcts           percentage points (99 99.9 99.99 99.999 99.9999)
%     seed           RandStream seed (20261120)
%     parallel_flag  1 to spread chunks over a pool with parfor (0)
%
%   out fields:
%     pct        table: n_draws, pct, mag_dBi, ci_lo_dBi, ci_hi_dBi, n_above
%                (one row per checkpoint and percentage)
%     exceed     table: mag_dBi (grid lower edges), exceed (P(MAG >= x))
%                at n_draws
%     n_draws    draws made
%
%   Draws in chunk c use RandStream("mrg32k3a") substream c, so
%   parallel_flag = 1 gives the same counts as parallel_flag = 0.
%
%   Example:
%     a = struct("g_max_dBi",30,"bw_deg",2,"gain_fn",@(az,t) 30*(az<=1), ...
%         "backlobe_dBi",0);
%     m = macd_monte_carlo_rev1(a, a, n_draws=1e6, chunk=1e5);
%
%   DRAFT//PRE-DECISIONAL
%
%   See also macd_rev1, report_figures_macd_rev2.

arguments
    ant1 (1,1) struct
    ant2 (1,1) struct
    opts.n_draws (1,1) double {mustBePositive} = 1e8
    opts.chunk (1,1) double {mustBePositive, mustBeInteger} = 1e7
    opts.checkpoints (1,:) double = []
    opts.bin_dB (1,1) double {mustBePositive} = 0.001
    opts.pcts (:,1) double {mustBeInRange(opts.pcts, 0, 100)} = [99 99.9 99.99 99.999 99.9999]'
    opts.seed (1,1) double {mustBeInteger, mustBeNonnegative} = 20261120
    opts.parallel_flag (1,1) double {mustBeMember(opts.parallel_flag, [0 1])} = 0
end
A = {fill_mc(ant1), fill_mc(ant2)};
A{2}.bearing_deg = mod(A{1}.bearing_deg + 180 + 180, 360) - 180;
n_chunks = ceil(opts.n_draws/opts.chunk);
cp = unique([ceil(opts.checkpoints/opts.chunk), n_chunks]);
cp = cp(cp >= 1 & cp <= n_chunks);

% Count grid: from the lowest possible MAG to the highest.
w = opts.bin_dB;
lo = floor((A{1}.backlobe_dBi + A{2}.backlobe_dBi - 1)/w)*w;
hi = A{1}.g_max_dBi + A{2}.g_max_dBi + 1;
nb = ceil((hi - lo)/w) + 1;
counts = zeros(nb, 1);
rows = {};
done = 0;
for s = cp
    idx = done + 1:s;
    part = zeros(nb, numel(idx));
    seed = opts.seed; ch = opts.chunk;
    if opts.parallel_flag == 1
        parfor j = 1:numel(idx)
            part(:, j) = chunk_counts(A, idx(j), seed, ch, lo, w, nb);
        end
    else
        for j = 1:numel(idx)
            part(:, j) = chunk_counts(A, idx(j), seed, ch, lo, w, nb);
        end
    end
    counts = counts + sum(part, 2);
    done = s;
    n = done*opts.chunk;
    rows = [rows; pct_rows(counts, n, lo, w, opts.pcts)]; %#ok<AGROW>
end
out.pct = cell2table(rows, 'VariableNames', ["n_draws", "pct", "mag_dBi", "ci_lo_dBi", "ci_hi_dBi", "n_above"]);
x = lo + (0:nb - 1)'*w;
ex = flip(cumsum(flip(counts)))/(done*opts.chunk);
k = find(counts > 0, 1):find(counts > 0, 1, "last");
out.exceed = table(x(k), ex(k), 'VariableNames', ["mag_dBi", "exceed"]);
out.n_draws = done*opts.chunk;
end

% -------------------------------------------------------------------------
function c = chunk_counts(A, j, seed, ch, lo, w, nb)
rs = RandStream("mrg32k3a", "Seed", seed);
rs.Substream = j;
g = zeros(ch, 1);
for a = 1:2
    az = A{a}.mainbeam_az_deg + A{a}.scan_width_deg*(rand(rs, ch, 1) - 0.5);
    rel = abs(mod(az - A{a}.bearing_deg + 180, 360) - 180);
    p = A{a}.modes.prob/sum(A{a}.modes.prob);
    if isscalar(p)
        tilt = A{a}.modes.tilt_deg*ones(ch, 1);
    else
        m = discretize(rand(rs, ch, 1), [0; cumsum(p(:))]);
        tilt = A{a}.modes.tilt_deg(m);
    end
    g = g + A{a}.gain_fn(rel, tilt);
end
k = min(max(floor((g - lo)/w) + 1, 1), nb);
c = accumarray(k, 1, [nb 1]);
end

function r = pct_rows(counts, n, lo, w, pcts)
cdf = cumsum(counts);
r = cell(numel(pcts), 6);
val = @(rank) lo + (find(cdf >= min(max(rank, 1), n), 1) - 1)*w;   % lower edge of the bin holding that rank
for i = 1:numel(pcts)
    q = pcts(i)/100;
    h = 1.96*sqrt(n*q*(1 - q));
    r(i, :) = {n, pcts(i), val(ceil(n*q)), val(floor(n*q - h)), val(ceil(n*q + h)), n - cdf(find(cdf >= ceil(n*q), 1))};
end
end

function A = fill_mc(A)
d = struct("scan_width_deg", 360, "mainbeam_az_deg", 0, "bearing_deg", 0, ...
    "modes", table(0, 1, 'VariableNames', ["tilt_deg", "prob"]));
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(A, f{i})
        A.(f{i}) = d.(f{i});
    end
end
for f = ["g_max_dBi", "gain_fn", "backlobe_dBi"]
    if ~isfield(A, f)
        error("macd_monte_carlo_rev1:missingField", "Antenna struct needs field %s.", f);
    end
end
end
