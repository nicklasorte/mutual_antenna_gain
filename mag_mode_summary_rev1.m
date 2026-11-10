function s = mag_mode_summary_rev1(g, pct)
%MAG_MODE_SUMMARY_REV1 Compact summary of one scan mode's MAG distribution.
%   s = mag_mode_summary_rev1(g, pct) keeps what mag_mode_envelope_rev1
%   needs from one scan mode's MAG, so a sweep over many modes need not
%   hold every mode's samples in memory. g is one of:
%     - MAG samples (dBi), as from radar_mutual_gain_cdf_rev1;
%     - a distribution struct with fields gain_dBi and pmf (out.dist21 of
%       mag_exact_cdf_rev1);
%     - a CDF output struct (out.cdf of nexrad_asr11_mag_rev1): its g21_dBi
%       samples when present, otherwise its dist21.
%   Fields:
%     .pct      percentiles requested
%     .q_dBi    MAG at each percentile (prctile of the samples, or the
%               smallest gain with P(MAG > gain) <= 1 - pct/100)
%     .edges    fixed 0.1 dB bin edges, -80.05 to 100.05 dBi
%     .exceed   percent of the time above each bin's upper edge (for a
%               distribution, at or above: gains are rounded up)
%     .n        number of samples (Inf for a distribution)
%   The envelope built from summaries equals the envelope built from the
%   samples (same percentiles, same exceedance grid).
%
%   Example:
%     s = mag_mode_summary_rev1(randn(1e5, 1), [99 99.9]);
%
%   See also mag_mode_envelope_rev1, example_nexrad_asr11_rev1,
%   mag_exact_cdf_rev1.

arguments
    g {mustBeNonempty}
    pct (1,:) double {mustBeGreaterThan(pct, 0), mustBeLessThan(pct, 100)}
end
if isstruct(g) && isfield(g, "g21_dBi")             % CDF output struct
    if ~isempty(g.g21_dBi)
        g = g.g21_dBi;
    else
        g = g.dist21;
    end
end
if isstruct(g)
    s = from_distribution(g, pct);
    return
end
if ~isnumeric(g)
    error("radar_coupling:summaryInput", "g must be samples, a distribution struct or a CDF struct.");
end
g = double(g(:));
s.pct = pct;
s.q_dBi = prctile(g, pct(:));
s.edges = (-80.05:0.1:100.05)';
if min(g) < s.edges(1) || max(g) > s.edges(end)
    error("radar_coupling:summaryRange", "MAG samples must lie in [-80, 100] dBi.");
end
c = histcounts(g, s.edges);
s.exceed = 100*(1 - cumsum(c(:))/numel(g));
s.n = numel(g);
end

% -------------------------------------------------------------------------
function s = from_distribution(d, pct)
if ~all(isfield(d, ["gain_dBi", "pmf"]))
    error("radar_coupling:summaryInput", "A distribution needs fields gain_dBi and pmf.");
end
x = d.gain_dBi(:);
p = d.pmf(:)/sum(d.pmf);
E = cumsum(p, "reverse") - p;                    % P(MAG > x)
s.pct = pct;
s.q_dBi = zeros(numel(pct), 1);
for j = 1:numel(pct)
    s.q_dBi(j) = x(find(E <= 1 - pct(j)/100 + 1e-13, 1));
end
s.edges = (-80.05:0.1:100.05)';
nz = p > 0;
if min(x(nz)) < s.edges(1) || max(x(nz)) > s.edges(end)
    error("radar_coupling:summaryRange", "MAG distribution must lie in [-80, 100] dBi.");
end
up = s.edges(2:end);
s.exceed = zeros(numel(up), 1);
for i = 1:numel(up)
    s.exceed(i) = 100*sum(p(x >= up(i) - 1e-9));
end
s.n = Inf;
end
