function s = mag_mode_summary_rev1(g, pct)
%MAG_MODE_SUMMARY_REV1 Compact summary of one scan mode's MAG samples.
%   s = mag_mode_summary_rev1(g, pct) keeps what mag_mode_envelope_rev1
%   needs from the MAG samples g (dBi) of one scan mode, so a sweep over
%   many modes need not hold every mode's samples in memory:
%     .pct      percentiles requested
%     .q_dBi    MAG at each percentile (prctile of g)
%     .edges    fixed 0.1 dB bin edges, -80.05 to 100.05 dBi
%     .exceed   percent of samples above each bin's upper edge
%     .n        number of samples
%   The envelope built from summaries equals the envelope built from the
%   samples (same percentiles, same exceedance grid).
%
%   Example:
%     s = mag_mode_summary_rev1(randn(1e5, 1), [99 99.9]);
%
%   See also mag_mode_envelope_rev1, example_nexrad_asr11_rev1.

arguments
    g {mustBeNumeric, mustBeNonempty}
    pct (1,:) double {mustBeGreaterThan(pct, 0), mustBeLessThan(pct, 100)}
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
