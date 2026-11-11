%MAKE_REV6_REFERENCE_REV1 Pin the revision 6 outputs before the speed work (run once).
%   Writes, under the revision 6 defaults (ROC specification pattern, 216
%   scan modes, ASR-11):
%     test_data/rev6_pair_reference_rev1.csv  mag_pair_cdf_rev1 envelope
%         percentiles and worst modes at six takeoff-angle pairs
%     test_data/rev6_ts_reference_rev1.csv    mag_pair_timeseries_rev1
%         per-scan maxima (every ATC revolution, WSR-88D cut and volume) of
%         four 1 h trials at the worst mode, takeoff angles 0.5/0.3 deg,
%         seed 20261105, and the harm counts per 0.5 deg cut
%   test_regression_rev6 checks later code against both files: the speed
%   work (docs/rev6_redteam_plan_rev1.md, V1 to V8) must reproduce them
%   exactly. Rerun only to re-pin on purpose, and say why in the commit
%   message.
%
%   See also rev6_reference_rev1, test_regression_rev6, mag_pair_cdf_rev1, mag_pair_timeseries_rev1.

here = fileparts(which("make_rev6_reference_rev1"));
out_dir = fullfile(here, "test_data");
if ~isfolder(out_dir), mkdir(out_dir); end
[ref_pair, ref_ts] = rev6_reference_rev1();
writetable(ref_pair, fullfile(out_dir, "rev6_pair_reference_rev1.csv"));
writetable(ref_ts, fullfile(out_dir, "rev6_ts_reference_rev1.csv"));
fprintf("Wrote %d pair rows and %d time-series rows to %s\n", height(ref_pair), height(ref_ts), out_dir);
