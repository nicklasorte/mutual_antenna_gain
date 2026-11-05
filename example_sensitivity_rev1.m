%EXAMPLE_SENSITIVITY_REV1 One-at-a-time sensitivity of the WSR-88D/ASR-11 MAG.
%   Changes one antenna or scan assumption at a time from the worst NEXRAD
%   mode and reports the MAG change in dB at each percentile
%   (mag_sensitivity_rev1). Requested at the 2.7 GHz Retuning LOE meeting
%   of November 4, 2026 (sensitivity section). The steps are study
%   choices, not measured alternatives:
%     NEXRAD pattern     ROC specification envelope; F.1245 (working paper)
%     NEXRAD back lobe   -45 and -65 dB (baseline -55 dB, NSSL far level
%                        held to 180 deg; no measured back lobe)
%     ASR-11 pattern     F.1245 (working paper)
%     ASR-11 remote SL   +0.5 and +6.5 dBi (baseline +3.5 dBi, M.1464-2)
%     ASR-11 tilt        -2.5 and +2.5 deg (M.1464-2 Radar C range)
%     Path elevation     0.5, 1 and 2 deg at both ends (takeoff angles of a
%                        transhorizon path; baseline 0)
%     Rate tolerance     0 and 2 deg/s (baseline 0.5); time series only
%     NEXRAD sector      180, 90, 30 and 10 deg centered on the ASR-11
%                        (bounding case for a scan aimed at the ATC radar)
%   Set sens_base before running to use another base mode; the default is
%   VCP 212 SAILSx3 AVSET 6.4, the worst mode of example_nexrad_asr11_rev1.

if ~exist("sens_base", "var")
    sens_base = struct("seed", 20261105, "vcp", 212, "sails", 3, "avset_top_deg", 6.4, ...
        "n_samples", 5e6, "duration_s", 3600, "n_trials", 4);
end

cs = struct("group", {}, "label", {}, "opts", {});
add = @(cs, g, l, o) [cs, struct("group", g, "label", l, "opts", o)];
cs = add(cs, "NEXRAD pattern", "ROC specification envelope", struct("nexrad_pattern", "roc_spec"));
cs = add(cs, "NEXRAD pattern", "F.1245, 46 dBi", struct("nexrad_pattern", "f1245"));
cs = add(cs, "NEXRAD back lobe", "-45 dB", struct("nexrad_back_dB", -45));
cs = add(cs, "NEXRAD back lobe", "-65 dB", struct("nexrad_back_dB", -65));
cs = add(cs, "ASR-11 pattern", "F.1245, 34 dBi", struct("asr_pattern", "f1245"));
cs = add(cs, "ASR-11 remote sidelobe", "+0.5 dBi", struct("asr_remote_sl_dBi", 0.5));
cs = add(cs, "ASR-11 remote sidelobe", "+6.5 dBi", struct("asr_remote_sl_dBi", 6.5));
cs = add(cs, "ASR-11 tilt", "-2.5 deg", struct("asr_tilt_deg", -2.5));
cs = add(cs, "ASR-11 tilt", "+2.5 deg", struct("asr_tilt_deg", 2.5));
cs = add(cs, "Path elevation", "0.5 deg", struct("el_path_deg", 0.5));
cs = add(cs, "Path elevation", "1 deg", struct("el_path_deg", 1));
cs = add(cs, "Path elevation", "2 deg", struct("el_path_deg", 2));
cs = add(cs, "ASR-11 rate tolerance", "0 deg/s", struct("asr_rate_jitter_deg_s", 0));
cs = add(cs, "ASR-11 rate tolerance", "2 deg/s", struct("asr_rate_jitter_deg_s", 2));
for sw_deg = [180 90 30 10]
    cs = add(cs, "NEXRAD sector", sprintf("%d deg", sw_deg), struct("nexrad_sector_deg", sw_deg));
end

Tsens = mag_sensitivity_rev1(sens_base, cs);
disp(Tsens(:, ["group", "label", "p99_9_dBi", "d99_9_dB", "p99_99_dBi", "d99_99_dB", ...
    "p99_999_dBi", "d99_999_dB", "mbmb_path_dBi", "rev_pct_ge_45_dBi", "max_run_45_revs"]))
