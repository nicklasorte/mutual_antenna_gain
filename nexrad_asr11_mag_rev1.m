function out = nexrad_asr11_mag_rev1(opts)
%NEXRAD_ASR11_MAG_REV1 Mutual antenna gain between a WSR-88D and an ASR-11.
%   out = nexrad_asr11_mag_rev1(opts) computes the mutual antenna gain
%   (MAG) between a WSR-88D (NEXRAD) weather radar running a volume
%   coverage pattern (VCP) and an ASR-11 airport surveillance radar:
%     g12 = G_tx,NEXRAD(toward ASR-11) + G_rx,ASR-11(toward NEXRAD)  (dBi)
%     g21 = G_tx,ASR-11(toward NEXRAD) + G_rx,NEXRAD(toward ASR-11)  (dBi)
%   Each radar uses one pattern to transmit and receive, so g12 = g21.
%
%   It returns:
%     1. The time-weighted MAG CDF (radar_mutual_gain_cdf_rev1, random
%        mode, Rec. ITU-R M.2069-0 Annex 1 method).
%     2. The maximum MAG per scan iteration, from a time series in which
%        both antennas rotate from random start phases: per ASR-11
%        revolution, per NEXRAD elevation cut, and per NEXRAD volume.
%
%   Every parameter value is set in assumptions() below with a source key;
%   citations() lists each source with its URL for the working paper.
%   out.assumptions and out.citations return both as tables.
%
%   opts struct (defaults in parentheses):
%     .seed             required RNG seed
%     .nexrad_scan      "vcp" (default) | "constant". "constant" is the
%                       2.7 GHz link budget working paper model: one
%                       horizon cut at 3 rpm (Table 4.4-1).
%     .vcp              12 | 31 | 32 | 34 | 35 | 112 | 212 (default) | 215
%     .sails, .mrle, .avset_top_deg, .cd_rate_deg_s
%                       VCP options, see nexrad_vcp_rev1 (0, 0, Inf, [])
%     .nexrad_pattern   "nssl_measured" (default) | "roc_spec" | "f1245"
%     .asr_pattern      "m1851_csc2" (default) | "f1245"
%     .asr_tilt_deg     ASR-11 beam tilt (0, TODO: site value)
%     .asr_rate_jitter_deg_s  ASR-11 rotation-rate variance, uniform in
%                       +/- this value about 75 deg/s (0.5, M.2069-0 Annex 1
%                       drive accuracy; the FAA tolerance is in Order JO
%                       6310.30B, FAA-internal, TODO)
%     .asr_jitter_mode  "per_trial" (default, one offset per trial) |
%                       "per_rev" (new offset every revolution)
%     .el_path_deg      path elevation at both ends (0: horizontal
%                       line-of-sight path, location-agnostic geometry)
%     .n_samples        random-mode samples for the CDF (5e6)
%     .duration_s       time series length per trial (3600)
%     .n_trials         time series trials (2)
%     .dt_s             time series step (0.002 s: about one tenth of the
%                       ASR-11 beam transit time, 1.45/75 = 0.019 s)
%     .pct              percentiles reported ([50 90 99 99.9 99.99 99.999])
%     .parallel_flag    0 (default) | 1, passed to the CDF run
%     .cdf_pct          percent grid of out.cdf.cdf_array12/21, 0 to 100
%                       (default 0:1:100 plus 99.5 to 99.999)
%     .keep_ts1         true: return trial 1 of the time series in
%                       out.scan_max.ts1 (false; about 70 MB per hour)
%     .plot_flag        true: plot the exceedance curve with the reference
%                       levels (false)
%
%   out struct:
%     .cdf              radar_mutual_gain_cdf_rev1 output (random mode),
%                       including cdf_array21: [percent 0-100, MAG dBi]
%     .ref_levels       table of main-beam and sidelobe combinations on the
%                       link path: combination, mag_dBi, nexrad_gain_dBi,
%                       asr11_gain_dBi, time_exceeded_pct (time-weighted CDF),
%                       asr_rev_reached_pct (% of ASR-11 revolutions whose
%                       maximum reaches the level, within the 0.2 dB
%                       time-step loss). Rows: main beam to main
%                       beam; NEXRAD main beam + ASR-11 near sidelobe and back
%                       lobe; ASR-11 main beam + NEXRAD near sidelobe and
%                       back lobe. Main beam = beam azimuth on the path (NEXRAD
%                       at its lowest cut); near sidelobe = largest gain at
%                       >= 2 beamwidths off the path; back lobe = 180 deg.
%     .sensitivity      table of MAG choices for the sensitivity analysis
%                       requested at the November 4, 2026 Retuning LOE
%                       meeting: peak-gain main beam to main beam (Gmax +
%                       Gmax), main beam to main beam on the path, largest
%                       time-series MAG, and the 99% to 99.999% MAG. Columns:
%                       mag_choice, mag_dBi, delta_vs_peak_mbmb_dB,
%                       time_exceeded_pct, asr_rev_reached_pct
%     .scan_max.asr_rev      table, one row per ASR-11 revolution: trial,
%                            rev, t_start_s, rate_deg_s, max_g12_dBi, max_g21_dBi,
%                            nexrad_el_at_max_deg, complete
%     .scan_max.nexrad_cut   table, one row per NEXRAD elevation cut:
%                            trial, volume, cut, el_deg, t_start_s,
%                            max_g12_dBi, max_g21_dBi, complete
%     .scan_max.nexrad_volume table, one row per NEXRAD volume: trial,
%                            volume, t_start_s, max_g12_dBi, max_g21_dBi,
%                            el_at_max_deg, complete
%                            The asr_rev table also gives, at each
%                            revolution maximum, each radar's gain and
%                            beam azimuth off the path (explains each peak).
%     .scan_max.plateaus     hand check of the per-revolution plateaus:
%                            NEXRAD main beam on the path at each tilt +
%                            ASR-11 floor, with predicted and simulated
%                            share of revolutions reaching each level
%     .scan_max.ts1          trial 1 samples when opts.keep_ts1 is true
%     .scan_max.by_elevation  per NEXRAD elevation: number of complete
%                            cuts, median and maximum of the cut maxima
%     .scan_max.asr_rev_pct   percentiles of the per-revolution maxima
%     .ts_pct           MAG percentiles from the time series samples (a
%                       check on the random-mode CDF)
%     .sys_nexrad, .sys_asr11   systems passed to the model
%     .assumptions      table: parameter, value, unit, source, note
%     .citations        table: key, reference, url
%     .opts             options used
%
%   Partial scans at the start and end of a trial have complete = false;
%   leave them out of scan statistics (by_elevation and asr_rev_pct do).
%
%   Example:
%     out = nexrad_asr11_mag_rev1(struct("seed", 1, "vcp", 212));
%     disp(out.cdf.pct21); disp(out.scan_max.by_elevation)
%
%   See also radar_mutual_gain_cdf_rev1, nexrad_vcp_rev1,
%   antenna_prepare_rev1, antenna_gain_rev1.

arguments
    opts (1,1) struct
end

opts = fill_opts(opts);
A = assumptions();
[sys_n, sys_a] = build_systems(A, opts);

% ---- 1. Time-weighted CDF (M.2069-0 random mode). The constant-rate
%      case has rates in a 6:25 ratio; random mode assumes the rates drift
%      (drive tolerance), so the commensurate-rate warning is expected.
o = struct("seed", opts.seed, "mode", "random", "n_samples", opts.n_samples, ...
    "el12_deg", opts.el_path_deg, "el21_deg", opts.el_path_deg, ...
    "pct", opts.pct, "parallel_flag", opts.parallel_flag, "cdf_pct", opts.cdf_pct);
w = warning("off", "radar_coupling:commensurateRates");
cleanup = onCleanup(@() warning(w));
cdf = radar_mutual_gain_cdf_rev1(sys_n, sys_a, o);
clear cleanup

% ---- 2. Time series and maximum MAG per scan iteration.
[scan_max, ts_pct] = scan_maxima(sys_n, sys_a, opts);

% ---- 3. Main-beam / sidelobe reference levels for the CDF plot.
ref = ref_levels(sys_n, sys_a, opts, cdf, scan_max);

% ---- 4. Sensitivity of the MAG choice against main beam to main beam.
sens = sensitivity_table(sys_n, sys_a, cdf, ref, scan_max, opts);

out = struct();
out.cdf = cdf;
out.ref_levels = ref;
out.sensitivity = sens;
out.scan_max = scan_max;
out.ts_pct = ts_pct;
out.sys_nexrad = sys_n;
out.sys_asr11 = sys_a;
out.assumptions = A.table;
out.citations = citations();
out.opts = opts;
if opts.plot_flag
    plot_cdf(out);
end
end

% =========================================================================
function R = ref_levels(sys_n, sys_a, opts, cdf, S)
% Main-beam and sidelobe combinations along the link path. Each antenna's
% azimuth cut toward the other station is evaluated at its own tilt (the
% NEXRAD at its lowest cut, where its main beam is closest to the path):
%   main beam    gain with the beam azimuth on the path
%   near sidelobe  largest gain with the beam azimuth >= 2 beamwidths
%                  off the path (first-sidelobe region)
%   back lobe    gain with the beam azimuth 180 deg off the path
pn = antenna_prepare_rev1(sys_n.tx);
pa = antenna_prepare_rev1(sys_a.tx);
e = opts.el_path_deg;
x = (0:0.01:180)';
[~, el_n] = nexrad_cuts(sys_n.scan);
gn = antenna_gain_rev1(pn, 0, e, state(-x, min(el_n) + 0*x));
ga = antenna_gain_rev1(pa, 180, e, state(180 - x, sys_a.scan.tilt_deg + 0*x));
% [main, near SL, back]; 2 beamwidths = 4 x the 3 dB half-width of the cut
lev = @(g) [g(1), max(g(x >= 4*x(find(g <= g(1) - 3, 1)))), g(end)];
Ln = lev(gn);
La = lev(ga);
name = ["Main beam to main beam"
    "NEXRAD main beam + ASR-11 near sidelobe"
    "NEXRAD main beam + ASR-11 back lobe"
    "ASR-11 main beam + NEXRAD near sidelobe"
    "ASR-11 main beam + NEXRAD back lobe"];
g = [Ln(1) + La(1); Ln(1) + La(2); Ln(1) + La(3); La(1) + Ln(2); La(1) + Ln(3)];
g_n = [Ln(1); Ln(1); Ln(1); Ln(2); Ln(3)];
g_a = [La(1); La(2); La(3); La(1); La(1)];
exceed = arrayfun(@(v) 100*mean(cdf.g21_dBi >= v), g);
Rv = S.asr_rev(S.asr_rev.complete, :);
% The time-series peak can fall up to 0.2 dB short of the true peak (ASR-11
% beam moves 0.15 deg per 0.002 s step); count a level reached within that.
rev_pct = arrayfun(@(v) 100*mean(Rv.max_g21_dBi >= v - 0.2), g);
R = table(name, g, g_n, g_a, exceed, rev_pct, 'VariableNames', ["combination", ...
    "mag_dBi", "nexrad_gain_dBi", "asr11_gain_dBi", "time_exceeded_pct", "asr_rev_reached_pct"]);
end

% -------------------------------------------------------------------------
function T = sensitivity_table(sys_n, sys_a, cdf, ref, S, opts)
% MAG choices from most to least conservative, each with its difference
% from the peak-gain main beam to main beam value (the FAA practice, 2.7
% GHz Retuning LOE meeting of November 4, 2026). A negative delta is the
% relief, in dB, that the choice takes relative to that case; the same dB
% applies one for one to the required path loss.
g_pk = antenna_prepare_rev1(sys_n.tx).g_max_dBi + antenna_prepare_rev1(sys_a.tx).g_max_dBi;
q = cdf.pct21(cdf.pct21.pct >= 99, :);
name = ["Peak gains, main beam to main beam"; "Main beam to main beam on the path"; ...
    "Largest MAG in the time series"; compose("%g%% time-weighted MAG", q.pct)];
g = [g_pk; ref.mag_dBi(1); max(S.nexrad_volume.max_g21_dBi); q.gain_dBi];
exceed = arrayfun(@(v) 100*mean(cdf.g21_dBi >= v), g);
Rv = S.asr_rev(S.asr_rev.complete, :);
rev_pct = arrayfun(@(v) 100*mean(Rv.max_g21_dBi >= v - 0.2), g);
T = table(name, g, g - g_pk, exceed, rev_pct, 'VariableNames', ["mag_choice", ...
    "mag_dBi", "delta_vs_peak_mbmb_dB", "time_exceeded_pct", "asr_rev_reached_pct"]);
T.Properties.Description = sprintf("el_path_deg = %g; NEXRAD lowest cut %.1f deg", ...
    opts.el_path_deg, min(sys_n.scan.tilt_deg));
end

% -------------------------------------------------------------------------
function plot_cdf(out)
% MAG exceedance (ASR-11 -> NEXRAD) with the main-beam/sidelobe levels.
c = out.cdf.cdf21;
figure("Name", "WSR-88D / ASR-11 MAG");
semilogy(c.gain_dBi, max(c.exceed_pct, 1e-6), "LineWidth", 1.5);
hold on; grid on;
ylim([1e-4 100]);
R = out.ref_levels;
sty = ["-", "--", ":", "--", ":"];
for i = 1:height(R)
    xline(R.mag_dBi(i), sty(i), sprintf("%s (%.1f dBi)", R.combination(i), R.mag_dBi(i)), ...
        "LabelOrientation", "aligned", "LabelVerticalAlignment", "bottom", "FontSize", 8);
end
xlabel("Mutual antenna gain (dBi)");
ylabel("Time exceeded (%)");
title("WSR-88D and ASR-11 mutual antenna gain");
end

% =========================================================================
function A = assumptions()
% Every model input, with its source key (see citations()).
% NEXRAD = WSR-88D; ASR-11 = M.1464-2 Table 1 Radar C (25 kW solid state,
% 1.0/89 us pulses, 89:1 compression: the ASR-11 transmitter of FI_ASR11).

% ---- WSR-88D antenna
A.nex.g_max_dBi = 45.5;     % ROC_TECH: 45.5 dB at 2850 MHz incl. radome (M1849: 45.0)
A.nex.bw_deg = 0.925;       % ROC_TECH: 0.925 deg at 2850 MHz (M1849: 0.92)
% Sidelobe envelopes (dB below peak) versus off-axis angle (deg):
% NSSL_1998: KOUN measurements, first sidelobe about -35 dB at 2 deg,
% falling linearly to about -48 dB at 6 deg and -55 dB at 20 deg;
% sidelobes beyond 20 deg practically all below -55 dB. No back-lobe data,
% so -55 dB is held to 180 deg (TODO: measured back lobe).
A.nex.nssl_env_deg = [0 2 6 20 180];
A.nex.nssl_env_dB = [-35 -35 -48 -55 -55];
% ROC_TECH specification: sidelobes 29 dB below peak beyond the main beam,
% 40 dB below peak beyond 10 deg (an upper bound, held flat).
A.nex.roc_env_deg = [0 10 10.001 180];   % 0.001 deg ramp stands in for the step at 10 deg
A.nex.roc_env_dB = [-29 -29 -40 -40];
A.nex.f1245_g_max_dBi = 46; % WP_REV3 Table 4.4-1 (F.1245 model)
A.nex.paper_rate_rpm = 3;   % WP_REV3 Table 4.4-1
A.nex.paper_tilt_deg = 0;   % WP_REV3 gives no elevation; horizon assumed

% ---- ASR-11 antenna (M1464 Table 1, Radar C)
A.asr.g_max_dBi = 34;       % M1464 Radar C (FI_ASR11: > 33.5 dB; WEBER_1999: ~34 dB)
A.asr.az_bw_deg = 1.45;     % M1464 Radar C (FI_ASR11 and WIKI_ASR11: 1.40)
A.asr.el_bw_deg = 4.8;      % M1464 Radar C (WEBER_1999: 5 deg)
A.asr.csc2_end_deg = 30;    % M1464 Radar C: cosecant-squared to +30 deg
A.asr.first_sl_dBi = 9.5;   % M1464 Radar C first sidelobe level
A.asr.remote_sl_dBi = 3.5;  % M1464 Radar C remote sidelobe level
A.asr.az_sll_dB = A.asr.first_sl_dBi - A.asr.g_max_dBi;    % -24.5 dB -> M1851 Table 9 cosine taper
A.asr.floor_dB = A.asr.remote_sl_dBi - A.asr.g_max_dBi;    % -30.5 dB below peak
A.asr.tilt_range_deg = [-2.5 2.5];  % M1464 Radar C vertical scan type; default tilt 0 (TODO)
A.asr.rate_rpm = 12.5;      % M1464 Radar C 75 deg/s; FI_ASR11 12.5 rpm; WP_REV3 Table 4.4-1
% Rotation-rate variance. The FAA standard and tolerance for the ASR-11
% antenna rotation rate are in FAA_JO6310_30B (maintenance handbook,
% FAA-network only; not public). Until that value is in hand, use the
% M2069 Annex 1 drive accuracy of +/-0.5 deg/s (TODO: JO 6310.30B value).
A.asr.rate_jitter_deg_s = 0.5;
A.asr.f1245_g_max_dBi = 34; % WP_REV3 Table 4.4-1 (F.1245 model)

% ---- Table for out.assumptions
p = {
"NEXRAD peak gain",              A.nex.g_max_dBi, "dBi", "ROC_TECH", "2850 MHz, includes radome; M1849 Radar 1 gives 45.0"
"NEXRAD 3 dB beamwidth",         A.nex.bw_deg, "deg", "ROC_TECH", "M1849 Radar 1 gives 0.92; NSSL_1998 0.85-0.95 across the band"
"NEXRAD main beam shape",        NaN, "", "M1849", "Gmax - 12(theta/bw)^2 until it meets the sidelobe envelope (fills the first null)"
"NEXRAD sidelobes, nssl_measured", NaN, "dB", "NSSL_1998", "-35 dB to 2 deg, -48 at 6, -55 at 20, -55 to 180 (back lobe TODO)"
"NEXRAD sidelobes, roc_spec",    NaN, "dB", "ROC_TECH", "-29 dB beyond main beam, -40 dB beyond 10 deg (specification bound)"
"NEXRAD F.1245 peak gain",       A.nex.f1245_g_max_dBi, "dBi", "WP_REV3", "Table 4.4-1; F.1245 overestimates pencil-beam width (M1849 Sec. 7.3)"
"NEXRAD scan",                   NaN, "", "ROC_VCP_SHEET", "VCP cut tables decoded from Level II (LEVEL2_AWS, METPY); see nexrad_vcp_rev1"
"NEXRAD paper scan rate",        A.nex.paper_rate_rpm, "rpm", "WP_REV3", "Table 4.4-1, constant-rate model"
"NEXRAD polarization",           NaN, "", "ROC_TECH", "Dual linear H+V; polarization mismatch not modeled (study assumption)"
"ASR-11 peak gain",              A.asr.g_max_dBi, "dBi", "M1464", "Radar C; FI_ASR11 > 33.5 dB; WEBER_1999 ~34 dB"
"ASR-11 azimuth beamwidth",      A.asr.az_bw_deg, "deg", "M1464", "Radar C; FI_ASR11 1.40"
"ASR-11 elevation beamwidth",    A.asr.el_bw_deg, "deg", "M1464", "Radar C; WEBER_1999 5 deg"
"ASR-11 elevation pattern",      A.asr.csc2_end_deg, "deg", "M1464", "Cosecant-squared to +30 deg, M1851 Eqs (25)-(27)"
"ASR-11 azimuth first sidelobe", A.asr.az_sll_dB, "dB", "M1464", "+9.5 dBi first SL; selects M1851 Table 9 cosine taper"
"ASR-11 floor (remote sidelobes)", A.asr.floor_dB, "dB", "M1464", "+3.5 dBi remote SL relative to 34 dBi"
"ASR-11 sidelobe envelope",      NaN, "", "M1851", "Peak envelope (single-entry interference), Section 2.1.3"
"ASR-11 tilt",                   0, "deg", "M1464", "Radar C range +2.5 to -2.5; 0 puts the beam nose on the horizon (TODO: site value)"
"ASR-11 beams",                  NaN, "", "WEBER_1999", "Low beam assumed for tx and rx (used beyond about 15 nmi); high beam not modeled"
"ASR-11 rotation rate",          A.asr.rate_rpm, "rpm", "M1464", "Radar C 75 deg/s; FI_ASR11 12.5 rpm; WP_REV3 12.5 rpm"
"ASR-11 rotation rate variance", A.asr.rate_jitter_deg_s, "deg/s", "M2069", "Uniform +/-0.5 deg/s drive accuracy; FAA tolerance is in FAA_JO6310_30B (FAA-internal, TODO)"
"ASR-11 F.1245 peak gain",       A.asr.f1245_g_max_dBi, "dBi", "WP_REV3", "Table 4.4-1"
"Path elevation",                0, "deg", "", "Location-agnostic horizontal line-of-sight path (study assumption)"
"MAG method",                    NaN, "", "M2069", "Annex 1: time-weighted independent pointing (random mode) and time series"
};
A.table = cell2table(p, "VariableNames", ["parameter", "value", "unit", "source", "note"]);
end

% =========================================================================
function C = citations()
% Sources for the assumptions, with URLs, for the working paper.
c = {
"ROC_TECH", "NOAA Radar Operations Center, ""NEXRAD Technical Information"" (WSR-88D antenna, transmitter, pedestal and VCP summary; RDA/RPG Build 23).", ...
    "https://roc.noaa.gov/public-documents/wsr88d/NEXRAD-Technical-Information.pdf"
"ROC_VCP_SHEET", "NOAA Radar Operations Center, ""Convective VCP Decision Matrix"" (VCP Info Sheet), last updated February 2026.", ...
    "https://www.roc.noaa.gov/public-documents/operations-branch/VCP_Info_Sheet.pdf"
"ROC_VCP112", "NOAA Radar Operations Center, ""Theory and Concept of Operations for VCP 112,"" 19 March 2019.", ...
    "https://roc.noaa.gov/WSR88D/PublicDocs/NewTechnology/Theory_ConOps_VCP112.pdf"
"ROC_VCP_PARADIGM", "NOAA Radar Operations Center, ""New VCP Paradigm,"" October 2015.", ...
    "https://www.weather.gov/media/roc/NewScience/New_VCP_Paradigm_Public_Oct_2015.pdf"
"LEVEL2_AWS", "NOAA NEXRAD Level II data on AWS (unidata-nexrad-level2 bucket); VCP messages (Message Type 5) from KDIX, KLOT, KFWS, KTBW, KOAX, KHGX, KLWX files.", ...
    "https://www.ncei.noaa.gov/products/radar/next-generation-weather-radar"
"METPY", "Unidata MetPy 1.7, metpy.io.Level2File (Level II decoder used to read the VCP messages).", ...
    "https://unidata.github.io/MetPy/"
"NSSL_1998", "NSSL, errata and supplement to the 1998 WSR-88D polarimetric upgrade report (D3a): KOUN antenna beamwidth and sidelobe measurements, Section II.", ...
    "https://www.nssl.noaa.gov/publications/wsr88d_reports/ERRATA1998Rpt_D3a(WC).pdf"
"M1849", "Recommendation ITU-R M.1849-3 (02/2023), Technical and operational aspects of ground-based meteorological radars; Annex 2 Table 6 (Radar 1) and Section 7.3 (antenna patterns).", ...
    "https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1849-3-202302-I!!PDF-E.pdf"
"M1464", "Recommendation ITU-R M.1464-2 (02/2015), Characteristics of radiolocation radars ... in the frequency band 2 700-2 900 MHz; Table 1, Radar C.", ...
    "https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1464-2-201502-I!!PDF-E.pdf"
"M1851", "Recommendation ITU-R M.1851-2, Mathematical models for radiodetermination radar systems antenna patterns for use in interference analyses.", ...
    "https://www.itu.int/rec/R-REC-M.1851/en"
"M2069", "Recommendation ITU-R M.2069-0, Technical characteristics and protection criteria for radiodetermination ... (Annex 1 mutual antenna gain method).", ...
    "https://www.itu.int/rec/R-REC-M.2069/en"
"FI_ASR11", "Forecast International, ""ASR-11 (DASR),"" Radar Forecast, August 2014 (gain > 33.5 dB, 1.40 deg azimuth beamwidth, 12.5 rpm, 1 and 89 us pulses).", ...
    "https://www.forecastinternational.com/archive/disp_pdf.cfm?DACH_RECNO=1032"
"WEBER_1999", "M. E. Weber, ""FAA Surveillance Radar Data as a Complement to the WSR-88D Network,"" MIT Lincoln Laboratory ATC-266, 1999, Section 3.1 (ASR-11 beams).", ...
    "https://archive.ll.mit.edu/mission/aviation/publications/publication-files/atc-reports/Weber_1999_ATC-266_WW-15318.pdf"
"WIKI_ASR11", "Wikipedia, ""ASR-11"" (12.5 rpm, 1.4 deg horizontal and 5 deg vertical beamwidth, 25 kW).", ...
    "https://en.wikipedia.org/wiki/ASR-11"
"FAA_ASR11", "FAA, ""Airport Surveillance Radar (ASR-11)"" (2700-2900 MHz, 25 kW, 60 nmi).", ...
    "https://www.faa.gov/air_traffic/technology/asr-11/"
"FAA_JO6310_30B", "FAA Order JO 6310.30B, ""Maintenance Handbook for Airport Surveillance Radar, ASR-11 Facilities,"" 20 July 2015 (standards and tolerances, including antenna rotation rate; accessible only from the FAA network).", ...
    "https://www.faa.gov/regulations_policies/orders_notices/index.cfm/go/document.information/documentID/1028067"
"WP_REV3", "DRAFT 2.7 GHz Link Budget Analysis, Rev3 redlines, 27 Oct 2026, Table 4.4-1 (F.1245 model, 34 and 46 dBi, 12.5 and 3 rpm). Project document, no public URL.", ...
    ""
};
C = cell2table(c, "VariableNames", ["key", "reference", "url"]);
end

% =========================================================================
function [sys_n, sys_a] = build_systems(A, opts)
% ---- NEXRAD antenna
switch opts.nexrad_pattern
    case "nssl_measured"
        ant_n = pencil_table(A.nex.g_max_dBi, A.nex.bw_deg, A.nex.nssl_env_deg, A.nex.nssl_env_dB);
    case "roc_spec"
        ant_n = pencil_table(A.nex.g_max_dBi, A.nex.bw_deg, A.nex.roc_env_deg, A.nex.roc_env_dB);
    case "f1245"
        ant_n = struct("type", "f1245", "g_max_dBi", A.nex.f1245_g_max_dBi);
end
% ---- NEXRAD scan
switch opts.nexrad_scan
    case "vcp"
        scan_n = nexrad_vcp_rev1(opts.vcp, sails=opts.sails, mrle=opts.mrle, ...
            avset_top_deg=opts.avset_top_deg, cd_rate_deg_s=opts.cd_rate_deg_s);
    case "constant"
        scan_n = struct("type", "mech_rotating", "rate_rpm", A.nex.paper_rate_rpm, ...
            "tilt_deg", A.nex.paper_tilt_deg);
end
sys_n = struct("name", "WSR-88D", "tx", ant_n, "rx", ant_n, "scan", scan_n);

% ---- ASR-11 antenna (low beam, tx = rx)
switch opts.asr_pattern
    case "m1851_csc2"
        ant_a = struct("type", "rect", "g_max_dBi", A.asr.g_max_dBi, ...
            "az_bw_deg", A.asr.az_bw_deg, "az_sll_dB", A.asr.az_sll_dB, ...
            "el_bw_deg", A.asr.el_bw_deg, "el_pattern", "csc2", ...
            "csc2_end_deg", A.asr.csc2_end_deg, "csc2_floor_dB", A.asr.floor_dB, ...
            "envelope", "peak", "floor_dB", A.asr.floor_dB);
    case "f1245"
        ant_a = struct("type", "f1245", "g_max_dBi", A.asr.f1245_g_max_dBi);
end
if opts.asr_tilt_deg < A.asr.tilt_range_deg(1) || opts.asr_tilt_deg > A.asr.tilt_range_deg(2)
    error("radar_coupling:asrTilt", "asr_tilt_deg must lie in the M.1464-2 Radar C range [%g %g] deg.", ...
        A.asr.tilt_range_deg(1), A.asr.tilt_range_deg(2));
end
scan_a = struct("type", "mech_rotating", "rate_rpm", A.asr.rate_rpm, "tilt_deg", opts.asr_tilt_deg, ...
    "rate_jitter_deg_s", opts.asr_rate_jitter_deg_s, "jitter_mode", opts.asr_jitter_mode);
sys_a = struct("name", "ASR-11", "tx", ant_a, "rx", ant_a, "scan", scan_a);
end

% -------------------------------------------------------------------------
function ant = pencil_table(g_max, bw, env_deg, env_dB)
% Pencil beam: Gaussian main beam, Gmax - 12*(theta/bw)^2 (the M.1851 /
% M.1849 parabolic main-lobe form), then the sidelobe envelope where the
% main beam falls below it. Filling the first null is conservative.
th = unique([(0:0.005:5)'; (5:0.05:180)'; env_deg(:)]);
env = interp1(env_deg, env_dB, th, "linear");
main = -12*(th/bw).^2;
x = find(main < env, 1);                        % main beam meets the envelope
g = main;
g(x:end) = env(x:end);
ant = struct("type", "offaxis_table", "g_max_dBi", g_max, "offaxis_deg", th, ...
    "gain_dB", g, "az_bw_deg", bw);
end

% =========================================================================
function [S, ts_pct] = scan_maxima(sys_n, sys_a, opts)
% Time series of both antennas from random start phases; maximum MAG per
% ASR-11 revolution, NEXRAD cut and NEXRAD volume.
pn = struct("tx", antenna_prepare_rev1(sys_n.tx), "rx", antenna_prepare_rev1(sys_n.rx));
pa = struct("tx", antenna_prepare_rev1(sys_a.tx), "rx", antenna_prepare_rev1(sys_a.rx));
[rate_n, el_n] = nexrad_cuts(sys_n.scan);
dur_n = 360./rate_n;
t_vol = sum(dur_n);
cut_start = [0; cumsum(dur_n(1:end-1))];
rate_a = 6*sys_a.scan.rate_rpm;
t_rev = 360/rate_a;

t = (0:opts.dt_s:opts.duration_s)';
n_t = numel(t);
[rev_tab, cut_tab, vol_tab] = deal(cell(opts.n_trials, 1));
% Time-series percentiles from a 0.01 dB histogram, so long runs (many
% trials) need not keep every sample in memory.
h_edges = (-80:0.01:100)';
h_counts = zeros(numel(h_edges) - 1, 1);
S.ts1 = [];
for j = 1:opts.n_trials
    rs = RandStream("mrg32k3a", "Seed", opts.seed);
    rs.Substream = 1000 + j;                    % apart from the CDF substreams
    az0_n = 360*rand(rs);
    t0_n = t_vol*rand(rs);
    az0_a = 360*rand(rs);

    % NEXRAD: volume, cut, azimuth (one revolution per cut, as in the model)
    tau = t + t0_n;
    vol = floor(tau/t_vol) + 1;
    w = tau - (vol - 1)*t_vol;
    cut = discretize(w, [cut_start; t_vol]);
    cut(isnan(cut)) = numel(dur_n);
    st_n = state(mod(az0_n + rate_n(cut).*(w - cut_start(cut)), 360), el_n(cut));
    % ASR-11: revolution counted from azimuth 0 (pointing along the link
    % direction 0 deg); one pass of the NEXRAD bearing per revolution.
    % Rate variance, uniform in +/- jitter about 75 deg/s:
    %   "per_trial"  one rate offset per trial (M2069 Annex 1 treatment)
    %   "per_rev"    a new rate offset each revolution
    jit = sys_a.scan.rate_jitter_deg_s;
    switch sys_a.scan.jitter_mode
        case "per_trial"
            r = rate_a + jit*(2*rand(rs) - 1);
            ang = az0_a + r*t;
            rev = floor(ang/360) + 1;
            r_rev = r + zeros(max(rev), 1);
        case "per_rev"
            n_rev = ceil((az0_a + (rate_a + jit)*opts.duration_s)/360) + 1;
            r_rev = rate_a + jit*(2*rand(rs, n_rev, 1) - 1);
            tb = [0; cumsum(360./r_rev)] - az0_a/r_rev(1);   % revolution start times
            rev = discretize(t, tb);
            ang = 360*(rev - 1) + (t - tb(rev)).*r_rev(rev);
    end
    st_a = state(mod(ang, 360), sys_a.scan.tilt_deg + zeros(n_t, 1));

    % Geometry as in radar_mutual_gain_cdf_rev1: radar 1 = NEXRAD, the
    % ASR-11 at 0 deg from it; the NEXRAD at 180 deg from the ASR-11.
    e = opts.el_path_deg;
    g12 = antenna_gain_rev1(pn.tx, 0, e, st_n) + antenna_gain_rev1(pa.rx, 180, e, st_a);
    gA = antenna_gain_rev1(pa.tx, 180, e, st_a);     % ASR-11 transmit gain toward the NEXRAD
    gN = antenna_gain_rev1(pn.rx, 0, e, st_n);       % NEXRAD receive gain toward the ASR-11
    g21 = gA + gN;
    h_counts = h_counts + histcounts(g21, h_edges)';
    offA = wrapTo180(180 - st_a.mech_az);            % ASR-11 beam azimuth off the path
    offN = wrapTo180(0 - st_n.mech_az);              % NEXRAD beam azimuth off the path
    if j == 1 && opts.keep_ts1
        S.ts1 = struct("t_s", t, "g21_dBi", g21, "asr11_gain_dBi", gA, "nexrad_gain_dBi", gN, ...
            "asr11_az_off_deg", offA, "nexrad_az_off_deg", offN, "nexrad_el_deg", el_n(cut));
    end

    % ---- per ASR-11 revolution (rev, cut key and vol are nondecreasing in
    %      t, so each scan is one contiguous run of samples)
    [first, last] = runs(rev);
    [m21, k21] = run_max(g21, first, last);
    m12 = run_max(g12, first, last);
    n_r = numel(first);
    rev_tab{j} = table(j + zeros(n_r, 1), (1:n_r)', t(first), r_rev(rev(first)), m12, m21, ...
        el_n(cut(k21)), gN(k21), gA(k21), offN(k21), offA(k21), t(k21), edge_flag(n_r), ...
        'VariableNames', ["trial", "rev", "t_start_s", "rate_deg_s", "max_g12_dBi", ...
        "max_g21_dBi", "nexrad_el_at_max_deg", "nexrad_gain_at_max_dBi", "asr11_gain_at_max_dBi", ...
        "nexrad_az_off_at_max_deg", "asr11_az_off_at_max_deg", "t_at_max_s", "complete"]);

    % ---- per NEXRAD cut
    [first, last] = runs((vol - vol(1))*numel(dur_n) + cut);
    c21 = run_max(g21, first, last);
    c12 = run_max(g12, first, last);
    n_c = numel(first);
    cut_tab{j} = table(j + zeros(n_c, 1), vol(first) - vol(1) + 1, cut(first), el_n(cut(first)), ...
        t(first), c12, c21, edge_flag(n_c), 'VariableNames', ["trial", "volume", "cut", ...
        "el_deg", "t_start_s", "max_g12_dBi", "max_g21_dBi", "complete"]);

    % ---- per NEXRAD volume
    [first, last] = runs(vol);
    [v21, kv] = run_max(g21, first, last);
    v12 = run_max(g12, first, last);
    n_v = numel(first);
    vol_tab{j} = table(j + zeros(n_v, 1), (1:n_v)', t(first), v12, v21, el_n(cut(kv)), ...
        edge_flag(n_v), 'VariableNames', ["trial", "volume", "t_start_s", "max_g12_dBi", ...
        "max_g21_dBi", "el_at_max_deg", "complete"]);
end
S.asr_rev = vertcat(rev_tab{:});
S.nexrad_cut = vertcat(cut_tab{:});
S.nexrad_volume = vertcat(vol_tab{:});
S.t_rev_s = t_rev;
S.t_volume_s = t_vol;

% ---- summaries over complete scans
C = S.nexrad_cut(S.nexrad_cut.complete, :);
vn = ["el_deg", "n_cuts", "median_max_g21_dBi", "max_max_g21_dBi"];
if isempty(C)
    % Run shorter than one complete cut: no per-elevation statistics.
    S.by_elevation = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), 'VariableNames', vn);
else
    [G, el] = findgroups(C.el_deg);
    S.by_elevation = table(el, splitapply(@numel, C.max_g21_dBi, G), ...
        splitapply(@median, C.max_g21_dBi, G), splitapply(@max, C.max_g21_dBi, G), ...
        'VariableNames', vn);
end
R = S.asr_rev(S.asr_rev.complete, :);
S.asr_rev_pct = table(opts.pct(:), prctile(R.max_g21_dBi, opts.pct(:)), ...
    'VariableNames', ["pct", "max_g21_dBi"]);

F = cumsum(h_counts)/sum(h_counts);
q = zeros(numel(opts.pct), 1);
for i = 1:numel(opts.pct)
    q(i) = h_edges(find(F >= opts.pct(i)/100, 1) + 1);   % upper edge: within 0.01 dB
end
ts_pct = table(opts.pct(:), q, 'VariableNames', ["pct", "gain_dBi"]);
S.plateaus = plateau_check(pn.rx, pa.tx, el_n, t_vol, t_rev, S.asr_rev, opts);
end

% -------------------------------------------------------------------------
function P = plateau_check(pn, pa, el_n, t_vol, t_rev, Rv, opts)
% Hand calculation of the per-revolution plateaus. Each NEXRAD cut makes
% one revolution, so its main beam crosses the ASR-11 bearing once per cut.
% During an ASR-11 revolution that contains that crossing, the MAG reaches
% at least G_NEXRAD(on path, cut tilt) + G_ASR-11(floor). The share of
% ASR-11 revolutions that contain a crossing from cuts at tilt e or lower
% is about (number of those cuts per volume) x T_rev / T_volume. The
% simulated share counts revolutions within 0.2 dB of each level. Below
% about 30 dBi, the ASR-11 main beam on NEXRAD sidelobes also reaches the
% level, so the hand count is a lower bound there.
z = zeros(1, 1);
g_floor = antenna_gain_rev1(pa, 180, opts.el_path_deg, ...
    struct("mech_az", z, "tilt", z, "steer_az", z, "steer_el", z));   % ASR-11 beam pointing away
el = unique(el_n);
n_el = numel(el);
[g_n, n_cut, lev, pred, sim] = deal(zeros(n_el, 1));
Rc = Rv(Rv.complete, :);
for i = 1:n_el
    st = struct("mech_az", 0, "tilt", el(i), "steer_az", 0, "steer_el", 0);
    g_n(i) = antenna_gain_rev1(pn, 0, opts.el_path_deg, st);
    lev(i) = g_n(i) + g_floor;
    n_cut(i) = sum(el_n == el(i));
end
for i = 1:n_el
    k = lev >= lev(i) - 1e-9;                 % cuts whose crossing reaches this level
    pred(i) = 100*min(1, sum(n_cut(k))*t_rev/t_vol);
    sim(i) = 100*mean(Rc.max_g21_dBi >= lev(i) - 0.2);
end
% Keep only levels above ASR-11 main beam + NEXRAD back lobe, which every
% revolution reaches; at or below it the hand count does not apply.
z1 = struct("mech_az", 180, "tilt", el_n(1), "steer_az", 0, "steer_el", 0);
g_base = antenna_gain_rev1(pa, 180, opts.el_path_deg, struct("mech_az", 180, "tilt", 0, ...
    "steer_az", 0, "steer_el", 0)) + antenna_gain_rev1(pn, 0, opts.el_path_deg, z1);
keep = lev > g_base + 1e-9;
[el, n_cut, g_n, lev, pred, sim] = deal(el(keep), n_cut(keep), g_n(keep), lev(keep), pred(keep), sim(keep));
n_el = numel(el);
P = table(el, n_cut, g_n, g_floor + zeros(n_el, 1), lev, pred, sim, 'VariableNames', ...
    ["el_deg", "cuts_per_volume", "nexrad_gain_on_path_dBi", "asr11_floor_dBi", "level_dBi", ...
    "predicted_rev_pct", "simulated_rev_pct"]);
P.Properties.Description = sprintf("T_rev = %.3f s, T_volume = %.2f s", t_rev, t_vol);
end

% -------------------------------------------------------------------------
function [rate, el] = nexrad_cuts(scan)
if isfield(scan, "rate_rpm") && ~isempty(scan.rate_rpm)
    rate = 6*scan.rate_rpm;
else
    rate = scan.rate_deg_s(:);
end
el = scan.tilt_deg(:) + zeros(size(rate));
end

function st = state(az, tilt)
z = zeros(numel(az), 1);
st = struct("mech_az", az(:), "tilt", tilt(:), "steer_az", z, "steer_el", z);
end

function [first, last] = runs(idx)
% Start and end samples of each run of equal values in idx.
b = find(diff(idx) ~= 0);
first = [1; b + 1];
last = [b; numel(idx)];
end

function [m, k] = run_max(g, first, last)
% Maximum of g over each run and the sample index where it occurs.
m = zeros(numel(first), 1);
k = zeros(numel(first), 1);
for i = 1:numel(first)
    [m(i), j] = max(g(first(i):last(i)));
    k(i) = first(i) + j - 1;
end
end

function c = edge_flag(n)
% First and last scans of a trial are cut by the trial boundaries.
c = true(n, 1);
c([1 n]) = false;
end

% =========================================================================
function opts = fill_opts(opts)
if ~isfield(opts, "seed") || isempty(opts.seed)
    error("radar_coupling:missingSeed", "opts.seed is required for reproducibility.");
end
d = struct("nexrad_scan", "vcp", "vcp", 212, "sails", 0, "mrle", 0, "avset_top_deg", Inf, ...
    "cd_rate_deg_s", [], "nexrad_pattern", "nssl_measured", "asr_pattern", "m1851_csc2", ...
    "asr_tilt_deg", 0, "asr_rate_jitter_deg_s", 0.5, "asr_jitter_mode", "per_trial", ...
    "el_path_deg", 0, "n_samples", 5e6, "duration_s", 3600, ...
    "n_trials", 2, "dt_s", 0.002, "pct", [50 90 99 99.9 99.99 99.999], "parallel_flag", 0, ...
    "plot_flag", false, "keep_ts1", false, "cdf_pct", unique([0:1:100, 99.5 99.9 99.95 99.99 99.995 99.999]));
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || (isempty(opts.(f{i})) && f{i} ~= "cd_rate_deg_s")
        opts.(f{i}) = d.(f{i});
    end
end
opts.nexrad_scan = string(opts.nexrad_scan);
opts.nexrad_pattern = string(opts.nexrad_pattern);
opts.asr_pattern = string(opts.asr_pattern);
if ~ismember(opts.nexrad_scan, ["vcp", "constant"])
    error("radar_coupling:badNexradScan", "opts.nexrad_scan must be ""vcp"" or ""constant"".");
end
if ~ismember(opts.nexrad_pattern, ["nssl_measured", "roc_spec", "f1245"])
    error("radar_coupling:badNexradPattern", "opts.nexrad_pattern must be ""nssl_measured"", ""roc_spec"" or ""f1245"".");
end
if ~ismember(opts.asr_pattern, ["m1851_csc2", "f1245"])
    error("radar_coupling:badAsrPattern", "opts.asr_pattern must be ""m1851_csc2"" or ""f1245"".");
end
opts.asr_jitter_mode = string(opts.asr_jitter_mode);
if ~ismember(opts.asr_jitter_mode, ["per_trial", "per_rev"])
    error("radar_coupling:badJitterMode", "opts.asr_jitter_mode must be ""per_trial"" or ""per_rev"".");
end
if opts.asr_rate_jitter_deg_s < 0 || opts.asr_rate_jitter_deg_s >= 75
    error("radar_coupling:badJitter", "opts.asr_rate_jitter_deg_s must be in [0, 75) deg/s.");
end
if opts.dt_s <= 0 || opts.dt_s > 0.01
    error("radar_coupling:badDt", "opts.dt_s must be in (0, 0.01] s to resolve the 1.45 deg ASR-11 beam at 75 deg/s.");
end
if opts.duration_s/opts.dt_s*opts.n_trials > 2e8
    error("radar_coupling:tooManySteps", "duration_s/dt_s*n_trials exceeds 2e8 samples.");
end
end
