# mutual_antenna_gain

MATLAB model of the mutual antenna gain (MAG, antenna coupling) between scanning radars in 2.7 GHz to 3.1 GHz, for the NTIA 2.7 GHz Retuning LOE. It follows the Rec. ITU-R M.2069-0 Annex 1 method and includes a worked example between a NOAA WSR-88D (NEXRAD) weather radar and an FAA ASR-11 airport surveillance radar.

**DRAFT//PRE-DECISIONAL.** The methodology report is `docs/mutual_antenna_gain_methodology_report_rev4.md` (and `.docx`).

## Requirements

MATLAB R2024b or later, with the Signal Processing Toolbox (`taylorwin`, for array tapers). The Parallel Computing Toolbox is optional (`opts.parallel_flag`). This is the product list that `matlab.codetools.requiredFilesAndProducts` reports for the repository in R2025b.

## Run

| To | Run |
| --- | --- |
| Get the MAG CDF for two radars | `radar_mutual_gain_cdf_rev1(sys1, sys2, opts)` |
| Run the WSR-88D and ASR-11 model for one VCP mode | `nexrad_asr11_mag_rev1(struct("seed", 1, "vcp", 212, "sails", 3))` |
| Sweep every VCP mode and option combination (74 modes) | `example_nexrad_asr11_rev1` |
| Run the sensitivity sweep | `example_sensitivity_rev1` |
| Find how long equal or multiple scan rates stay aligned | `example_asr11_drift_rev1` |
| Make every figure and table in the report | `report_figures_rev4` (writes `docs/figures`, about 5 minutes) |

Every model input is in `nexrad_asr11_mag_rev1` (`assumptions` and `citations`) or `nexrad_vcp_rev1`, with a source key and URL. Values without a public source are marked TODO.

## Test

```matlab
runtests(["test_radar_mutual_gain_cdf_rev1", "test_nexrad_asr11_rev1", ...
    "test_scripts_rev1", "test_validation_rev1"])
```

The four test files hold 149 tests and cover 2,409 of 2,411 lines (MATLAB code coverage). `test_scripts_rev1` runs every example script, the C# crosscheck, and `report_figures_rev4` end to end with small runs, writing to a temporary folder.

## Files

| File | Purpose |
| --- | --- |
| `radar_mutual_gain_cdf_rev1.m` | Sampling engine: random mode, time series, synchronized rotation, sector scans |
| `antenna_prepare_rev1.m`, `antenna_gain_rev1.m`, `m1851_rect_field_rev1.m` | Antenna patterns (Statgain, F.1245, M.1851-2, measured, off-axis table, planar array) |
| `nexrad_vcp_rev1.m` | WSR-88D VCP cut tables decoded from Level II data, with SAILS, MRLE, and AVSET |
| `nexrad_asr11_mag_rev1.m` | WSR-88D and ASR-11 model: per-scan maxima, reference levels, pattern-region pairs, persistence |
| `mag_mode_envelope_rev1.m`, `mag_mode_summary_rev1.m` | Bound on the MAG for any mix of scan modes |
| `mag_sensitivity_rev1.m` | One-at-a-time sensitivity |
| `alignment_dwell_rev1.m` | Alignment of equal, multiple, or volume-locked scan rates |
| `peak_label_rev1.m`, `merge_overrides_rev1.m` | Figure 5 labels; small-run overrides for the script tests |
| `example_*.m`, `report_figures_rev4.m` | Entry points and report figures |
| `crosscheck_csharp_rev1.m` | Comparison with the C# reference code of the 2.7 GHz project |
| `test_*.m` | Tests |
