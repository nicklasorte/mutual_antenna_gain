# mutual_antenna_gain

MATLAB model of the mutual antenna gain (MAG, antenna coupling) between scanning radars in 2.7 GHz to 3.1 GHz, for the NTIA 2.7 GHz Retuning LOE. It follows the Rec. ITU-R M.2069-0 Annex 1 method and includes a worked example between a NOAA WSR-88D (NEXRAD) weather radar and an FAA ASR-11 airport surveillance radar, with the FAA ASR-9 as a second ATC radar.

**DRAFT//PRE-DECISIONAL.** The methodology report is `docs/mutual_antenna_gain_methodology_report_rev5.md` (and `.docx`).

## Requirements

MATLAB R2024b or later, with the Signal Processing Toolbox (`taylorwin`, for array tapers). The Parallel Computing Toolbox is optional (`opts.parallel_flag`). This is the product list that `matlab.codetools.requiredFilesAndProducts` reports for the repository in R2025b.

## Run

| To | Run |
| --- | --- |
| Get the MAG CDF for two radars | `radar_mutual_gain_cdf_rev1(sys1, sys2, opts)` |
| Run the WSR-88D and ASR-11 model for one VCP mode | `nexrad_asr11_mag_rev1(struct("seed", 1, "vcp", 212, "sails", 3))` |
| Run the same model with an exact CDF (no sampling noise) and a takeoff angle at each end | `nexrad_asr11_mag_rev1(struct("seed", 1, "cdf_method", "exact", "el_nexrad_deg", 1.7, "el_atc_deg", 2.0))` |
| Get the exact random-mode MAG distribution of any two rotating radars, including sectors | `mag_exact_cdf_rev1(sys1, sys2, opts)` |
| Aim a NEXRAD sector off the ATC radar, sweep it back and forth, or give the ATC radar a sector | `nexrad_asr11_mag_rev1(struct("seed", 1, "nexrad_sector_deg", 30, "nexrad_sector_center_deg", 10, "nexrad_sweep_mode", "bidirectional", "nexrad_turnaround_s", 0.5))`; `"atc_sector_deg"`, `"atc_sector_center_deg"` |
| Get a sector radar's pair CDF at its bearing to the other radar | `mag_pair_cdf_rev1(T_tx, T_rx_sector, el_tx, el_rx, struct("brg_rx_deg", b))` (sector tables from `gain_dist_table_rev1`) |
| Run the time series on event windows only (exact maxima above the floor), or in parfor | `nexrad_asr11_mag_rev1(struct("seed", 1, "ts_method", "events", "parallel_flag", 1))` |
| Build a time-series event library once per pair of scan modes, then get any pair's per-scan maxima, persistence and corrupted-radial counts | `lib = scan_event_library_rev1(out.scan_max.clock_spec, struct("seed", 1, "n_trials", 120))`; `mag_pair_timeseries_rev1(lib, ant_n, ant_a, el_n_deg, el_a_deg)` |
| Run the same model with the ASR-9 | `nexrad_asr11_mag_rev1(struct("seed", 1, "vcp", 212, "sails", 3, "atc_radar", "asr9"))` |
| Sweep every VCP mode for the ASR-9 and compare it with the ASR-11 | `example_asr9_rev1` |
| Sweep the 74 VCP modes of the Rev 5 study | `example_nexrad_asr11_rev1` |
| Run the sensitivity sweep | `example_sensitivity_rev1` |
| Find how long equal or multiple scan rates stay aligned | `example_asr11_drift_rev1` |
| Make every figure and table in the report | `report_figures_rev5` (writes `docs/figures`) |
| Build a gain table once per radar (all 216 WSR-88D modes: every VCP option and AVSET termination; `nexrad_mode_scans_rev1("lowest")` gives the 74 Rev 5 modes) | `[s, n] = nexrad_mode_scans_rev1(); T = gain_dist_table_rev1(ant, s, struct("mode_names", n))` |
| Get the MAG CDF of many pairs from their two takeoff angles, worst-mode envelope | `R = mag_pair_cdf_rev1(T_tx, T_rx, el_tx_deg, el_rx_deg)` (about 1.1 ms per pair for thousands of pairs) |
| Make a large pair study faster by sharing slices between nearby angles (each 0.01 deg cell gets its upper envelope, so results are never below the unrounded ones) | `mag_pair_cdf_rev1(T_tx, T_rx, el_tx, el_rx, struct("el_round_deg", 0.01))` |
| Build both tables and run a worked pair and 1,000 random pairs | `example_mag_lookup_rev1` |
| Check the two-step method against random-mode Monte Carlo | `validate_mag_lookup_rev1` (writes `docs/figures/mag_lookup_validation_rev1.csv`) |
| Get one antenna's exact gain distribution at a tilt and takeoff angle, full rotation or a sector | `h = gain_slice_rev1(gain_slice_prep_rev1(ant), tilt_deg, el_deg, az_lim_deg)` |
| Check that the tests catch planted bugs (red team F2) | `batch(parcluster("Processes"), @redteam_fault_seed_rev1, 1, {pwd}, "AutoAddClientPath", false)` |
| Run the red-team attacks and rebuild the findings register (`docs/redteam_register_rev1.csv`) | `batch(parcluster("Processes"), @redteam_runs_rev1, 1, {"all"}, "AutoAddClientPath", false)`, then `redteam_runs_rev1("register")` (about an hour; split the attacks across several batch jobs) |

Every model input is in `nexrad_asr11_mag_rev1` (`assumptions` and `citations`) or `nexrad_vcp_rev1`, with a source key and URL. Values without a public source are marked TODO.

## Test

```matlab
runtests(["test_radar_mutual_gain_cdf_rev1", "test_nexrad_asr11_rev1", ...
    "test_scripts_rev1", "test_validation_rev1"])
```

Run `test_timeseries_rev1` for the shared scan clock, event windows and pair time series, `test_sectors_rev1` for sector scans in both models, `test_exact_cdf_rev1` for the exact CDF in the original model, and `test_mag_lookup_rev1` and `test_gain_slice_rev1` for the two-step lookup (gain tables, gain slices and pair CDFs), `test_redteam_rev1` for the red-team runner's code paths, and `test_regression_rev5` to check both models against pinned Rev 5 and two-part reference outputs (`docs/figures/vcp_sweep_rev5.csv`, `test_data/mag_lookup_reference_rev1.csv`).

The four test files hold 152 tests and cover 2,562 of 2,564 lines (MATLAB code coverage). `test_scripts_rev1` runs every example script, the C# crosscheck, and `report_figures_rev5` end to end with small runs, writing to a temporary folder.

## Files

| File | Purpose |
| --- | --- |
| `radar_mutual_gain_cdf_rev1.m` | Sampling engine: random mode, time series, synchronized rotation, sector scans |
| `antenna_prepare_rev1.m`, `antenna_gain_rev1.m`, `m1851_rect_field_rev1.m` | Antenna patterns (Statgain, F.1245, M.1851-2, measured, off-axis table, planar array) |
| `nexrad_vcp_rev1.m` | WSR-88D VCP cut tables decoded from Level II data, with SAILS, MRLE, and AVSET |
| `nexrad_asr11_mag_rev1.m` | WSR-88D and ASR-11 or ASR-9 model (`opts.atc_radar`): per-scan maxima, reference levels, pattern-region pairs, persistence |
| `vcp_mode_cases_rev1.m` | The WSR-88D modes (216 by default, every AVSET termination; 74 with "lowest") and the paper baseline, for the sweeps |
| `mag_mode_envelope_rev1.m`, `mag_mode_summary_rev1.m` | Bound on the MAG for any mix of scan modes |
| `mag_sensitivity_rev1.m` | One-at-a-time sensitivity |
| `alignment_dwell_rev1.m` | Alignment of equal, multiple, or volume-locked scan rates |
| `peak_label_rev1.m`, `merge_overrides_rev1.m` | Figure 5 labels; small-run overrides for the script tests |
| `example_*.m`, `report_figures_rev5.m` | Entry points and report figures |
| `scan_clock_rev1.m`, `scan_clock_state_rev1.m` | Shared WSR-88D / ATC scan clock (random phases, per-trial, per-revolution or random-walk rate drift) for every time series |
| `scan_event_windows_rev1.m`, `scan_event_library_rev1.m`, `ts_floor_rev1.m` | Time steps where a beam is near the path, the library of them for many trials, and the MAG bound outside them |
| `mag_pair_timeseries_rev1.m` | Per-scan maxima, persistence and corrupted 0.5 deg radials for one pair from an event library |
| `sector_sweep_rev1.m` | Shared sector sweep law (restart, or back and forth with a turnaround dwell) for both models' time series |
| `mag_exact_cdf_rev1.m` | Exact counterpart of the random mode: same outputs, no samples; sector scans supported |
| `gain_slice_prep_rev1.m`, `gain_slice_rev1.m` | Exact gain distribution of one rotating antenna at one tilt and elevation: closed form for patterns of the off-axis angle, shifted azimuth cut for separable apertures, azimuth grid otherwise; full rotation or sector |
| `make_mag_lookup_reference_rev1.m`, `redteam_fault_seed_rev1.m` | Pin the two-part reference outputs; plant known bugs and check the tests catch them |
| `redteam_runs_rev1.m` | Red-team attacks A to F and the NOAA case (plan section 7); evidence in `docs/figures/redteam`, register in `docs/redteam_register_rev1.csv` |
| `gain_dist_table_rev1.m`, `mag_pair_cdf_rev1.m`, `nexrad_mode_scans_rev1.m` | Two-step MAG: exact gain distribution of each antenna by takeoff angle and scan mode (run once), then any pair's CDF by convolution, without a Monte Carlo run per pair |
| `crosscheck_csharp_rev1.m` | Comparison with the C# reference code of the 2.7 GHz project |
| `test_*.m` | Tests |
