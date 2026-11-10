# Plan: optimize both MAG models, add time series and sector scans to both, and red-team both

**DRAFT//PRE-DECISIONAL**

Revision 1, 2026-11-10. Repository: `nicklasorte/mutual_antenna_gain`, main at `8c980b7`.

This plan covers the two mutual antenna gain (MAG) models in the repository:

- **Original model:** `radar_mutual_gain_cdf_rev1` (the engine) and `nexrad_asr11_mag_rev1` (the WSR-88D and ATC radar wrapper). It samples pointing states, either at random or along a clock.
- **Two-part model:** `gain_dist_table_rev1` (one gain table per radar) and `mag_pair_cdf_rev1` (any pair's CDF by convolution).

Nothing in this plan is built yet, except one prototype that is labeled as such. Every baseline number was measured on 2026-11-10 in MATLAB R2025b on the Apple silicon Mac. Every target is an estimate until its work package is done.

## Status, 2026-11-10 (steps 1 to 8)

| Package | Status | Result |
| --- | --- | --- |
| K3 Pin references | Done | `test_regression_rev5`: three `vcp_sweep_rev5.csv` rows reproduced to 1e-10 relative; two-part outputs pinned in `test_data/mag_lookup_reference_rev1.csv` (5 pairs, 74 modes, 6 percentiles) |
| K2 Gain-slice kernel | Done | `gain_slice_prep_rev1`, `gain_slice_rev1`: closed form, separable, grid; full rotation or sector (azimuth interval). Covers every antenna type |
| T1 Closed-form slice | Done | Within one bin of the azimuth grid for every pattern, geometry and sector tested. Each table cell is split until it spans at most one bin, so the bound stays within one bin |
| T2 Shift rule (separable) | Done in the kernel | Identical to the azimuth grid (same samples), binary search per slice |
| T3 Part 2 | Done (with step 7) | `opts.interp = "exact"` added (33 ms per pair, within 0.1 dB of interpolation). Interpolation and percentile search vectorized in step 3 (4.5 ms per pair against 6.4 ms). Step 7 meets the 2 ms target for large studies; see the T4 and T5 row |
| T4, T5 Batching (step 7) | Done | Pairs grouped by the takeoff angle (and bearing) of the side with more slice columns; that side's slices and FFT made once per group, the other side's once per distinct angle (T5). One FFT per slice column, one inverse FFT for all tilt pairs (T4, `conv_method = "fft"`, default; `"direct"` kept for checks). The interpolation's sort is made once per table interval and reused for every angle in it, and one sort replaces three. Optional `el_round_deg` rounds angles before grouping. Percentiles, per-mode percentiles, worst modes and maximum path gains equal the step 6 code exactly for 200 random pairs in both directions; envelope curves within 7e-14 percentage points. Pair time, WSR-88D and ASR-11 tables, random angles from -5 to 10 deg: 5,000 pairs, 1.12 ms against 2.88 ms (2.6 times faster; 0.84 ms with `el_round_deg = 0.01`); 200 pairs, 2.18 ms against 3.28 ms (few pairs share an interval, so the target holds only for larger studies or rounded angles: 1.90 ms at 0.01 deg). FFT against direct convolution at 5,000 pairs: 1.12 ms against 1.47 ms |
| K1 Shared scan clock | Done (step 6) | `scan_clock_rev1`, `scan_clock_state_rev1`: one clock for the wrapper, the event library and the pair time series; the wrapper reproduces the Rev 5 rows bit for bit on it |
| O1 Exact CDF in the original model | Done | `mag_exact_cdf_rev1` and `nexrad_asr11_mag_rev1` option `cdf_method = "exact"`. CDF step 0.087 s against 1.09 s for 5 million samples (13 times faster); percentiles inside the Monte Carlo intervals or within the 0.1 dB round-up above them. Monte Carlo stays the default, so Rev 5 reproduces unchanged |
| O2 Mode and sensitivity sweeps | Done | `mag_mode_summary_rev1` takes a distribution; `example_overrides.cdf_method = "exact"` runs the 74-mode sweep exactly. Exact envelope matches `mode_envelope_rev5.csv` within the Rev 5 intervals |
| O5 Takeoff angle at each end | Done | `el_nexrad_deg`, `el_atc_deg` (default `el_path_deg`), used in the CDF, the time series, the reference levels, the plateau check and the lobe pairs |
| Exceedance near plateaus | Finding | Rounding gains up moves mass that sits just below a level to above it. At the 45.49 dBi plateau the exact time share is 0.027% against 0.007% by Monte Carlo: an upper bound equal to the share at the level minus at most 0.1 dB. Where a time share at a plateau matters, call `mag_exact_cdf_rev1` with a finer `bin_dB` (the wrapper uses 0.05 dB) |
| S1 to S3 Sector scans | Done | `sector_sweep_rev1` (restart, or bidirectional with a turnaround dwell) shared by the engine and the wrapper; wrapper options `nexrad_sector_center_deg`, `nexrad_sweep_mode`, `nexrad_turnaround_s`, `atc_sector_deg`, `atc_sector_center_deg`; sector tables and `brg_tx_deg`, `brg_rx_deg` in the two-part model; `gain_slice_rev1` edge dwell. Rev 5 "NEXRAD sector 10 deg" row reproduced bit for bit; engine, exact and two-part agree (`test_sectors_rev1`). A bearing sweep of 26 bearings takes about 0.1 s |
| E1 (red team, run early) | Finding | "Centered on the path is the worst case" holds for a restart sweep (99.99%: 69.95 dBi centered, 64.40 dBi at the edge) but fails for a bidirectional sweep with a 0.5 s turnaround: with the path at the sector edge the dwell parks the main beam on it, 75.40 against 68.05 dBi at 99.99% (+7.4 dB) and +7.6 dB at 99.9% (30 deg sector, VCP 212 SAILSx3 AVSET 6.4, `docs/figures/sector_bearing_sweep_rev1.csv`). The dwell is a bounding assumption (the beam radiates while parked); no public sweep law is known (TODO). Sector studies must sweep the bearing, not assume the center |
| Step 6 Time series | Done | Event windows (`scan_event_windows_rev1`: a beam within 4 deg of the path, a subset of the fixed grid, about 5% of it) with a per-pair floor (`ts_floor_rev1`, 15.4 dBi at 0/0). Library once per mode pair (`scan_event_library_rev1`, 0.29 s for 20 h); pair time series (`mag_pair_timeseries_rev1`) 18 ms per pair-hour against 312 ms for the fixed grid (17 times faster, 2.2 s per pair and mode for 120 h; the plan estimated 0.6 s). Reproduces `persistence_120h_rev5.csv` and `persistence_volume_120h_rev5.csv` exactly; maxima above the floor equal the fixed grid. Wrapper: `ts_method = "events"` (O4), parfor trials (O3, identical results), `asr_jitter_mode = "random_walk"`. Harm metric (D4): corrupted 0.5 deg radials per cut, checked against a brute-force count |
| Step 8 Red team | Done, with three items open | `redteam_runs_rev1` ran attacks A1 to F1 and the NOAA case; findings register `docs/redteam_register_rev1.csv` (33 rows), evidence in `docs/figures/redteam`. Section 7.4 has the results. The runs found and fixed three defects (X1 to X3). Open: F3 and parts of C3 and C5 need an RFS workbook (decision 5); layer 2 review done (next row) |
| Layer 2 independent review | Done; three fixes | A separate agent broke four of seven claims (register rows R1 to R6). Fixed: R1, bidirectional sector time series left part of the sector unvisited (per-volume maxima 13 to 39 dB low; CDFs unaffected); R3, AVSET terminations above 6.4 deg were missing from the envelope (up to 4.35 dB low at 99.9% for WSR-88D takeoff angles above about 5.5 deg; default mode set now 216 modes); R4, "nearest" and rounded angles could be up to 0.85 dB low (now "upper" and cell envelopes, never below quantile interpolation). Open: R5 (low). R6 fixed afterwards |
| D1 (red team, run early) | Finding | Rate-drift model and persistence, 120 h, worst mode: volumes reaching 70 dBi 14.2% per-trial offset (Rev 5 baseline, fully correlated), 13.9% random walk (0.02 deg/s per revolution), 10.6% independent per revolution; longest run 4, 3 and 2 volumes. The Rev 5 baseline is the conservative end; an independent per-revolution draw would understate persistence. The random-walk step has no public source (TODO). `docs/figures/drift_persistence_rev1.csv` |
| F2 Fault seeding | Done (first set) | `redteam_fault_seed_rev1`: 10 planted faults across both models. First run caught 9; the miss (time-series trials all starting at the same point of the scan cycle) got a new test, `timeSeriesTrialsStartAnywhereInScanCycle`. Rerun: 10 of 10 caught (`docs/figures/fault_seed_rev1.csv`). Step 7: faults 5, 6 and 8 pointed at code that steps 5 and 7 had rewritten (the harness stops on a missing fault text, so none was silently skipped); they now plant the same bugs in the new code, and the rerun catches 10 of 10 |

Measured on 2026-11-10: WSR-88D table (74 modes, 21 tilts, 301 elevations) 2.0 s, down from 127 s; ASR-11 table 0.4 s, down from 15 s. Worked pair (2.0°/1.7°): every percentile unchanged except 99.9%, up one bin (43.55 to 43.60 dBi), as the conservative closed form allows.

## 1. Summary

1. **Neither model is slow per sample.** The original model costs 132 ns per time step and 189 ns per random sample. The cost comes from how many samples a study needs: thousands of pairs, 74 scan modes, and 120 hours of clock time.
2. **The largest saving is to stop sampling where an exact answer exists.** Use the exact convolution for every CDF, in both models. Keep the clock simulation only for what needs a clock: persistence, volume locks, and per-volume statistics.
3. **The time series can be split in two, the same way the CDF was.** The scan clocks don't depend on the takeoff angles. Build a library of beam-passage events once per scan-mode pair, then apply each pair's gains to it.
4. **Sector scans fit the two-part model if the gain slice is computed in closed form.** A sector makes the gain distribution depend on the bearing, so a fixed table no longer works. A closed-form slice does, and the prototype shows it is accurate and fast.
5. **The red team is planned as a register of attacks.** Each attack names the claim it tests, the test, and a pass or fail criterion.

## 2. Where the two models stand

| Capability | Original model | Two-part model |
| --- | --- | --- |
| Random-mode CDF | Yes, Monte Carlo (5 million samples by default) | Yes, exact |
| Separate takeoff angle at each end | Engine: yes (`el12_deg`, `el21_deg`). Wrapper: no, one `el_path_deg` for both ends | Yes |
| All 74 WSR-88D modes | One mode per run; sweep by `example_nexrad_asr11_rev1` | Yes, envelope in one call |
| Time series (per-revolution and per-volume maxima, persistence) | Yes | **No** |
| Sector scans | Engine: yes (`sector_deg`, `sector_center_deg`). Wrapper: WSR-88D only, centered on the ATC radar | **No** (raises an error) |
| Lobe-pair breakdown, plateau check | Yes | No |
| Parallel option | Random mode: yes. Wrapper time series: no | Yes, both parts |

The original model already has a time series and sector scans. For the original model this plan therefore covers speed, the wrapper gaps, and realism. For the two-part model it covers adding both features.

## 3. Measured baselines

| Item | Measurement |
| --- | --- |
| Original, default run (VCP 212 SAILSx3 AVSET 6.4, 5 million samples, two 1-hour trials) | 1.6 s warm |
| Original, random mode | 189 ns per sample |
| Original, time series | 132 ns per step; 0.24 s per 1-hour trial (1.8 million steps at `dt_s` = 0.002) |
| Original, 120-hour time series, one pair, one mode | About 29 s (120 trials at 0.24 s; computed from the line above) |
| Original, fixed overhead per call | 0.15 s |
| Two-part, WSR-88D table (21 tilts, 301 elevations, 360,000 azimuths each) | 127 s; 29 MB |
| Two-part, ASR-11 table | 15 s; 0.7 MB |
| Two-part, pair CDF with the 74-mode envelope | 5.3 ms per pair (1,000 pairs in 5.3 s) |
| Two-part, share of pair time spent in elevation interpolation (`quantile_interp`) | About 60% (profile of 10 pairs) |
| Monte Carlo reference, 20 million samples, one pair, one mode | 12.7 s |

**Scaling problem.** A 120-hour time series for 1,000 pairs and 74 modes is 74,000 runs at 29 s, about 25 days on one core. That number, not the per-sample cost, is what the optimization must remove.

**Hardware limit.** MATLAB `gpuArray` needs an NVIDIA GPU, so it is not an option on the Apple silicon Mac. The plan uses vectorization and `parfor` only.

## 4. Optimization plan

Each work package lists the change, the target, and the acceptance test. Targets marked "estimate" are computed from the baselines, not measured.

### 4.1 Shared kernels (do first)

| ID | Change | Why | Acceptance test |
| --- | --- | --- | --- |
| K1 | Move the scan clock (volume, cut, azimuth against time, with jitter) into one function, `scan_kinematics_rev1`, used by the engine, the wrapper, and the new event library | Today the engine (`state_time_series`) and the wrapper (`scan_maxima`) each have their own clock code. Three copies would drift apart | Same seeds give sample-identical states before and after, in both callers |
| K2 | One gain-slice function, `gain_slice_rev1(ant, tilt, el, az_interval)`, that returns the gain distribution of one antenna at one tilt and one elevation | Parts 1 and 2 and the sector work all need it | Matches the current azimuth-grid histogram within one bin for every antenna type in the repository |
| K3 | Pin the Rev 5 outputs as regression references before any change | Every later package is judged against them | `docs/figures/*_rev5.csv` reproduced exactly by the unchanged code |

### 4.2 Original model

| ID | Change | Target | Acceptance test |
| --- | --- | --- | --- |
| O1 | Add `opts.cdf_method = "exact"` to the wrapper: take the CDF from the convolution instead of 5 million random draws | Removes the random-mode draws and the tail noise (the 99.999% value has a 95% interval of 66.8 to 71.1 dBi today) | Exact percentiles inside the Monte Carlo interval at every level; Rev 5 summary reproduced within the stated bias |
| O2 | Mode sweep and sensitivity sweep call the exact CDF | 74 modes cost one table build at one elevation plus 74 convolutions, in place of 74 Monte Carlo runs | `mode_envelope_rev5.csv` reproduced within the Monte Carlo intervals |
| O3 | `parfor` over trials in `scan_maxima`, with the histogram and lobe counts as reductions | Near-linear in workers for long runs | Parallel result equals serial result for the same seed |
| O4 | Event-based time series (section 5) in place of fixed 2 ms steps | About 47 times fewer gain evaluations (estimate: about 38,000 per simulated hour, down from 1.8 million) | Section 5 acceptance tests |
| O5 | Separate `el_tx_deg` and `el_rx_deg` in the wrapper | Closes the red-team gap from 2026-11-07: the wrapper can't run the asymmetric geometry | Equals the engine with `el12_deg` and `el21_deg` set |

Not planned: a faster per-sample pattern evaluation. At 132 ns per step, with `antenna_gain_rev1` already vectorized, the gain there is small next to O1 and O4.

### 4.3 Two-part model

| ID | Change | Target | Acceptance test |
| --- | --- | --- | --- |
| T1 | Closed-form slice for rotationally symmetric patterns (WSR-88D). For a uniform azimuth, the fraction of time with off-axis angle below γ is `acos((cos γ − sin(el)·sin(tilt)) / (cos(el)·cos(tilt))) / π`. Push that through the one-dimensional pattern table; no azimuth grid | **Prototype measured:** 0.65 ms per slice against 5.2 ms (8 times faster); percentiles within one 0.05 dB bin, on the high side, in six test geometries. WSR-88D table estimate: about 4 s, down from 127 s | All 6,321 slices of the current table reproduced within one bin; test suite passes |
| T2 | Shift rule for separable patterns (ASR-9, ASR-11). The gain is `g_max + max(g_az(az) + g_el(el, tilt), floor)`, so the distribution at any elevation is one stored azimuth histogram shifted by `g_el` and clamped at the floor | The ATC table becomes one histogram plus a function. No elevation grid, no interpolation error, any site tilt for free | Equals the current table at every grid elevation |
| T3 | Remove elevation interpolation where T1 and T2 make an exact slice cheap; otherwise precompute the cumulative sums that `quantile_interp` rebuilds on every call | Pair time under 2 ms (estimate; interpolation is about 60% of 5.3 ms today) | Results unchanged within one bin |
| T4 | Batch the convolutions: one FFT per slice, products for all tilt pairs, one inverse FFT; vectorize the percentile search over modes | Pair time under 1 ms with T3 (estimate) | Results equal `conv` within 1e-12 in probability; negative round-off clipped and reported |
| T5 | Group pairs that share takeoff angles (rounded to the grid) and cache slices | Large studies reuse most slices; the saving depends on the study | Same results as the ungrouped run |

**Decision for T1 to T3.** With exact slices this cheap, the fixed −5° to 10° grid becomes a cache, not a requirement. I recommend keeping the table file as the cache and adding exact evaluation for any elevation outside it or between nodes.

## 5. Time series in both models

### 5.1 What the time series must deliver

The CDF comes from the exact convolution. The time series is needed only for clock-based outputs:

- maximum MAG in each ATC revolution, each WSR-88D cut, and each volume;
- persistence (runs of consecutive revolutions or volumes above a threshold);
- volume locks and how long they last;
- a harm metric that counts corrupted radials in each volume and cut (red-team gap 1 from 2026-11-07).

### 5.2 Design: an event library, then gains for each pair

The scan clocks don't depend on the takeoff angles. Only the gains do. That gives the same two-part split as the CDF.

**Part A, once per scan-mode pair: `scan_event_library_rev1`.** Run the clocks (K1) for N trials and record two kinds of events:

1. **ATC passages:** each time the ATC beam crosses the bearing to the WSR-88D (once per revolution, 750 per hour at 12.5 rpm). Store the time, the WSR-88D azimuth offset from the path, its rate, the tilt, and the cut and volume numbers.
2. **WSR-88D passages:** each time the WSR-88D beam crosses the bearing to the ATC radar (once per cut). Store the time, the ATC azimuth offset, the tilt, and the cut and volume numbers.

The library holds a few numbers per event, so 120 hours is a few megabytes.

**Part B, for each pair: `mag_pair_timeseries_rev1`.** For each event, evaluate both gains on a short window of times around it (for example 41 points), at the pair's two takeoff angles, as one matrix of events by window points. The maximum in each revolution is the largest of its event windows.

**Censoring floor.** Outside the windows, MAG can't exceed the largest sidelobe-plus-sidelobe sum. The function computes that bound and reports any maximum below it as "at or below the floor." It raises an error if a requested persistence threshold is below the floor.

**Cost estimate.** About 38,000 gain evaluations per simulated hour, against 1.8 million: about 0.6 s per pair and mode for 120 hours at today's 132 ns per evaluation. The estimate needs a measurement before it goes in any report.

**Mode screening.** Part 2 of the CDF model already names the worst mode at each percentile. Run the time series for each pair's worst few modes, not all 74.

### 5.3 Original model

The original time series stays as the reference implementation. Changes: K1, O3, O4 as an option (`opts.ts_method = "events"`), O5, and two realism options the red team asked for:

- **Correlated rate drift:** a random-walk option beside the current independent draws per trial or per revolution.
- **Radial-count harm metric:** corrupted radials in each volume at each cut, above a threshold the user sets.

### 5.4 Acceptance tests

1. Same seeds, same clocks: per-revolution maxima from the event method match the fixed-step method within 0.1 dB for at least 99.9% of revolutions above the floor. The fixed-step method samples every 2 ms, so the event method should be equal or slightly higher.
2. `persistence_120h_rev5.csv` and `persistence_volume_120h_rev5.csv` reproduced: shares within the binomial interval, longest runs equal.
3. The time-percentile CDF built from a long event-method run agrees with the exact convolution (checks that the two parts describe the same system).
4. Parallel equals serial.

## 6. Sector scans in both models

### 6.1 What changes with a sector

With full rotation, the bearing between the radars doesn't matter. With a sector, it does: the gain distribution depends on where the path sits in the sector. New input for each sector-scanning radar: the bearing of the other station relative to the sector center.

### 6.2 Original model

The engine already supports sectors. Gaps to close:

| ID | Change |
| --- | --- |
| S1 | Wrapper: `opts.nexrad_sector_center_deg` (bearing of the sector center from the path) in place of the fixed "centered on the ATC radar" case, which stays as the default worst case |
| S2 | Wrapper: optional sector for the ATC radar, through the engine's existing fields |
| S3 | Engine: sweep law option. Today the sweep restarts from the leading edge (`mod` in `state_time_series`). Add a back-and-forth sweep with a turnaround time. TODO: confirm the sweep law of any real sector-scanning radar in scope from a public source before using it |

### 6.3 Two-part model

A fixed table would need a bearing dimension at azimuth resolution, which is thousands of times the storage. Use the closed-form slice (K2, T1) on demand:

- **Symmetric patterns:** restrict the closed-form azimuth fraction to the sector's azimuth interval.
- **Separable patterns:** histogram the stored azimuth cut over the sector's interval, then apply the T2 shift.
- **Non-uniform sweep:** weight azimuth by time spent (one over the sweep rate) when S3's turnaround is on.

`mag_pair_cdf_rev1` takes two more optional inputs per pair, `brg_tx_deg` and `brg_rx_deg`. When the sector center is unknown, the default is the worst case: sector centered on the path.

The event library (section 5) records passages only when the path is inside the sector.

### 6.4 Acceptance tests

1. A 360° sector equals the full-rotation result exactly.
2. Engine random mode with a sector (existing code) matches the two-part result inside the Monte Carlo interval, for the path at the sector center, at the edge, and outside.
3. The Rev 5 sensitivity row "sector 10°" (+13.9 dB at 99.9%, +28.4 dB at 99.99%) is reproduced.

## 7. Red-team plan for both models

### 7.1 Method

- **Stance:** a pessimistic NOAA radar engineer for harm to NEXRAD, and a skeptical numerical analyst for the code.
- **Rule:** every attack names a claim, a test, and a pass or fail criterion decided before the test runs.
- **Independence:** the model can't grade itself. Three layers: (1) scripted tests in the repository; (2) a separate review agent that hasn't seen this plan or the code history, which needs your go-ahead; (3) human review by NOAA and FAA subject-matter experts before any result leaves NTIA.
- **Sources:** public and cleared material only. Where a test needs restricted data, the register names the data and stops.
- **Output:** a findings register (`docs/redteam_register_rev1.csv`) with one row per attack: result, evidence file, severity, and fix. A "NOAA case" run sits next to the baseline in the Rev 6 report.

### 7.2 Attack register

**A. Is the convolution valid? (two-part model, and random mode in the original)**

| ID | Claim attacked | Test | Pass criterion |
| --- | --- | --- | --- |
| A1 | Pointing states are independent | For all 74 modes, compare the exact CDF with the time-percentile CDF of a long clock run | Agreement inside the run's confidence interval for every mode; list every mode that fails |
| A2 | The ATC rate tolerance doesn't matter | Search the tolerance band on a fine grid of rates, for each mode, for ratios between the ATC revolution time and the WSR-88D cut and volume times that repeat | No rate in the band where clock and exact CDFs differ by more than 1 dB at 99.9%; otherwise report the rates and durations |
| A3 | The envelope bounds any mix of modes | Check the proof for mixtures; test mode changes in the middle of a volume | Mixture bound holds; transitions add no exceedance above the envelope |
| A4 | The 74 modes are the complete set | Compare against public NWS and ROC VCP documentation, including site-specific lowest elevation angles below 0.5°. TODO: confirm from a public source | Every operational VCP and lowest angle is in the table, or the gap is listed with its effect |

**B. Numerical accuracy (both models)**

| ID | Claim attacked | Test | Pass criterion |
| --- | --- | --- | --- |
| B1 | Percentiles are conservative by at most `bias_dB` | Rebuild with 0.005 dB bins and a 0.0002° azimuth step | Reported value is between the refined value and the refined value plus the bias, everywhere |
| B2 | Interpolation error is at most 0.1 dB | Search all elevations in the grid at 0.005° spacing, all tilts, not only 24 midpoints | Maximum error stated and inside the bias; otherwise tighten the grid or use exact slices (T3) |
| B3 | Single-precision storage doesn't hurt the tail | Compare with a double-precision table at 99.999% and 99.9999% | Same bin |
| B4 | The closed-form slice (T1) is exact | Compare with the azimuth sweep for every slice, and at the edge cases: elevation equal to tilt, negative elevation, tilt 19.5° | Within one bin, never low |
| B5 | Monte Carlo intervals in the original are honest | Repeat the baseline with 20 seeds | The stated 95% interval covers about 95% of seeds |

**C. Inputs (both models; carried over from the 2026-11-07 red team)**

| ID | Claim attacked | Test | Pass criterion |
| --- | --- | --- | --- |
| C1 | The NSSL KOUN measured pattern represents the fleet | Build the ROC-specification table and the F.1245 table; rerun the pair set | Difference reported for every percentile; the "NOAA case" uses the specification |
| C2 | The ASR-9 remote sidelobe (+3.5 dBi, borrowed) is safe | Run +0.5, +3.5, and +6.5 dBi | One-for-one sensitivity confirmed; value stays TODO until a public source is found |
| C3 | Takeoff angles from the RFS tool mean what the model assumes | Check sign, reference (local horizontal), and earth-curvature handling in `CouplingCalculations.cs` and `SpreadsheetIO.cs` against `el_tx_deg` and `el_rx_deg`, on three hand-computed paths | Same angle within 0.01° |
| C4 | Site tilt doesn't matter | Run the ATC tilt from −1° to +3° | Effect stated; tilt becomes a per-site input (free with T2) |
| C5 | Far-field patterns apply | Compute 2D²/λ for each antenna and compare with the shortest pair distance in the study | Every pair beyond the far-field distance, or flagged |

**D. Time series (both models)**

| ID | Claim attacked | Test | Pass criterion |
| --- | --- | --- | --- |
| D1 | Independent rate jitter is realistic | Rerun the persistence tables with correlated drift | Lock probability and longest run reported for both; the larger one is used |
| D2 | 120 hours is enough | Bootstrap intervals on every persistence statistic | Intervals stated; extend the run if the lock probability interval is wider than a factor of two |
| D3 | The event method misses nothing | Section 5.4 tests, and a search for maxima outside every window | No maximum above the floor outside a window |
| D4 | Time percentiles describe harm | Compute corrupted radials in each volume at the 0.5° cut for the pair set | Harm metric reported beside the percentiles |

**E. Sector scans (both models)**

| ID | Claim attacked | Test | Pass criterion |
| --- | --- | --- | --- |
| E1 | "Centered on the path" is the worst case | Sweep the bearing across and beyond the sector | Maximum at the center, or the true worst bearing reported |
| E2 | Restarting sweep equals back-and-forth sweep | Compare the two sweep laws, with turnaround time | Difference stated; a path at the sector edge is tested |

**F. Software (both models)**

| ID | Claim attacked | Test | Pass criterion |
| --- | --- | --- | --- |
| F1 | The two models agree | Property test: random patterns, tilts, scan rates, and takeoff angles; the exact result must fall inside the Monte Carlo band | 500 random cases, no failure |
| F2 | The tests would catch a bug | Seed ten deliberate faults (wrong sign on an angle, swapped tx and rx, weight not normalized) | Every fault fails at least one test |
| F3 | The model agrees with the C# tool | Extend `crosscheck_csharp_rev1` to the tool's gain-distribution and MACD outputs for one pair with nonzero takeoff angles. Needs an RFS workbook with those columns | Same percentiles within the tool's bin width |

**G. Scope (stated limits, to be written into the report)**

No test can close these; the report must state them: no aggregate interference from several ATC radars, no polarization or dual-polarization product effects, no pulse-level coupling or frequency-dependent rejection, no validation against observed interference. ROC interference records are restricted and are only named here.

### 7.4 Results of the step 8 runs

All runs used public inputs, the worst mode (VCP 212 SAILSx3 AVSET 6.4) unless stated, and the 0/0 deg geometry unless stated. Severity is the red team's call, for NOAA and FAA reviewers to confirm.

| ID | Result | What it shows |
| --- | --- | --- |
| A1 | Pass | Exact and time-series percentiles agree for all 74 modes: 2 of 296 checks fall outside the 95% interval, fewer than the 15 expected by chance |
| A2 | Pass | No ATC rate in 74.5 to 75.5 deg/s moves the 99.9% value by more than 0.40 dB (limit 1 dB) |
| A3 | Finding (low) | The envelope bounds every mixture of whole volumes (exceedance is a weighted mean). A mode change in the middle of a volume repeats the low cuts: that one window can sit 1.35 dB above the envelope at 99% and 0.65 dB at 99.9% (0 at 99.99%) |
| A4 | Finding (medium) | The mode table starts at 0.5 deg; public NWS material lists sites with a supplemental cut at -0.2 to +0.4 deg. One such cut raises the worst mode's MAG by up to 3.1 dB |
| B1, B3, B4, B5 | Pass | Grid bias inside the stated bound; single precision costs nothing in the tail; closed form never low; Monte Carlo intervals cover 19 or 20 of 20 seeds |
| B2 | Fail (low) | Interpolation between grid angles: 9 of 32,400 checks over 0.1 dB, largest 0.25 dB; one bin low in 0.25%. Within the corrected bias; use `interp = "exact"` where 0.1 dB matters |
| C1 | Reported (medium) | The ROC specification envelope gives MAG 10 to 15 dB higher than the NSSL measured pattern at 50 to 99.9% (median), about 1 dB at 99.99%. The NOAA case uses the specification |
| C2 | Reported (medium) | The ASR-9 remote sidelobe moves the MAG one-for-one at 50, 90, 99.9 and 99.99% |
| C3 | Partial | Angle convention matches the C# tool (takeoff - tilt, local horizontal, positive up). New C# defect X4. Path-by-path check needs the RFS workbook |
| C4 | Finding (medium) | ATC tilt -1 to +2.5 deg (the model's range) moves the 99% value by 2.6 dB and 99.999% by 3.3 dB; 99.99% does not move |
| C5 | Partial | Far field: WSR-88D 1.3 to 1.5 km, ASR about 0.6 km. Pairs under 1.5 km need flagging once the pair list is available |
| D2 | Pass | 120 h is enough: volume share reaching 70 dBi 12.3% (95% interval 10.6 to 14.2%); no lock-probability interval wider than a factor of two |
| F1 | Pass after fixes | 500 random cases, 0 outside the Monte Carlo band. The first run found X2 and X3 |
| NOAA | Reported | ROC specification pattern, ASR remote sidelobe +6.5 dBi, worst tilt: 99.99% MAG 50.3 dBi against 45.5 baseline; volumes reaching 70 dBi unchanged at 12.3% (main beam to main beam does not depend on sidelobes) |
| X1 | Fixed | The wrapper crashed when the ATC gain toward the path is flat (ASR tilted 2.5 deg over a +6.5 dBi floor) |
| X2 | Fixed (high) | Closed-form slices crashed when a sector fell inside one segment of a coarse table (10 of 500 F1 cases) |
| X3 | Fixed | Closed-form slices round up twice, so the stated bias was half the true bound: now 0.2 dB for `mag_exact_cdf_rev1` and 0.3 dB for `mag_pair_cdf_rev1`. Values are unchanged |

### 7.3 Status of the 2026-11-07 findings

| Finding | Status |
| --- | --- |
| Tail percentiles carry Monte Carlo noise | Fixed in the two-part model; fixed in the original by O1 |
| Wrapper can't run asymmetric takeoff angles | Fixed in the two-part model; O5 fixes the wrapper |
| Headline metric hides per-volume harm | Harm metric built and tested (D4, step 6); to be reported per pair in Rev 6 |
| Percentile choice rests on meeting statements | Open; needs the source reports |
| Non-conservative pattern defaults | Quantified (C1, C2, C4); the NOAA case carries the conservative choices |
| Independent rate jitter | Closed: the per-trial baseline is the most persistent drift model (D1); 120 h is enough (D2) |
| One pair at a time (no aggregate) | Open; out of scope for this plan, listed in G |

## 8. Order of work

| Step | Packages | Depends on | Size |
| --- | --- | --- | --- |
| 1 | K3 (pin references), F2 scaffold | | Small |
| 2 | K1, K2 (shared kernels) | 1 | Medium |
| 3 | T1, T2, T3 (exact slices) | 2 | Medium |
| 4 | O1, O2, O5 (exact CDF in the original) | 3 | Small |
| 5 | Sector scans, two-part and wrapper (section 6) | 3 | Medium |
| 6 | Event library and pair time series (section 5), O3, O4 | 2 | Large |
| 7 | T4, T5 (batching), only if step 3 leaves the pair time above target (done) | 3 | Small |
| 8 | Red-team runs A to F, findings register, NOAA case (done) | 3 to 6 | Large |
| 9 | Rev 6 methodology report | 8 | Medium |

Steps 3 and 6 can run side by side after step 2.

## 9. Decisions needed

1. **Grid or exact slices?** Keep the fixed elevation grid as a cache and add exact slices (recommended), or keep the grid only.
2. **Sector scope.** Which radars in the study actually sector-scan? If none do today, sector work can move after the time series.
3. **Time-series modes.** Run each pair's worst few modes by the CDF envelope (recommended), or all 74.
4. **Independent review agent.** Done (approved 2026-11-10); see the status table.
5. **RFS workbook.** F3, C3 and C5 need a workbook with takeoff-angle and ground-distance columns. The project workbook is over the 20 MB in-session limit. Can a smaller, cleared extract be provided?
6. **Supplemental low cut (A4).** Add a per-site lowest-angle option to the WSR-88D scan model and run the affected sites with it? (Recommended.)
7. **Pattern for the report (C1).** Report the NSSL measured pattern as the baseline and the ROC specification as the NOAA case side by side (recommended), or pick one.

## 10. Numbers ledger

| Number | Value | Source |
| --- | --- | --- |
| Time-series cost | 132 ns per step | Profile, 2026-11-10 |
| Random-mode cost | 189 ns per sample | Profile, 2026-11-10 |
| WSR-88D table build | 127 s, 29 MB | `example_mag_lookup_rev1`, 2026-11-09 |
| ASR-11 table build | 15 s | Same run |
| Pair CDF | 5.3 ms per pair | Same run, 1,000 pairs |
| Closed-form slice | 0.65 ms against 5.2 ms; within 0.05 dB | Prototype, 2026-11-10, six geometries |
| 99.999% Monte Carlo interval | 66.8 to 71.1 dBi | `convergence_rev5.csv` |
| Sector 10° sensitivity | +13.9 dB at 99.9%, +28.4 dB at 99.99% | `sensitivity_sweep_rev5.csv` |
| ATC revolutions | 750 per hour | 12.5 rpm |
| Pair CDF after step 7 | 1.12 ms per pair (5,000 pairs) | Step 7 benchmark, 2026-11-10 |
| Exact CDF bias (corrected) | 0.2 dB (`mag_exact_cdf_rev1`), 0.3 dB (`mag_pair_cdf_rev1`) | Red team X3 |
| ROC spec minus NSSL measured | Median 15.0, 13.2, 10.6, 11.3, 1.0, 0.0 dB at 50 to 99.999% | `docs/figures/redteam/C1_rev1.csv` |
| NOAA case 99.99% MAG | 50.3 dBi (baseline 45.5) | `docs/figures/redteam/NOAA_rev1.csv` |
| Volumes reaching 70 dBi, 120 h | 12.3% (95% 10.6 to 14.2%) | `docs/figures/redteam/D2_rev1.csv` |
| Supplemental low cut effect | Up to +3.1 dB | `docs/figures/redteam/A4_rev1.csv` |
