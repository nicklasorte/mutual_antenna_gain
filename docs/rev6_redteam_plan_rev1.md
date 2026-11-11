**DRAFT//PRE-DECISIONAL**

# Red team of the revision 6 working paper, and the way forward

Draft of November 11, 2026. Subject: `docs/mutual_antenna_gain_methodology_report_rev6.md` (commit `dab4a44`). New evidence is in `docs/figures/redteam_rev6/`. This plan covers two things: what is wrong with or missing from the paper, and how to make the code fast enough for an RFS-scale study.

## Status, 2026-11-11

| Item | Status |
| --- | --- |
| P1 to P14 | Answered in revision 7 of the report (`docs/mutual_antenna_gain_methodology_report_rev7.md`, Section 11.6 lists each) |
| Step 1, pin revision 6 references | Done: `make_rev6_reference_rev1`, `test_regression_rev6` (pair percentiles, worst modes, every per-scan maximum and harm count) |
| Step 2, P6, P10, P11, P12 | Done: exact CDF is the wrapper default; bootstrap, bias wording and timings corrected in revision 7 |
| Step 3, V1, V2, V4 | Done, differently from the plan: one batched `mag_pair_timeseries_rev1` computes each trial's states once for a batch of pairs (V1 and V3 together), builds tables once (V2), and evaluates the ATC elevation pattern once per angle (V4). 120 h, one pair: 3.5 s to 1.4 s; batch of 64: 0.42 s per pair. Outputs identical |
| Step 4, analysis runs | Done: takeoff surface, pattern bracket, drift and run length (`report_numbers_rev7`, Figures 10 to 14) |
| Step 5, register under revision 6 inputs | Done for A1, A2, A3, A4, C4, D2 and the NOAA case (`redteam_runs_rev1` option `inputs = "rev6"`) |
| Step 6, V5 to V8 | V5 tried and dropped: 155 of the 216 mode weight vectors are hull vertices, so pruning saved under 30%. V6 done as a bisection percentile search (1.85 to about 1 ms per pair, same results). V7 not needed after V6. V8: the batched time series takes `parallel_flag` over chunks of pairs |
| Step 7, figures | Done: Figures 10 to 14 |
| Step 8, revision 7 | Drafted |
| Step 9, independent review | Open: needs a reviewer who didn't write revision 7 |

## Summary

1. **The headline numbers are at the wrong takeoff angle.** The summary states the 99.9% MAG as "at most 44.8 dBi," but that value is for a horizontal path only. With both takeoff angles at 0.5°, the same envelope gives 48.1 dBi at 99.9% and 73.2 dBi at 99.999%: 3.3 dB and 4.5 dB higher. Revision 5 showed this in Section 8, but the revision 6 summary dropped it.
2. **The persistence results are capped by the trial length.** Each trial runs for one hour, so no run can last longer than 10 volumes. With 12-hour trials and the same model, the longest run at 60 dBi grows from 6 volumes (35 min) to 14 volumes (83 min), and at 70 dBi from 4 to 11 volumes. The "35 minutes" in the summary is an artifact of how the run is set up.
3. **Some metrics sit on a knife edge.** With the specification pattern, 99.5% of volumes reach 50 dBi with the ASR-11 but only 31.9% with the ASR-9. The only difference is a 0.5 dB lower peak gain, which leaves the ATC main beam plus the specification's flat near sidelobe at exactly 50.0 dBi. A result that swings by 68 points when an input moves 0.5 dB shouldn't be in a summary.
4. **The paper says the code defaults to methods it doesn't use.** Section 12.1 lists the exact CDF and the event time series as revision 6 defaults. The wrapper still defaults to Monte Carlo and the 2 ms grid. Anyone who reruns `nexrad_asr11_mag_rev1` with default options gets different methods from the ones the paper describes.
5. **The new defaults weren't red-teamed.** All 39 register attacks ran on the revision 5 inputs. The paper also cites a bootstrap interval computed on a different run (NSSL pattern, seed 1) beside the headline estimate.
6. **The code is fast enough per pair, but not yet for a study.** Measured today: 1.85 ms per pair for the 216-mode CDF, and 3.5 s per pair and mode for a 120-hour time series. Of the 216 modes, only 33 ever set the envelope across 2,000 random pairs, and about 40% of the time-series cost is table building and recomputing scan states that don't depend on the pair. The vectorization plan in Section 5 targets about 0.9 ms per pair for the CDF and about 1.4 s per pair and mode for the time series, and far less per pair once pairs are batched (estimates).

## 1. How the review ran

The review tested the claims a NOAA or FAA reviewer, or an opposing filer, would attack first. It used four methods:

- **Read against the code.** Every method and default claimed in Sections 11 to 16 was compared with the code at `dab4a44`.
- **Rerun the claims most likely to be challenged.** These were the takeoff-angle dependence, the bootstrap interval on the headline run, the effect of trial length on persistence, and the set of modes that set the envelope. The new runs use the revision 6 defaults.
- **Check against sources.** MIT Lincoln Laboratory ATC-138 was reread for the ASR-9 sidelobe value; it says "maximum sidelobes −28 dB down" for the azimuth pattern.
- **Profile the code** with the MATLAB profiler, using the revision 6 tables and the worst mode.

**Limit.** The author of the paper also ran this review, so the model is grading its own work. Before the paper goes to NOAA or FAA reviewers, a separate reviewer that hasn't seen the drafting should repeat it (Section 4.4).

## 2. Findings

| ID | Severity | Claim in the paper | What the review found | Evidence | Fix |
| --- | --- | --- | --- | --- | --- |
| P1 | High | The 99.9% MAG is at most 44.8 dBi (ASR-11) over all 216 modes | True only at 0°/0°. At 0.5°/0.5°: 48.05 dBi at 99.9% and 73.15 dBi at 99.999%. With only the WSR-88D takeoff angle varied (ATC at 0°), 99.999% reaches 73.35 dBi at 0.75° | `takeoff_sweep_symmetric_rev6.csv`, `takeoff_sweep_wsr88d_rev6.csv` | State headline values as a curve against takeoff angle, give the worst value over a stated angle range, and say which angle each number uses |
| P2 | High | Longest run of volumes at 60 dBi or more: six volumes, 35 min | Capped by the one-hour trials (10 volumes per trial). With 12-hour trials: 14 volumes (83 min) at 60 dBi, 11 volumes (65 min) at 70 dBi, 16 volumes (94 min) at 55 dBi. The share of volumes reaching each level stays within 1.9 points | `trial_length_rev6.csv` | Run trials long enough that the longest run ends inside a trial, or report runs as censored; tie lock length to the rate-drift model, not to the trial length |
| P3 | High | 99.5% of volumes reach 50 dBi (specification pattern) | Knife edge: the ATC main beam plus the specification's flat −29 dB sidelobe gives 50.5 dBi for the ASR-11 and 50.0 dBi for the ASR-9, so 99.5% against 31.9% | `timeseries_rev6.csv` | Report per-volume exceedance curves, not fixed thresholds that equal sums of pattern constants; add a sidelobe-margin sensitivity (P4) |
| P4 | Medium | The ROC specification envelope is the WSR-88D pattern | A specification envelope is a limit on peak sidelobes. Used as the gain at every angle, it overstates the time spent near the sidelobe level: the 50% to 99% values, the 50 dBi persistence, and the harm counts at 50 dBi | Section 12.2, ref levels | Bracket the result: specification envelope (upper), the envelope lowered by a stated average-to-peak margin, and the measured pattern (lower). Ask NOAA for a fleet pattern sample |
| P5 | Medium | ASR-9 remote sidelobe +5.5 dBi from ATC-138 | ATC-138 gives "maximum sidelobes −28 dB down" for the azimuth pattern. The model's first sidelobe for the ASR-9 (M.1464-2, −26.2 dB) already exceeds that maximum, and −28 dB is used as the floor at every angle away from the main beam. The two sources disagree, and the floor is a bound, not a typical value | ATC-138 Section II.A; `nexrad_asr11_mag_rev1.m` assumptions | Say so in Section 12.1 and Open item 2; run the ASR-9 at +3.5 and +5.5 dBi as a bracket until an FAA pattern is available |
| P6 | Medium | Section 12.1 table: exact convolution and event windows are the revision 6 defaults | `fill_opts` in `nexrad_asr11_mag_rev1` still defaults `cdf_method` to `"monte_carlo"` and `ts_method` to `"fixed"` | `nexrad_asr11_mag_rev1.m`, `fill_opts` | Either change the wrapper defaults (and re-pin the revision 5 reproductions with explicit options) or correct the table |
| P7 | Medium | The red team tested the results | All 39 register attacks used the revision 5 inputs. A1 (independence), A3 (mode change), D1 (drift), D2 (bootstrap), C4 (tilt) and E1 (sector) weren't rerun with the specification pattern and 216 modes | `redteam_runs_rev1.m`, `base_context` | Add an input switch to `base_context`, rerun the register under the revision 6 defaults, and keep both columns |
| P8 | Medium | Sections 7 to 10 keep the revision 5 inputs | Several of their conclusions now contradict Section 12: Section 8 says "the NEXRAD pattern matters less," Section 10.2 says the ASR-9 99.9% value is unchanged, and Figures 1 to 9 all show the NSSL pattern. A reader who stops at Section 10 gets the old answer | Report text | Move Sections 7 to 10 to an appendix titled for the revision 5 inputs, or rerun their figures under the revision 6 defaults |
| P9 | Medium | No figures for the revision 6 results | Every figure is from revision 5. The new results are tables only; there is no exceedance curve, takeoff-angle curve or harm chart for a non-technical reader | Report | Add a figure script for revision 6 (Section 4.2) |
| P10 | Low | Bootstrap interval 10.6% to 14.2% for volumes reaching 70 dBi | That interval is from the NSSL run with seed 1 (point estimate 12.3%), printed beside the headline 14.2%. On the headline run the interval is 11.9% to 16.4% | `bootstrap_headline_rev6.csv` | Replace the interval |
| P11 | Low | Each value is high by at most 0.3 dB | At grid nodes the red team found values up to 0.005 dB below a refined calculation (B1), and between nodes interpolation can be one bin low in 0.2% of checks (B2) | Register B1, B2 | Write "within 0.3 dB, and conservative except for rare one-bin cases" |
| P12 | Low | Event time series: 18 ms per pair-hour, 2.2 s per 120 h | Measured today: 29 ms per pair-hour for a 120-hour run with the specification pattern (3.5 s), and 35 to 54 ms per pair-hour for 10-hour runs | Section 5.1 | State the machine, pattern and run length with each timing, or rerun under stated conditions |
| P13 | Low | NOAA case | Since the specification became the baseline, the NOAA case differs only in the ATC remote level, so it no longer bounds anything the baseline doesn't already cover except the sidelobe. A NOAA reviewer's real worst case would also include P1 (takeoff angle) and A4 (supplemental low cut) | Section 15.4 | Redefine the NOAA case as the worst defensible combination: specification pattern, +6.5 dBi ATC remote level, worst takeoff angle in the pair's range, and the site's supplemental cut |
| P14 | Info | The summary gives single values | Decision makers get no range. Every headline number moves 3 to 15 dB across the pattern, takeoff angle and sidelobe inputs | Sections 12 and 15 | Give each headline value as a range with the input that moves it most |

## 3. The three high findings in detail

### 3.1 P1: takeoff angle

ASR-11, specification pattern, envelope over 216 modes, same takeoff angle at both ends (`takeoff_sweep_symmetric_rev6.csv`):

| Takeoff angle | 99% | 99.9% | 99.99% | 99.999% | Largest MAG on the path |
| --- | --- | --- | --- | --- | --- |
| −0.5° | 26.8 | 39.4 | 50.0 | 58.1 | 65.4 |
| 0° (paper) | 31.1 | 44.8 | 50.3 | 68.7 | 76.0 |
| 0.5° | 33.2 | 48.1 | 50.3 | 73.2 | 79.4 |
| 1° | 33.3 | 47.3 | 50.0 | 72.6 | 78.9 |
| 2° | 28.5 | 44.6 | 48.5 | 66.4 | 76.9 |

All values in dBi. The worst angle is near 0.5°, where the 0.5° cut puts the WSR-88D main beam on the path. A horizontal path isn't a bounding case; it sits on the steep side of the curve. For real paths, the takeoff angles from the propagation model decide the answer, and the two-part model already takes them per pair. The paper should lead with the curve and state the worst value over the angle range of the actual site pairs.

### 3.2 P2: trial length

Worst mode, specification pattern, the same per-trial rate-offset model (`trial_length_rev6.csv`):

| Level | Volumes reaching it, 120 × 1 h | Volumes reaching it, 10 × 12 h | Longest run, 1 h trials | Longest run, 12 h trials |
| --- | --- | --- | --- | --- |
| 55 dBi | 28.9% | 28.1% | 8 volumes (47 min) | 16 volumes (94 min) |
| 60 dBi | 23.7% | 23.3% | 6 volumes (35 min) | 14 volumes (83 min) |
| 70 dBi | 14.2% | 12.3% | 4 volumes (24 min) | 11 volumes (65 min) |

The shares agree, so the time-weighted results stand. The run lengths don't. Under the per-trial model, the ATC rate stays fixed for the whole trial, so a volume lock lasts until the trial ends or the phase walks out of the window. Neither 1 hour nor 12 hours is a physical value: how long a lock lasts depends on how the ATC rate and the WSR-88D volume time drift, which is open item 3. The fix is to make the drift model, not the trial length, set the run length: long trials, a random-walk drift whose step is a stated input, and runs reported against that step (D1 already has the random-walk option).

### 3.3 P3: knife-edge thresholds

The specification envelope is flat at −29 dB from the edge of the main beam to 10°. So the WSR-88D gain is exactly 16.5 dBi over a wide range of angles, and the ATC main beam adds exactly its peak gain. A threshold at 50 dBi then sits between the ASR-11 (50.5 dBi) and the ASR-9 (50.0 dBi). Any threshold within about 1 dB of a sum of pattern constants will behave this way. Two changes fix it:

- Report the per-volume maximum as an exceedance curve (share of volumes whose maximum exceeds x, for every x). The curve shows the step at 50.5 dBi as a step instead of hiding it in one threshold.
- Add the sidelobe-margin bracket of P4, which moves the step and shows how much the result depends on the flat envelope.

## 4. Way forward

### 4.1 Revision 7 of the paper (text and tables)

1. Rewrite the summary around ranges: each headline value with the worst takeoff angle, the pattern bracket, and the ATC sidelobe bracket (P1, P4, P5, P14).
2. Replace fixed-threshold persistence with per-volume exceedance curves, and report run lengths only with a stated drift model and uncensored trials (P2, P3).
3. Fix the defaults table or the code (P6), the bootstrap interval (P10), the bias wording (P11), and the timings (P12).
4. Move Sections 7 to 10 to an appendix for the revision 5 inputs, keeping the method walk-through in the body (P8).
5. Redefine the NOAA case (P13).

### 4.2 Analysis runs

| Run | What it answers | Inputs | Cost estimate |
| --- | --- | --- | --- |
| Takeoff-angle surface | Worst MAG over the angle range of real pairs | Two-part model, el_n and el_a from −1° to 3° in 0.05° steps, both ATC radars, both patterns | Under 1 min |
| Pattern bracket | How much the flat envelope drives each result | Specification, specification lowered by 3 and 6 dB beyond the main beam (TODO: a public average-to-peak sidelobe margin), measured | 1 min for CDFs; 4 × 3.5 s per time series |
| Drift and run length | How long a volume lock lasts | Random-walk step from 0.005 to 0.05 degrees per second per revolution, 10 trials of 24 h each | About 1 min per step |
| Register under revision 6 | Whether the attacks still pass with the new defaults | `redteam_runs_rev1("all")` with a pattern and mode-set switch | About 30 min |
| Revision 6 figures | Visuals for leadership | Exceedance curves, takeoff-angle curve, harm chart, NTIA style | Small |

### 4.3 Decisions for NOAA and FAA reviewers

1. **Which WSR-88D pattern stands:** the specification envelope (bounding, and it creates the P3 step), the measured KOUN envelope, or a fleet sample that only NOAA has.
2. **Which takeoff angles to use:** per-pair angles from the propagation model (recommended; the two-part model takes them), or a bounding 0.5° for every pair.
3. **Which persistence metric feeds the link budget:** a time percentile, a per-volume exceedance, or the harm metric, and at what drift model.
4. **Which ASR-9 sidelobe value:** the M.1464-2 first sidelobe, the ATC-138 maximum, or an FAA pattern.

### 4.4 Independent review

Before the paper leaves the working group, a reviewer that hasn't seen the drafting should repeat this red team against revision 7 (layer 2), and NOAA and FAA subject-matter experts should review the inputs (layer 3). The checks in Section 7 don't replace either.

## 5. Speed: vectorization plan

### 5.1 Measured baselines (November 11, 2026, Apple silicon Mac, 8 cores, MATLAB R2025b, no pool)

| Item | Time | Profile |
| --- | --- | --- |
| WSR-88D table, 216 modes, 21 tilts, 301 angles | 1.7 s | One-time |
| ATC table | 0.13 s | One-time |
| Pair CDF, 216-mode envelope, 2,000 random pairs | 1.85 ms per pair | `finish_pair` 45% (mode weighting, cumulative sums, percentile search), slices 35% (of which interpolation 27%), grouping and FFT the rest |
| Pair time series, 120 h, one mode | 3.5 s (29 ms per pair-hour) | For a 10-hour profile: ATC and WSR-88D gains 20%, table building 30%, window index expansion (`arrayfun`) 10%, scan-state recomputation 9%, floor 5% |
| Scan-event library, 10 h | 0.20 s | Window search 88% |

**Study scale.** An RFS study of 1,000 pairs needs one CDF per pair (about 2 s in all) and a time series per pair for each mode that sets its envelope. With the 33 modes that set the envelope anywhere, that is 33,000 runs at 3.5 s: 32 hours on one core, or about 4 hours on 8 workers. Restricting each pair to its own worst modes (one to four) cuts that to 1 to 4 hours on one core. The time series is the bottleneck; the CDF isn't.

### 5.2 Work packages

| ID | Change | Target (estimate unless measured) | Acceptance test |
| --- | --- | --- | --- |
| V1 | **Cache pointing states in the event library.** Store, for each window sample, the WSR-88D azimuth, tilt, cut and volume keys, and the ATC azimuth. They don't depend on the pair, but every pair recomputes them today (`scan_clock_state_rev1`, and the `arrayfun` index expansion) | −19% per pair; library memory about 5% of the 2 ms grid | Outputs identical to the current code on `test_timeseries_rev1` and the revision 5 persistence tables |
| V2 | **Numeric arrays inside the trial loop, tables once at the end.** Replace the per-trial `table`, `movevars` and `array2table` calls with preallocated arrays and one `accumarray` per statistic | −30% per pair | Same tables, `isequal` after sorting |
| V3 | **Batch pairs that share a mode pair.** Evaluate gains for P pairs at once (window samples × P, implicit expansion over the takeoff angles), per-scan maxima with `accumarray` on two-column subscripts, and one floor call for all P | 5 to 10 times faster per pair for P ≥ 100 | Each pair's result equals its single-pair run |
| V4 | **Separable ATC gain lookup.** Precompute the ATC azimuth cut on a fine grid with `griddedInterpolant` and add the elevation gain per takeoff angle, in place of evaluating the M.1851-2 field per sample | −10% per pair | Within 1e-9 dB of `antenna_gain_rev1` |
| V5 | **Prune modes that can't set the envelope.** At table build, find for each cell of the angle grid the modes that set the envelope at any percentile or gain at its corners; keep that set per cell (33 of 216 over 2,000 random pairs today) and keep the full set as a check option | `finish_pair` time cut by about 6 times; pair CDF about 1.2 ms | Envelope identical to the full set at every grid node and at 10,000 random off-grid pairs; any miss fails the build |
| V6 | **Batched percentile search.** Replace `max` over a logical array with `sum(E > lim) + 1` on the monotone exceedance, and run `finish_pair` over a page of pairs with `pagemtimes` | −30% of `finish_pair` | Identical percentiles |
| V7 | **Store quantile functions at table build** so the interpolation is a weighted sum of two stored vectors instead of a sort per interval | Interpolation 27% → under 10% | Within one bin of the current interpolation; B2 rerun |
| V8 | **`parfor` over pairs** in the study driver, with each worker taking whole mode-pair groups so V3 batches stay intact | Near-linear in workers | Parallel equals serial for the same seed (existing pattern) |

Not planned: GPU arrays (no NVIDIA GPU on the Mac), and rewriting `antenna_gain_rev1`, which is already vectorized.

**Combined target (estimate).** V1, V2 and V4 remove about 59% of the time-series cost: about 1.4 s per pair and mode for 120 hours, from 3.5 s. V3 then shares the gain and maxima work across a batch; at 5 to 10 times for batches of 100 or more pairs, and with V8 on 8 workers, the 1,000-pair study with each pair's worst modes would take a few minutes. Pair CDF with V5 to V7: about 0.9 ms per pair, from 1.85 ms.

### 5.3 Rules for every speed change

- Pin the current outputs first (the existing regression references plus a new pinned time series under the revision 6 defaults), and accept a change only if they reproduce exactly or within a stated tolerance.
- Rerun fault seeding after each package; a speed change must not let a planted fault through.
- Profile before and after with `timeit` and the profiler, and record the conditions.

## 6. Order of work

| Step | Work | Depends on | Size |
| --- | --- | --- | --- |
| 1 | Pin revision 6 references (CDF at five pairs, 120 h time series) | | Small |
| 2 | P6 (defaults), P10, P11, P12 corrections | 1 | Small |
| 3 | V1, V2, V4 (time series speed) | 1 | Medium |
| 4 | Analysis runs: takeoff surface, pattern bracket, drift and run length (Section 4.2) | 3 | Medium |
| 5 | Register rerun under revision 6 inputs (P7) | 1 | Medium |
| 6 | V5, V6, V7 (CDF speed) and V3, V8 (batching) | 3 | Medium |
| 7 | Revision 6 figure script (P9) | 4 | Small |
| 8 | Revision 7 of the paper (Section 4.1) | 4, 5, 7 | Medium |
| 9 | Independent review of revision 7 (Section 4.4) | 8 | Small |

## 7. Checks behind this plan

- Every rerun number in this plan comes from a file in `docs/figures/redteam_rev6/`, run on November 11, 2026 at commit `dab4a44` plus no code changes.
- The 1-hour rows of `trial_length_rev6.csv` equal `timeseries_rev6.csv` (the paper's run), which ties the new runs to the paper.
- The timings are single runs on a shared machine; expect ±30%.
- Not verified: the vectorization targets marked as estimates, and whether a public average-to-peak sidelobe margin for the WSR-88D exists (TODO).

## Numbers ledger

| Quantity | Exact value | Source | Form used |
| --- | --- | --- | --- |
| 99.9% and 99.999% MAG at 0.5°/0.5° | 48.05, 73.15 dBi | `takeoff_sweep_symmetric_rev6.csv` | 48.1, 73.2 dBi |
| Increase over 0°/0° at 99.9% and 99.999% | 3.30, 4.45 dB | same | 3.3, 4.5 dB |
| 99.999% MAG, WSR-88D at 0.75°, ATC at 0° | 73.35 dBi | `takeoff_sweep_wsr88d_rev6.csv` | 73.35 dBi |
| Longest run at 60 dBi, 1 h and 12 h trials | 6 (2,122 s), 14 (4,952 s) | `trial_length_rev6.csv` | 35 min, 83 min |
| Longest run at 70 dBi, 1 h and 12 h trials | 4 (1,415 s), 11 (3,891 s) | same | 24 min, 65 min |
| Longest run at 55 dBi, 12 h trials | 16 (5,659 s) | same | 94 min |
| Volumes reaching 50 dBi, ASR-11 and ASR-9 | 99.549%, 31.859% | `timeseries_rev6.csv` | 99.5%, 31.9% |
| Volumes reaching 70 dBi, 95% interval, headline run | 11.86% to 16.38% | `bootstrap_headline_rev6.csv` | 11.9% to 16.4% |
| Modes that set the envelope, 2,000 random pairs | 33 of 216 | `worst_mode_set_2000pairs_rev6.csv` | 33 |
| Pair CDF, 216 modes | 1.85 ms per pair | profiler run | 1.85 ms |
| Time series, 120 h, one pair and mode | 3.52 s | profiler run | 3.5 s |

**DRAFT//PRE-DECISIONAL**
