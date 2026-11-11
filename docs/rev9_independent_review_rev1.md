**DRAFT//PRE-DECISIONAL**

# Independent recheck of the revision 9 working paper

Subject: `docs/mutual_antenna_gain_methodology_report_rev9.md`, commit `0f81bb5` of `nicklasorte/mutual_antenna_gain`. Rechecked November 11, 2026, by a separate agent that didn't write revision 9 (red-team layer 2).

What I did:

- Read the prior recheck (`docs/rev8_independent_review_rev1.md`), `report_numbers_rev9.m`, `report_figures_rev9.m`, the `nexrad_asr11_mag_rev1.m` diff and the two new tests.
- Recomputed every number in scope from `docs/figures` in Python (pandas). The sources were `*_rev9.csv`, `ts_all_modes_rev8.csv`, `seeds_summary_rev8.csv`, `seeds_runs_rev8.csv`, `timing_rev8.csv` and `takeoff_worst_rev8.csv`.
- Ran short checks in MATLAB R2025b Update 4 on the Mac copy (`~/github-runs/mag_lookup`). Its `report_numbers_rev9.m` (14,777 bytes) matches the repo file in size.

I edited and committed nothing in the repository.

## Summary

- **14 findings:** 1 High, 2 Medium, 7 Low, 4 Info.
- **S1 to S11, Q2 and Q7:** 11 fixed, 2 partly fixed (S2 and Q2).

**Verdict.** Revision 9 does what the recheck asked, and its arithmetic is sound:

- Every per-hour value, rank, ratio and duration in Tables 9-2, 9-4, 9-5 and 13-3, and in Appendix B.1, matches its CSV.
- The 20-seed rerun reproduces every revision 8 column exactly (largest difference 0).
- `lock_fit_rev9.csv` recomputes exactly from `lock_trials_rev9.csv`, and the phase slip recomputes from the rate to 6e-11°.

Two run-length conclusions are wrong:

- **T1 (High).** The "longest uncensored" VCP 212 run (43 volumes, 4.2 h) isn't a separate event. It sits in the same 24-hour trial as the censored 54-volume run (seed 20261207, trial 2, slip 0.21° per volume). In MATLAB, that trial has runs at volumes 66–105, 142–184 and 190–243 of 243 at 0.5°. Volumes 142 to 243 never fall below 56.9 dBi, which is 102 volumes or 10.0 h, and the trial end cuts them. So the summary's "longest seen … 4.2 hours for VCP 212" and open item 3's "runs of 4 to 7 hours" understate what one lock did in this model. The longest independent uncensored VCP 212 run is 25 volumes (2.5 h).
- **T2 (Medium).** The VCP 34 mechanism isn't "open." Every VCP 34 trial with a run of 14 volumes or more has |k × slip| within 7.8° of a multiple of 360° for some k from 1 to 4. That is a higher-order lock, where k volumes span close to a whole number of ATC revolutions. No trial outside that condition has a run of 14 or more.

T3 (Medium) corrects a summary percentage. The rest is wording, ranking convention and test scope.

## Status of S1 to S11, Q2, and Q7

| ID | Status | Evidence |
| --- | --- | --- |
| S1 | Fixed | Table 9-2 matches `per_hour_modes_rev9.csv` and `per_hour_envelope_rev9.csv` in all 40 cells. Per-hour values recompute to within 5e-15 from `vol_pct_60 × 3600 / t_volume_s` and `cuts_hit_60_pct × cuts05_per_volume × 3600 / t_volume_s`. The revolution rate uses a 4.8 s ATC revolution. "Hit most often" and "worst for persistence" are reworded. Caveats: T4, T5, T6 and T11. |
| S2 | Partly fixed | Censoring is flagged (`run_stats`, `c = b == 1 \| e == numel(hit)`), the 54- and 51-volume runs are flagged, and "use as the bound" is gone. But the reported "longest uncensored run" of 43 volumes for VCP 212 is part of the same censored lock episode (T1), and the summary misstates how often a long run occurs (T3). |
| S3 | Fixed | Section 9.3 limits the claims to two modes and gives the spread over 216 modes from `run_modes_rev9.csv`: 95th-percentile run 1–3 and 2–4; share in runs of two or more 7.38%–55.02% (median 28.96%) and 12.09%–63.64% (median 41.08%). The summary bullet still blends the sources (T7). |
| S4 | Fixed | Rank and ratio are given: 0.36% (rank 4, ratio 0.939) and 0.47% (rank 18 dense, ratio 0.913). The ranks are dense ranks (T6). |
| S5 | Fixed | ASR-9 99.9% ATC range is −0.80° to 0.80° (94 nodes), and the 134-of-162 box note is added. NOAA 99.99% is 0.45° to 2.5°. The other NOAA range is loose (T12). |
| S6 | Fixed | The text now says "55, 60, and 70 dBi". The tie with VCP 35 SAILS×1 is stated: 34.82% against 33.83%, with a seed SD of 1.34 pp (per_trial, 0°). See T13. |
| S7 | Fixed | The Figure 13 alt text matches `seeds_runs_rev9.csv`: fixed offset gives 12 seeds at 5–8, 7 at 12–24, one at 40, and the 17 is hollow; walk 0.01 gives 7–13 with one 17, and the 9 is hollow; walk 0.05 gives 7–12, and the 9 is hollow; the medians of the 95th-percentile runs are 4, 4 and 4.4. The Gaussian step SD and clamp are stated. |
| S8 | Fixed | "3.3 dB, nearly the same". |
| S9 | Fixed | The warning and the `Description` flag work. In MATLAB, VCP 212 AVSET 6.4 warns at path angles of −2°, −1.5°, −1°, 7.8° and above, and doesn't warn at 0°, 6.4°, 7° or 7.5°. The plateau base takes the minimum over every cut. Caveats: T9. |
| S10 | Fixed | Timed at six modes: 0.603–0.721 s per pair in a batch of 32 (mean 0.646 s), so 216,000 × 0.646 s = 38.8 core-hours. 41.9/38.8 is a 7.4% difference. 216 libraries × 2.38 s = 8.6 min. Ratios to the revision 8 run are 1.394, 1.134 and 1.032. Caveats: T8. |
| S11 | Fixed | The SAILS assumption is stated in Section 11.4. |
| Q2 | Partly fixed | The S3 part is closed. The S2 part is partly fixed (T1, T2). |
| Q7 | Fixed | Table 13-1 is marked historical. Tables 13-2 and 13-3 have evidence files (`timing_rev8.csv`, `timing_modes_rev9.csv`, and the machine files). |

## New findings

| ID | Severity | Claim (section) | What I did to test it | Result | Evidence | Suggested fix |
| --- | --- | --- | --- | --- | --- | --- |
| T1 | High | Summary bullet 6: "the longest seen are 4.2 hours for VCP 212 … plus a 5.3-hour run that the trial end cut." Table 9-4: longest uncensored run 43 (4.2 h). Section 9.3: "a rate drawn closer to an exact lock holds the lock longer." Open item 3: "runs of 4 to 7 hours occur." | Matched `seeds_runs_rev9.csv` with `lock_trials_rev9.csv`. In MATLAB, reran seed 20261207 (VCP 212 SAILS×3 AVSET 6.4, fixed offset, ten 24-h trials) and listed the runs and gap values in trial 2. | Broken | The seed 20261207 row has `max_uncensored_run_vols` = 43 and `max_censored` = 54. Every other trial of that seed has a longest run of 6 or less, and only one VCP 212 trial reaches 40 (`lock_fit_rev9.csv`, `n_trials_ge_N` = 1). So the 43 is in the same trial as the 54. MATLAB, 0.5°, trial 2 of 243 volumes: runs at 66–105 (40), 142–184 (43) and 190–243 (54). The 5-volume gap (185–189) peaks at 56.9–59.6 dBi; the 36-volume gap (106–141) falls to the 50.4 dBi plateau. At 0°: 72–99 (28), 144–182 (39) and 193–243 (51), with the 183–192 gap at 53.6–59.7 dBi. One rate draw (slip 0.211° per volume) gives at least 102 volumes (10.0 h) of high coupling, all at 56.9 dBi or more and censored. It also gives 137 of 178 volumes at 60 dBi from volume 66 to the trial end (17.5 h). Run length isn't monotonic in slip either: the same lock breaks into three runs. The longest VCP 212 run from another trial is 25 volumes (2.5 h, slip 0.58°). | Report "longest uncensored run from a trial whose longest run is uncensored" (25 volumes, 2.5 h, for VCP 212), or report lock episodes per trial. State that the 54-volume trial holds 102 consecutive volumes at 56.9 dBi or more (10 h, censored). Revise open item 3 and summary bullet 6. |
| T2 | Medium | Section 9.3: "for VCP 34, the mechanism is open … don't come from the simple lock." Figure 14 caption: "for VCP 34, long runs don't follow the slip." Open item 11. | For each fixed-offset trial, computed m4 = min over k = 1 to 4 of \|mod(k × slip + 180, 360) − 180\|. In MATLAB, listed the run positions in the long VCP 34 trials. | Broken | VCP 34, 0.5°, 60 dBi: the 8 trials with a run of 14 or more have m4 of 2.6° (40 volumes, k = 3, slip −119.1°), 5.0° (24, k = 4, slip 91.2°), 3.7° (19, k = 1), 5.3° (18, k = 4), 0.6° (17, k = 2), 4.3° (16, k = 3), 3.9° (15, k = 2) and 7.8° (14, k = 3). P(run ≥ 14 \| m4 < 8°) = 38%, against 0% for m4 ≥ 8° (179 trials). The recurrence fits as well. Seed 20261202, trial 3 (4 × slip ≈ 5.0°, expected period 360/5.0 = 72 volumes) has runs starting at volumes 39 and 110, 71 apart. Seed 20261205, trial 6 (slip 3.74°, period 96 volumes) has runs at 17, 64 and 113. Figure 14 shows the clusters near 90°, 120° and 180°. The "48% independent volumes" comparison is true (P ≈ 1e-11) but beside the point. | Replace "open" with the higher-order lock: k volumes ≈ an integer number of ATC revolutions, for k = 2 to 4. Plot run length against m4 (or the nearest-rational distance) instead of \|slip\|. Close open item 11. Count rate draws near every low-order resonance, not only k = 1, when you relate lock length to the FAA tolerance. |
| T3 | Medium | Summary bullet 6: "about 1% of rate draws give a run of 2 hours or more." | Converted each trial's longest run to hours with `t_volume_s`. | Broken for VCP 34 | At 0.5°, VCP 34 has 9 of 200 trials (4.5%) at 2 h or more (12 volumes or more at 638.2 s); VCP 212 has 1.0%. At 0°: 1.5% and 1.0%. Section 9.3 states it correctly as 1% for 20 volumes or more (2 h for VCP 212, 3.5 h for VCP 34). | Say "about 1% give a run of 20 volumes or more (2 h for VCP 212, 3.5 h for VCP 34); 4.5% of VCP 34 draws give 2 h or more." |
| T4 | Low | Section 9.2: "The modes with three SAILS rescans and early AVSET termination revisit the 0.5° cut most often, so they give the most high-coupling volumes, cuts, and revolutions per hour." | Computed 0.5° revolutions per hour (`cuts05_per_volume × 3600 / t_volume_s`) and ranked them against the per-hour counts. | Holds, with caveats | The leaders are VCP 12 SAILS×3 AVSET 10.0 with a CD cut (cuts per hour; 76.0 revisits per hour; late termination, longer volume) and VCP 112 SAILS×1 AVSET 10.0 (volumes per hour at 0.5°; one rescan; 62.7 revisits per hour). The mode with the most revisits, VCP 12 SAILS×3 AVSET 6.4 (107.6 per hour), ranks 6th to 7th for cuts per hour. The Spearman correlation of cuts per hour with revisits is 0.71 (0°) and 0.72 (0.5°). For volumes per hour against volumes per hour of the mode, r = −0.04. | Say the fast modes as a group lead and the revisit rate explains most of the order, but not the top. Drop "early AVSET termination". |
| T5 | Low | Summary bullet 7 and the changes list: the CDF-worst mode is "within 10% of the largest" for cuts per hour. | Ratio | Weakened (rounding) | 2.523/2.807 = 0.8989, which is 10.1% below the largest (0.5°). | Say "within about 10%" or "90% of the largest". |
| T6 | Low | Table 9-2 ranks ("rank of 216, 1 = largest"); Section 9.2 "rank is 3 to 19 of 216". | Recomputed dense ranks and competition ranks. | Weakened (convention) | `unique(-x)` gives dense ranks. Ties lift some ranks: share of cuts 151 → 152 (0°) and 145 → 148 (0.5°); revolutions per hour at 0.5° 18 → 20, because two modes tie above. The standard rank range is 3 to 20. | Use competition rank (`tiedrank`) or say "dense rank". |
| T7 | Low | Summary bullet 6: "the 95th-percentile run is one to five volumes over all modes." | Compared `run_modes_rev9.csv` with `seeds_summary_rev9.csv`. | Weakened (wording) | Over the 216 modes (one-hour trials) it is 1–3 (0°) and 2–4 (0.5°). The 5 comes only from the 24-hour seed study of two modes. | Say "1 to 4 over all 216 modes in one-hour trials; up to 5 in the 24-hour trials of the two modes studied". |
| T8 | Low | Section 13: "The cost hardly depends on the mode … VCP 212 is the most expensive of the six." | Compared the mode spread with the run-to-run spread, and read `timing_rev9_machine.txt` and the open-pool logic. | Holds, with caveats | The mode spread is 0.603–0.721 s (±9%). The same VCP 212 code ran 3% to 39% slower than in revision 8, so "most expensive" can't be resolved. The machine file says the parallel pool for seeds was open during the timings, because the pool stays open after the `parfor`. The six modes leave out the modes with the most low-cut revisits (VCP 12 SAILS×3 AVSET 6.4: 108 0.5° revolutions per hour, against 81 for VCP 212). Single runs. | Time with no pool, repeat 3 to 5 times, and include the VCP 12 SAILS×3 modes. Drop "most expensive". |
| T9 | Info | S9 warning: "best on-path gain more than 20 dB below peak"; Section 8.2: warns for "a path above the top cut". | MATLAB sweep of `el_path_deg` (VCP 212 SAILS×3 AVSET 6.4, ROC). | Holds, with caveats | The flag turns on between 7.5° and 7.8°, which is 1.1° to 1.4° above the 6.4° top cut. At 7.5°, the unflagged "main beam" row uses 28.5 dBi (17 dB below peak); main beam to main beam reads 49.6 dBi. When the flag is on, `Description` reports "NEXRAD cut 0.50 deg", because the first cut wins the 16.5 dBi tie, even with the path near the 6.4° cut. It also warns for paths 1.5° or more below the lowest cut (−1° and −2° tested), which the text doesn't mention. `testRefLevelsFlagNoMainBeamOnPath` tests only 0° and 10°. | Report the gain gap in `Description`, pick the nearest cut on ties, and say "above the top cut or below the lowest cut by more than about 1.3°". |
| T10 | Info | Section 15.1: "`testReportNumbersRev9` checks that the longest run of every seed is either its longest uncensored or its longest censored run"; the test comment says "a run touching a trial end is flagged censored". | Read `test_scripts_rev1.m`. | Holds, with caveats | `max(max_uncensored, max_censored) == max_run` is true by construction (`s(4) = max(len)`, `mu` and `mcn` split `len`), so the test doesn't exercise the censoring rule. The data are right: I checked the 54- and 51-volume runs in MATLAB. | Add a unit test of `run_stats` with a synthetic hit vector whose run ends at the last volume. |
| T11 | Low | Section 11.8 and the changes list: "Counts per hour beside every share." | Read Tables 9-1, 9-3 and 9-6. | Weakened (wording) | Per-hour counts are given at 60 dBi only (Table 9-2). The 55 and 70 dBi shares, the 50 and 70 dBi cut shares, and Tables 9-3 and 9-6 have none. | Say "at 60 dBi", or add the other levels to `per_hour_modes_rev9.csv`. |
| T12 | Low | Section 11.4: the largest NOAA values occur "with the WSR-88D takeoff angle between 0.45° and 1° at 99%, 99.9%, and 99.999%". | Read `takeoff_worst_rev8.csv`. | Weakened (wording) | NOAA: 99% at 0.90°–1.00°, 99.9% at 0.60°–0.65°, 99.999% at 0.70°. The union is 0.60°–1.00°; 0.45° comes only from the 99.99% row. | Say "0.6° to 1°". |
| T13 | Info | Section 9.2: the near tie between VCP 34 and VCP 35 SAILS×1. | Checked the other levels. | Holds, with caveats | At 0.5° and 70 dBi, VCP 34 leads by 1.0 pp (29.29% against 28.26%), against a seed SD of 1.74 pp (per_trial, 0.5°), so it also ties there. The 1.3 pp SD comes from sets of ten 24-hour trials (10 rate draws), not from the 120 one-hour trials of Table 9-1. | Note the second tie and the source of the SD. |
| T14 | Info | `run_stats`: the per-trial longest run and its censored flag. | Read the code. | Holds, with caveats | `[tr.longest(jj), im] = max(l)` takes the first maximum. If a censored run ties with an uncensored one, `longest_censored` depends on order, while Figure 13 marks such ties hollow (`max_censored >= y`). This can't be checked from the CSVs, and no reported number seems to depend on it. | Flag censored when any longest run is censored, as Figure 13 does. |

## Claims tested that held

- **Table 9-2.** All 40 cells match. Ratios to the largest: 0.966 and 0.931 (volumes per hour, "within 7%"), 0.909 and 0.899 (cuts per hour), 0.939 and 0.913 (revolutions per hour, "within 9%"). VCP 31 and VCP 32 rank 209 to 216 for volumes per hour (209, 211, 214 and 216) and 96 to 111 for cuts per hour ("near the middle").
- **Summary per-hour bullets and Section 9.4.** 2.5 and 3.4; 2.4 and 3.1; 2.0 and 2.7; 1.2 and 1.3 against 2.3 and 2.5; largest 2.5 and 2.8. A VCP 32 volume is 569.0 s (9.5 min); the VCP 212 volume is 353.7 s (5.9 min) with 8 0.5° revolutions. Figure 15 caption: 9.8% against a median of 3.0% is 3.3 times.
- **Table 9-4.** Every cell matches `seeds_summary_rev8.csv` and `seeds_summary_rev9.csv`. Durations: 43 × 353.69 s = 4.22 h; 40 × 638.22 s = 7.09 h; 10 → 58.9 min; 9 → 53.1 min; 17 → 3.01 h; 12 → 2.13 h; 54 → 5.31 h. The censored maxima are 54, 4, 7, 17, 9 and 9.
- **Table 9-5.** All 8 cells match `lock_fit_rev9.csv`, which I recomputed exactly from `lock_trials_rev9.csv`: 6.5%, 2.5%, 1.0%, 0.5% and 16.0%, 5.5%, 1.0%, 0.5%. Hours: 0.59, 0.98, 1.96, 3.93, 1.06, 1.77, 3.55 and 7.09.
- **Figure 14 alt text.** Longest VCP 212 runs per trial: 54 (hollow), 25 and 12, at slips 0.211°, 0.580° and 1.465°. 94.5% of VCP 212 trials have \|slip\| > 10°, and 90.5% have runs of 2 to 5. VCP 34 runs of 14 to 40 sit at slips of 3.7° to 179.7°.
- **Phase slip.** `mod(r × T_vol + 180, 360) − 180` recomputes from `rate_deg_s` and `t_volume_s` to 6e-11°. The rates span 74.52 to 75.50 deg/s, which fits a uniform ±0.5 deg/s offset. It is the right quantity for the k = 1 lock. T2 adds the k = 2 to 4 locks.
- **Code logic.** The censoring rule (a run touching the first or last complete volume) is right, and I confirmed it in MATLAB for the 54- and 51-volume runs. `modesC` picks VCP 34 SAILS×1 (the largest `vol_pct_60` at 0.5°, with no tie) and the CDF-worst mode. The `seeds_summary` grouping maps NaN steps to 0 correctly. The `reuse_seeds` path reads `longest_censored` back as 0 or 1, and `nnz(k & L.longest_censored)` works on doubles. The Figure 13 hollow rule (`max_censored >= max_run`) is conservative at ties. `mode_p = modes(1)` is VCP 34, because `ndgrid` varies the mode fastest.
- **Reproduction of revision 8.** All 480 seed rows agree on every revision 8 column (largest difference 0), and `max_run_vols == max(uncensored, censored)` in every row.
- **"No upper bound in the fixed-offset model."** This holds as a model property, since the rate is continuous and the run length grows as the slip goes to 0. It is a statement about the limit; in the data, a run of 20 volumes or more comes from 1% of draws.
- **"Every random-walk run ends within about three hours."** It holds over 4,800 h per mode and model: the longest is 17 volumes (3.01 h), and no censored random-walk run exceeds 9 volumes. It is an empirical maximum, not a bound.
- **Section 8.2 and 11.4 wording.** S5 and S8 are fixed as stated. The NSSL 99.9% box (2 × 81 = 162 nodes) holds 159 nodes.
- **Section 13.** 38.8 core-hours, within 10% of 41.9; libraries about 9 min; VCP 212 3% to 39% above revision 8.
- **Summary count of findings.** 3 + 6 + 14 + 14 + 11 = 48.
- **S9 code.** The plateau base now takes the minimum NEXRAD back-lobe gain over every cut tilt, the right conservative level for "every revolution reaches". The warning fires at 10° and not at 0°, as the test asserts.

## Checks that needed MATLAB

I ran these on the Mac copy (R2025b Update 4), each in under 10 s:

- VCP 212 SAILS×3 AVSET 6.4, fixed offset, seed 20261207, ten 24-hour trials at 0° and 0.5°: run positions and gap values in trial 2 (T1).
- VCP 34 SAILS×1, seeds 20261205 and 20261202, 0.5°: run positions in every trial with a run of 8 or more (T2).
- A sweep of the S9 warning and `Description` over path angles from −2° to 10° (T9).

Not run: the full `report_numbers_rev9.m`, the test suite (the 229-test claim isn't verified), the timing runs, and a rerun of `per_hour_modes_rev9.csv` from scratch. Those numbers were checked only against their CSVs.

## Spot check by the author

The author reran T1's trial in MATLAB (VCP 212 SAILS×3 AVSET 6.4, fixed offset, seed 20261207, trial 2, 24 hours, 0.5° at both ends). It reproduces: 243 complete volumes; runs at 60 dBi at volumes 66–105, 142–184, and 190–243; every volume from 142 to the end at 57.0 dBi or more; 137 of the 178 volumes from 66 to the end reach 60 dBi. The 43-volume and 54-volume runs are parts of one near-lock that lasts most of the trial.

T2 wasn't rechecked by the author.

## Way forward

| Step | Fix | Findings |
| --- | --- | --- |
| 1 | Report lock events, not only runs at one threshold: group runs in the same trial separated by short dips (or report the longest stretch above a lower level), give the 10-hour near-lock as one event cut by the trial end, and drop "4.2 hours" and "4 to 7 hours" as the longest seen | T1, S2 |
| 2 | Replace "mechanism open" with the higher-order lock: a whole number of ATC revolutions over 2 to 4 volumes; compute the generalized slip for every trial and redraw Figure 14 against it; update open item 11 | T2 |
| 3 | Correct the summary's run-length share for VCP 34 (4.5% of draws give 2 hours or more), the "over all modes" 95th-percentile range (1 to 4), the "within 10%" cut figure, and the per-hour mechanism wording | T3 to T5, T7 |
| 4 | Use standard ranks or say ranks are dense; state that the six-mode timing spread is within run-to-run noise and that a pool was open; time the modes that revisit the 0.5° cut most | T6, T8 |
| 5 | Make the S9 warning name the cut it used and fire on the actual main-beam edge; document that it also fires below the lowest cut; make the censoring test independent of the code under test; fix T11 to T14 wording | T9 to T14 |

## Status

As of November 11, 2026, revision 10 (`docs/mutual_antenna_gain_methodology_report_rev10.md`) answers every finding; its Section 11.9 maps each to its fix. Lock episodes are defined by a 5 dB margin below the level rather than by a gap length, because a gap rule joins nearly every run when 30% to 50% of volumes are high.

**DRAFT//PRE-DECISIONAL**
