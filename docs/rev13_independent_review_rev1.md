**DRAFT//PRE-DECISIONAL**

# Independent recheck of revision 13

Subject: `docs/mutual_antenna_gain_methodology_report_rev13.md`, commit `ddbc4da` of `nicklasorte/mutual_antenna_gain` (diff against `4e3f76c`). Rechecked November 12, 2026, by a separate agent that didn't write revision 13. I read the changes list, the summary, Section 9.3 (the W1 to W14 bullets, Table 9-6, Figure 15), Sections 11.12 and 15.1, Appendixes A.1 and B.1, and open items 1, 3, and 11. I read the new or changed code: `periodic_runs_rev1.m`, `report_numbers_rev13.m`, `report_figures_rev13.m`, `testPeriodicRunsLockOverK`, `testReportNumbersRev12`, and `testReportNumbersRev13`, plus `lock_drift_multiple_rev1.m`, `lock_drift_m4_rev1.m`, and the pool check in `report_numbers_rev10.m`. I recomputed every changed number in Python (pandas 3.0.5, SciPy 1.18.1) from the `*rev13*` CSV files, wrote an independent Python null, and ran short checks in MATLAB R2025b on the Mac copy. I edited and committed nothing in the repository.

## Verdict

The arithmetic of revision 13 is sound: every cell of Table 9-6, every interval, the ratio ends (11.64 to 92.69, 9.58 to 31.64, 19.66 to 39,825), the counts 8 of 10, 14 of 22, 10 of 37, 22 of 108, the seven trials, the walk chances, and the reproduction of revision 12 all recompute. The code changes do what they say, and the rev 12 review's own W3 label (117 at "≤ 9.62°") was wrong, as revision 13 notes: 117 is the count at 10.29°. But three of the four Medium answers trade one biased statistic for another. The W1 span rate is above average by construction, and the volumes between the lattice hits are below average. The W2 ratio range ignores the sampling of the 200 trial shares, which dominates the VCP 212 null. The W3 "one in four or five" depends on a threshold read from the data. Section 15.1 also describes a period-7 unit test that doesn't exist.

## Status of the revision 12 findings

| ID | Status | Evidence |
| --- | --- | --- |
| W1 | Partly fixed | "Regular, not heavier" is withdrawn and `per_span_high` is computed correctly (first to last high volume of the period-4 run, `report_numbers_rev13.m` lines 291 to 295). Medians 4.286 and 2.848 per hour, averages 3.238 and 2.659, 8 of 10 and 14 of 22, and 45 (2.99 per hour) and 72 (3.12 per hour) all recompute. The new conclusion "at least as often as on average" is biased upward and doesn't hold for 2 of 10 and 8 of 22 runs (Y1). |
| W2 | Partly fixed | The null now draws model trial 1 + mod(i − 1, 200) for null trial i, so each of the 200 shares and lengths is used 100 times (lines 176 to 186, 225 to 239). The substreams don't overlap (base 1e6 × cell, i ≤ 20,000, `mrg32k3a`). My Python null gives 0.19% (VCP 212, 41 volumes) and 0.68% (VCP 34, 46 volumes), consistent with 0.145% and 0.615% within Monte Carlo error. The Clopper-Pearson interval of the null still covers only Monte Carlo error, not the sampling of the 200 shares (Y2). |
| W3 | Partly fixed | The slip-alone counts are given and correct. `lock_any16` matches min over k = 1 to 16 of abs(mod(k·slip + 180, 360) − 180) to 8e-12, with every k identical, and handles negative slips. The threshold is the largest drift among the long runs, so the share depends on a data-chosen cut (Y4). |
| W4 | Partly fixed | `periodic_runs_rev1` with `k_max` = 16 works on hand inputs (MATLAB, below). The seven trials, the periods (equal to the lock or a divisor: 16 → 8), 121 to 241 volumes, and 13.4 to 23.7 hours recompute. The claim that chance can't be separated at periods up to 16 is contradicted by the table's own 8-hour row (Y6), and "13 to 24 hours" hides that these runs hold 9 to 49 high volumes (Y7). No test covers a period above 4 (Y3). |
| W5 | Fixed | The last paragraph of Section 9.3 and open item 3 now name which statistic counts locks over 1 to 4 volumes and which over up to 16. |
| W6 | Fixed | Rewritten: the null "isn't a model of the scans". |
| W7 | Fixed | The summary names 0.5° and 60 dBi and splits the ranges by mode. |
| W8 | Partly fixed | Table 9-6 now gives intervals. The VCP 34 1-hour chance is still printed as "100%" for 99.98% (upper end 99.995%), and open item 3 still says "within a few degrees of a lock" where VCP 34 needs 9.6° (Y11). |
| W9 | Partly fixed | The legend says circles and squares. But the x axis stops at 20°, so 10 of 400 points (20.02° to 21.13°) aren't drawn, and the alt text's "0.1 to 21 degrees" describes one of them (Y5). The y axis still mixes 354 s and 638 s volumes. |
| W10 | Fixed | `lock_trials_rev13.csv` has no `lock16_rate_deg_per_vol` column. |
| W11 | Fixed | Section 13 gives 15% and 12% of the median, and 17% and 13% of the fastest repeat. |
| W12 | Fixed | `testReportNumbersRev12` (line 345) checks `binofit(5, 200)` against [0.0081662 0.0573744]. SciPy gives 0.00816617 and 0.05737435, and MATLAB prints the same values. |
| W13 | Fixed | "At least 90%" is right: Bonferroni needs no independence, and Clopper-Pearson intervals are conservative. It holds only for Monte Carlo error (Y2). |
| W14 | Fixed | 0.922, 0.992, 0.996, 0.873, and the lowest value 0.221 recompute. The CSV column is mislabeled (Y10). |
| V2 | Partly fixed | Mechanism and counts are in. The slip-alone share depends on the threshold (Y4), and periods above 4 are untested (Y3). |
| V3 | Partly fixed | The rate now counts every high volume, but the comparison is biased upward (Y1). |
| V4 | Partly fixed | Per-trial shares are used. The ratio interval still leaves out share-sampling uncertainty (Y2). |

## New findings

| ID | Severity | Finding | Evidence | Fix |
| --- | --- | --- | --- | --- |
| Y1 | Medium | "A long periodic run carries average or above-average exposure … at least as often as on average" (Section 9.3, line 392; changes list line 9; summary) is above average by construction, and it doesn't hold for every run. | The span runs from one high volume to another, and every k-th volume in it is high by definition, so even independent volumes would put it above the average. With each trial's own share, independence predicts a median span rate of 5.36 per hour (VCP 212) and 3.50 per hour (VCP 34), above the observed 4.29 and 2.85. The volumes between the lattice hits are high less often than average: median share 0.167 against a mode share of 0.318 (VCP 212), and 0.318 against 0.471 (VCP 34). 2 of 10 VCP 212 runs (2.87 and 2.99 per hour) and 8 of 22 VCP 34 runs (2.12 to 2.61 per hour) are below the average. The VCP 34 median is only 7% above it. The denominator and time base are the same (`per_vols` and `t_v`, lines 111 to 112), so the bias comes from the selection, not the units. | Say that a long run concentrates high coupling on a fixed lattice, and that the volumes between are high less often than average (give the off-lattice share). Report the span rate as "about the average for VCP 34 and about 1.3 times it for VCP 212, partly by construction". Drop "at least as often as on average". |
| Y2 | Medium | The ratio range "12 to 93" (summary, line 396, Table 9-6) covers Monte Carlo error only. The VCP 212 null at 4 hours depends on a few of the 200 shares. | I fitted the null chance at 41 volumes as a function of share (40,000 Python trials per grid point) and bootstrapped the 200 shares. That gives a null of 0.18%, with a share-sampling 95% range of 0.073% to 0.387%, against the reported 0.097% to 0.208%. The 10 highest-share trials give 67% of the null chance. The highest share (0.568) is the 54-volume lock trial itself. Leaving out the 10 long-run trials halves the null (0.086%). With the share range, the lower ratio end becomes 2.42 / 0.387 = 6.3, not 11.6. VCP 34 is stable: 0.55% to 0.88%, and 0.71% without the long-run trials. | Give the VCP 212 range as about 6 to 93. Say that the null interval doesn't include the sampling of the trial shares, and that the null includes the lock trials' own shares. Consider a leave-long-runs-out null. |
| Y3 | Medium | Section 15.1 (line 765) says `testPeriodicRunsLockOverK` checks "a run at period 7 found with k_max = 16 and missed with 4". The test has no such check. A.1 cites it as evidence for periods up to 16. | `test_nexrad_asr11_rev1.m` lines 317 to 337: the only `k_max` = 16 check uses the period-3 sequence and expects the same [7 3 3]. No test exercises detection or censoring at k > 4. In MATLAB I ran a period-7 sequence by hand: `k_max` = 4 gives [1 1 1 1 0] and 16 gives [43 7 7 43 43], so the code works. The report misstates the evidence, though. | Add the period-7 test that Section 15.1 describes, with censoring at k > 4, or correct Section 15.1 and A.1. |
| Y4 | Medium | "About one tight lock in four or five becomes a long periodic run" (line 394; open item 11) depends on a "tight" threshold set to the loosest lock among the long runs (`report_numbers_rev13.m` lines 129 to 134). Open item 11 passes the share on without any threshold. | The share falls steadily with the threshold. VCP 212: 5 of 11 (45%) at 1°, 8 of 24 (33%) at 2°, 10 of 37 (27%) at 3.26°, 10 of 68 (15%) at 5°, and 10 of 93 (11%) at 8°. VCP 34: 5 of 8 (62%) at 1°, 12 of 21 (57%) at 2°, 19 of 54 (35%) at 5°, 22 of 108 (20%) at 9.62°. At 9.62°, 54% of all VCP 34 rate draws count as "tight". Every slip has a lock within 360/17 = 21.2° over 1 to 16 volumes (the largest observed is 21.13°), so a 9.6° cut selects weakly. "One in four or five" is a fair summary of 27% and 20% only at these cuts. | Give the share as a function of the drift (a small table: 1°, 2°, 5°, 10°). In open item 11, say that a study should take the chance of a long run from that curve and the rate tolerance, not one share. |
| Y5 | Low | Figure 15 clips 10 of 400 points, and the alt text describes a point that isn't drawn. The caption keeps "the tighter the lock …, the longer". | `report_figures_rev13.m` line 178: `xlim([1e-2 20])`. `lock_any16_drift_deg` reaches 21.13°. In the PNG, no marker is right of 20°. The alt text's "runs of 20 to 60 volumes spread over drifts of 0.1 to 21 degrees" uses the hidden 21-volume run at 21.04°. Locks of 0.041°, 0.16°, and 0.26° have runs of only 4 to 13 volumes (the W4 trials), which W4 suggested the caption admit. | Set `xlim` to [1e-2 25] (or 360/17). Correct the alt text. Use "long runs need a tight lock" as the caption. |
| Y6 | Low | "At periods of 5 to 16 volumes, the comparison with chance can't separate locks from chance" (line 395) is checked only at 4 hours. Table 9-6's 8-hour row shows that VCP 212 exceeds chance. | `lock_fit_rev13.csv`, VCP 212, 0.5°/60, 8 h, period up to 16: model 29.5% (23.3% to 36.3%) against chance 16.24% (15.73% to 16.76%), ratio 1.39 to 2.31. At 0.5°/70, 8 h: 13% against 1.28%, ratio 6.0 to 16.4. | State the 8-hour result. Say that the period-16 statistic saturates at 4 hours but not at 8. |
| Y7 | Low | "Periodic runs of 13 to 24 hours at their lock period" (changes list line 12, line 395, B.1) reads as hours of recurring high coupling. These runs hold 9 to 49 high volumes, recurring every 5 to 16 volumes. | `lock_long_k_rev13.csv`: 9 hits over 121 volumes at k = 15 (VCP 34, one high volume every 2.7 hours), and 11, 15, 16, 17, 19, and 49 hits for the others. Duration (h − 1)k + 1 favors large k: two high volumes 16 apart (17 volumes) outrank a 16-volume consecutive run. Censoring also gets broad at k ≤ 16, because a tie with a multiple period ORs in the wider window. In MATLAB, an uncensored k = 7 run (volumes 10 to 52 of 60) is flagged censored through its tied k = 14 sub-run. All seven runs are "cut". | Give hits and the recurrence interval beside the hours. Say that at k > 4 the censoring flag means "within up to 16 volumes of an end". |
| Y8 | Low | "All 237 tests … pass" (Section 15.1, line 766) combines two environments. No single run passes all 237. | `testReportNumbersRev10` fails whenever a pool is open, because `report_numbers_rev10.m` lines 67 to 68 raise `poolOpen`. In the batch run, an earlier test leaves a pool open. The report discloses it, but the headline overstates it. `testReportNumbersRev13` passes `parallel_flag` = 0, but line 58 still uses any open pool. | Make the test close or skip on an open pool (`assumeEmpty(gcp("nocreate"))` or a fixture). Then state the single-run result. |
| Y9 | Low | `testReportNumbersRev13`'s three checks (lines 363 to 365) hold by construction. They don't test the W2 change. | `per_span_high` ≥ `per_hits` (the lattice hits lie inside the span). `per16_vols` ≥ `per_vols` (k ≤ 16 includes k ≤ 4). `lock_any16` ≤ `lock_run` (minimum over a superset). Only an indexing bug could fail them. Nothing checks that null trial i uses model trial 1 + mod(i − 1, 200), or any hand-worked span count. | Add a hand-made trial with known off-lattice hits, and a check of the null's share and length mapping. |
| Y10 | Info | `lock_walk_null_rev13.csv` column `null_trial_ge_pct` holds a fraction, not a percent. | 0.0127 means 1.27%, since 1 − (1 − 0.0127)^200 = 0.9224 as tabulated. `report_numbers_rev13.m` lines 207 and 209 don't multiply by 100. The same is true for rev 12. | Multiply by 100 or rename the column. |
| Y11 | Info | Wording in the changed text. | (1) Line 396: "null trial i draws model trial i's own share". It is model trial 1 + mod(i − 1, 200), each used 100 times. (2) "Within 9.6°" is below the cut of 9.6228°: at 9.6°, 21 of 107. (3) Table 9-6 prints 99.98% as "100%". (4) Open item 3 says "within a few degrees of a lock", but VCP 34 needs 9.6°. (5) The changes list says "average or above-average" although 2 and 8 runs are below. | Fix each. |
| Y12 | Info | The summary's "with a random walk of the rate, every run and episode ends within about 3 hours" leaves out that periodic recurrence under the walk reaches 8.0 hours. | `lock_walk_null_rev13.csv`: VCP 34, 0.5°/60, step 0.01, 45 volumes (7.98 hours), one volume below the 46-volume 8-hour cut. Section 9.3 gives it, but the summary doesn't. | Add "(periodic recurrence up to 8 hours, within chance)". |

## Checks that needed MATLAB

Ran on the Mac copy (R2025b), each call under 2 minutes. The sizes of `periodic_runs_rev1.m` (3,786 bytes), `report_numbers_rev13.m` (15,440), both test files (47,679 and 20,102), and `lock_trials_rev13.csv` (428,681) match the repository.

- `testPeriodicRunsLockOverK`: passes.
- `periodic_runs_rev1` on hand inputs. Period 7 (volumes 10 to 52 of 60): `k_max` = 4 gives [1 1 1 1 0], and 16 gives [43 7 7 43 43] with first 10, last 52, and the censored flag 1 through the k = 14 tie (Y7). Period 16 (volumes 5 to 53): [49 16 4 1 49], censored because 5 ≤ 16, which is correct.
- `binofit(5, 200)`: 0.0081662 and 0.0573744.
- `lock_drift_multiple_rev1(-159.997257, 1, 16)`: 0.0247 at k = 9. ±10.3° gives 10.3 at k = 1 (sign handled).

Not run: `testReportNumbersRev13`, `testReportNumbersRev10`, the full suite, `report_numbers_rev13.m`, and the 20,000-trial MATLAB null. So I haven't verified the 237 and 236 test counts. My null checks (Y2) are Python only.

## Spot check by the author

The author confirmed two findings from the repository before committing this review:

- **Y3:** `testPeriodicRunsLockOverK` (`test_nexrad_asr11_rev1.m`, lines 327 to 328) calls `periodic_runs_rev1` with k_max = 16 only on a period-3 sequence and checks that it finds the same run. No test checks a period above 4, so Section 15.1 and Appendix A.1 of revision 13 overstate the test.
- **Y1:** from `lock_trials_rev13.csv` (0.5°, 60 dBi), the median share of high volumes among the volumes between the periodic hits is 0.19 for the 10 VCP 212 runs of 4 hours or more and 0.32 for the 22 VCP 34 runs of 8 hours or more, against median trial shares of 0.33 and 0.48. The between volumes are high less often than average.

Y2, Y4, and Y5 to Y12 weren't rechecked by the author.

## Way forward

| Step | Fix | Findings |
| --- | --- | --- |
| 1 | Replace "at least as often as on average" with the off-lattice share against the mode share. Note that the span rate is partly above average by construction. | Y1 |
| 2 | Give the VCP 212 ratio with share-sampling uncertainty (about 6 to 93). Show the leave-long-runs-out null. | Y2 |
| 3 | Add the period-7 and k > 4 censoring tests that Section 15.1 describes, plus a span-count and null-mapping test. Make `testReportNumbersRev10` robust to an open pool. | Y3, Y8, Y9 |
| 4 | Report the long-run share as a function of the slip-alone drift, and use it in open item 11. | Y4 |
| 5 | Widen Figure 15's x axis, fix the alt text, and soften the caption. | Y5 |
| 6 | State the period-16 8-hour result, and give hits and recurrence intervals for the seven trials. | Y6, Y7 |
| 7 | Fix the CSV column label and the small wording items. | Y10, Y11, Y12 |

## Status after revision 14

Revision 14 (`docs/mutual_antenna_gain_methodology_report_rev14.md`, Section 11.13) answers Y1 to Y12:

- **Y1:** the share of high volumes between the repeating ones is given (17% against 34% for VCP 212, 32% against 47% for VCP 34), and "at least as often as on average" is withdrawn.
- **Y2:** the chance comparison resamples the 200 trial shares; the ratios become 6.3 to 112 (VCP 212, 4 hours) and 8.2 to 28.5 (VCP 34, 8 hours), close to this review's estimate of about 6 for the lower VCP 212 end.
- **Y3, Y8, Y9:** `testPeriodicRunsLongPeriod`, `testNullTrialMapping`, and `testReportNumbersRev14` are added, and `testReportNumbersRev10` passes with a pool open.
- **Y4:** Table 9-7 gives the share of tight locks that reach a long run at thresholds from 1° to 21.2°.
- **Y5 to Y7, Y10 to Y12:** fixed as listed in Section 11.13.


**DRAFT//PRE-DECISIONAL**
