**DRAFT//PRE-DECISIONAL**

# Independent recheck of the revision 10 working paper

Subject: `docs/mutual_antenna_gain_methodology_report_rev10.md`, commit `81c2626` of `nicklasorte/mutual_antenna_gain`. Rechecked November 11, 2026, by a separate agent that didn't write revision 10 (red-team layer 2).

What I did:

- Read the prior recheck (`docs/rev9_independent_review_rev1.md`), the revision 10 sections in scope, `run_stats_rev1.m`, `report_numbers_rev10.m`, `report_figures_rev10.m`, the `ref_levels` function in `nexrad_asr11_mag_rev1.m`, and the three tests named in the task.
- Recomputed every number in scope in Python (pandas) from `docs/figures`: `lock_trials_rev10.csv`, `lock_fit_rev10.csv`, `seeds_runs_rev10.csv`, `seeds_summary_rev10.csv`, `per_hour_modes_rev10.csv`, `per_hour_envelope_rev10.csv`, `per_hour_corr_rev10.csv`, `timing_modes_rev10.csv`, `timing_rev10_machine.txt`, `seeds_runs_rev9.csv`, `per_hour_modes_rev9.csv`, `run_modes_rev9.csv`, `timing_rev8.csv`, and `ts_all_modes_rev8.csv`.
- Ran short checks in MATLAB R2025b Update 4 on the Mac copy (`~/github-runs/mag_lookup`). Its `run_stats_rev1.m` (3,683 bytes) and `report_numbers_rev10.m` (15,053 bytes) match the repository files in size. Each call took under 15 s.

I edited and committed nothing in the repository.

## Summary

- **11 findings:** 1 High, 2 Medium, 4 Low, 4 Info.
- **T1 to T14, S2, and Q2:** 13 fixed, 3 partly fixed (T1, S2, Q2).

**Verdict.** The arithmetic of revision 10 is sound. Every cell of Tables 9-2, 9-4, 9-5, and 13-3, every rank, correlation, count, and duration in Section 9.3 and Appendix B.1 matches its CSV. The 20-seed rerun reproduces every shared revision 9 column exactly. m4, the phase slip, and the competition ranks recompute exactly. "Every run of 14 or more has m4 below 8°; none with m4 of 8° or more does" holds at all four combinations of path (0°, 0.5°) and level (60, 70 dBi), not only at 0.5° and 60 dBi.

Three of the new lock claims don't hold as written:

- **U1 (High).** "A lock can hold the coupling low" is wrong. The two VCP 212 trials that Figure 14 cites as frozen at an alignment that "never reaches 60 dBi" reach it every k volumes for hours. Seed 20261201, trial 3 (k = 3, m4 = 0.014°) has 30 volumes at 60 dBi or more, every third volume from volume 155 to volume 242 of 244. That is 8.6 hours of high coupling every 17.7 minutes, cut by the trial end. Seed 20261211, trial 3 (k = 2) reaches 60 dBi on every second volume from volume 2 to volume 36, peaking at 69.3 dBi. Runs and episodes both need consecutive volumes, so they report these hours-long periodic exposures as length 1.
- **U2 (Medium).** The 10.3-hour headline episode depends on the 5 dB margin and on the 0.5° path. At 0.5° with a 3 dB margin, the same lock gives 56 volumes (5.5 h). At 0° with the 5 dB margin, it gives 55 volumes (5.4 h). The dips that decide the result are 57.0 and 56.96 dBi at 0.5°, and 53.6 dBi at 0°. They aren't "much lower than 55 dBi".
- **U3 (Medium).** The "longest uncensored VCP 212 episode", 46 volumes (4.5 h), repeats the T1 error. It is volumes 63 to 108 of the same trial (seed 20261207, trial 2) as the 105-volume episode. The longest uncensored VCP 212 episode from another trial is 28 volumes (2.75 h).

## Status of T1 to T14, S2, and Q2

| ID | Status | Evidence |
| --- | --- | --- |
| T1 | Partly fixed | Lock episodes are added, and the 105-volume episode (139–243, 97 volumes at 60 dBi, censored) reproduces in MATLAB. But the "longest uncensored episode" of 46 is the earlier episode of the same trial (U3), and the 10.3 h depends on the margin and path (U2). |
| T2 | Fixed | m4 recomputes to 2e-12°; k matches in all 1,600 rows. The 8 of 21 and 2 of 29 counts are right. The rule holds at 0°/60, 0°/70, 0.5°/60, and 0.5°/70. Caveats: U1, U4. |
| T3 | Fixed | Table 9-5 recomputes exactly by hours with N = ceil(H × 3600 / t_v): VCP 212 N = 11, 21, 41; VCP 34 N = 6, 12, 23. Values: 2.0/2.0, 1.0/1.0, 0.5/0.5 and 16.0/27.0, 4.5/5.5, 1.0/1.5. The VCP 212 4-hour share (0.5%) is the one censored trial. |
| T4 | Fixed | Spearman 0.721, 0.710, 0.659 (0°) and 0.587, 0.717, 0.449 (0.5°) recompute. Revisits: 76.0 and 62.7 per hour; VCP 12 SAILS×3 AVSET 6.4 is 107.6 per hour and ranks 6th (0°) and 7th (0.5°) for cuts per hour. Caveat: U5. |
| T5 | Fixed | 2.523/2.807 = 0.899, "at 90%"; at 0° it is 0.909. |
| T6 | Fixed | All 3,024 rank cells match competition ranks from pandas `rank(method="min")`. Table 9-2 ranks: 23, 6, 50, 19, 152, 4, 148, 3, 4, 20. |
| T7 | Fixed | `run_modes_rev9.csv` gives 1–3 (0°) and 2–4 (0.5°); seeds give up to 5. |
| T8 | Fixed | No pool during timing: the timing part runs first and errors if `gcp("nocreate")` isn't empty. The wrapper and the time series default to `parallel_flag` 0. Batch repeats agree within 0.0%–3.8%. Caveats: U7, U9. |
| T9 | Fixed | In MATLAB: flagged at −0.4° (11.4 dB down) and 7.3°, not at −0.3° and 7.2° (9.0 dB down); names the 0.50° and 6.40° cuts. Caveat: U6. |
| T10 | Fixed | `testRunStatsCensoringAndEpisodes` exercises censoring, ties, episodes at 5 and 3 dB, and incomplete volumes on hand-made data; I checked each expected value by hand. Caveat: U10. |
| T11 | Fixed | Section 9.2 and Table 9-2 say "at 60 dBi". The historical row S1 in Section 11.8 still says "beside every share" (U11). |
| T12 | Fixed | "between 0.6° and 1°" in Section 11.4. |
| T13 | Fixed | 29.3% against 28.3% and SD 1.7 pp are stated, with the source of the SD. |
| T14 | Fixed | `tr.censored(jj) = any(c(l == max(l)))` and the same for episodes; tested by the V2 case. |
| S2 | Partly fixed | The same error as T1 recurs at the episode level (U3). |
| Q2 | Partly fixed | Follows S2. |

## New findings

| ID | Severity | Claim and section | What I did | Result | Evidence | Suggested fix |
| --- | --- | --- | --- | --- | --- | --- |
| U1 | High | Changes list: "A lock can also hold the coupling low." Figure 14 caption: "such a rate can also hold the coupling low." Section 9.3: "A lock holds the coupling where it starts, high or low. Two VCP 212 trials with m4 below 0.15° … the lock froze the two scans at an alignment that never reaches 60 dBi." | In MATLAB, reran VCP 212 SAILS×3 AVSET 6.4, fixed offset, seeds 20261201 and 20261211 (ten 24-hour trials each, 0.5°) and listed the volumes at 60 dBi in trial 3 of each. Also read `ep_span` and `ep_high` in `lock_trials_rev10.csv`. | Broken | The CSV already shows that both trials reach 60 dBi: longest run 1, `ep_high` 1, not 0. Seed 20261201, trial 3 (slip −119.995°, k = 3, m4 = 0.014°): 30 of 244 volumes (12.3%) at 60 dBi or more, at volumes 155, 158, …, 242 (every 3rd, all gaps 3), max 63.1 dBi. That is 88 volumes, 8.6 h at one hit per 17.7 min, ending at the trial end. Seed 20261211, trial 3 (slip −179.94°, k = 2): 18 volumes (7.4%), at volumes 2, 4, …, 36 (every 2nd), max 69.3 dBi, starting at the trial start. A k ≥ 2 lock repeats the alignment every k volumes, so it gives periodic hits, not low coupling. Runs and episodes miss them by design. The T1 trial (k = 1, 0.21° per volume) doesn't hold "where it starts" either: it goes high (63–108), down to the 50.4 dBi plateau (109–138), and high again (139–243). | Withdraw the "hold low" wording. Add a periodic-lock statistic, such as the longest stretch in which every k-th volume reaches the level (k = 1 to 4), or hits per hour over a sliding window. Report the 8.6-hour k = 3 exposure beside the episodes. Ask the link budget (open item 11) whether a hit every k volumes for hours counts as one event. |
| U2 | Medium | Summary bullet 6: "such a lock held one VCP 212 episode at or near 60 dBi for at least 10.3 hours." Section 9.3: "Most dips between runs fall much lower than 55 dBi, so an episode joins only runs that a near-lock holds together." Open item 3: "lock episodes of 8 to more than 10 hours occur." | In MATLAB, ran `run_stats_rev1` on seed 20261207 (VCP 212) and seed 20261205 (VCP 34) with margins of 3, 5, 8, and 9.4 dB at 0.5°, and 3, 5, and 8 dB at 0°. Listed the dips in trial 2. Computed the share of volumes from 55 to 60 dBi. | Weakened | VCP 212, 0.5°, longest episode: 56 (3 dB), 105 (5 dB), 106 (8 dB), 107 (9.4 dB). The 3 dB split comes at volumes 186–187 (57.04 and 56.96 dBi). VCP 212, 0°: 54 (3 dB), 55 (5 dB, 5.4 h; split at 53.57 dBi in volumes 183–192), 105 (8 dB). VCP 34, 0.5°: 44, 47, 50 (stable). The 55–60 dBi band holds only 4.2% (VCP 212) and 5.5% (VCP 34) of volumes, so a margin of about 5 dB is defensible because the volumes split into two groups. But the dips that decide the headline sit 2 dB above and 1.4 dB below the cut. The summary doesn't say the 10.3 h is at 0.5°; at 0° the same rate draw reads 5.4 h. | Give the range: "5.4 to 10.3 hours, depending on the path and on a 3 to 8 dB margin; censored". State the path in the summary. Replace "much lower than 55 dBi" with the dip values. Better: define the episode by the rate (one lock per rate draw) and report the hours at 50.5 dBi + x for several x. |
| U3 | Medium | Table 9-4: "Longest uncensored lock episode, VCP 212 fixed offset: 46 (4.5 h)". Section 9.3: "The longest VCP 212 episode that the trial didn't cut is 46 volumes (4.5 hours)." Appendix B.1. | Matched `seeds_runs_rev10.csv` (seed 20261207, 0.5°/60: `max_uncensored_episode_vols` 46) with `lock_trials_rev10.csv` (the longest `ep_span` in any other VCP 212 trial is 28). In MATLAB, listed the 55 dBi stretches of trial 2. | Broken | Trial 2 has 55 dBi stretches at 63–108 (46 volumes, holding the 40-volume run 66–105) and 139–243 (105). So the "uncensored" 46 is the first half of the same lock, the error T1 found for runs. The same holds at 0°: 45 is 141–185 of the same trial. The longest uncensored VCP 212 episode from an independent rate draw is 28 volumes (2.75 h; seed 20261206, trial 5, slip 0.58°). | Report the longest uncensored episode from trials whose longest episode is uncensored (28, 2.75 h), or one value per rate draw. |
| U4 | Info | Section 9.3 and Figure 14: "every run of 14 volumes or more has m4 below 8°"; "Above 8 degrees no run exceeds 12 volumes." | Checked the rule at all four path and level combinations and looked at the m4 of the long runs. | Holds, with caveats | The rule holds at 0°/60 (3 of 21 VCP 34, 2 of 29 VCP 212), 0°/70 (0 of 21, 2 of 29), 0.5°/60 (8 of 21, 2 of 29), and 0.5°/70 (2 of 21, 2 of 29). No trial with m4 of 8° or more has a run of 14 or more at any combination. But 8° sits just above the largest m4 seen (7.76°), so the threshold was set after looking. Runs of 11 and 12 occur at m4 of 18.4° and 15.3°, and the episode maximum with m4 of 8° or more is 16 (VCP 34, 0.5°/60). For VCP 212, no k ≥ 2 lock gives a run above 10. m4 is the drift over k volumes; the drift per volume is m4/k. | Say the 8° split is descriptive, not a fitted threshold. Note the 16-volume episode above 8°. |
| U5 | Low | Section 9.2: the 0.5° revisit rate "explains most of the order" (Spearman 0.59 to 0.72). Changes list: "Spearman 0.59 to 0.72". | Recomputed the correlations. | Weakened (wording) | ρ = 0.45 to 0.72, so ρ² = 0.20 to 0.52. The revisit rate accounts for at most half of the rank variance, and 20% for revolutions per hour at 0.5°. The changes list leaves out the 0.45 to 0.66 for revolutions. | Say "explains part of the order (Spearman 0.45 to 0.72)". |
| U6 | Low | Section 8.2: the warning fires "for a path 0.9° or more above the top cut or below the lowest cut". Warning text: "the path, which is above the top cut or below the lowest". | MATLAB: VCP 212 SAILS×3 with no AVSET (cuts at 8.0°, 10.0°, and 12.5°), paths at 9.0° and 11.25°. | Weakened | Flagged at 9.0° (14.0 dB down, names the 8.00° cut) and 11.25° (21.9 dB down, names 10.00°), both between cuts. The flag is right; the message and Section 8.2 describe the wrong cause. Any mode with a gap of about 1.7° or more between cuts flags mid-gap, and the takeoff grid reaches 10°. | Say "0.9° or more from every cut (beyond the top or lowest cut, or in a gap of 1.7° or more between cuts)" and reword the warning. |
| U7 | Low | Section 13 and open item 10: "about 42 hours on one core … at the cost of VCP 212, the most expensive mode timed". | Compared with `timing_modes_rev10.csv`. | Weakened (stale) | 42 h uses the revision 8 time of 0.699 s. Revision 10 times VCP 212 at 0.735 s, so 216,000 × 0.735 s = 44.1 core-hours. The eight-mode mean is 0.644 s, or 38.6 core-hours (8% below 41.9). | Update to 44 h at the cost of VCP 212, or 39 h at the eight-mode mean, and keep open item 10 consistent. |
| U8 | Low | Summary bullet 6: "every run and episode in 4,800 hours of each mode ends within 3.0 hours". | Duration of 17 × 638.22 s. | Weakened (rounding) | 3.014 h, slightly over 3.0 h. Section 9.3 says "about three hours". | Say "about 3 hours". |
| U9 | Info | Table 13-3 and `timing_rev10_machine.txt`. | Read the file and the min/max columns. | Holds, with caveats | The machine file says both "3 repeats, median" and "single wall-clock runs (tic/toc)". The batch repeats agree within 4%, as stated. The library and one-pair repeats spread up to 17.2% (VCP 212 library: 2.23 to 2.61 s) and 13.0%. "parallel pool open: 0" is hard-coded, but the `gcp` guard makes it true. | Fix the machine-file text. Give the spread of the library and one-pair columns too. |
| U10 | Info | Section 15.1: `testReportNumbersRev10` checks competition ranks, m4, and episode span. | Read the test. | Holds, with caveats | The rank and m4 checks repeat the formulas in the code (`1 + sum(v' > v, 2)` and the same `min(abs(mod(...)))`), so they guard against edits, not errors in the rule. `ep_span >= longest_run` is true by construction. A 2 × 0.5-hour run has no ties or locks. `testRunStatsCensoringAndEpisodes` doesn't test that a stretch at thr − margin with no volume at thr is dropped, that an episode at volume 1 is censored, or that two trials stay separate. `testRefLevelsFlagNoMainBeamOnPath` doesn't test the 0.9° boundary or a gap between cuts. My independent recomputations agree with the code, so no result is affected. | Add the three `run_stats` cases and a test at 7.2°/7.3° and at 9.0° with no AVSET. |
| U11 | Info | Section 11.8, row S1: "Counts per hour beside every share". | Read the row. | Holds, with caveats | Historical row; still overstates (T11). | Append "at 60 dBi". |

## Claims tested that held

- **Table 9-4, all 42 cells.** The ranges of vol_pct, median run, 95th-percentile run, longest run, longest episode, longest uncensored episode, and share in runs of two or more match `seeds_runs_rev10.csv` and `seeds_summary_rev10.csv`. Durations: 105 × 353.69 s = 10.32 h; 46 × 353.69 s = 4.52 h; 47 × 638.22 s = 8.33 h; 10 → 58.9 min; 17 × 353.69 s = 1.67 h; 17 × 638.22 s = 3.01 h; 15 × 638.22 s = 2.66 h.
- **Medians.** The median run is 2 in exactly 5 of the 60 VCP 34 seed sets at 0.5°/60 (4 fixed offset, 1 walk 0.05) and 1 everywhere else.
- **Episode of T1.** In MATLAB, runs at 66–105, 142–184, and 190–243; episode 139–243 with 97 volumes at 60 dBi, censored; the 55 dBi stretch starts at 139 (55.7 dBi) on a ramp of about 2 dB per volume.
- **Reproduction of revision 9.** All 480 rows; every shared column has a largest difference of 0.
- **Table 9-2.** All values recompute from `per_hour_modes_rev10.csv` and equal the revision 9 values (difference 0). VCP 31 and VCP 32 rank 209 to 216 for volumes per hour and 96 to 111 for cuts per hour.
- **Rate-draw counts.** 21 and 29 trials with m4 below 8°; 179 and 171 at 8° or more; 7 of the 8 long VCP 34 runs have k from 2 to 4; both long VCP 212 runs have k = 1. Figure 14 alt text values (54 at 0.21°, 25 at 0.58°, VCP 34 at 0.6° to 7.8°, nothing above 12 beyond 8°) match.
- **Code logic.** Episode definition, censoring of runs and episodes (a stretch touching the first or last complete volume), tie flags, the per-trial loop, and the `lock_fit` by-hours rows are correct. m4 measures what Section 9.3 says it measures. The Figure 14 hollow rule reads `longest_censored`; exact locks are clipped to 0.05° on the log axis. The `part` logic runs timing before any pool and returns when the part is `"timing"`.
- **Table 13-3.** All 32 cells round correctly. Mean 0.644 s, 38.6 core-hours, 216 libraries × 2.11 s = 7.6 min.
- **Section 8.2 wording.** The 10 dB criterion matches the ROC main beam (45.5 − 12(θ/0.925)² reaches 10 dB down at θ = 0.84°); the −0.4°/−0.3° and 7.3°/7.2° boundaries reproduce in MATLAB, and the tie rule names the nearest cut.
- **Section 11.4.** "between 0.6° and 1°" at 99%, 99.9%, and 99.999% (T12).

## Checks that needed MATLAB

Run on the Mac copy (R2025b Update 4), each call under 15 s:

- VCP 212 SAILS×3 AVSET 6.4, fixed offset, seed 20261207, ten 24-hour trials, 0° and 0.5°: episodes with margins 3, 5, 8, and 9.4 dB, the stretch boundaries, and the dip values (U2, U3).
- Seeds 20261201 and 20261211 (VCP 212): the volumes at 60 dBi in trial 3 (U1).
- VCP 34 SAILS×1, seed 20261205, 0.5°: episodes at margins 3, 5, and 8 dB (U2).
- `ref_levels` description at −0.4°, −0.3°, 7.2°, and 7.3° (AVSET 6.4), and at 9.0° and 11.25° (no AVSET) (T9, U6).

Not run: the full `report_numbers_rev10.m`, the test suite (the 231-test claim isn't verified), the timing runs, and a full rerun of the 216-mode tables. Those numbers were checked only against their CSVs. The margin sensitivity was run for one seed of each mode, not all 20.

## Spot check by the author

The author reran U1's first trial in MATLAB (VCP 212 SAILS×3 AVSET 6.4, fixed offset, seed 20261201, trial 3, 24 hours, 0.5° at both ends). It reproduces: 244 complete volumes; 30 reach 60 dBi, at volumes 155, 158, ..., 242, every third volume without exception. The lock repeats high coupling every 17.7 minutes for about 8.6 hours until the trial ends; runs and episodes, which need consecutive volumes, count it as length 1.

U2 and U3 weren't rechecked by the author.

## Way forward

| Step | Fix | Findings |
| --- | --- | --- |
| 1 | Add a lock statistic that doesn't need consecutive volumes: within a trial, the longest stretch in which a volume at the level recurs at least every k volumes (k = 1 to 4), with its censoring; report it beside runs and episodes; withdraw "a lock can hold the coupling low" | U1 |
| 2 | Report the episode's sensitivity to the margin (3, 5, and 8 dB) and path angle; correct "dips fall much lower than 55 dBi"; give the headline as a range with the margin stated | U2 |
| 3 | Report the longest episode from a separate trial (28 volumes, 2.75 h for VCP 212) instead of a second episode from the 105-volume trial | U3, T1 |
| 4 | Fix the wording and small numbers: the 8° cut (U4), "explains most of the order" (U5), the warning's firing cases (U6), the 42 core-hour figure (U7), 3.0 h rounding (U8), the timing machine file (U9), and the historical S1 row (U11) | U4 to U9, U11 |
| 5 | Replace self-referential tests with checks against independent hand calculations | U10 |

## Status

As of November 11, 2026, revision 11 (`docs/mutual_antenna_gain_methodology_report_rev11.md`) answers every finding; its Section 11.10 maps each to its fix. Step 1's periodic runs use an exact period k (volumes v, v + k, v + 2k, ...), not "at least every k volumes", because a looser rule joins chance recurrences when 30% to 50% of volumes are high; a comparison with independent volumes shows which periodic runs exceed chance.

**DRAFT//PRE-DECISIONAL**
