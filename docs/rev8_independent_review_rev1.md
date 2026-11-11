**DRAFT//PRE-DECISIONAL**

# Independent recheck of the revision 8 working paper

Subject: `docs/mutual_antenna_gain_methodology_report_rev8.md`, commit `0b6cad5` of `nicklasorte/mutual_antenna_gain`. Rechecked November 11, 2026, by a separate agent that didn't write revision 8 (red-team layer 2). I read the prior review (`docs/rev7_independent_review_rev1.md`), `report_numbers_rev8.m`, `report_figures_rev8.m`, the Q6 code change and its test. I recomputed every number in the scope from `docs/figures/*_rev8.csv` in Python. I also ran short checks in MATLAB R2025b Update 4 on the Mac copy (`~/github-runs/mag_lookup`). That copy is not a git checkout, but its `report_numbers_rev8.m` (20,229 bytes) matches the repo file and contains the Q6 change. I edited and committed nothing.

## Summary

- **11 findings:** 0 High, 2 Medium, 6 Low, 3 Info.
- **Q1 to Q14:** 12 fixed, 2 partly fixed (Q2 and Q7).

**Verdict.** Revision 8 fixes what the prior review asked for. The time series now runs at all 216 modes. Run lengths use 20 seeds for each of three drift models. The reference levels use the cut nearest the path. Every revision 8 table comes from `report_numbers_rev8.m`. Every number I checked in the summary, Sections 8, 9, 11.4 and 13, and Appendix B.1 matches its CSV, with consistent rounding. That covers ranges over seeds, durations (volumes × `t_volume_s`), differences, ratios and core-hour estimates. In MATLAB I reproduced these exactly:

- the VCP 34 SAILS×1 and VCP 32 rows of `ts_all_modes_rev8.csv`;
- six takeoff-surface values for the baseline and the NOAA case;
- the Q6 fix, swept over takeoff angles and both patterns.

Two conclusions need changing:

- **S1.** "Worst for persistence" and "hit most often" depend on the unit counted. Per hour, the mode that is worst for the CDF gives more high-coupling volumes than VCP 34 SAILS×1, and about 1.9 times as many exposed 0.5° cuts as VCP 31 or VCP 32.
- **S2.** The headline fixed-offset run of 54 volumes (5.3 h) is cut off by the end of its 24-hour trial. It is a lower bound, so it can't serve as "the bound" that Section 9.3 tells studies to use.

The rest is wording and scope: claims stated "for every mode" that were tested on two modes, location ranges, and figure alt text.

## Status of Q1 to Q14

| ID | Status | Evidence |
| --- | --- | --- |
| Q1 | Fixed | `report_numbers_rev8.m` lines 116–160 run `mode_ts` for every mode in `vcp_mode_cases_rev1` except the working-paper baseline (`mc(2:end)`), at 0/0 and 0.5/0.5, seed 20261105: 432 rows, 216 distinct modes. In MATLAB, VCP 34 SAILS×1 and VCP 32 reproduce exactly (34.8214/48.5714% of volumes at 60 dBi; 9.8459/10.5827% of cuts). The per-hour view is missing (S1). |
| Q2 | Partly fixed | 20 seeds × 3 drift models × 2 modes; shares in runs of two or more and three or more; Figure 13 plots every seed. `seeds_summary_rev8.csv` recomputes exactly from `seeds_runs_rev8.csv` (largest difference 0). But `run_stats` (lines 385–402) doesn't flag runs cut off at the end of a trial, and the 54-volume headline run is one of them (S2). The "every mode" claims rest on two modes (S3). |
| Q3 | Fixed | Section 9.2 gives the bootstrap interval (21.0–26.5% at 60 dBi; 11.9–16.4% at 70 dBi, which matches `bootstrap_rev8.csv`). It also gives the spread over seeds, 3.5–5.5 pp (computed: 3.55, 3.70, 4.77, 5.10, 5.39 and 5.48 pp). |
| Q4 | Fixed | 0.5/0.5: +3.0, +0.75, +2.45, +1.7 dB, and 0 at 99.999%. "0.8 to 3 dB" is right. |
| Q5 | Fixed | The largest value is taken over the 81 × 81 surface. In MATLAB, NOAA 99.9% = 50.65 at (0.6°, −1°), 99% = 34.45 at (0.95°, −0.25°), and baseline 99.999% = 73.40 at (0.7°, −0.2°). All match `noaa_case_rev8.csv`. |
| Q6 | Fixed | `ref_levels` now takes the cut tilt with the largest on-path gain (`nexrad_asr11_mag_rev1.m` lines 282–285). I swept the path from −1° to 10°, with and without a −0.2° cut, for the ROC and NSSL patterns. The reported NEXRAD main-beam gain equals the largest on-path gain over all cuts in all 40 cases. At ties (for example 0.7°, between the 0.5° and 0.9° cuts) it is 44.939 dBi either way. The test fails on the old code (72.5 dBi against 79.38 dBi) and passes on the new code. The labels mislead at angles above the top cut (S9). |
| Q7 | Partly fixed | Every `*_rev8` file is written by `report_numbers_rev8.m`, and Figures 12 to 14 come from `report_figures_rev8.m`. The other cited files have producers: the `bracket_*_rev7` and `takeoff_surface_rev7` files from `report_numbers_rev7.m`, `harm_rev6`, `ref_levels_rev6` and `lookup_pair_rev6` from `report_numbers_rev6.m`, and the sector sweep from `redteam_runs_rev1.m`. Still open: Table 13-1 (the revision 6 and 7 timings, "Measured on November 11") has no evidence file, and the prior review asked to "save the timing runs". |
| Q8 | Fixed | Section 9.2 now gives 45.5 + 5.5 = 51.0 dBi, rests it on the disputed +5.5 dBi, and gives 49.0 dBi at +3.5 dBi. |
| Q9 | Fixed | `worst_modes_rev8.csv` supports every bullet in Section 8.3: 99.9% worst mode with 2 tied modes in all bracket cases; the 99.999% exception for ASR-9 at 0.5/0 (VCP 12 SAILS×3 AVSET 6.4, CD, 1 mode); the NSSL 99% exceptions; 2 to 202 modes tied at 99.99%; 8 to 33 tied for NSSL at 90%. One NOAA row (0.5/0, 99.9%) has 4 tied modes, but the NOAA case isn't in "the bracket". |
| Q10 | Fixed | 79.5 dBi at (0.5°, 0°) in `worked_example_rev8.csv` and `ref_levels_rev8.csv`. |
| Q11 | Fixed | 0.05° grid, ranges given. Against `takeoff_worst_rev7.csv`, the largest increase is 0.10 dB (NSSL 99.9%, 47.55 to 47.65), as stated. Some range statements in the text are wrong (S5). |
| Q12 | Fixed | One computational thread (`timing_rev8_machine.txt`). 3.136/1.270 = 2.47; 2.281/0.699 = 3.26; 1,000 × 0.699 s = 11.7 min; 1,000 × 2.281 s = 38 min; 216,000 × 0.699 s = 41.9 h; ÷ 6 = 7.0 h. The study estimate assumes every mode costs what VCP 212 costs (S10). |
| Q13 | Fixed | `vol_exceed_modes_rev8.csv`, mode that is worst for the CDF, 0°: 99.549% at 50.4 dBi and 32.040% at 50.5 dBi. |
| Q14 | Fixed | Documented in the `antenna_gain_rev1` help (the review offered "assert or document"). There is no assertion. |

## New findings

| ID | Severity | Claim (section) | What I did to test it | Result | Evidence | Suggested fix |
| --- | --- | --- | --- | --- | --- | --- |
| S1 | Medium | Summary bullets 5 and 7, Section 9.2, Table 8-2 caption ("the mode that is worst for persistence"), Section 9.4, Figure 14 caption: VCP 34 SAILS×1 is worst for the share of volumes, and the slow clear-air cuts "are hit most often" | Turned the per-volume and per-cut shares into rates per hour: high volumes per hour = `vol_pct_60` × 3600 / `t_volume_s`. In MATLAB I counted the 0.5° cuts per volume: VCP 31 and VCP 32 have 2, VCP 34 SAILS×1 has 4, and VCP 212 SAILS×3 AVSET 6.4 has 8. | Weakened | **Volumes per hour at 60 dBi, 0.5° path:** VCP 34 SAILS×1 gives 2.74 per hour; VCP 212 SAILS×3 AVSET 6.4 gives 3.14 (rank 19 of 216); the largest, 3.37, is VCP 112 SAILS×1 AVSET 10.0. At 0°: 1.96 against 2.42; the largest is 2.50. **Exposed 0.5° cuts per hour at 60 dBi:** VCP 32 gives 12.65 × 10.58% = 1.34 and VCP 31 gives 1.30, against 81.43 × 3.10% = 2.52 for VCP 212. At 0°: 1.25 against 2.25. A longer volume, or a cut that comes round less often, raises the share per unit only because there are fewer units. The explanations hold per unit, but "hit most often" is false per hour, and per hour the mode that is worst for the CDF is within 4% to 7% of the largest. | Report per-hour counts beside the shares, and say which unit the link budget needs: the share of products affected, or how often events happen. Reword "hit most often" and "worst for persistence" as "largest share of volumes or cuts". |
| S2 | Medium | Summary bullet 6, Section 9.3, Table 9-3, open item 3: fixed offset gives "up to 54 volumes (5.3 hours) for VCP 212 … a study … should use the fixed-offset result as the bound"; "runs of 5 to 7 hours occur" | Reran seed 20261207 (VCP 212 SAILS×3 AVSET 6.4, fixed offset, ten 24-hour trials, 0.5/0.5) and three other long-run seeds in MATLAB, and found where each longest run sits in its trial | Weakened | The 54-volume run reproduces: trial 2, volumes 190 to 243 **of 243**. It ends at the last complete volume, so the trial end cut it. The ATC rate drawn for that trial is 75.3197 deg/s, which gives 73.9995 revolutions per volume: a phase slip of about 0.19° per volume, close to an exact lock. `run_stats` (`report_numbers_rev8.m` lines 385–402) has no censoring flag, unlike `persistence_volume.censored_runs`. The other runs I checked are not censored: VCP 34, seed 20261205, 40 volumes at 51–90 of 135; seed 20261202, 24 volumes at 110–133 of 134; VCP 212, seed 20261206, 25 volumes at 59–83 of 244. So the largest fixed-offset run is a lower bound, limited by the 24-hour trial, and can't be "the bound". The 0° value of 51 volumes comes from the same seed (20261207) and is probably cut off too; I didn't confirm it. This could be read as High, because a study would apply it, but the qualitative point ("hours with a fixed offset") stands. | Add a censoring column to `seeds_runs_rev8.csv`. Report the fixed-offset maximum as "at least 54 volumes; the trial end cut it". Give the lock length as a function of the rate offset, from the phase slip per volume. Drop "use as the bound". |
| S3 | Low | Summary bullet 6 and Section 9.3: "the median run is one volume and the 95th percentile two to five volumes for every mode, seed, and drift model tested", "40% to 76% of the volumes reaching 60 dBi at 0.5°"; summary preamble: "time-series results are for 120 hours at every one of the 216 modes" | Checked `ts_all_modes_rev8.csv` | Weakened | 6 modes have a 95th-percentile run of 1 volume at 0°, 60 dBi (for example VCP 212 AVSET 6.4 and VCP 12). No CSV has a median run per mode. Across the 216 modes, the share of 60 dBi volumes in runs of two or more at 0.5° (one-hour trials) is 12.1% to 63.6%, median 41.1%. The 40% to 76% range covers only the two modes in the seed study. The run-length results come from 2 modes × 4,800 h, not from 216 modes. | Limit these statements to the two modes tested, or add the median and the share of runs per mode to the table of all modes. |
| S4 | Low | Summary bullet 5 and Section 9.2: the share of ATC revolutions reaching 60 dBi "is near its largest at the mode that is worst for the CDF … as expected" | Ranked `rev_pct_60` | Holds, with caveats | At 0°: 0.359% against the largest 0.383% (rank 4). At 0.5°: 0.468% against 0.513% (rank 20 of 216; VCP 35 SAILS×1, 0.505%, ranks higher). Within 9%, but not "near the largest" in rank at 0.5°. | Give the rank or the ratio. |
| S5 | Low | Section 8.2: "for … the ASR-9, the largest 99.9% and 99.99% values hold over a plateau that spans most ATC takeoff angles from −1° to 3°". Section 11.4: "the largest NOAA values occur with the WSR-88D takeoff angle between 0.45° and 1°" | Read `takeoff_worst_rev8.csv` | Broken (wording) | ASR-9 99.9%: ATC −0.80° to 0.80° only (94 nodes), not a plateau. NOAA 99.99%: WSR-88D 0.45° to 2.5° (803 nodes). The bounding box also overstates where a value holds; for example, ROC 99.99% has 134 nodes inside an 18 × 9 = 162-node box. | Correct both ranges, and say the box isn't filled. |
| S6 | Low | Section 9.2: VCP 34 SAILS×1 is largest "at every level from 55 to 70 dBi and both takeoff angles"; Table 9-1 names one mode per statistic | Checked which levels are tabulated, and the margins against the spread over seeds | Holds, with caveats | Only 55, 60 and 70 dBi are tabulated; 65 dBi isn't. At 0°, 60 dBi, VCP 34 SAILS×1 (34.8%) leads VCP 35 SAILS×1 (33.8%) by 1.0 pp. The seed-to-seed standard deviation for VCP 34 is 1.3 pp, so the named mode isn't robust at 0°, though the value is: the mean over 20 seeds is 35.0%. | Say "at 55, 60, and 70 dBi". Note the near-tie with VCP 35 SAILS×1. |
| S7 | Low | Figure 13 alt text: "a few 12 to 24"; Section 9.3: a random walk "with a step of 0.01 or 0.05 degrees per second per revolution" | Read the seed rows and `scan_clock_rev1.m` lines 59–68 | Weakened (wording) | 7 of the 20 fixed-offset seeds fall from 12 to 24 volumes, which is more than "a few". The walk step is Gaussian with that standard deviation, clamped at ±0.5 deg/s, not a fixed step. | Say "7 of 20 seeds" and "a random walk with a step standard deviation of …". |
| S8 | Low | Text after Table 8-2: main-beam levels rise 3.5 dB, and the 99.9% value "rises by the same 3.3 dB" | Arithmetic | Broken (wording) | 79.5 − 76.0 = 3.5 dB; 48.05 − 44.75 = 3.3 dB. They aren't the same. | Say "by 3.3 dB, nearly the same". |
| S9 | Info | Q6 fix: labels at large takeoff angles; related code | Swept the path to 10° in MATLAB | Holds, with caveats | Above the top cut of the mode (VCP 212 AVSET 6.4 at 9° and 10°), "Main beam to main beam" reads 36.0 or 35.1 dBi, with a NEXRAD "main beam" of 16.5 dBi (a sidelobe) and a "near sidelobe" of 5.5 dBi (the far level). `plateau_check` still takes its base level at `el_n(1)` (line 878), which is the supplemental cut when one is added. `ref_levels` ignores a NEXRAD sector. | Flag when no cut has its main beam on the path, and note the sector limit. |
| S10 | Info | Section 13: 42 core-hours per 1,000 pairs at every mode; 7 hours on six workers | Read the timing code (lines 293–318) | Holds, with caveats | Timed only at VCP 212 SAILS×3 AVSET 6.4. The cost per pair at other modes, and the 216 event libraries (216 × 1.9 s for each ATC radar), aren't measured. "7 hours" assumes perfect scaling over six workers. | State that the cost was measured at one mode. |
| S11 | Info | Section 11.4: a cut at −0.2° lowers the largest 99.999% value by 0.65 dB, "because it lengthens the volume" | Read `nexrad_vcp_rev1.m` lines 189–206 | Holds | In "add" mode one CS and CD pair at −0.2° goes first (+45.8 s), and the SAILS rescans stay at 0.5°. The 0.5° share of the volume falls from 51.7% to 45.8% (−0.53 dB of time), which fits the 0.65 dB drop (73.40 to 72.75 dBi). Whether real sites insert SAILS rescans at the supplemental cut instead is a NOAA question. | None. Optionally note the SAILS assumption. |

## Claims tested that held

- **Summary and Table 8-1.** Every value matches `takeoff_worst_rev8.csv`: 33.65, 48.05, 50.45 and 73.40 dBi; NSSL 22.05, 47.65, 49.00 and 73.40 dBi; ASR-9 32.80, 49.70, 51.00 and 72.35 dBi; and each location range. The summary's "0.4° to 1.25°" covers the ROC ranges.
- **Tables 8-2 and 8-3.** All 24 cells match `worked_example_rev8.csv`, and all 20 cells match `ref_levels_rev8.csv`. "One revolution in 2,800" = 1 / 0.000356 = 2,809.
- **Table 8-4, the summary pattern bullet, and the Section 8.3 bullets** match `bracket_cdf_rev7.csv` and `worst_modes_rev8.csv`. The 99.9% value moves by 44.75 − 44.00 = 0.75, which rounds to 0.8 dB.
- **Table 9-1.** All 40 cells match `ts_envelope_rev8.csv`. Revolution share ranges: 0.148% to 0.383% at 0° and 0.244% to 0.513% at 0.5°. The share of volumes at 60 dBi beats the mode that is worst for the CDF at 22 modes (0°) and 49 modes (0.5°). Median share of cuts hit: 3.03% and 3.32%.
- **Table 9-3.** All 42 cells match `seeds_summary_rev8.csv`. Durations hold: 54 × 353.69 s = 5.31 h; 40 × 638.22 s = 7.09 h; 17 × 638.22 s = 3.01 h; 10 × 353.69 s = 58.9 min; 9 × 353.69 s = 53.1 min; 12 × 638.22 s = 2.13 h. 20 × 10 × 24 h = 4,800 h.
- **Trial length.** 6 volumes (35.4 min) against 14 (82.5 min); 9 against 14. The largest change in the share of volumes is 1.88 pp, under 2 pp.
- **Volume and 0.5° times** (`nexrad_vcp_rev1.m`). VCP 34 SAILS×1: 4 × 0.5° revolutions = 2 × (360/5.043 + 360/8.492) = 227.6 s, volume 638.2 s. VCP 212 SAILS×3: 4 × (360/21.149 + 360/12.524) = 183.1 s, volume 353.7 s. VCP 31 and 32 rotate the 0.5° cut at 4.5 to 5.1 deg/s; VCP 212 at 12.5 to 21.1 deg/s.
- **Figure 12** alt text and caption: 59.3% and 36.6% at 50.4 dBi; 48.6% and 30.9% at 60 dBi; 29.3% and 20.0% at 70 dBi; last nonzero bin 79.3 dBi; ratio 1.57. **Figure 14:** 9.85% and 10.58% at the top, falling to about 4.7% by the ninth mode and 1.9% to 2.0% at the bottom. **Figure 13:** walk 0.01 gives 7 to 13 volumes and one 17; walk 0.05 gives 7 to 12; the median of the 95th-percentile runs is 4, 4 and 4.4.
- **Section 11.4 and Table 11-1.** All 35 cells match `noaa_case_rev8.csv`. Differences over the surface: 0.8 dB (99%) to 3.0 dB (50%); 50.65 − 44.75 = 5.9 dB. A 0.2° cut changes the largest values by −0.05 to +0.10 dB.
- **Appendix B.1.** All 17 rows match their sources and rounding.
- **Script logic.** `takeoff_worst` uses the right grid and mask. Ties in `worst_modes` use an exact 1e-9 match on per-mode percentiles. The `seeds_summary` grouping maps NaN to 0 correctly, and I recomputed it exactly. The bootstrap resamples whole trials (the right unit) with a fixed stream. The NOAA rows with a cut pair the NOAA ATC table (`Ta{6}`) with a rebuilt NEXRAD table. The "largest" rows of `noaa_case_rev8.csv` take the maximum of each column, so each percentile can occur at a different angle pair, as the Table 11-1 caption says.
- **The 15-volume run in one-hour trials is a cap.** VCP 12 AVSET 8.0, CD has 229 s volumes, so at most 15 complete volumes fit in one hour.
- **Q6 test.** `testRefLevelsUseCutNearestThePath` tells the old code from the new. `el_path_deg` sets both takeoff angles (`fill_opts`, lines 1126–1127).

## Checks that needed MATLAB

I ran all of them on the Mac copy, as reported above. I didn't rerun the full `report_numbers_rev8.m` or the test suite, and I didn't confirm that the 51-volume fixed-offset run at 0° is cut off by the trial end.

## Spot check by the author

The author recomputed S1 from the revision 8 tables and the cut counts (eight 0.5° revolutions per volume for VCP 212 SAILS×3 AVSET 6.4, two for VCP 32), at a 0.5° path. It reproduces:

| Mode | Volumes at 60 dBi per hour | 0.5° cuts with a radial at 60 dBi per hour |
| --- | --- | --- |
| VCP 212 SAILS×3 AVSET 6.4 | 30.87% × 3,600 / 353.7 s = 3.14 | 3,600 / 353.7 × 8 × 3.10% = 2.52 |
| VCP 34 SAILS×1 | 48.57% × 3,600 / 638.2 s = 2.74 | Not computed |
| VCP 32 | Not computed | 3,600 / 569.0 × 2 × 10.58% = 1.34 |

S2 wasn't rechecked by the author.

## Way forward

| Step | Fix | Findings |
| --- | --- | --- |
| 1 | Add a per-hour view beside every share: high-coupling volumes per hour and exposed 0.5° cuts per hour, with the worst mode for each; reword "worst for persistence" and "hit most often" by the unit counted | S1, S4 |
| 2 | Flag runs cut off at a trial end in the seed runs (as `censored_runs` does for the persistence table); report the longest uncensored run and the censored runs as lower bounds; rerun long fixed-offset trials (or many more rate draws) so the lock length has a distribution; change "use the fixed-offset result as the bound" | S2, Q2 |
| 3 | Limit "every mode" run-length claims to the two modes tested, or report them over all 216 modes from `ts_all_modes_rev8.csv` | S3 |
| 4 | Correct the location ranges, the 55 to 70 dBi claim, the Figure 13 alt text, the random-walk step description and the 3.3 dB wording | S5 to S8 |
| 5 | Label reference levels above a mode's top cut; time the all-mode estimate over several modes; give Table 13-1 an evidence file or mark it historical | S9, S10, Q7 |

## Status

As of November 11, 2026, revision 9 (`docs/mutual_antenna_gain_methodology_report_rev9.md`) answers every finding; its Section 11.8 maps each to its fix. Steps 1 to 5 are done. Step 2's longer fixed-offset runs weren't added: revision 9 instead states that a near-exact lock has no upper bound in the fixed-offset model and counts how often a rate draw gives a long run.

**DRAFT//PRE-DECISIONAL**
