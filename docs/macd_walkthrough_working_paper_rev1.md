**DRAFT//PRE-DECISIONAL**

# How MACD computes the mutual antenna gain: a step-by-step walkthrough

Working paper, revision 1, draft of November 20, 2026. Repository `nicklasorte/mutual_antenna_gain`, `main` at commit `c909041` plus this revision's changes. Every number and picture in this paper comes from the MATLAB code in the repository, run in MATLAB R2025b. Appendix A lists the code.

This paper grew out of Section 7.4 of the mutual antenna gain methodology report, revision 17 (`docs/mutual_antenna_gain_methodology_report_rev17.md`). Revision 18 of that report now points here.

## Summary

- **What MACD is.** The 2.7 GHz RFS tool gets the mutual antenna gain (MAG) of each radar pair from a C# routine called the Mutual Antenna Coupling Distribution (MACD). MACD builds each antenna's distribution of gain toward the other radar, combines the two, and reads the MAG that's exceeded only a small share of the time.
- **How this paper explains it.** It walks through MACD's five steps one at a time, with a picture and a worked number example for each. The pictures use a MATLAB copy of the C# steps (`macd_rev1`) with steps 100 times finer than the C# code, so each picture shows what the method computes, not the rounding.
- **A Monte Carlo check confirms the method.** Ten billion random, independent pointings per case give the same MAG distribution as the step-by-step calculation, within 0.025 dB at every percentile from 99% to 99.9999% (`macd_monte_carlo_rev1`).
- **The C# settings round the tail up.** With the C# step (a quarter of the beamwidth) and 1 dB bins, the 99.999% MAG comes out 1.1 to 1.7 dB high, and it moves by up to 3 dB with the bearing between the radars.
- **Three defects in the C# code (X5 to X7).** One reads a pattern under the wrong condition, one counts the scan start twice, and one sets two different step sizes.
- **What's missing.** The C# library the RFS tool calls (`AntennaCoupling`) and the file that does the binning and convolution (`DistributionAnalysis.cs`) aren't in the 2.7 GHz project, so three conventions here are assumptions.

## 1. Purpose and audience

This paper is for readers who need to understand, check, or explain the MAG numbers that the RFS tool produces: engineers reviewing the 2.7 GHz link budget, and staff briefing its results. It assumes no knowledge of the C# code. It does assume the reader knows that a radar antenna has a narrow main beam and weaker sidelobes, and that gains add in decibels.

## 2. The idea in one paragraph

MACD answers one question: if you look at the two radars at a random moment, what's the chance that the MAG is at or above a given level? Each antenna turns at its own steady rate, and nothing ties one antenna's pointing to the other's, so the pointing of each antenna is random and independent of the other. MACD therefore does three things:

1. It works out how often each antenna points each way, and so how often it has each gain toward the other radar. This is the antenna's *gain distribution*.
2. It combines the two gain distributions into one distribution of their sum, the MAG. The math for that is a *convolution*.
3. It reads the MAG that's exceeded only a small share of the time, such as 1% (the 99% value).

## 3. The setup

Figure 1 shows the setup from above. The two radars sit at the ends of a path. Each beam points some angle away from the path: φ<sub>1</sub> for the WSR-88D weather radar and φ<sub>2</sub> for the ASR-11 airport surveillance radar. Each antenna's gain toward the other radar depends only on its own angle, so

MAG = G<sub>WSR-88D</sub>(φ<sub>1</sub>) + G<sub>ASR-11</sub>(φ<sub>2</sub>), in dBi.

The MAG is high only when both angles are small, so that both main beams point along the path at the same moment.

![Top view of two radars at the ends of a dashed path. The WSR-88D on the left points its beam about 28 degrees above the path, marked as angle phi 1. The ASR-11 on the right points its beam about 140 degrees away from the direction to the WSR-88D, marked as angle phi 2. Above them, the equation MAG equals G of the WSR-88D at phi 1 plus G of the ASR-11 at phi 2.](figures/fig_macd_setup_rev1.png)

*Figure 1. The MAG adds each antenna's gain at its own pointing angle from the path. Top view; beam widths are exaggerated. Source: `report_figures_macd_rev1.m` (2026).*

The walkthrough uses the two antennas of the methodology report's simple version (Section 1.1 of that report), on a horizontal path:

- **WSR-88D (ROC specification envelope).** G = 45.5 − 12(θ/0.925)² dBi until it meets the sidelobe limit: 16.5 dBi out to 10° off axis and 5.5 dBi beyond. θ is the angle between the beam axis and the path.
- **ASR-11 (M.1464-2 Radar C levels).** G = 34 − 12(φ/1.45)² dBi until it meets the first sidelobe, 9.5 dBi, held out to 3.55° from the beam; 3.5 dBi beyond.

It runs two cases:

- **One cut.** The WSR-88D scans only its 0.5° cut.
- **VCP 212.** The WSR-88D runs volume coverage pattern 212 with three SAILS rescans and AVSET ending the volume at 6.4°. Each cut counts in proportion to its share of the volume time.

## 4. What the C# code holds

The C# files in the 2.7 GHz project hold the single-mode MACD and a spreadsheet driver for it. They don't hold the library that the RFS tool calls.

- **In the project:**
  - `AntennaGainDistributions.cs` holds the MACD calculation: `DetermineMutualCouplingDistribution`, `DetermineAntennaGainDistribution`, and `DetermineAntennaGainAtPct`. Its driver, `GetParameters`, reads the `Inputs` sheet of a workbook, and the result goes to a `MAC Distribution` sheet.
  - `AntennaData.cs` holds the antenna object: the 3 dB beamwidth, the backlobe gain (the gain at 180°), and the step divisor `ANTENNA_SCAN_STEP_SIZE_DIVISOR = 4`.
  - `AntennaGainCalculations.cs` holds the antenna gain models.
  - `CouplingCalculations.cs` (main beam to main beam) and `ScanRateAnalysis.cs` (the scan rate analysis of the link budget working paper, Appendix C) hold the other two MAG methods.
- **Not in the project:**
  - The RFS tool calls `AntennaCoupling.AntennaCoupling.DetermineAntennaCoupling`, from `RFS4SMR.DetermineAntennaCoupling` and `RFSAnalysisDriver.DetermineCoupling`. It passes the two takeoff angles, the coupling percentage, and options for synchronized radars and the inverse scan rate analysis. That library holds the multi-mode loop that the link budget working paper describes in Appendix B.
  - `DistributionAnalysis.cs` holds the histogram, the convolution, and the percentage-point lookup.

So Steps 1, 2, and 4 below follow code in the project, Step 3 follows Appendix B, and three small rules inside `DistributionAnalysis` are assumptions. Steps 2, 4, and 5 say which.

*Table 1. MACD settings in the C# code and their MATLAB equivalents (`macd_rev1`).*

| Setting | C# | C# value here | `macd_rev1` |
| --- | --- | --- | --- |
| Pointing step | `bearingIncrement = BeamWidth_deg / 4.0` | 0.231° (WSR-88D), 0.3625° (ASR-11) | `step_div` |
| Number of pointings | `(int)(AzimuthScanwidth_deg / bearingIncrement) + 1` | 1,557 (WSR-88D), 994 (ASR-11) | same rule |
| First pointing | `MainbeamAzimuth_deg − AzimuthScanwidth_deg / 2` | −180° for a full scan centered on 0° | same rule |
| Gain toward the other antenna | `CalculateAntGain_dBi(model, gain, |relative azimuth|, pattern)` | the pattern of each antenna | `gain_fn(rel_az, tilt)` |
| Histogram bins | `binWidth_dB = 1.0`, starting at `(int)BacklobeGain_dBi` | 1 dB | `bin_dB` |
| Weight per pointing | `UpdateHistogram(gain, 1.0, false)` | 1 | 1 |
| Weight per mode and cut | Appendix B, Eq. (B-1): P(mode, el) = P<sub>mode</sub> × P<sub>el</sub> | share of volume time per cut | `modes.prob` |
| MAG histogram | `DistributionAnalysis.ConvolveHistograms` | | `conv` |
| MAG at a percentage | `ComputeValueAtPercentagePoint` | | bin upper edge, or interpolated |

## 5. Walkthrough settings: finer than the C# code

The walkthrough runs MACD's steps with finer settings than the C# code, so the pictures show the method rather than its rounding. Table 2 lists the four calculations in this paper.

*Table 2. The four calculations and their settings.*

| Calculation | Pointing step | Bins | Pointings per turn (WSR-88D, ASR-11) | Used for |
| --- | --- | --- | --- | --- |
| C# settings | beamwidth / 4 | 1 dB | 1,557; 994 | Section 9 |
| Walkthrough | beamwidth / 100 | 0.01 dB | 38,919; 24,828 | Steps 1 to 5 |
| Fine reference | beamwidth / 400 | 0.001 dB | 155,676; 99,311 | Sections 8 and 9 |
| Monte Carlo | none: random, continuous pointing | 0.001 dB counts | 10<sup>10</sup> random draws per case | Section 8 |

The walkthrough's steps are 25 times finer than the C# step and its bins 100 times finer. Its percentiles land within 0.03 dB of the fine reference at every level (Table 3), so its pictures show the answer the method converges to.

## 6. The five steps

### 6.1 Step 1: Sample the pointings

**What MACD does.** For each antenna, MACD steps the beam all the way around its scan in small, equal steps. At each step, or *pointing*, it works out the angle between the beam and the direction to the other radar, and looks up the antenna's gain at that angle (Figure 2).

**Why a fixed step works.** Each beam turns at a constant rate, so it spends the same time at every pointing. Counting pointings is then the same as measuring time.

**How big the step is.** The C# step is a quarter of the antenna's 3 dB beamwidth: 0.925° / 4 = 0.231° for the WSR-88D and 1.45° / 4 = 0.3625° for the ASR-11. That puts 12 pointings in each main beam above its sidelobe level. The walkthrough's step is a hundredth of the beamwidth, 0.0093° and 0.0145°, which puts 292 WSR-88D pointings and 286 ASR-11 pointings in the main beams. In Figure 2 the walkthrough's pointings form a continuous curve, and the C# pointings are the circles on it.

**Example.** A WSR-88D pointing 0.144° from the path in azimuth, with the beam tilted up 0.5°, is 0.520° off the path, so its gain toward the ASR-11 is 45.5 − 12 × (0.520 / 0.925)² = 41.7 dBi. An ASR-11 pointing 0.3625° from the path has a gain toward the WSR-88D of 34 − 12 × (0.3625 / 1.45)² = 33.25 dBi.

![Two panels of gain toward the other radar against pointing from the path, from minus 2.5 to 2.5 degrees. Dense blue dots, one every 0.0093 degree for the WSR-88D and every 0.0145 degree for the ASR-11, trace each pattern as a continuous curve; dark circles mark the C# pointings, one every 0.231 and 0.362 degree. Left: the WSR-88D at 0.5 degree tilt, peaking near 42 dBi with a 16.5 dBi sidelobe. Right: the ASR-11, peaking at 34 dBi with a 9.5 dBi first sidelobe; two C# circles sit together on the path.](figures/fig_macd_sampling_rev2.png)

*Figure 2. Step 1: the walkthrough samples each pattern 25 times as often as the C# code, which puts about a dozen samples in each main beam. WSR-88D at 0.5° and ASR-11, horizontal path. Source: `macd_rev1.m`, `report_figures_macd_rev2.m` (2026).*

### 6.2 Step 2: Count the gains in bins

**What MACD does.** MACD sorts the gains from Step 1 into narrow bins and counts the pointings in each bin. Dividing each count by the number of pointings turns it into a share of time (Figure 3). The C# bins are 1 dB wide; the walkthrough's are 0.01 dB. The bins start at the backlobe gain and end at the peak gain.

**Example.** The 41.7 dBi pointing of Step 1 goes in the walkthrough bin from 41.70 to 41.71 dBi, which holds 2 of the 38,919 WSR-88D pointings, 0.005% of the time. In the C# version it goes in the 1 dB bin from 41 to 42 dBi, which holds 2 of 1,557 pointings, 0.13% of the time.

**What the distributions show.** Almost all the time, each antenna points away from the other radar. The far sidelobe holds 94.45% of the WSR-88D's time and 98.03% of the ASR-11's. The WSR-88D's near sidelobe, within 10° of the path, holds about 5%. The main beams hold well under 1%: 0.75% for the WSR-88D above its 16.5 dBi sidelobe. Figure 3 plots *share of time per dB*, so bins of different widths compare directly; the walkthrough's bins are drawn summed to 0.25 dB, and bins with no pointings drop to the bottom of the plot.

**Assumption.** `macd_rev1` puts each gain in the bin whose lower edge is at or below it, and puts any gain outside the bins into the end bin. `DistributionAnalysis.UpdateHistogram` isn't in the project, so the C# may differ.

![Two plots of share of time per dB against gain toward the other radar, log scale. Blue lines show the walkthrough distributions; gray bars show the C# 1 dB histograms. Left, WSR-88D with 38919 pointings: a tall spike at the 5.5 dBi far sidelobe, a second at the 16.5 dBi near sidelobe, and a low, rising plateau from 17 to 42 dBi for the main beam. Right, ASR-11 with 24828 pointings: a spike at 3.5 dBi, one at 9.5 dBi, and a plateau from 10 to 34 dBi rising to a peak at 34 dBi.](figures/fig_macd_hist_rev2.png)

*Figure 3. Step 2: each antenna's gain distribution is almost all far sidelobe, with a thin spread of main-beam gains. Share of time per dB, WSR-88D at 0.5° and ASR-11. Source: `macd_rev1.m`, `report_figures_macd_rev2.m` (2026).*

### 6.3 Step 3: Average over modes and cuts

**What MACD does.** Some antennas scan in more than one mode, or at more than one elevation. MACD makes a distribution for each mode and elevation, then averages them. Each one counts in proportion to its probability, P(mode, el) = P<sub>mode</sub> × P<sub>el</sub> (link budget working paper, Appendix B, Eq. (B-1) and Section B.7).

**Example.** In VCP 212 with AVSET ending the volume at 6.4°, the WSR-88D scans nine elevations from 0.5° to 6.4°, 18 cuts in all. The 0.5° elevation has eight of them: its two base cuts and two more for each of the three SAILS rescans. Each cut counts in proportion to its share of the volume time, and the 0.5° cuts take 52% of it. On a horizontal path only the 0.5°, 0.9°, and 1.3° cuts reach above the 16.5 dBi sidelobe, up to 42.0, 34.1, and 21.8 dBi (Figure 4). The higher cuts see the ASR-11 only through their sidelobes.

**What averaging does.** Near the peak, the average holds about half the 0.5° cut's share, because the 0.5° cuts take about half the time. Above the sidelobe overall, it holds two-thirds of the 0.5° cut's share, 0.51% against 0.75%, because the 0.9° and 1.3° cuts add some.

**One more rule in Appendix B.** For its built-in gain models, Appendix B raises each mode's peak gain by 10 log<sub>10</sub>(reference beamwidth / smaller beamwidth) (Eq. (B-2)). The walkthrough gives the patterns directly, so this rule doesn't apply.

![Stair plot of share of time per dB against WSR-88D gain toward the ASR-11, log scale. Gray lines for the 0.5, 0.9, and 1.3 degree cuts end at about 42, 34, and 22 dBi. The blue time-weighted average of all VCP 212 cuts has a spike at 16.5 dBi and runs about half as high as the 0.5 degree cut near 42 dBi.](figures/fig_macd_modes_rev2.png)

*Figure 4. Step 3: averaging the cuts by time spreads the WSR-88D main beam over several cuts and halves its share near the peak. VCP 212 SAILS×3 AVSET 6.4, horizontal path. Source: `macd_rev1.m`, `report_figures_macd_rev2.m` (2026).*

### 6.4 Step 4: Combine the two distributions

**What MACD does.** MACD pairs every bin of one antenna with every bin of the other. A pair's chance is the product of the two shares, because the antennas point independently. The pair's MAG is the sum of the two gains. MACD then adds up the chances of all pairs with the same sum. This is the convolution of the two distributions (`DistributionAnalysis.ConvolveHistograms`; `conv` in MATLAB).

**How to read Figure 5.** Panel a is a grid with one square for each pair of bins: WSR-88D bins across, ASR-11 bins up. A darker square has a larger chance. Squares on a diagonal line share the same sum, so every square on the 60 dBi line gives a MAG near 60 dBi. Panel b adds up each diagonal, which gives the MAG distribution. The grid is drawn with 0.25 dB bins so the squares are visible.

**Example.**

- The bottom-left square pairs the two far sidelobes: 0.9445 × 0.9803 = 92.6% of the time, at a MAG of 5.5 + 3.5 = 9.0 dBi.
- The top-right squares pair the two main beams. Their chances are tiny, a few millionths of the time for a whole 1 dB band of MAG, because each main beam takes well under 1% of its antenna's time.

**How the bins add.** The MAG bin gets the sum of the two lower edges, but each gain can sit anywhere in its bin, so the true sum can be up to two bin widths above that edge. With the C# 1 dB bins, a pair of 41.9 and 34.0 dBi (75.9 dBi) lands in the 74 dBi MAG bin. With the walkthrough's 0.01 dB bins the shift is at most 0.02 dB.

![Two panels. Panel a: a grid of squares, WSR-88D gain bin across from 5 to 42 dBi and ASR-11 gain bin up from 3 to 34 dBi, shaded by the log of the share of time; the darkest square, bottom left, is both far sidelobes at 92.6 percent; purple diagonal lines mark sums of 50, 60, and 70 dBi. Panel b: share of time per dB against MAG from 9 to 76 dBi on a log scale, with a tall spike at 9 dBi, steps near 20, 40, 45, and 50 dBi, and a flat tail near 3e-6 per dB from 51 to 76 dBi.](figures/fig_macd_convgrid_rev2.png)

*Figure 5. Step 4: every pair of bins gives one MAG, and adding the pairs along each diagonal gives the MAG distribution. WSR-88D at 0.5° and ASR-11. Source: `macd_rev1.m`, `report_figures_macd_rev2.m` (2026).*

### 6.5 Step 5: Read the MAG at a percentage

**What MACD does.** MACD returns the MAG that's not exceeded for a given share of time. The link budget working paper's example uses 99%. To find it, MACD adds up the MAG distribution from the bottom until the running total reaches the percentage. The bin where that happens holds the answer (Figure 6).

**Example at 99%.** With the walkthrough's bins, the bins below 34.94 dBi hold 98.99726% of the time, and the 34.94 dBi bin holds 0.0038%. The running total passes 99% inside that bin, so the 99% MAG is 34.95 dBi. A straight-line interpolation inside the bin gives 34.94 + 0.01 × (0.99 − 0.9899726) / 0.0000381 = 34.947 dBi.

**Example at 99.999%.** The running total passes 99.999% inside the 72.31 dBi bin, so the 99.999% MAG is 72.31 dBi, or 72.314 dBi interpolated.

**What changes with the C# bins.** With 1 dB bins, the same lookup lands in a 1 dB bin, and the answer can be the bin's upper edge, which rounds up by as much as 1 dB, or an interpolated point inside it (Section 9). `ComputeValueAtPercentagePoint` isn't in the project, so `macd_rev1` reports both.

![Two zoomed plots of share of time at or above each MAG, in percent, on a linear scale. Left, around 1 percent: the blue stair line, with steps 0.01 dB wide, crosses the purple 1 percent line inside one bin; a dot marks the 99 percent value, 34.95 dBi. Right, around 0.001 percent: the stair line crosses the purple 0.001 percent line at 72.31 dBi, the 99.999 percent value.](figures/fig_macd_lookup_rev2.png)

*Figure 6. Step 5: the percentage falls inside one 0.01 dB bin, which fixes the MAG to within 0.01 dB. WSR-88D at 0.5° and ASR-11. Source: `macd_rev1.m`, `report_figures_macd_rev2.m` (2026).*

## 7. Results of the walkthrough

Table 3 gives the MAG percentiles from the walkthrough, the fine reference, and the Monte Carlo of Section 8.

*Table 3. MAG percentiles, horizontal path: walkthrough (step beamwidth/100, 0.01 dB bins, interpolated), fine reference (beamwidth/400, 0.001 dB), and Monte Carlo (10<sup>10</sup> draws, 95% interval). Source: `docs/figures/macd_pct_rev2.csv`.*

| Share of time below | One cut: walkthrough | One cut: reference | One cut: Monte Carlo | VCP 212: walkthrough | VCP 212: reference | VCP 212: Monte Carlo |
| --- | --- | --- | --- | --- | --- | --- |
| 99% | 34.95 | 34.94 | 34.925 (34.922 to 34.927) | 31.70 | 31.68 | 31.668 (31.665 to 31.672) |
| 99.9% | 45.30 | 45.30 | 45.301 (45.301 to 45.301) | 44.72 | 44.73 | 44.729 (44.728 to 44.730) |
| 99.99% | 50.48 | 50.47 | 50.471 (50.471 to 50.472) | 50.32 | 50.32 | 50.315 (50.313 to 50.316) |
| 99.999% | 72.31 | 72.31 | 72.30 (72.27 to 72.32) | 68.88 | 68.87 | 68.85 (68.80 to 68.90) |
| 99.9999% | 75.62 | 75.63 | 75.62 (75.62 to 75.63) | 75.28 | 75.28 | 75.28 (75.27 to 75.30) |

What the numbers say:

- **Half the time, and 90% of the time, the MAG is 9.0 dBi:** both antennas point away from each other.
- **At 99% and 99.9%,** one main beam crosses the path while the other antenna points away: 31.7 to 45.3 dBi.
- **At 99.99%,** one main beam meets the WSR-88D's near sidelobe: 50.3 to 50.5 dBi.
- **Only at 99.999% and beyond do both main beams meet:** 68.9 to 75.6 dBi. Time-weighting the VCP 212 cuts lowers the 99.999% value by 3.4 dB, from 72.3 to 68.9 dBi, because the WSR-88D spends only half its time at 0.5°.

The one-cut values match the exact calculation of the methodology report's Section 1.1 (34.9, 45.3, 50.5, 72.3, and 75.6 dBi), and the VCP 212 99.999% value is close to the full model's 68.7 dBi on a horizontal path, which uses the M.1851-2 ASR-11 pattern.

## 8. Monte Carlo check

**Why a second method.** The convolution of Section 6 rests on one assumption: the two antennas point independently. A Monte Carlo tests the whole chain without the steps and bins. It doesn't sample a grid; it draws random pointings directly.

**How it works.** For each draw, `macd_monte_carlo_rev1` takes:

1. A random azimuth for each antenna, uniform over its full turn, with no step.
2. For the WSR-88D in VCP 212, a random cut, chosen with probability equal to the cut's share of the volume time.
3. Each antenna's gain toward the other at that pointing, from the same pattern equations, and their sum, the MAG.

It repeats this 10 billion times per case, in chunks of 10 million draws spread over six MATLAB workers, and counts the MAGs in 0.001 dB bins, so the memory needed doesn't grow with the number of draws. Each chunk uses its own random-number substream, so the parallel run gives exactly the same counts as a serial run (`testParallelEqualsSerial`). A 95% interval for each percentile comes from the binomial spread of the number of draws below it.

**Why 10 billion draws.** At 99.999%, only one draw in 100,000 lands above the answer. With 10 billion draws, about 100,000 draws land there, enough to fix the value to a few hundredths of a decibel. Figure 8 shows how the estimate settles as the draws grow: at 100,000 draws a single draw sets the 99.999% value, and its interval spans 7 dB; at 10 million draws the interval is still 1.2 dB wide; at 10 billion it's 0.05 dB wide for the one-cut case.

**What it shows.** The Monte Carlo and the convolution give the same distribution, down to a share of time of one in a billion (Figure 7). Their percentiles agree within 0.025 dB at every level in Table 3. At 99.999% and 99.9999%, the fine reference falls inside the Monte Carlo's 95% interval. From 99% to 99.99% the Monte Carlo's interval is only a few thousandths of a decibel wide, and the reference sits up to 0.013 dB above it: that's the reference's own step and bin rounding, not a disagreement in the method. The exact calculation of the methodology report's Section 1.1, on a 0.0001° grid, gives 34.926 dBi at 99%, against the Monte Carlo's 34.925 dBi. The Monte Carlo's 0.001 dB counts read up to 0.001 dB low.

![Two panels of share of time at or above each MAG on a log scale from 1 to 1e-9. A thick gray line, the Monte Carlo with 1e10 draws, and a thin blue line, the convolution, lie on top of each other from 9 to 76 dBi in both panels: left, WSR-88D at 0.5 degree only; right, VCP 212 with cuts weighted by time.](figures/fig_macd_mc_rev2.png)

*Figure 7. Ten billion random pointings give the same MAG distribution as the convolution, down to one in a billion. Horizontal path. Source: `macd_monte_carlo_rev1.m`, `report_figures_macd_rev2.m` (2026).*

![Two panels of the Monte Carlo 99.999 percent MAG, with its 95 percent interval as a shaded band, against the number of draws from 1e5 to 1e10 on a log scale. Left, WSR-88D at 0.5 degree: the estimate jumps between 70.6 and 74.3 dBi at 1e5 and 1e6 draws, then settles on the convolution line at 72.3 dBi, while the band narrows from 7 dB to 0.05 dB; a dashed line marks the C# settings value, 74 dBi. Right, VCP 212: the estimate settles on 68.9 dBi; the C# settings value is 70 dBi.](figures/fig_macd_mc_conv_rev2.png)

*Figure 8. The Monte Carlo's 99.999% value settles on the convolution's as the draws grow, and the C# settings sit above both. Source: `docs/figures/macd_mc_convergence_rev2.csv` (2026).*

## 9. How close the C# settings get

With the C# step and 1 dB bins, MACD's percentiles land within about 1 dB of the fine values from 50% to 99.99%, and 1.1 to 1.7 dB high at 99.999% (Figure 9, Table 4).

![Two panels of share of time at or above each MAG on a log scale from 1 to 1e-7, the C# settings as a blue stair line and the fine reference as a gray line. Left, WSR-88D at 0.5 degree only: marked bin upper edges of 36, 45, 51, and 74 dBi at 99, 99.9, 99.99, and 99.999 percent. Right, VCP 212: 33, 45, 51, and 70 dBi. The stair lines follow the gray lines and sit slightly to the right near the tail.](figures/fig_macd_csharp_rev2.png)

*Figure 9. The C# settings track the fine distribution, and their tail percentiles come out up to 1.7 dB high. Horizontal path. Source: `docs/figures/macd_pct_rev2.csv` (2026).*

*Table 4. MAG percentiles with the C# settings (bin upper edge; interpolated in parentheses) against the fine reference, horizontal path (`docs/figures/macd_pct_rev2.csv`).*

| Share of time below | One cut: C# settings | One cut: fine | VCP 212: C# settings | VCP 212: fine |
| --- | --- | --- | --- | --- |
| 50% | 9 (8.5) | 9.0 | 9 (8.5) | 9.0 |
| 90% | 9 (9.0) | 9.0 | 9 (9.0) | 9.0 |
| 99% | 36 (35.9) | 34.9 | 33 (32.3) | 31.7 |
| 99.9% | 45 (44.5) | 45.3 | 45 (44.0) | 44.7 |
| 99.99% | 51 (50.6) | 50.5 | 51 (50.1) | 50.3 |
| 99.999% | 74 (73.1) | 72.3 | 70 (69.7) | 68.9 |

Three things move the C# answer away from the fine one:

- **Rounding to the bin.** The bin's upper edge rounds up by as much as 1 dB. The interpolated value can land on either side of the fine value: 44.5 against 45.3 dBi at 99.9%.
- **The step.** In the tail, the step matters more than the bins. With 0.1 dB bins, the 99.999% value is still 73.1 dBi at a quarter-beamwidth step; it falls to 72.4 dBi at a step of 1/40 of the beamwidth (`docs/figures/macd_resolution_rev1.csv`).
- **Where the samples land.** Only a dozen pointings fall in each main beam, so the answer depends on where they fall, and the bearing between the radars sets that. Over bearings from 0° to 2°, the one-cut 99.999% value (bin upper edge) ranges from 71 to 74 dBi, and the 99% value from 34 to 36 dBi. The spreadsheet driver's default bearing, 0°, gives the 74 dBi extreme, because there the ASR-11's twice-counted scan start lands on the path (X6). With the walkthrough's step, the bearing moves the answer by less than 0.1 dB.

## 10. Findings on the C# code

Reading the code for this walkthrough found three defects. They're numbered after X4 in the methodology report's defect register (its Section 11.1):

- **X5. Antenna 2's two-dimensional pattern depends on antenna 1's model.** Before `GetParameters` reads antenna 2's pattern sheet (`Rx2Dptrn`), it checks `antennaModel1 == AntennaModel.PatternProvided2D`, not `antennaModel2`. If only antenna 2 has a two-dimensional pattern, the code doesn't read it.
- **X6. The scan start is counted twice.** For a full turn, `numBearings = (int)(360 / step) + 1` puts the last pointing less than one step after the first: 0.04° for the ASR-11. So the direction where the scan starts carries two samples. It matters only when that direction is near the bearing to the other radar. With the spreadsheet defaults (scan centered on 0°, bearing 0°) it is, for antenna 2, and it raises the 99.999% value by 2 dB. Counting `(int)(360 / step)` pointings for a full turn would fix it.
- **X7. Two step sizes.** `BuildAntenna4MACD` sets the scan step to the smaller of 0.5° and half the beamwidth, but `DetermineAntennaGainDistribution` ignores that value and uses a quarter of the beamwidth. Which step the RFS tool's MACD uses depends on the missing `AntennaCoupling` library.

Two more points are limits, not defects:

- `DetermineAntennaGainDistribution` takes the gain from the azimuth offset alone. The takeoff angle and tilt enter only through the modes and two-dimensional patterns of the missing library.
- The lowest bin starts at the backlobe gain rounded toward zero, so a negative backlobe gain such as −5.5 dBi falls below the first bin, which starts at −5 dB. What happens then depends on `UpdateHistogram`.

## 11. Next steps

1. **Get the missing C# files.** `DistributionAnalysis.cs` and the `AntennaCoupling` library aren't in the 2.7 GHz project. With them, `macd_rev1` can be checked against the C# output for one workbook, and the binning, convolution indexing, percentage-point lookup, and multi-mode loop stop being assumptions.
2. **Send X5 to X7 to the code's owner.** X6 is a one-line fix; X7 needs a decision on which step the tool should use.
3. **Use a finer step in the RFS tool's MACD.** The walkthrough's settings, a step of 1/100 of the beamwidth with 0.01 dB bins, bring every percentile within 0.03 dB of the fine answer and cut the bearing's effect from 3 dB to under 0.1 dB. For a full turn that's about 39,000 WSR-88D pointings and 25,000 ASR-11 pointings, and `macd_rev1` runs it in well under a second. A step of 1/40 of the beamwidth with 0.1 dB bins isn't enough: the bearing still moves the 99% value by 0.3 dB.

## 12. Verification

- **Hand values.** `test_macd_rev1` checks the C# pointing count and first pointing; a box pattern's histograms, MAG histogram, and percentage points worked by hand, including the twice-counted scan start (X6); mode averaging; and the lowest bin from a negative backlobe gain. `test_macd_monte_carlo_rev1` checks the Monte Carlo's main-beam shares for a box pattern against 1/180 and 2(1/180)(179/180), and a two-mode antenna's against a quarter of that.
- **Agreement with the exact calculation.** At the fine reference settings, `macd_rev1` reproduces the exact percentiles of the methodology report's Section 1.1 within 0.05 dB (`testFineSettingsMatchExact`); with 10<sup>8</sup> draws, the Monte Carlo's 95% intervals hold the exact 99.9% and 99.99% values (`testMatchesExactDistribution`).
- **Parallel runs.** The Monte Carlo gives identical counts with and without a parallel pool (`testParallelEqualsSerial`).
- **Tests.** All 255 tests in the 15 test files pass in one batch job with a six-worker pool in MATLAB R2025b on November 20, 2026, including the four new tests in `test_macd_monte_carlo_rev1.m`.
- **Numbers.** Every number in this paper was checked against `docs/figures/macd_pct_rev2.csv`, `docs/figures/macd_mc_convergence_rev2.csv`, or a rerun of `macd_rev1`, as listed in Appendix B.

## Appendix A. Code

All code is MATLAB, in the repository root. `report_figures_macd_rev2` reproduces every picture and table of Sections 6 to 9 in about 20 minutes on six workers.

| File | What it does |
| --- | --- |
| `macd_rev1.m` | MATLAB copy of the C# MACD steps: pointings at a set fraction of the beamwidth, binned gains, mode averaging, convolution, and percentage points |
| `macd_monte_carlo_rev1.m` | Monte Carlo of random, independent pointings, counted in 0.001 dB bins, parallel over chunks with one random-number substream each |
| `report_figures_macd_rev2.m` | Runs both for the two cases and writes Figures 2 to 9 and both tables |
| `report_figures_macd_rev1.m` | Writes Figure 1 and the C#-settings tables of revision 17 of the methodology report, including the step-size study (`macd_resolution_rev1.csv`) |
| `test_macd_rev1.m`, `test_macd_monte_carlo_rev1.m` | The tests of Section 12 |

## Appendix B. Numbers ledger

| Quantity | Exact value | Source | Form used |
| --- | --- | --- | --- |
| Pointings per turn: C#; walkthrough; reference (WSR-88D, ASR-11) | 1,557, 994; 38,919, 24,828; 155,676, 99,311 | `macd_rev1` | Same |
| Main-beam pointings above the sidelobe: C#; walkthrough | 12, 12; 292, 286 | `macd_rev1` | Same |
| Far sidelobe share: WSR-88D; ASR-11; both | 0.94450; 0.98026; 0.92586 | `macd_rev1` (walkthrough) | 94.45%; 98.03%; 92.6% |
| Share above the WSR-88D sidelobe: one cut; VCP 212 | 0.00750; 0.00511 | `macd_rev1` (walkthrough) | 0.75%; 0.51% |
| Step 2 example: pointings in the 41.70 dBi bin; in the C# 41 dBi bin | 2 of 38,919; 2 of 1,557 | `macd_rev1` | 0.005%; 0.13% |
| Step 5 example, 99%: share below the 34.94 dBi bin; its share; interpolated | 0.9899726; 3.805 × 10<sup>−5</sup>; 34.9472 | `macd_rev1` (walkthrough) | 98.99726%; 0.0038%; 34.947 dBi |
| Step 5 example, 99.999%: bin; interpolated | 72.31; 72.3144 | `macd_rev1` (walkthrough) | 72.31; 72.314 dBi |
| Peak gains on the path: 0.5°, 0.9°, 1.3° cuts | 41.994, 34.136, 21.796 dBi | 45.5 − 12(tilt/0.925)² | 42.0, 34.1, 21.8 dBi |
| Table 3 values | as tabulated | `macd_pct_rev2.csv` | Table 3 |
| Monte Carlo 99.999% interval width: 10<sup>5</sup>; 10<sup>7</sup>; 10<sup>10</sup> draws (one cut) | 7.67; 1.24; 0.045 dB | `macd_mc_convergence_rev2.csv` | 7 dB; 1.2 dB; 0.05 dB |
| Largest Monte Carlo to reference gap (VCP 212, 99.999%, inside the interval); largest outside the interval (one cut, 99%) | 0.023; 0.013 dB | `macd_pct_rev2.csv` | 0.025 dB; 0.013 dB |
| Bearing spread, 0° to 2°, interpolated: walkthrough settings (99%, 99.999%); beamwidth/40 with 0.1 dB bins (99%, 99.999%) | 0.098, 0.023; 0.302, 0.156 dB | `macd_rev1` rerun | under 0.1 dB; 0.3 dB |
| Exact 99% value, one cut (methodology report, Section 1.1) | 34.926 dBi | `simple_mag_nexrad_asr11_rev1` | 34.926 dBi |
| Table 4 values | as tabulated | `macd_pct_rev2.csv` | Table 4 |
| Bearing spread, 0° to 2°, C# settings: 99.999%; 99% (upper edges) | 71 to 74; 34 to 36 | `macd_rev1` rerun | Same |

**DRAFT//PRE-DECISIONAL**
