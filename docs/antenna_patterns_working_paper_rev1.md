**DRAFT//PRE-DECISIONAL**

# Worst-case horizontal antenna patterns for the WSR-88D, ASR-11, and ASR-9 radars

Working paper, revision 1, draft of November 21, 2026. Repository `nicklasorte/mutual_antenna_gain`, `main` at commit `b3aaa4e` plus this revision's changes. Every number and figure in this paper comes from the MATLAB code in the repository, run in MATLAB R2025b. The [Code](#code) section lists the files.

## Summary

- **Purpose.** Coordination zones that protect radars in the 2700-2900 MHz band from 5G base stations need each radar's antenna gain toward a station on the horizon. This paper defines that gain for three radars: the WSR-88D weather radar (NEXRAD), the ASR-11 airport surveillance radar, and the ASR-9 airport surveillance radar.
- **One function per radar.** `nexrad_worst_case_slice_rev1`, `asr11_worst_case_slice_rev1`, and `asr9_worst_case_slice_rev1` each return an N-by-2 array: beam azimuth from 0° to 360° in the first column, and gain toward the station in dBi in the second. Each function is self-contained, so you can copy it into another project.
- **Worst case at the horizon.** Each slice puts the peak of the beam on the horizon, at 0° elevation. For these models, the gain toward a station at any elevation is at or below the slice, so the slice is an upper bound.
- **Public sources only.** The patterns come from the NOAA Radar Operations Center (ROC) specification, NSSL measurements, and Recommendations ITU-R M.1464-2, M.1849-3, and M.1851-2. The [Antenna model sources](#antenna-model-sources) section gives each citation.
- **Generic patterns don't bound the far sidelobes.** A fourth function, `statgain_worst_case_slice_rev1`, gives the NTIA statistical gain (statgain) pattern. Statgain and Rec. ITU-R F.1245 put the far sidelobes 9.5 to 18.1 dB below the radar-specific levels. Those sidelobes face a station for 94% to 98.5% of each turn.

## Who this paper is for

This paper is for engineers who size coordination zones or check link budgets for the 2700-2900 MHz band, and for staff who brief those results. It assumes that you know a radar antenna has a narrow main beam and weaker sidelobes, and that gains in decibels add.

## The worst-case horizontal slice

A coordination zone limits the interference that a 5G base station can put into a radar receiver. The radar's antenna gain toward the station is one term of that link budget. As the radar antenna turns, its gain toward a fixed station traces a curve over one turn. That curve is the *horizontal slice*.

Each function returns the slice for a station on the horizon, at 0° elevation:

- The first column is the beam azimuth relative to the bearing of the station, from 0° to 360°. At 0° and 360°, the beam points at the station. At 180°, it points away.
- The second column is the antenna gain toward the station, in dBi.

The default step is 0.1°, which gives 3,601 rows. You can set a finer or coarser step.

### Why the slice is the worst case

The slice is worst case in two ways:

- **The beam peak is on the horizon.** The default beam elevation (`tilt_deg = 0`) puts the peak of the beam at 0°, so the elevation pattern adds 0 dB toward the station. Any other tilt moves the station off the peak and lowers the gain at every azimuth (Figure 3b). The slice also bounds the gain toward a station above or below the horizon, because the elevation term never exceeds 0 dB in these models.
- **Each sidelobe level is an upper bound.** The WSR-88D slice uses the ROC sidelobe specification, which sits above the NSSL measurements. The ASR slices use the M.1851-2 *peak* envelope, which M.1851-2 gives for a single interferer. That envelope runs above the first sidelobe levels in M.1464-2 (see [ASR-11](#asr-11)).

Real radars don't sit at 0° tilt all the time. The WSR-88D's lowest operational cut is 0.5°, which puts 42.0 dBi on the horizon instead of 45.5 dBi. The `tilt_deg` option gives the slice at any tilt.

Figure 1 shows the three slices over a full turn. The main beam takes a small share of each turn. The WSR-88D slice is above its far sidelobe level for 5.6% of the turn, the ASR-11 slice for 2.0%, and the ASR-9 slice for 1.5%.

![Three stacked line plots of antenna gain in dBi against beam azimuth from 0 to 360 degrees. Top: WSR-88D, a peak of 45.5 dBi at 0 and 360 degrees, a shelf at 16.5 dBi out to 10 degrees on each side, and 5.5 dBi elsewhere. Middle: ASR-11, a peak of 34 dBi falling to a floor of 3.5 dBi within about 3.6 degrees. Bottom: ASR-9, a peak of 33.5 dBi falling to a floor of 5.5 dBi within about 2.8 degrees.](figures/fig_patterns_slices_rev1.png)

*Figure 1. Each radar's main beam faces the station for a few degrees of each turn; the far sidelobes face it for the rest. Worst-case horizontal slices at 0° tilt. Source: `report_figures_patterns_rev1.m` (2026).*

Figure 2 zooms in on the main beams.

![Line plot of gain in dBi against beam azimuth from -10 to 10 degrees for the three radars. The WSR-88D peaks at 45.5 dBi and drops to a flat 16.5 dBi shelf at about 1.4 degrees. The ASR-11 peaks at 34 dBi and falls smoothly to 3.5 dBi at about 3.6 degrees. The ASR-9, dashed, peaks at 33.5 dBi and falls to 5.5 dBi at about 2.8 degrees. A dotted line at 9.5 dBi marks the ASR-11 first sidelobe in M.1464-2; the ASR-11 curve stays above it out to about 2.5 degrees.](figures/fig_patterns_zoom_rev1.png)

*Figure 2. The WSR-88D main beam is the narrowest and highest; the ASR envelopes run above the ASR-11 first sidelobe level out to about 2.5°. Source: `report_figures_patterns_rev1.m` (2026).*

## WSR-88D (NEXRAD)

The WSR-88D uses a pencil beam from an 8.53 m parabolic reflector. The pattern is the same in every plane through the beam axis, so the gain depends only on the off-axis angle θ.

### Parameters

*Table 1. WSR-88D antenna parameters.*

| Parameter | Value | Source |
| --- | --- | --- |
| Peak gain | 45.5 dBi, including the radome | ROC_TECH (M1849 Radar 1: 45 dBi) |
| 3 dB beamwidth | 0.925° | ROC_TECH (M1849 Radar 1: 0.92°) |
| Main beam | Parabolic in dB, 3 dB down at half the beamwidth | M1851; M1849 Section 7.3 |
| Sidelobes, ROC specification | 29 dB below the peak out to 10°, 40 dB below beyond | ROC_TECH |
| Sidelobes, NSSL measurements | 35 dB below the peak to 2°, then linear to 48 dB below at 6° and 55 dB below at 20°; 55 dB below beyond | NSSL_1998 |
| Reflector diameter | 8.53 m (28 ft) | ROC_TECH |

Gain and beamwidth are quoted at 2850 MHz.

### Pattern

The main beam is

G(θ) = 45.5 − 12 (θ / 0.925)<sup>2</sup> dBi,

until it meets the sidelobe envelope. From that angle on, the gain follows the envelope. With the ROC specification, the main beam meets the −29 dB envelope at θ = 0.925 × (29/12)<sup>1/2</sup> = 1.438°. This step fills the first null, which is conservative.

On the horizon, the off-axis angle comes from the beam azimuth φ and the tilt t:

cos θ = cos t × cos φ.

At 0° tilt, θ equals the azimuth from the station.

### Options

`pattern="roc_spec"` (the default) uses the ROC specification, the upper bound. `pattern="nssl_measured"` uses the KOUN measurements, which are 2.6 to 21 dB lower in the sidelobes (Figure 6, top, and Figure 7). The NSSL report has no back-lobe data, so the function holds the 20° level to 180°.

### Range of validity

The pattern holds beyond the far-field distance 2D<sup>2</sup>/λ, about 1.38 km at 2850 MHz.

## ASR-11

The ASR-11 antenna is a reflector with a fan beam, narrow in azimuth and broad in elevation. A cosecant-squared skirt above the beam covers aircraft at high elevation. This paper models the low beam, which transmits and receives.

### Parameters

*Table 2. ASR-11 antenna parameters (low beam).*

| Parameter | Value | Source |
| --- | --- | --- |
| Peak gain | 34 dBi | M1464 Radar C (FI_ASR11: more than 33.5 dB; WEBER_1999: about 34 dB) |
| Azimuth 3 dB beamwidth | 1.45° | M1464 Radar C (FI_ASR11: 1.40°) |
| Elevation 3 dB beamwidth | 4.8° | M1464 Radar C (WEBER_1999: 5°) |
| First sidelobe | +9.5 dBi (24.5 dB below the peak) | M1464 Radar C |
| Remote sidelobes (floor) | +3.5 dBi (30.5 dB below the peak) | M1464 Radar C |
| Elevation pattern | Cosecant-squared to +30° | M1464 Radar C |
| Vertical tilt range | −2.5° to +2.5° | M1464 Radar C |

### Pattern

The gain combines an azimuth pattern and an elevation pattern in decibels (M.1851-2 Eq. (15)):

G(φ) = 34 + max[g<sub>az</sub>(φ) + g<sub>el</sub>(0°), −30.5] dBi.

**Azimuth pattern.** The −24.5 dB first sidelobe selects the cosine aperture taper (M.1851-2 Table 9, taper power n = 1). Near the peak, g<sub>az</sub> is the taper's pattern,

g<sub>az</sub> = 20 log<sub>10</sub> |F(μ)|, μ = π K sin φ / 1.45°, K = 68.8 (M.1851-2 Table 4),

where F(μ) = [S(μ + π/2) + S(μ − π/2)] / (2 F(0)) and S(b) = sin(b)/b. Past the point where the main beam is 14.4 dB down (1.416°), g<sub>az</sub> follows the M.1851-2 peak envelope (Table 6):

g<sub>az</sub> = −17.51 ln(2.33 |φ| / 1.45°).

The envelope reaches the −30.5 dB floor at 3.552°. With K = 68.8, the computed 3 dB beamwidth is 1.433°, 1.1% narrower than the 1.45° input. That difference comes from M.1851-2 and matches the repository model.

**Elevation pattern.** The beam peaks at the tilt angle. Below the peak, and up to the tilt plus half the elevation beamwidth (2.4°), the pattern is that of a uniform aperture. Above that, it follows the cosecant-squared law to 30° (M.1851-2 Eqs. (22) and (25) to (27)). At 0° tilt, g<sub>el</sub>(0°) = 0 dB. Figure 3a shows the pattern.

### The envelope sits above the first sidelobe

The peak envelope is smooth: it has no nulls, and it doesn't drop to the first sidelobe level between the main beam and the first sidelobe. It stays above the +9.5 dBi first sidelobe in M.1464-2 out to 2.52°. At 1.5°, the slice gives 18.6 dBi, 9.1 dB above that level. This is conservative.

## ASR-9

The ASR-9 antenna has the same fan-beam form as the ASR-11. WEBER_1999 describes the two beams as essentially identical. The function uses the same equations as the ASR-11 function, with Radar B values.

### Parameters

*Table 3. ASR-9 antenna parameters (low beam).*

| Parameter | Value | Source |
| --- | --- | --- |
| Peak gain | 33.5 dBi | M1464 Radar B (WEBER_1999: about 34 dB) |
| Azimuth 3 dB beamwidth | 1.3° | M1464 Radar B (WIKI_ASR9 and WEBER_1999: 1.4°) |
| Elevation 3 dB beamwidth | 4.8° | M1464 Radar B (WIKI_ASR9 and WEBER_1999: 5°) |
| First sidelobe | +7.3 dBi (26.2 dB below the peak) | M1464 Radar B |
| Remote sidelobes (floor) | +5.5 dBi (28 dB below the peak) | WEBER_1986 maximum azimuth sidelobe |
| Elevation pattern | Cosecant-squared to +30° | M1464 Radar B |

M.1464-2 gives no remote sidelobe level for Radar B. WEBER_1986 gives a maximum azimuth sidelobe of 28 dB below the peak, so the function uses that bound as the floor: 33.5 − 28 = 5.5 dBi. M.1464-2 doesn't name its radars. Radar B is the ASR-9 because M.1464-2 Annex 2 Section 2 describes it with the ASR-9's features, such as a klystron transmitter and two coherent processing intervals (WIKI_ASR9).

### Pattern

The −26.2 dB first sidelobe also selects the cosine taper, so the constants match the ASR-11's (K = 68.8, A = 17.51, B = 2.33). The envelope starts at 1.269° and reaches the −28 dB floor at 2.761°. As with the ASR-11, the envelope runs above the first sidelobe level, here out to 2.49°.

## Elevation and tilt

Figure 3a shows the ASR-11 gain at the beam azimuth against elevation. At 0° tilt, the beam peak is on the horizon. At 2.5° tilt, the horizon falls on the lower edge of the beam, and the horizon gain drops by 3.3 dB, to 30.7 dBi.

Figure 3b shows the peak gain toward the horizon against tilt for all three radars. The gain is highest at 0° tilt for each one. The WSR-88D's narrow beam makes it the most sensitive to tilt. Its horizon gain falls from 45.5 dBi at 0° to 42.0 dBi at its 0.5° cut, and to the 16.5 dBi sidelobe level by 1.44°.

![Two line plots. Left: ASR-11 gain in dBi against elevation from -6 to 40 degrees. At 0 degrees tilt, the curve peaks at 34 dBi at 0 degrees elevation and falls along the cosecant-squared skirt to about 9.5 dBi at 30 degrees, then drops to the 3.5 dBi floor. At 2.5 degrees tilt, dashed, the curve peaks at 2.5 degrees elevation and is about 30.7 dBi on the horizon. Right: peak gain toward the horizon against beam elevation from 0 to 2.5 degrees. The WSR-88D falls from 45.5 dBi to 16.5 dBi by about 1.44 degrees, with a marker at 42.0 dBi for the 0.5 degree cut. The ASR-11 and ASR-9 curves fall slowly, from 34 and 33.5 dBi to about 30.7 and 30.2 dBi.](figures/fig_patterns_elevation_rev1.png)

*Figure 3. Every radar has its highest horizon gain at 0° tilt, so the 0° slice bounds every tilt. a: ASR-11 elevation pattern at 0° and 2.5° tilt. b: Peak horizon gain against beam elevation. Source: `report_figures_patterns_rev1.m` (2026).*

Figures 4 and 5 show the whole main beam at several tilts. For the WSR-88D, Figure 4 uses the lowest cuts of the operational scan patterns. Each tilt lowers the peak on the horizon and narrows the part of the turn above the sidelobe level. At 1.3°, only a 21.8 dBi bump is left.

![Line plot of WSR-88D gain toward the station in dBi against beam azimuth from -5 to 5 degrees, for beam tilts of 0, 0.5, 0.9, and 1.3 degrees. The peaks on the horizon are 45.5, 42.0, 34.1, and 21.8 dBi. Every curve falls to the 16.5 dBi sidelobe shelf; the higher the tilt, the narrower the beam above the shelf.](figures/fig_patterns_nexrad_tilts_rev1.png)

*Figure 4. Each step up in WSR-88D tilt lowers the horizon peak, from 45.5 dBi at 0° to 21.8 dBi at 1.3°. Slices at 0°, 0.5°, 0.9°, and 1.3° tilt. Source: `report_figures_patterns_rev1.m` (2026).*

The ASR fan beams are broad in elevation, so tilt changes them much less (Figure 5). From 0° to 2.5° tilt, the ASR-11 and ASR-9 peaks each drop 3.3 dB. The floors don't change.

![Line plot of gain toward the station in dBi against beam azimuth from -6 to 6 degrees for the ASR-11 (blue) and ASR-9 (gray) at tilts of 0, 1, and 2.5 degrees. ASR-11 peaks are 34.0, 33.5, and 30.7 dBi; ASR-9 peaks are 33.5, 33.0, and 30.2 dBi. The ASR-11 curves fall to a 3.5 dBi floor and the ASR-9 curves to a 5.5 dBi floor.](figures/fig_patterns_asr_tilts_rev1.png)

*Figure 5. ASR tilt lowers the horizon peak by at most 3.3 dB, so the 0° slice is close to every operational case. Slices at 0°, 1°, and 2.5° tilt. Source: `report_figures_patterns_rev1.m` (2026).*

## Comparison with generic patterns

Two generic patterns are common in spectrum studies. Neither needs a beamwidth or sidelobe level:

- **Statgain** (NTIA TM-13-489, Section 4.4.2 and Table 4-1) needs only the peak gain. The function `statgain_worst_case_slice_rev1(g_max_dBi)` returns its slice in the same N-by-2 form.
- **Rec. ITU-R F.1245** is the pattern that the 2.7 GHz link budget working paper uses (WP_REV3 Table 4.4-1: 46 dBi for the WSR-88D and 34 dBi for the ASRs). The comparison uses the repository's F.1245 model.

Both patterns are rotationally symmetric, so they treat every antenna as a pencil beam. Their beamwidth follows from the peak gain alone: statgain implies 0.920° at 45.5 dBi, which matches the WSR-88D, but 3.456° at 34 dBi, 2.4 times the ASR-11's 1.45° azimuth beamwidth.

![Three line plots of gain in dBi against beam azimuth from 0.05 to 180 degrees on a log scale, one per radar. Each shows the radar-specific pattern in blue, statgain in gray, and F.1245 dashed in purple. WSR-88D: all three agree in the main beam; statgain and F.1245 sit above the ROC shelf from about 1.4 to 3.5 degrees, then fall below it, ending near -12 dBi against the 5.5 dBi specification; the NSSL measurements, dash-dot, sit below the specification. ASR-11 and ASR-9: statgain and F.1245 have main beams more than twice as wide, so they are far above the radar patterns from about 1 to 20 degrees, then fall below the radar floors, ending near -6 and -10 dBi.](figures/fig_patterns_compare_rev1.png)

*Figure 6. Generic patterns overstate the ASR beam near the peak and understate every radar's far sidelobes. Radar-specific slices against statgain and F.1245 at each radar's peak gain. Source: `report_figures_patterns_rev1.m` (2026).*

Figures 7 to 9 show the same comparison on a linear azimuth scale, one radar at a time, over the main beam. For the WSR-88D (Figure 7), all three patterns agree on the main beam. Outside it, statgain and F.1245 cross the ROC shelf near 3.5° and fall below it.

![Line plot of gain in dBi against beam azimuth from -10 to 10 degrees. The WSR-88D worst-case slice (blue) peaks at 45.5 dBi and drops to a 16.5 dBi shelf at about 1.4 degrees. The NSSL measurements (blue, dash-dot) fall to about 10.5 dBi at 2 degrees and -4.5 dBi at 10 degrees. Statgain at 45.5 dBi (gray) and F.1245 at 46 dBi (purple, dashed) match the main beam, run above the shelf from about 1.4 to 3.5 degrees, and fall to about 5 dBi at 10 degrees.](figures/fig_patterns_nexrad_compare_rev1.png)

*Figure 7. Statgain and F.1245 match the WSR-88D main beam, sit above the ROC shelf out to about 3.5°, then fall below it. Source: `report_figures_patterns_rev1.m` (2026).*

For the ASR-11 (Figure 8), the generic main beams are more than twice as wide as the radar's. The M.1851-2 average envelope, which M.1851-2 gives for aggregate interference, reaches the floor at 2.78° instead of 3.55°.

![Line plot of gain in dBi against beam azimuth from -10 to 10 degrees. The ASR-11 worst-case slice (blue) peaks at 34 dBi and falls to a 3.5 dBi floor at about 3.6 degrees. The M.1851-2 average envelope (dotted) follows the same main beam but reaches the floor at about 2.8 degrees. A dotted purple line marks the +9.5 dBi first sidelobe in M.1464-2. Statgain (gray) and F.1245 (purple, dashed) at 34 dBi have much wider main beams, about 18.5 dBi at 4 to 5 degrees and 7 to 11 dBi at 10 degrees.](figures/fig_patterns_asr11_compare_rev1.png)

*Figure 8. The generic patterns are far above the ASR-11 slice from 1° to 10°; the M.1851-2 average envelope is below it. Source: `report_figures_patterns_rev1.m` (2026).*

![Line plot of gain in dBi against beam azimuth from -10 to 10 degrees. The ASR-9 worst-case slice (blue) peaks at 33.5 dBi and falls to a 5.5 dBi floor at about 2.8 degrees. Statgain at 33.5 dBi (gray) and F.1245 at 34 dBi (purple, dashed) have much wider main beams and stay above the slice out to 10 degrees.](figures/fig_patterns_asr9_compare_rev1.png)

*Figure 9. The ASR-9 comparison matches the ASR-11's: the generic main beams are more than twice as wide. Source: `report_figures_patterns_rev1.m` (2026).*

Figure 10 shows the full turn. The generic patterns fall below the radar-specific levels beyond 14° to 20° for the ASRs and beyond 3.3° to 3.6° for the WSR-88D. They stay below them for the rest of the turn.

![Line plot of gain in dBi against beam azimuth from -180 to 180 degrees. The WSR-88D and ASR-9 slices sit at 5.5 dBi and the ASR-11 slice at 3.5 dBi outside their main beams, which are narrow spikes at 0 degrees. Statgain at 45.5 dBi (gray) falls to -11.75 dBi beyond 48 degrees, and statgain at 34 dBi (gray, dotted) to -6 dBi.](figures/fig_patterns_fullturn_compare_rev1.png)

*Figure 10. Over most of the turn, statgain is 9.5 to 17.3 dB below the radar-specific far sidelobes. The WSR-88D and ASR-9 lines overlap at 5.5 dBi. Source: `report_figures_patterns_rev1.m` (2026).*

Table 4 gives the gain at selected azimuths.

*Table 4. Gain toward a station on the horizon at 0° tilt, in dBi.*

| Azimuth | WSR-88D, ROC spec | WSR-88D, NSSL | Statgain 45.5 | F.1245 46 | ASR-11 | ASR-9 | Statgain 34 | F.1245 34 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0° | 45.5 | 45.5 | 45.5 | 46.0 | 34.0 | 33.5 | 34.0 | 34.0 |
| 0.5° | 42.0 | 42.0 | 42.0 | 41.8 | 32.6 | 31.7 | 33.7 | 33.7 |
| 1° | 31.5 | 31.5 | 31.3 | 29.4 | 27.8 | 25.6 | 33.0 | 32.9 |
| 1.5° | 16.5 | 13.9 | 25.8 | 25.0 | 18.6 | 16.2 | 31.7 | 31.6 |
| 2° | 16.5 | 10.5 | 22.7 | 21.9 | 13.6 | 11.1 | 30.0 | 29.7 |
| 3° | 16.5 | 7.3 | 18.3 | 17.5 | 6.5 | 5.5 | 25.0 | 24.4 |
| 5° | 16.5 | 0.8 | 12.8 | 12.0 | 3.5 | 5.5 | 18.5 | 15.0 |
| 10° | 16.5 | −4.5 | 5.3 | 4.4 | 3.5 | 5.5 | 11.0 | 7.4 |
| 20° | 5.5 | −9.5 | −2.3 | −3.1 | 3.5 | 5.5 | 3.5 | −0.1 |
| 90° | 5.5 | −9.5 | −11.8 | −12.6 | 3.5 | 5.5 | −6.0 | −9.6 |
| 180° | 5.5 | −9.5 | −11.8 | −12.6 | 3.5 | 5.5 | −6.0 | −9.6 |

Source: `patterns_compare_rev1.csv`. Statgain 34 is the comparison for the ASR-11; the ASR-9 comparison at 33.5 dBi is 0.25 dB higher in the far sidelobes.

### What the comparison means for coordination zones

- **Near the main beam, the generic patterns are higher for the ASRs.** Between 1° and about 20°, statgain is up to 18.5 dB above the ASR-11 slice (near 3°), because its main beam is more than twice as wide.
- **In the far sidelobes, the generic patterns are lower for every radar.** Beyond 48°, statgain is 17.3 dB below the WSR-88D specification, 9.5 dB below the ASR-11 floor, and 11.3 dB below the ASR-9 bound. F.1245 is 18.1 dB and 13.1 dB below the WSR-88D and ASR-11 levels.
- **The far sidelobes matter most for zone size.** They face a given station for most of each turn (94.4% for the WSR-88D, 98.0% for the ASR-11, and 98.5% for the ASR-9). A zone sized with a generic pattern's far sidelobes comes out smaller than one sized with the radar-specific bounds in this paper.

## Other antenna patterns for 2700-2900 MHz radars

You can model other radars in this band, or other forms of these radars, with the following public sources:

- **M.1851-2 average envelope.** M.1851-2 gives a peak envelope for a single interferer and an average envelope for aggregate interference from many sources. For the cosine taper, the average envelope is 4.32 dB below the peak envelope past its breakpoint. Aggregate interference from many 5G stations might call for the average envelope; that choice belongs to the study team. The repository's `antenna_prepare_rev1` supports both envelopes (`envelope="peak"` or `"average"`).
- **M.1851-2 circular aperture.** M.1851-2 Section 4 models a parabolic reflector with a circular aperture, which fits the WSR-88D dish. The repository supports it as antenna type `"parabolic"`.
- **Other radars in M.1464-2 Table 1.** The table lists further ATC radars (Radars A, D, E, and F) and two combined ATC and weather radars (Radars F1 and F2). The ASR functions take the same equations; each radar needs its own gain, beamwidths, and sidelobe levels read from the table.
- **Meteorological radars in M.1849-3.** M.1849-3 Annex 2 lists ground-based meteorological radars in this band and discusses their antenna patterns in Section 7.3. Radar 1 is the WSR-88D.
- **M.1851-2 planar array.** M.1851-2 Section 6 models a phased-array antenna. A phased-array radar in this band needs it. The repository supports it as antenna type `"array_planar"`.
- **Measured patterns.** Measured antenna patterns from the radar operators can replace the specification and envelope bounds. They aren't public. The repository accepts digitized cuts as antenna type `"measured"`.

## Use the functions

Each function takes name-value options and returns the slice and an `info` struct with every parameter, its source, and the citations.

```matlab
s_nex = nexrad_worst_case_slice_rev1();                    % ROC specification, 0.1 deg step
s_a11 = asr11_worst_case_slice_rev1(step_deg=0.01);        % 36,001 rows
s_a9  = asr9_worst_case_slice_rev1();
s_sg  = statgain_worst_case_slice_rev1(45.5);              % generic pattern at 45.5 dBi

[s, info] = nexrad_worst_case_slice_rev1(pattern="nssl_measured", tilt_deg=0.5);
disp(info.params)                                          % every value with its source
disp(info.sources)                                         % citations with URLs
```

The options are as follows:

- **`step_deg`**: azimuth step in degrees (default 0.1). 360 must be a whole number of steps.
- **`tilt_deg`**: beam elevation in degrees (default 0, the worst case). For the ASRs, the range is −2.5° to +2.5°, and the start of the cosecant-squared skirt (tilt + 2.4°) must stay above 0°.
- **`pattern`** (WSR-88D only): `"roc_spec"` (default) or `"nssl_measured"`.

To get the gain toward a station at bearing β when the beam points at azimuth α, look up the slice at mod(α − β, 360).

The same slices at a 0.1° step are in `docs/figures` as CSV files: `nexrad_worst_case_slice_rev1.csv`, `asr11_worst_case_slice_rev1.csv`, and `asr9_worst_case_slice_rev1.csv`.

## Limitations and open items

- **Horizon only.** The slices assume the station is on the horizon. Because the elevation term never exceeds 0 dB, the slices also bound stations above or below the horizon, but they overstate the gain there.
- **Bounds, not measurements.** The ROC specification and the M.1851-2 peak envelope are upper bounds. Real antennas have nulls and lower sidelobes. The NSSL option shows how far below the bound one measured WSR-88D antenna sits.
- **ASR-9 floor.** The +5.5 dBi floor is the WEBER_1986 maximum sidelobe, not a remote sidelobe level. An FAA value can replace it (open item).
- **ASR site tilt.** The ASR tilt at each site isn't public. The default of 0° bounds every tilt (open item).
- **WSR-88D back lobe.** Neither source gives a measured back lobe. The ROC specification holds −40 dB to 180°; the NSSL option holds −55 dB.
- **Near field.** The WSR-88D pattern holds beyond about 1.38 km. This paper doesn't give a far-field distance for the ASRs, because no public source gives their aperture size (open item).
- **Polarization.** The slices don't include polarization mismatch between the radar and the 5G antenna.

## Verification

- **Tests.** `test_worst_case_slices_rev1.m` has seven tests, and all pass. They check the output shape and symmetry, hand-worked values (peaks, the 3 dB points, envelope levels, breakpoints, floors, and the cosecant-squared law), and agreement with the repository model:
  - The ASR and statgain slices match `antenna_prepare_rev1` and `antenna_gain_rev1` within 10<sup>−9</sup> dB at 0°, 1°, and 2.5° tilt.
  - The WSR-88D slices match the repository's off-axis table within 10<sup>−3</sup> dB, except within 0.006° of the main-beam corner, where the table's 0.005° grid rounds the corner by up to 0.25 dB. The function computes the corner exactly.
- **Code analysis.** MATLAB `checkcode` reports no messages for the four functions and the figure script.
- **Figures and tables.** `report_figures_patterns_rev1.m` writes all 10 figures and every CSV file in this paper. Every value in Table 4 comes from `patterns_compare_rev1.csv`.

## Code

*Table 5. MATLAB files for this paper.*

| File | What it does |
| --- | --- |
| `nexrad_worst_case_slice_rev1.m` | WSR-88D slice, ROC specification or NSSL measurements, any tilt |
| `asr11_worst_case_slice_rev1.m` | ASR-11 slice, M.1851-2 cosine taper and cosecant-squared elevation pattern |
| `asr9_worst_case_slice_rev1.m` | ASR-9 slice, same model with Radar B values and the WEBER_1986 floor |
| `statgain_worst_case_slice_rev1.m` | Generic statgain slice for any peak gain of 10 dBi or more |
| `report_figures_patterns_rev1.m` | Writes Figures 1 to 10, the slice CSV files, `patterns_compare_rev1.csv`, and `patterns_params_rev1.csv` |
| `test_worst_case_slices_rev1.m` | The tests in [Verification](#verification) |

## Antenna model sources

- **FI_ASR11**: Forecast International, "[ASR-11 (DASR)](https://www.forecastinternational.com/archive/disp_pdf.cfm?DACH_RECNO=1032)," Radar Forecast, August 2014. Gain more than 33.5 dB, 1.40° azimuth beamwidth.
- **M1464**: Recommendation ITU-R M.1464-2 (02/2015), "[Characteristics of radiolocation radars, and characteristics and protection criteria for sharing studies for aeronautical radionavigation and meteorological radars in the radiodetermination service operating in the frequency band 2 700-2 900 MHz](https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1464-2-201502-I!!PDF-E.pdf)." Table 1, Radar B (ASR-9) and Radar C (ASR-11); Annex 2 Section 2.
- **M1849**: Recommendation ITU-R M.1849-3 (02/2023), "[Technical and operational aspects of ground-based meteorological radars](https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1849-3-202302-I!!PDF-E.pdf)." Annex 2 Table 6 (Radar 1) and Section 7.3.
- **M1851**: Recommendation ITU-R M.1851-2, "[Mathematical models for radiodetermination radar systems antenna patterns for use in interference analyses](https://www.itu.int/rec/R-REC-M.1851/en)." Tables 4, 6, 7, and 9; Eqs. (15), (22), and (25) to (27).
- **NSSL_1998**: NSSL, "[Errata and supplement to the 1998 WSR-88D polarimetric upgrade report (D3a)](https://www.nssl.noaa.gov/publications/wsr88d_reports/ERRATA1998Rpt_D3a(WC).pdf)," Section II. KOUN antenna beamwidth and sidelobe measurements.
- **ROC_TECH**: NOAA Radar Operations Center, "[NEXRAD Technical Information](https://roc.noaa.gov/public-documents/wsr88d/NEXRAD-Technical-Information.pdf)." WSR-88D antenna gain, beamwidth, sidelobe specification, and reflector size.
- **TM_13_489**: NTIA Technical Memorandum TM-13-489, "[Antenna Models for Electromagnetic Compatibility Analyses](https://www.ntia.gov/report/2012/antenna-models-electromagnetic-compatibility-analyses)," October 2012. Section 4.4.2 and Table 4-1 (statgain).
- **WEBER_1986**: M. E. Weber, "[Assessment of ASR-9 Weather Channel Performance: Analysis and Simulation](https://archive.ll.mit.edu/mission/aviation/publications/publication-files/atc-reports/Weber_1986_ATC-138_WW-15318.pdf)," MIT Lincoln Laboratory ATC-138, 1986, Section II.A. Maximum azimuth sidelobes 28 dB below the peak.
- **WEBER_1999**: M. E. Weber, "[FAA Surveillance Radar Data as a Complement to the WSR-88D Network](https://archive.ll.mit.edu/mission/aviation/publications/publication-files/atc-reports/Weber_1999_ATC-266_WW-15318.pdf)," MIT Lincoln Laboratory ATC-266, 1999, Section 3.1. ASR-9 and ASR-11 beams: 1.4° azimuth, 5° elevation, about 34 dB.
- **WIKI_ASR9**: Wikipedia, "[ASR-9](https://en.wikipedia.org/wiki/ASR-9)." 1.4° azimuth beamwidth, 5° elevation beamwidth, cosecant-squared high and low beams.
- **WP_REV3**: "DRAFT 2.7 GHz Link Budget Analysis," Rev3 redlines, October 27, 2026, Table 4.4-1 (F.1245 model). Project document, no public URL.

## Numbers ledger

*Table 6. Numbers ledger: each value in this paper, its exact value, and its source.*

| Quantity | Exact value | Source | Form used |
| --- | --- | --- | --- |
| WSR-88D main beam meets the ROC envelope | 1.4380° | `info.x_deg` | 1.438° |
| WSR-88D far-field distance | 1,383 m | `info.far_field_m` | about 1.38 km |
| WSR-88D horizon gain at the 0.5° cut | 41.994 dBi | 45.5 − 12(0.5/0.925)<sup>2</sup> | 42.0 dBi |
| ASR-11 breakpoint; floor reached | 1.4158°; 3.5522° | `info.brk_deg`; `info.floor_az_deg` | 1.416°; 3.552° |
| ASR-9 breakpoint; floor reached | 1.2693°; 2.7610° | `info.brk_deg`; `info.floor_az_deg` | 1.269°; 2.761° |
| ASR-11 computed 3 dB beamwidth | 1.4334° | `test_worst_case_slices_rev1` | 1.433°, 1.1% narrow |
| Envelope above the M.1464-2 first sidelobe to: ASR-11; ASR-9 | 2.5215°; 2.4913° | −17.51 ln(2.33φ/bw) = SLL | 2.52°; 2.49° |
| ASR-11 slice at 1.5° | 18.595 dBi | `patterns_compare_rev1.csv` | 18.6 dBi, 9.1 dB above +9.5 |
| ASR-11 horizon gain at 2.5° tilt | 30.71 dBi | `info.el_fn` at tilt 2.5° | 30.7 dBi, 3.3 dB lower |
| Share of the turn above the far level: WSR-88D; ASR-11; ASR-9 | 5.56%; 1.98%; 1.54% | `report_figures_patterns_rev1` | 5.6%; 2.0%; 1.5% |
| Statgain implied beamwidth: 45.5; 34; 33.5 dBi | 0.920°; 3.456°; 3.661° | `info.bw_implied_deg` | Same |
| Far sidelobe gap, statgain: WSR-88D; ASR-11; ASR-9 | 17.25; 9.5; 11.25 dB | Table 4 | 17.3; 9.5; 11.3 dB |
| Far sidelobe gap, F.1245: WSR-88D; ASR-11 | 18.075; 13.075 dB | Table 4 | 18.1; 13.1 dB |
| Largest statgain excess over the ASR-11 slice | 18.50 dB at 3° | Table 4 | 18.5 dB |
| WSR-88D horizon peak at 0.9° and 1.3° tilt | 34.140; 21.798 dBi | `nexrad_worst_case_slice_rev1` | 34.1; 21.8 dBi |
| ASR peak at 1° and 2.5° tilt: ASR-11; ASR-9 | 33.507, 30.708; 33.007, 30.208 dBi | `asr11_`, `asr9_worst_case_slice_rev1` | 33.5, 30.7; 33.0, 30.2 dBi |
| Generic pattern crossings: statgain WSR-88D, ASR-11, ASR-9; F.1245 WSR-88D, ASR | 3.55°, 19.95°, 17.0°; 3.28°, 14.4° | Solve each Table 4-1 or F.1245 sidelobe law for the radar level | 3.6°, 20°, 17°; 3.3°, 14° |
| ASR-11 average envelope reaches the floor | 2.776° | 1.45/2.33 exp((30.5 − 4.32)/17.51) | 2.78° |
| Table 4 values | as tabulated | `patterns_compare_rev1.csv` | Table 4 |

**DRAFT//PRE-DECISIONAL**
