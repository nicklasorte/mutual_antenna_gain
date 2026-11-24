**DRAFT//PRE-DECISIONAL**

# Radar antenna patterns for 2700-2900 MHz: source research

Research report, revision 1, November 24, 2026. This report answers the research prompt in `antenna_patterns_next_steps_plan_rev1.md` (Steps 1 to 4). It uses public sources only. Every table was read from the PDF page image, not from extracted text. Values marked *derived* show their arithmetic. Values marked *graph reading* come from a plotted figure and carry lower confidence.

## Summary

- **The core values verify.** Gain, beamwidth, and sidelobe values for the WSR-88D, ASR-11 (Radar C), and ASR-9 (Radar B) match their primary sources, and so do the M.1851-2 cosine-taper constants and the statgain formulas.
- **The WSR-88D sidelobe specification has a second, more detailed public form.** Federal Meteorological Handbook No. 11 Part B gives −27 dB for the first sidelobe, a line from −29 dB at 2° to −34 dB at 10°, and −40 dB from 10° to 180°. The model's flat −29 dB is 2 dB below this bound between 1.44° and 2°, and up to 5 dB above it near 10°.
- **The WSR-88D back lobe is no longer open.** The same handbook states the −40 dB bound out to 180° and plots a typical measured pattern over ±180°, with far sidelobes near −55 to −70 dB and spillover lobes near −50 dB at about ±105° (graph reading).
- **ASR tilt has public values.** MIT Lincoln Laboratory ATC-138 gives a typical ASR-9 tilt of 2.0° and a low-beam two-way gain 6 dB down on the horizon, which is 3 dB one way. ATC-266 Figure 2 shows the same for the ASR-11. The model's 0° tilt is about 3 dB above that operational case.
- **The ASR-9 far sidelobe level stays open.** No public source gives a measured remote sidelobe level. The NTIA Manual standard for fixed radars in this band is a median gain of −10 dBi or less in the horizontal plane, as quoted in NTIA TM-13-489.
- **ASR aperture size stays open.** No public primary source gives it. A derived estimate from the beamwidth gives a far-field distance of about 480 to 600 m.
- **NTIA has published representative patterns for this band.** TM-13-489 Table 6-123 (meteorological radar, 46 dBi) and Table 6-95 (aeronautical radionavigation, 34 dBi) use the −12(θ/bw)² main beam and statgain sidelobes.
- **Precedent favors M.1851.** NTIA Report 15-517 (3.5 GHz) used M.1851 cosine patterns with 0° tilt in an aggregate analysis. NTIA Report 13-490 (this band) used main-beam gains of 45 dBi and 32 dBi. Report ITU-R M.2112 used the handbook's worst-case envelope for the meteorological radar.
- **M.1851-2 assigns the peak mask to a single interferer and the average mask to multiple interferers.** It also says to use a real, digitized pattern when one is available.
- **Four citation details need correcting** in the working paper: the titles of M.1464-2, M.1851-2, and ATC-266, and the source of the parabolic main-beam form.

## Table A: verification of the current model

Confidence is high unless the note says otherwise.

| # | Value in the model | Source, edition | Table or section, page | Short quote | Result |
| --- | --- | --- | --- | --- | --- |
| A1 | WSR-88D gain 45.5 dBi at 2850 MHz, with radome | ROC, NEXRAD Technical Information | Antenna characteristics, p. 1 | "45.5 dB (including radome loss)" | Match |
| A2 | WSR-88D beamwidth 0.925° | ROC, NEXRAD Technical Information | Antenna characteristics, p. 1 | "0.925 deg (Pencil Beam)" | Match. FMH-11 Part B gives 0.93°; M.1849-3 gives 0.92°; NTIA Report 13-490 gives 0.95° |
| A3 | WSR-88D sidelobes −29 dB, −40 dB beyond 10° | ROC, NEXRAD Technical Information | Antenna characteristics, p. 1 | "29 dBc beyond main beam" | Match to this source. FMH-11 Part B gives a different envelope; see B3 |
| A4 | WSR-88D reflector 8.53 m | ROC, NEXRAD Technical Information | Antenna characteristics, p. 1 | "8.53 m/28 ft in diameter" | Match |
| A5 | NSSL measured sidelobes: −35 dB at 2°, −48 dB at 6°, −55 dB at 20° | NSSL, errata to the April 1998 polarimetric upgrade report (updated March 13, 2017) | Section II, pp. 6-8 | "about -35 dB at 2°" | Match. Measured by Andrew Canada on an antenna without a radome. The original URL returned HTTP 410 to a direct download on November 24, 2026; the text was still readable through a web fetch. Medium confidence on page numbers |
| A6 | ASR-11 is Radar C: 34 dBi, 4.8° elevation, 1.45° azimuth | Rec. ITU-R M.1464-2 (02/2015) | Table 1, p. 5 | "34"; "4.8"; "1.45" | Match |
| A7 | Radar C sidelobes +9.5 dBi first, +3.5 dBi remote | Rec. ITU-R M.1464-2 | Table 1, p. 5 | "+9.5 3.5" | Match |
| A8 | Radar C cosecant-squared to 30° | Rec. ITU-R M.1464-2 | Table 1, p. 4 | "Cosecant-squared 6 to +30" | Partial. The table says 6° to 30° for Radars C, D, and E. The model starts the skirt at tilt + 2.4° |
| A9 | Radar C tilt range −2.5° to +2.5° | Rec. ITU-R M.1464-2 | Table 1, p. 5 | "+2.5 to –2.5" | Match |
| A10 | ASR-9 is Radar B: 33.5 dBi, 4.8°, 1.3° | Rec. ITU-R M.1464-2 | Table 1, p. 5 | "33.5"; "4.8"; "1.3" | Match |
| A11 | Radar B first sidelobe +7.3 dBi, no remote level | Rec. ITU-R M.1464-2 | Table 1, p. 5 | "+7.3" | Match. The remote entry is blank |
| A12 | Radar B is the ASR-9 | Rec. ITU-R M.1464-2 | Annex 2, Sections 2.2 and 2.3, pp. 12-13 | "high- and low-beam horns" | Consistent, not named. The description (60 NM range, 256 azimuth intervals of 1.4°, ten pulses then eight) matches ATC-138. NTIA Report 13-490 Table C-1 gives ASR-9 values of 1.32 MW and 1.05 µs against Radar B's 1 320 kW and 1.03 µs |
| A13 | Radar C is the ASR-11 | Rec. ITU-R M.1464-2; NTIA Report 13-490 | Table 1, p. 4; Table C-1, p. 81 | "25 kW"; "89 µs" | Consistent, not named. Radar C (25 kW, 1.0 and 89 µs, solid state) matches the FAA-provided ASR-11 values |
| A14 | ASR-9 maximum azimuth sidelobes −28 dB | MIT Lincoln Laboratory ATC-138 (Weber, 1986) | Section II.A, p. 3 | "maximum sidelobes -28 dB down" | Match. The same page gives 1.4° azimuth beamwidth; Table II-1 (p. 4) gives 34 dB gain |
| A15 | ASR-11 beams like the ASR-9: 1.4°, 5°, about 34 dB | MIT Lincoln Laboratory ATC-266 (Weber, July 12, 1999) | Section 3.1, p. 9 | "essentially identical" | Match. The report's title differs from the working paper's citation; see [Citation corrections](#citation-corrections) |
| A16 | M.1851-2 cosine taper K = 68.8 | Rec. ITU-R M.1851-2 (12/2023) | Table 4, Eq. (10), p. 8 | "68.8" | Match. The table gives −23 dB as this taper's first sidelobe |
| A17 | Peak mask −17.51 ln(2.33 θ/θ3), breakpoint −14.4 dB | Rec. ITU-R M.1851-2 | Table 6, Eq. (18), p. 13 | "−14.4" | Match |
| A18 | Average mask offset −4.32 dB, breakpoint −20.6 dB, floor −50 dB | Rec. ITU-R M.1851-2 | Table 6, Eq. (18), p. 13 | "−4.32"; "−20.6"; "−50" | Match. Section 2.1.3 also says in general terms that the average mask is the peak mask "minus 4 dB" |
| A19 | Taper chosen from the first sidelobe: 20 to under 30 dB gives n = 1 | Rec. ITU-R M.1851-2 | Section 2.4, Table 9, p. 22 | "20 to < 30" | Match |
| A20 | Cosecant-squared pattern, Eqs. (22) and (25) to (27) | Rec. ITU-R M.1851-2 | Section 2.2, Table 7, pp. 15-16 | "Cosecant floor level" | Match on Eq. (25) (floor, example −55 dB). The other equation numbers were not re-read line by line; medium confidence |
| A21 | Statgain formulas | NTIA TM-13-489 (Wang and Keech, October 2012) | Section 4.4.2, Table 4-1, p. 4-3 | "0.75 × Gmax - 7" | Match, all three gain categories |
| A22 | F.1245 | Rec. ITU-R F.1245-3 (01/2019) | recommends 2 and 3 | "average and related radiation patterns" | Match to the Recommendation. Its title says it models average patterns for point-to-point fixed links |
| A23 | Parabolic main beam −12(θ/bw)², cited to M.1849-3 Section 7.3 and M.1851-2 | Rec. ITU-R M.1849-3 (02/2023) | Section 7.3, p. 20 | "no defined radar antenna" | Mismatch in the citation. Section 7.3 gives no main-beam equation. The form appears in NTIA TM-13-489 Table 6-123 (p. 6-185) and NTIA Report 15-517 Eq. (4) (p. 7) |
| A24 | M.1849-3 Radar 1: 45 dBi, 0.92° | Rec. ITU-R M.1849-3 | Annex 2, Table 6, p. 42 | "45.0"; "0.92" | Match. The sidelobe entry reads "20" under a dB heading; TM-13-489 Table 5-67 reads it as a 20 dBi first sidelobe |

## Table B: values for the gaps

| Gap | Finding | Source | Confidence |
| --- | --- | --- | --- |
| B1. ASR-9 remote sidelobe level | Not public as a measured value. Three public bounds exist: (1) −28 dB maximum azimuth sidelobe (ATC-138, p. 3, "maximum sidelobes -28 dB down"); (2) "At least 25 dB below" the main beam for the ASR-8, ASR-9, and ASR-11, as provided by the FAA (NTIA Report 13-490, Table C-1, p. 81); (3) the NTIA Manual Section 5.5 standard for fixed radars in 2700-2900 MHz with 360° rotation: "median gain" of −10 dBi or less in the principal horizontal plane (quoted in NTIA TM-13-489, Table 5-41, p. 5-59). The FAA holds the measured patterns | ATC-138; NTIA Report 13-490; NTIA TM-13-489 | High for the bounds. The current NTIA Manual text was not re-read |
| B2. ASR tilt | ASR-9: "typical antenna tilt of 2.0 degrees," which puts the low-beam maximum at 2.0° and the high-beam maximum at 6.5° (ATC-138, Section II.A, p. 3). At that tilt, the low-beam two-way gain is "about" 6 dB down on the horizon and the high beam 21 dB down (ATC-138, ground clutter section). *Derived:* one-way low-beam horizon gain is −6/2 = −3 dB from the peak. ASR-11: ATC-266 Figure 2 (p. 9) shows the low-beam peak near 2° and about −3 dB at 0° (graph reading). The FAA siting handbook measures tilt "to nose of beam" and its worked example sets 2.9° as the maximum acceptable tilt for an ASR-8 site (FAA Order 6310.6, 1976, with Change 1). NTIA says ASR antennas "are slightly angled upward" (Spectrum Compendium, 2700-2900 MHz, February 2021, p. 6). Site-by-site tilts aren't public; the FAA and DoD hold them | ATC-138; ATC-266; FAA Order 6310.6; NTIA Compendium | High for the ATC-138 text; medium for graph readings |
| B3. WSR-88D far sidelobes and back lobe | Specification: first sidelobe −27 dB or lower; between 2° and 10°, below a line from −29 dB at 2° to −34 dB at 10°; between 10° and 180°, −40 dB or lower (FMH-11 Part B, December 2005, Section 3.6, p. 3-28, "less than or equal to -27 dB"). The handbook adds that the actual pattern is "about 5 dB below" the envelope beyond 2°. Typical measured pattern over ±180° (Figure 3-18, p. 3-30): far sidelobes about −55 to −70 dB, spillover lobes about −50 dB near ±105°, about −60 dB at 180° (graph reading). NSSL: nearly all measured sidelobes beyond 20° are below −55 dB | FMH-11 Part B; NSSL errata | High for the specification; medium for graph readings |
| B4. ASR aperture size | Not found in a public primary source. M.1464-2 gives only "Parabolic reflector" and an 8 m antenna height. *Derived estimate* from M.1851-2 Table 4 (θ3 = 68.8 λ/l, λ = 0.1071 m at 2800 MHz): l = 68.8 × 0.1071/1.45 = 5.08 m for the ASR-11 and 68.8 × 0.1071/1.3 = 5.67 m for the ASR-9. Far field 2l²/λ = 482 m and 600 m. The FAA and the manufacturers hold the dimensions | Derived | Low |

One more value bears on B2 for the WSR-88D. The ROC document says the scan runs from 0.5° to 19.5°, and "0.2 to 19.5 degrees at a single coastal location." NTIA Report 13-490 (footnote 16) says NEXRAD antennas can be tilted to 0° and below, "but such tilt angles are not commonly used operationally."

## Table C: radar systems in 2700-2900 MHz

NTIA's Spectrum Compendium for this band (February 2021) lists ASR-8, ASR-9, and ASR-11 airport surveillance radars and the NEXRAD network. It counts ASR systems at more than 250 FAA airports, about 150 DoD ASR systems, and 159 NEXRAD sites. Footnote G15 also allows military air defense radiolocation in the band.

| System | Where listed | Antenna type | Gain (dBi) | Beamwidth, az / el | Sidelobes | Elevation pattern | Scan | Best-fit M.1851-2 model |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Radar A (consistent with ASR-8: klystron, 0.6 µs) | M.1464-2 Table 1 | Parabolic reflector | 33.5 | 1.35° / 4.8° | Not given | Cosecant-squared to +30° | 75°/s | Cosine taper plus cosecant-squared |
| Radar B (consistent with ASR-9) | M.1464-2 Table 1 | Parabolic reflector | 33.5 | 1.3° / 4.8° | +7.3 dBi first | Cosecant-squared to +30° | 75°/s | Cosine taper plus cosecant-squared |
| Radar C (consistent with ASR-11; the DoD version is the AN/GPN-30, not verified here) | M.1464-2 Table 1 | Parabolic reflector | 34 | 1.45° / 4.8° | +9.5 dBi first, +3.5 dBi remote | Cosecant-squared 6° to +30° | 75°/s | Cosine taper plus cosecant-squared |
| Radar D (magnetron, 450 kW) | M.1464-2 Table 1 | Parabolic reflector | 32.8 | 1.6° / 4° | Not given | Cosecant-squared 6° to +30° | 90°/s | Cosine taper plus cosecant-squared; sidelobe level needed |
| Radar E (solid state, 22 kW) | M.1464-2 Table 1 | Parabolic reflector | 34.3 low beam, 33 high beam | 1.4° / 4.8° | Not given | Cosecant-squared 6° to +30° | 75°/s | Cosine taper plus cosecant-squared; sidelobe level needed |
| Radar F (TWT, 70 kW) | M.1464-2 Table 1 | Parabolic reflector | 33.5 | 1.5° / 5.0° | +7.5 dBi first, 0 to −3 dBi remote | Cosecant-squared, enhanced to +40° | 90°/s or 60°/s | Cosine taper plus cosecant-squared |
| Radar F1 (same as M.1849-3 Radar 3) | M.1464-2 Table 1; M.1849-3 Table 6 | Phased array, four faces, 4 m diameter | 41 | 1.6° to 2.7° both planes | 17 dB transmit, 25 dB receive | Pencil beam to 70 000 ft | Irregular, covers 360° | Planar array (Section 6) |
| Radar F2 (same as M.1849-3 Radar 4) | M.1464-2 Table 1; M.1849-3 Table 6 | Phased array, four faces, 8 m diameter | 46 | 0.9° to 1.4° / 0.9° to 1.5° | 17 dB transmit, 25 dB receive | Pencil beam to 100 000 ft | Irregular, covers 360° | Planar array (Section 6) |
| Radars I to M (radiolocation, 2700-3400 MHz) | M.1464-2 Table 2 | Shaped reflector, planar array, or active array | 33.5 to more than 40 | 1.1° to 5° / 1° to 30° | Radar I: 26 dB first, 35 dB remote | Cosecant-squared or pencil | 30 to 360°/s | By type; Radar I fits the ASR model |
| Radar 1 (WSR-88D) | M.1849-3 Table 6 | Parabolic reflector | 45.0 | 0.92° / 0.92° | "20" (see A24) | Pencil; steps from 0.5° to 20° | 4 to 30°/s | Circular aperture (Section 4) or the ROC envelope |
| Radar 2 (non-Doppler, 2700-2900 MHz) | M.1849-3 Table 6 | Parabolic reflector | 38.0 | 2.0° / 2.0° | +15 (estimated) | Pencil; −2° to +60° | 18°/s | Circular aperture |
| Radars 5 and 6 | M.1849-3 Table 6 | Parabolic reflector | Up to 46 | 0.9° to 1.25° both planes | −26 dB | Pencil; volume scan | 0 to 36°/s | Circular aperture |

The mapping of Radars D, E, and F to named products is not in the Recommendation, and this report doesn't assign one. The Recommendation describes Radars F1 and F2 as multipurpose systems that can replace a meteorological radar (Table 1, footnote 7). Planned U.S. programs were not verified in this pass; see [Claims not verified](#claims-not-verified).

## Table D: precedent studies

| Study | Band | Radar | Radar antenna pattern used | Single entry or aggregate | Reason the report gives |
| --- | --- | --- | --- | --- | --- |
| NTIA Report 13-490, Sanders and others, 2012 | 2700-2900 MHz (adjacent-band WiMAX) | WSR-88D, ASR-8, ASR-9, ASR-11 | Main-beam gain only: 45 dBi for NEXRAD and 32 dBi for the ASRs (Table D-2, p. 107) | Single entry, ITM run in reverse for separation distance | "mainbeam coupling for worst case" |
| NTIA Report 15-517 (revision 1), Drocella and others, 2015 | 3550-3650 MHz | Shipborne and ground-based radars | M.1851 cosine patterns in azimuth and elevation; 0° up-tilt; −12(θ/BW)² for the shipborne off-axis gain (pp. 6-7, Section 4.5, p. 11) | Aggregate Monte Carlo over access points | "ability to track down to the horizon." The report lists both the peak and average mask equations and doesn't say which one it applied |
| Report ITU-R M.2112 (2007) | 2700-2900 MHz, IMT | M.1464 Radars A, C, and G (meteorological) | Meteorological radar: FMH-11 Part B "worst case envelope" (p. 36; the Report writes "§ 3.28," which matches page 3-28 of the handbook). ASR patterns shown as plots (Figures 17 and 19) | Single and aggregate base stations; radar rotation searched for the worst angle | Worst-case analysis tied to the safety function; polarization loss set to zero |
| NTIA TM-13-489, Wang and Keech, 2012 | All federal bands | Meteorological radar and aeronautical radionavigation radar in this band | Representative patterns: 46 − 12(θ/bw)², 27.4 dBi to 1.25°, 30 − 25 log θ to 48°, −12 dBi beyond (Table 6-123, p. 6-185); 34 − 12(θ/bw)², 18.5 dBi to 5°, 36 − 25 log θ to 48°, −6 dBi beyond (Table 6-95, p. 6-164) | For spectrum sharing analyses in general | Built from Government Master File gain data: 160 assignments at 46 dBi and 155 at 34 dBi |
| Rec. ITU-R M.1849-3 (2023) | Meteorological radars | WSR-88D class | Says to use the actual pattern where possible; otherwise F.1245, F.699, M.1652, or M.1851 (Section 7.3, p. 20) | Not specified | Says the first three "tend to overestimate the beamwidth" of a pencil beam |

No FCC docket for 2700-2900 MHz, AWS-3 report, or 3450-3550 MHz report was read in this pass.

## M.1851-2 guidance on peak and average envelopes

- **Scope, p. 1:** the Recommendation defines both patterns, "peak, for single interferer," and average "for multiple interferers."
- **Section 2.1, p. 6:** "Peak or average masks are recommended" for studies "with single or multiple interferers respectively." The mask runs down to a floor that represents far sidelobes and back lobes.
- **Section 2.1, p. 6:** if real patterns are available, they "should be digitized and used."
- **Section 2.1.3, pp. 11-12:** the peak mask "rides over the side-lobe peaks" and is for a single-entry interferer. The average mask approximates the integral of the theoretical pattern and is for "aggregated multiple interferers."
- **Section 2.1.3, p. 12:** the peak-to-average difference in one principal plane is "approximately 4 dB." Table 6 gives the exact constant for each taper: −4.32 dB for the cosine taper.

M.1851-2 (12/2023) also adds a second model family "with pedestal" (Tables 2, 3, and 5), in which the beamwidth factor and mask constants vary continuously with the sidelobe level. For a first sidelobe of −24.5 dB or −26.2 dB, Table 3 places the antenna in its cosine-squared range (−40 to −22.7 dB). The working paper uses the model without pedestal (Tables 4, 6, and 9), which is also the model NTIA Report 15-517 printed.

## What these findings change in the model

These are observations for the study team, not recommendations:

1. **WSR-88D envelope.** The FMH-11 Part B envelope and the ROC summary differ. In dBi at 45.5 dBi peak, the handbook gives 18.5 dBi for the first sidelobe, 16.5 dBi at 2° falling to 11.5 dBi at 10°, and 5.5 dBi from 10° to 180°. The model holds 16.5 dBi from 1.44° to 10°.
2. **ASR horizon gain.** The public tilt data put the operational low-beam horizon gain about 3 dB below the peak. NTIA Report 13-490 used 32 dBi for ASR main-beam coupling.
3. **ASR floor.** The model's floors (+3.5 dBi and +5.5 dBi) are 13.5 dB and 15.5 dB above the NTIA Manual's −10 dBi median-gain standard. The floor is a peak-style bound; the standard is a median.
4. **M.1464-2 scope.** The Recommendation says it "should not be used for radar to radar analyses." That doesn't affect the 5G work, but it bears on the mutual antenna gain paper, which uses these values for radar-to-radar coupling.
5. **High beam.** The ASRs receive on the high beam at short range (inside about 15 nmi) and on the low beam beyond. M.1464-2 chose the low beam for its tests because interference from ground-based emitters most likely enters through it. The high beam is about 17 to 21 dB down on the horizon (graph reading and ATC-138 two-way value).

## Citation corrections

- **M.1464-2 title.** The 02/2015 edition is titled "Characteristics of non-meteorological radiolocation radars, and characteristics and protection criteria for sharing studies for aeronautical radionavigation and radars in the radiodetermination service operating in the frequency band 2 700-2 900 MHz." The working paper uses the earlier edition's title.
- **M.1851-2 title and date.** The edition is 12/2023, titled "Mathematical models for radiodetermination radar and aeronautical mobile systems antenna patterns for use in interference analyses."
- **ATC-266 title.** The report at the cited URL is "Wind Shear Detection Using the Next Generation Airport Surveillance Radar (ASR-11)," M. E. Weber, July 12, 1999. The Section 3.1 content cited is correct.
- **Main-beam form.** Cite NTIA TM-13-489 Table 6-123 and NTIA Report 15-517 Eq. (4) for −12(θ/bw)², not M.1849-3 Section 7.3.
- **Radar B description.** M.1464-2 Annex 2 Section 2.3.1 does describe the high and low beam horns, so the working paper can cite it directly.

## Claims not verified

- The names of real products for M.1464-2 Radars A, D, E, and F. Radar A's values are consistent with the ASR-8, but the Recommendation doesn't name it.
- That the AN/GPN-30 is the DoD designation of the ASR-11, and which other military radars (for example, AN/GPN-27) remain in service in this band.
- Planned systems: NOAA's next-generation weather radar program and FAA surveillance radar replacement plans, and whether either would use a phased array like Radars F1 and F2.
- The current text of NTIA Manual Section 5.5. This report relies on TM-13-489's 2012 quotation.
- Radar patterns used in AWS-3 studies, the 3450-3550 MHz studies, and any FCC docket on this band.
- How many WSR-88D sites scan below 0.5°. The ROC document says one coastal site scans at 0.2°.
- FAA Order 6310.6 page numbers and whether its tilt guidance applies to the ASR-9 and ASR-11. The scanned text was readable only in part.
- The Forecast International ASR-11 values (gain above 33.5 dB, 1.40°) and the Wikipedia ASR-9 values; neither was re-read against a primary source.

## Sources

- NOAA Radar Operations Center, "[NEXRAD Technical Information](https://roc.noaa.gov/public-documents/wsr88d/NEXRAD-Technical-Information.pdf)."
- Office of the Federal Coordinator for Meteorology, "[Federal Meteorological Handbook No. 11, Doppler Radar Meteorological Observations, Part B, Doppler Radar Theory and Meteorology](https://www.icams-portal.gov/resources/ofcm/fmh/FMH11/fmh-11B-2005.pdf)," FCM-H11B-2005, December 2005.
- NSSL, "[Errata for the report 'Polarimetric Upgrades to Improve Rainfall Measurements'](https://www.nssl.noaa.gov/publications/wsr88d_reports/ERRATA1998Rpt_D3a(WC).pdf)," April 1998, errata updated March 13, 2017.
- Recommendation ITU-R [M.1464-2](https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1464-2-201502-I!!PDF-E.pdf) (02/2015).
- Recommendation ITU-R [M.1849-3](https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1849-3-202302-I!!PDF-E.pdf) (02/2023).
- Recommendation ITU-R [M.1851-2](https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1851-2-202312-I!!PDF-E.pdf) (12/2023).
- Recommendation ITU-R [F.1245-3](https://www.itu.int/dms_pubrec/itu-r/rec/f/R-REC-F.1245-3-201901-I!!PDF-E.pdf) (01/2019).
- Report ITU-R [M.2112](https://www.itu.int/dms_pub/itu-r/opb/rep/R-REP-M.2112-2007-PDF-E.pdf) (2007), "Compatibility/sharing of airport surveillance radars and meteorological radar with IMT systems within the 2 700-2 900 MHz band."
- C. W. Wang and T. Keech, "[Antenna Models for Electromagnetic Compatibility Analyses](https://www.its.ntia.gov/publications/download/13-489.pdf)," NTIA TM-13-489, October 2012.
- F. H. Sanders, R. L. Sole, J. E. Carroll, G. S. Secrest, and T. L. Allmon, "[Analysis and Resolution of RF Interference to Radars Operating in the Band 2700-2900 MHz from Broadband Communication Transmitters](https://www.ntia.gov/files/ntia/publications/13-490_1_.pdf)," NTIA Report 13-490, 2012.
- E. Drocella, J. Richards, R. Sole, F. Najmy, A. Lundy, and P. McKenna, "[3.5 GHz Exclusion Zone Analyses and Methodology](https://its.ntia.gov/publications/download/TR-15-517r1.pdf)," NTIA Report 15-517, 2015.
- NTIA, "[Spectrum Compendium, 2700-2900 MHz](https://www.ntia.doc.gov/files/ntia/publications/compendium/2700.00-2900.00-02092021.pdf)," February 2021.
- M. E. Weber, "[Assessment of ASR-9 Weather Channel Performance: Analysis and Simulation](https://archive.ll.mit.edu/mission/aviation/publications/publication-files/atc-reports/Weber_1986_ATC-138_WW-15318.pdf)," MIT Lincoln Laboratory ATC-138, 1986.
- M. E. Weber, "[Wind Shear Detection Using the Next Generation Airport Surveillance Radar (ASR-11)](https://archive.ll.mit.edu/mission/aviation/publications/publication-files/atc-reports/Weber_1999_ATC-266_WW-15318.pdf)," MIT Lincoln Laboratory ATC-266, July 12, 1999.
- FAA, "[Primary/Secondary Terminal Radar Siting Handbook](https://www.faa.gov/documentLibrary/media/Order/6310.6.pdf)," Order 6310.6, July 20, 1976, with Change 1.
- Wikipedia, "[ASR-9](https://en.wikipedia.org/wiki/ASR-9)" (checked only for antenna dimensions; none given).

**DRAFT//PRE-DECISIONAL**
