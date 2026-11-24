**DRAFT//PRE-DECISIONAL**

# Antenna patterns: next steps plan

Plan, revision 1, November 24, 2026. Follows the antenna patterns working paper, revision 1 (`docs/antenna_patterns_working_paper_rev1.docx`, commit `091b6b0`). Work sources are limited to public and cleared material.

## Goal

The working paper gives a worst-case horizontal gain slice for the WSR-88D, ASR-11, and ASR-9. The slices feed coordination zones that protect these radars from 5G base stations. The next steps aim at one outcome: a radar antenna model that a coordination-zone study can defend, with each choice tied to a public source, a measured check, and a stated decision.

The plan has four parts:

- **Verify.** Confirm every value in the paper against the primary source, page by page.
- **Close the open items.** Find public values for the four gaps the paper lists.
- **Extend.** Add the other radars and patterns in 2700-2900 MHz that a zone study needs.
- **Decide.** Show how much each pattern choice moves a zone, so the choices can be made on evidence.

The research prompt at the end of this plan covers the source work in Steps 1 to 4. You can paste it into Claude Fable as is.

## What the paper leaves open

The working paper lists these gaps:

| # | Open item | Why it matters for a zone |
| --- | --- | --- |
| O1 | ASR-9 floor is the WEBER_1986 maximum sidelobe (−28 dB), not a remote sidelobe level | Sets the ASR-9 gain for about 98.5% of each turn |
| O2 | ASR tilt at each site isn't public; the paper uses 0° | Up to 3.3 dB on the horizon peak |
| O3 | No public WSR-88D back-lobe measurement; the ROC spec holds −40 dB to 180° | Sets the gain behind the radar |
| O4 | No public ASR aperture size, so no ASR far-field distance | Decides how close a station can be before the pattern fails |

Three questions came up in the paper that it doesn't answer:

- **Peak or average envelope.** M.1851-2 gives a peak envelope for one interferer and an average envelope for aggregate interference. 5G deployments are aggregate. The two differ by 4.32 dB past the breakpoint for the ASRs.
- **Specification or measurement for the WSR-88D.** The NSSL measurements sit 2.6 to 21 dB below the ROC specification in the sidelobes.
- **Whether a generic pattern is acceptable.** Statgain and F.1245 put the far sidelobes 9.5 to 18.1 dB below the radar-specific bounds.

## Steps

Each step names its deliverable and the outside check that shows it's right. Steps 1 to 4 are research. Steps 5 to 8 are analysis and code. Step 9 is review. Step 10 is the decision package.

### Step 1: Verify every source value

Re-read each primary source and record the exact value, table, and page for every number the four functions use.

- **Sources:** ROC NEXRAD Technical Information; NSSL 1998 errata, Section II; Recommendations ITU-R M.1464-2 Table 1 and Annex 2; M.1849-3 Annex 2 and Section 7.3; M.1851-2 Tables 4, 6, 7, and 9 and Eqs. (15), (22), and (25) to (27); NTIA TM-13-489 Table 4-1; Rec. ITU-R F.1245.
- **Deliverable:** a source ledger (`docs/antenna_patterns_source_ledger_rev1.csv`) with columns for value, unit, source, table or equation, page, quote of five words or fewer, and match status against the code.
- **Outside check:** reproduce one published figure from M.1851-2 (the cosine-taper peak and average envelopes) with `asr11_worst_case_slice_rev1` internals, and overlay it on the Recommendation's figure.

A text extraction of M.1464-2 Table 1 earlier in this work came back with scrambled columns, so this step reads the PDF table visually, not from extracted text.

### Step 2: Close the open items from public sources

Search for public values for O1 to O4:

- **O1, ASR-9 sidelobes:** FAA specifications and test reports (for example, the ASR-9 system specification and MIT Lincoln Laboratory ATC reports), and IEEE or AMS papers on the ASR-9 weather channel.
- **O2, ASR tilt:** FAA orders and siting criteria (for example, FAA Order 6310 series maintenance handbooks), and published ASR site surveys or environmental assessments that state antenna tilt.
- **O3, WSR-88D back lobe:** ROC and NSSL antenna measurement reports, Doviak and Zrnić, and NEXRAD radome and antenna studies in AMS journals and conference preprints.
- **O4, ASR aperture:** manufacturer data sheets (Raytheon or Northrop Grumman for ASR-11, Westinghouse for ASR-9), Forecast International, and FAA or DoD system descriptions (the ASR-11 is also the DoD AN/GPN-30 DASR).

**Deliverable:** for each item, the best public value with citation and confidence, or a clear statement that no public value exists and which agency holds it.

**Outside check:** each new value must have two independent public sources, or one primary source (a specification, test report, or Recommendation).

### Step 3: Find the other radars and patterns in the band

List the radar systems a 2700-2900 MHz coordination study must protect, and the public antenna data for each:

- The other radars in M.1464-2 Table 1 (Radars A, D, E, and F, and the combined ATC and weather Radars F1 and F2), each mapped to a real system where the Recommendation's description allows it.
- The meteorological radars in M.1849-3 Annex 2 that operate in this band.
- Federal radars in 2700-2900 MHz listed in public NTIA, FAA, DoD, or FCC documents, including military ATC radars (for example, AN/GPN-30 and AN/FPN-66 if still in service).
- Planned systems with public descriptions, such as NOAA's next-generation weather radar work and FAA surveillance radar replacements, and the antenna type each implies (for example, a planar phased array, M.1851-2 Section 6).

**Deliverable:** a system table with antenna type, gain, beamwidths, sidelobe levels, elevation pattern, scan, the source of each value, and which M.1851-2 model fits.

**Outside check:** cross-reference each system against the NTIA Government Master File summary data or the published band plan, where public.

### Step 4: Find the precedent for radar patterns in 5G studies

Find which radar antenna patterns past public studies used when protecting radars from commercial mobile systems, and why:

- NTIA technical reports for 3.5 GHz (for example, NTIA TR 15-517 and the CBRS exclusion zone reports), 3.45 GHz, and AWS-3.
- ITU-R studies for IMT and radars (for example, Report ITU-R M.2112 and studies under WRC-23 agenda items).
- FCC dockets on 2700-2900 MHz or adjacent bands, and FAA or NOAA comments in those dockets.

**Deliverable:** a precedent table: study, band, radar, pattern used (peak, average, statgain, F.1245, measured), single-entry or aggregate, and the stated reason.

**Outside check:** each row quotes the report's own section or table number.

### Step 5: Measure how much each pattern choice moves a zone

Size a simple single-radar coordination distance with each pattern option and report the difference. Use the slices, a 5G base station EIRP from the repository's IMT parameters, and a propagation model (free space as a bound, then ITM through the project's vectorized MATLAB port).

- **Options:** ROC spec or NSSL for the WSR-88D; peak or average envelope for the ASRs; radar-specific or statgain or F.1245; 0° or operational tilt.
- **Deliverable:** `zone_sensitivity_patterns_rev1.m` and a table and figure of zone distance by option, for main-beam, near-sidelobe, and far-sidelobe geometry.
- **Outside check:** the free-space result must match the closed-form distance from the link budget to within 0.1 dB.

Follow the project rule: reuse the most recent version of each function across repositories, such as the newest vectorized ITM and P.2108 code.

### Step 6: Add the gain distribution over a turn

A rotating radar spends most of each turn on its sidelobes. Aggregate Monte Carlo studies need the distribution of gain over a turn, not only the slice.

- **Deliverable:** a function that returns the exceedance curve of each slice (share of the turn at or above each gain), plus the time-averaged gain over a turn. Add a figure to the paper.
- **Outside check:** the share of the turn above the floor must match the paper's values (5.6%, 2.0%, and 1.5%).

### Step 7: Extend the functions

- Add `envelope="peak"|"average"` to the ASR functions.
- Add a station elevation input (`station_elev_deg`) to all four functions, so that stations above the horizon (on terrain or towers) get their real gain, not the horizon bound.
- Add one function per new radar from Step 3, with the same N-by-2 output.
- **Outside check:** the repository test pattern: match `antenna_prepare_rev1` and `antenna_gain_rev1` within 10<sup>−9</sup> dB, plus hand values.

### Step 8: Port the slices into the zone studies

Copy the slice functions into the repositories that size zones (`aggregate_explainer` and `2dot7_agent_retuning`) and replace any older radar antenna code there. Check each repository's current antenna code first and keep the most recent version.

- **Outside check:** each repository's test suite passes, and a reference case reproduces its previous result except where the pattern changed.

### Step 9: Red-team review

Have a reviewer that didn't build the work try to break it. The reviewer takes the position of a commercial filer arguing the zones are too large:

- Are the specification bounds unrealistic compared with measured antennas?
- Is a 0° tilt defensible when the lowest operational cut is 0.5°?
- Is the peak envelope right for aggregate interference?
- Does the analysis double-count conservatism (bound pattern, bound tilt, bound propagation)?

**Deliverable:** a findings table with a response or fix for each finding.

### Step 10: Decision package and revision 2

- **Deliverable:** the working paper, revision 2, with the verified ledger, closed open items, new radars, and the zone sensitivity. Add a short decision memo that lists the choices in [Decisions for Nick](#decisions-for-nick), with options and the zone effect of each.
- **Human review:** the decision memo and revision 2 need your review before they go outside the office.

## Decisions for Nick

These choices belong to you or the study team. Step 5 supplies the zone effect of each.

| Decision | Options | Evidence the plan supplies |
| --- | --- | --- |
| ASR envelope for aggregate 5G interference | Peak (M.1851-2 single entry) or average (M.1851-2 aggregate) | Step 4 precedent, Step 5 zone effect |
| WSR-88D sidelobes | ROC specification or NSSL measurement | Step 2 back-lobe data, Step 5 zone effect |
| Beam tilt | 0° bound or lowest operational cut | Step 2 ASR tilt data, Step 5 zone effect |
| Generic patterns | Allow statgain or F.1245, or require radar-specific patterns | Step 4 precedent, Step 5 zone effect |
| ASR-9 floor | WEBER_1986 bound or a value from Step 2 | Step 2 |

## Order and dependencies

1. Run Steps 1 to 4 together with the research prompt. They need no code.
2. Run Steps 5 to 7 after Step 1, because they build on verified values.
3. Run Step 8 after Step 7.
4. Run Step 9 on the Step 5 results and the extended functions.
5. Write Step 10 last.

## Who works on it

These roles describe methods, not voices:

- **RF engineer:** link budgets, antenna models, and propagation for Steps 5 to 7.
- **Radar systems engineer:** maps M.1464-2 and M.1849-3 radars to real systems, tilt, and scan for Steps 2 and 3.
- **Federal operations planner:** which radars and sites must be protected, and their operating modes, for Step 3.
- **Regulatory analyst:** precedent in NTIA, ITU-R, and FCC records for Step 4.
- **Data visualization reviewer:** one clear chart of zone distance by pattern choice for Step 10.
- **Red team:** the commercial filer in Step 9.

## Research prompt for Claude Fable

Paste the following into Claude Fable. It covers Steps 1 to 4.

```text
You are researching radar antenna patterns for a U.S. spectrum study. The
study sizes coordination zones that protect radars in the 2700-2900 MHz band
from 5G base stations. Use public sources only: ITU-R Recommendations and
Reports, NTIA technical reports and memoranda, NOAA Radar Operations Center
(ROC) and NSSL documents, FAA orders and specifications, MIT Lincoln
Laboratory reports, IEEE and AMS papers, FCC dockets, Federal Register
notices, and manufacturer data sheets. Don't use or infer controlled,
proprietary, or nonpublic data. If a value isn't public, say so and name the
agency that holds it.

For every value you report, give the source title, edition or date, table or
equation number, page, a quote of five words or fewer, the URL, and your
confidence (high, medium, or low). Read tables from the PDF page itself, not
from extracted text; one earlier extraction of M.1464-2 Table 1 scrambled
the columns. Never estimate a value and present it as sourced. Mark any
value you derived and show the arithmetic.

Our current model, to verify or correct:
- WSR-88D: 45.5 dBi peak, 0.925 deg 3 dB beamwidth at 2850 MHz, sidelobes
  -29 dB to 10 deg and -40 dB beyond (ROC NEXRAD Technical Information);
  NSSL 1998 errata Section II measured sidelobes -35 dB to 2 deg, -48 dB at
  6 deg, -55 dB at 20 deg; 8.53 m reflector.
- ASR-11 (M.1464-2 Radar C): 34 dBi, 1.45 deg azimuth and 4.8 deg elevation
  beamwidths, +9.5 dBi first sidelobe, +3.5 dBi remote sidelobes,
  cosecant-squared to 30 deg, tilt range -2.5 to +2.5 deg.
- ASR-9 (M.1464-2 Radar B): 33.5 dBi, 1.3 deg and 4.8 deg beamwidths,
  +7.3 dBi first sidelobe, no remote sidelobe level given; we use the
  -28 dB maximum azimuth sidelobe from MIT Lincoln Laboratory ATC-138
  (Weber, 1986).
- Pattern equations: Rec. ITU-R M.1851-2 cosine taper (n = 1, K = 68.8,
  peak envelope A = 17.51, B = 2.33, breakpoint -14.4 dB, average envelope
  offset -4.32 dB), cosecant-squared elevation Eqs. (22), (25)-(27);
  NTIA TM-13-489 Table 4-1 (statgain); Rec. ITU-R F.1245.

Research tasks:

1. Verify each value above against its primary source. Report any
   mismatch, and confirm that Radar B is the ASR-9 and Radar C is the
   ASR-11 from M.1464-2 Annex 2's descriptions.

2. Find public values for these gaps:
   a. ASR-9 remote (far) sidelobe level.
   b. ASR-9 and ASR-11 antenna tilt at typical sites, or FAA siting
      guidance that sets it.
   c. WSR-88D back-lobe and far-sidelobe measurements beyond 20 deg.
   d. ASR-9 and ASR-11 antenna aperture dimensions (for far-field
      distance).

3. List every radar system that operates in 2700-2900 MHz in the United
   States, from public sources: the other radars in M.1464-2 Table 1
   (Radars A, D, E, F, F1, F2) mapped to real systems where possible, the
   meteorological radars in M.1849-3 Annex 2, military ATC radars (for
   example AN/GPN-30), and planned systems (for example NOAA's
   next-generation weather radar and FAA surveillance radar replacements).
   For each, give antenna type, gain, beamwidths, sidelobe levels,
   elevation pattern, scan rate, and the best-fitting M.1851-2 model.

4. Find which radar antenna patterns past public studies used to protect
   radars from commercial mobile systems: NTIA reports for 3.5 GHz (for
   example TR 15-517 and the CBRS zone reports), 3.45 GHz, and AWS-3;
   Report ITU-R M.2112 and other IMT-radar studies; FCC dockets for
   2700-2900 MHz or adjacent bands. For each, report the band, radar,
   pattern used (peak, average, statgain, F.1245, measured), whether the
   analysis was single-entry or aggregate, and the reason the report gives.

5. Report what M.1851-2 says about when to use the peak envelope and when
   to use the average envelope, with the section number and a short quote.

Output, as one Markdown document:
- A summary of no more than 10 bullets: what changed in our model, and
  what stays open.
- Table A: verification of each current value (value, source, table,
  page, match or mismatch, confidence).
- Table B: values for the gaps in task 2, or "not public" with the
  holding agency.
- Table C: radar systems in the band (task 3).
- Table D: precedent studies (task 4).
- The M.1851-2 envelope guidance (task 5).
- A source list with full citations and URLs.
- A list of claims you couldn't verify, kept separate from verified
  findings.

Write in plain, direct sentences, present tense, and active voice. Don't
recommend a policy choice; report what the sources say.
```

## Checks on this plan

- Every open item and question in this plan comes from the working paper, revision 1, or from the tests and figures behind it.
- The numbers quoted here (3.3 dB, 4.32 dB, 2.6 to 21 dB, 9.5 to 18.1 dB, and 5.6%, 2.0%, and 1.5%) match the paper's numbers ledger.

**DRAFT//PRE-DECISIONAL**
