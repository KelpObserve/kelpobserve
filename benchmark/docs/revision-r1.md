# Sesimbra habitat benchmark: review check and exploratory re-analysis r1

Executed 2026-10-02, after an external review of v1 the same day. **The review's six
findings hold. Correcting them leaves v1's ranking in place. v1's negative result is
inconclusive: in an injection test, the registered models did not reliably recover
injected signals.** (Wording narrowed after a second review; see the end.)

- **Consistent QC changes the numbers, not the ranking.** Depth-only is still the best
  July candidate on Pléiades Neo, raw or glint-corrected, and on the one partly usable
  Sentinel-2 scene.
- **The visible bands show nothing NIR does not.** A near-infrared band, which cannot see
  a 7–17 m seabed, does as well wherever either appears to help.
- **Injection: the registered models did not reliably recover injected signals.** Even an
  optimistic kelp-on-bare-substrate signal with perfect positions was missed by the ridge
  in 54% of draws natively and 91% on the 10 m grid. Under the registered 2.5 m patch
  model, any position error of 5 m or more erased it.
- **Fixed June and July dates rarely get a clean Sentinel-2 overpass.** Within a week of
  them, 2 of 9 years had one likely to pass; September had one in 7 of 9.
- **ACOLITE's own glint corrections improve the water radiometry without reading labels.**
  On 3 September they give the first plausible PNeo QAA fit. 8 July stays implausible, so
  the corrected branch still cannot be trained.

r1 is exploratory. Every cover value had been seen in v1, so nothing here replaces a v1
result. The design, variants and thresholds were registered in
[`registration/revision-r1.json`](../registration/revision-r1.json) before any r1 result
was computed. v1's outputs, lock and reports are unchanged.

## The review's claims, checked

| Claim | Verdict | Evidence |
|---|---|---|
| Summer glint recurs at Sesimbra, worst on orbit R080 | Holds, 2017–2026 | 2026: every June–July R080 candidate reads 0.10–0.35 at 1610 nm over water, every R037 candidate 0.031–0.041 except 2 July (0.008) and 12 July (0.166). Across ten summers R080 never passes in June or July (Q4) |
| The July Sentinel-2 comparison was excluded on a coverage rule, 14 of 18 cells against 15 required | Holds | `scene-selection.md`; r1 reproduces the 14 cells |
| "No usable scenes matched the campaigns" is too categorical | Holds, with a correction | The phrase is not in these reports. The nearest statements are `coastal-validation/README.md` ("only winter and spring scenes are usable"), which 2 July and 31 August contradict, and `scene-selection.md` consequence 5, which calls recovery a field-campaign problem rather than a processing one |
| The winter availability figures do not answer the campaign question | Holds, and more so | Those 25–33 usable dates per season were screened over the whole 36 × 45 km site, not the 620 m quadrat area |
| 1. The PNeo optics rule could never pass | Holds | The as-run 20 September fits already carry `a_cdm_443_not_computable`, which fires for any sensor without a band within 15 nm of 411 nm, including Sentinel-2 (B1 is 32 nm away). Backscatter reproduced exactly: 0.126 and 0.120 m⁻¹ automatic, 0.062 and 0.069 forced 560 nm |
| … but the radiometric problem is real | Holds | No flag or mask variant brings backscatter under 0.02 m⁻¹ (Q1); only glint correction does, and only on 3 September (Q2). The "offshore NIR about 0.03" is an unmasked mean; medians are 0.023 (July) and 0.011 (September) |
| 2. The glint screen was not water-only | Holds; no verdict changes | Water-only values are 0.003–0.005 lower; 9 June falls from 0.0354 to 0.0312, and the nearest miss, 31 August R080, from 0.0251 to 0.0214 |
| 3. Pixel QC was not propagated into the predictors | Holds, and wider | Also the QAA reference region was unmasked (6% and 10% of its pixels). One July unit was affected as well (`july-shallow-T3-Q2`, scored in the July comparison) |
| 4. Glint correction options were not exhausted | Holds | `dsf_residual_glint_correction=False` in `registration/acolite/common.txt`. The Sesimbra glint tests in coastal-validation used winter scenes with too little glint to fit; the processor comparison ran at Ria Formosa |
| 5. The lidar attenuation regression cannot isolate attenuation at Sesimbra | Holds | E2c says so; r1's NIR control shows the same depth-correlated signal reaching the habitat predictors (Q3) |
| 6. The positional and support sensitivities were registered but never run | Holds | `plan.md` section 8, item 2. r1 runs them (Q3) |

One correction to the review: its reproduction script looked up `july-shallow-T3-Q2` in the
September predictions, so it reported that unit as unscored. It was scored, in the July
comparison.

## Q1: which observations survive consistent QC

Valid pixel, as v1's scene checks define it: water (lidar or DTM), finite blue, green and
red, and ACOLITE's cloud and glint mask passed (PNeo rhot at 826 nm ≤ 0.1; Sentinel-2 at
1610 nm ≤ 0.0215). A native unit needs all nine pixels valid; a 10 m cell needs 80% valid
area. Under v1's own finite-only rule, r1's code reproduces every v1 unit to 1e-15 and
every v1 July metric exactly.

| Sensor, date, scale | v1 units | r1 kept | Excluded |
|---|---:|---:|---|
| PNeo 8 Jul, native | 30 | 29 | shallow T3 Q2 (masked pixel) |
| PNeo 8 Jul, 10 m | 18 | 18 | – |
| PNeo 3 Sep, native | 22 | 19 | deep T1 Q5 and shallow T2 Q4 (masked), shallow T3 Q4 (land) |
| PNeo 3 Sep, 10 m | 13 | 12 | cell 9-10, shallow T3 Q4 (23% water) |
| S2 2 Jul, 10 m | not run | 14 of 18 | four shallow cells (mask) |
| S2 31 Aug, 10 m | not run | 12 of 13 | cell 9-10 (land) |

All but one exclusion are shallow, where July cover is highest; the exception is deep T1 Q5
in September.

QAA over the registered PNeo reference region, current library:

| Date | Variant | Pixels | Reference | bbp(555) m⁻¹ | Kd(490) m⁻¹ | Fatal flags |
|---|---|---:|---|---:|---:|---|
| 8 Jul | automatic, unmasked (v1) | 199,238 | 655 nm | 0.126 | 1.32 | red reference, bbp |
| | automatic, masked | 186,592 | 655 nm | 0.101 | 1.16 | red reference, bbp |
| | forced 560, unmasked | 199,238 | 563 nm | 0.062 | 0.66 | bbp, a below water at 655 nm |
| | forced 560, masked | 186,592 | 563 nm | 0.049 | 0.57 | bbp, a below water at 655 nm |
| 3 Sep | automatic, unmasked (v1) | 199,484 | 655 nm | 0.120 | 1.26 | red reference, bbp |
| | automatic, masked | 180,096 | 655 nm | 0.067 | 0.88 | red reference, bbp |
| | forced 560, unmasked | 199,484 | 563 nm | 0.069 | 0.73 | bbp, a below water at 655 nm |
| | forced 560, masked | 180,096 | 563 nm | 0.040 | 0.53 | bbp, a below water at 655 nm |

The plausible range for bbp(555) is 0 to 0.02 m⁻¹. Separating the informational flag and
masking the region both move the fit the right way, by up to a factor of three, and none
reaches plausibility. The corrected branch stays unavailable.

## Q2: can summer radiometry be improved without labels?

One registered variant per sensor, both settings of the pinned ACOLITE: residual glint
correction from SWIR for Sentinel-2, and the "alternative" correction from the 826 nm band
for PNeo. The PNeo settings are ACOLITE's own defaults for WorldView-2 and GeoEye-1. No
label was read.

| Scene | Offshore NIR V0 → V1 | NIR/green depth-slope ratio V0 → V1 | Contrast V0 → V1 | Land change | Verdict |
|---|---|---|---|---:|---|
| S2 28 May R037 | 0.0090 → 0.0029 | 0.52 → 0.07 | 3.05 → 2.90 | 0 | improves |
| S2 2 Jun R037 | no valid water | – | – | 0 | not assessable |
| S2 2 Jul R037 | 0.0111 → 0.0082 | 2.93 → 0.17 | 0.76 → 0.72 | 0 | improves |
| S2 9 Jul R037 | no valid water | – | – | 0 | not assessable |
| S2 5 Jul R080 | no valid water | – | – | 0 | not assessable |
| S2 31 Aug R037 | 0.0088 → 0.0030 | 2.88 → 0.21 | 0.35 → 0.90 | 0 | improves |
| PNeo 8 Jul | 0.0165 → 0 | 1.11 → 0.00 | 0.82 → 1.83 | 0 | improves, see limit |
| PNeo 3 Sep | 0.0060 → 0 | 1.86 → 0.01 | 0.50 → 1.42 | 0 | improves, see limit |

- **On Sentinel-2 the SWIR correction does what it should** wherever there is water to
  judge it on. Offshore NIR falls by 0.003–0.006. The near-shore NIR depth slope, the
  common-mode signature, almost vanishes while the green slope remains, and land is
  untouched.
- **It cannot rescue the glinted summer scenes under the registered mask.** That mask
  reads top-of-atmosphere 1610 nm, which correction does not change. On 2 June, 9 July and
  5 July it removes every water pixel. Using those scenes needs a glint-tolerant mask: a
  new registration, judged on the same label-free measures.
- **Limit for PNeo:** its correction takes glint from the NIR band, so it zeroes offshore
  NIR and the NIR slope by construction. Two of the four registered checks pass
  automatically (registration erratum). The informative evidence is QAA:

| Date | Variant, masked reference region | Reference | bbp(555) m⁻¹ | Kd(490) m⁻¹ | Available |
|---|---|---|---:|---:|---|
| PNeo 3 Sep | V0, forced 560 | 563 nm | 0.040 | 0.53 | no |
| PNeo 3 Sep | **V1, automatic and forced 560** | 563 nm | **0.0059** | **0.24** | **yes** |
| PNeo 8 Jul | V0, forced 560 | 563 nm | 0.049 | 0.57 | no |
| PNeo 8 Jul | V1, forced 560 | 563 nm | ≈ 0 | 0.12 | no: over-corrected |
| S2 2 Jul / 31 Aug | V1, forced 560 | 559 nm | 0.016 / 0.009 | 0.17 / 0.15 | no: a below pure water at 665 nm |

- **3 September gives the first plausible PNeo reference fit.** The automatic rule itself
  picks 563 nm. On 8 July the same correction leaves the offshore green at 0.001, too dark
  for clear water: the July NIR term is not all glint, and subtracting it as glint
  over-corrects.
- **The corrected branch still cannot run.** It needs a July fit to train on.
- **The Sentinel-2 fits fail only the red pure-water floor,** which r1 registered as fatal.
  E2c's plausibility criterion (reference, bbp, Kd) would pass all three.

**Q3 repeated on corrected PNeo** (the registered follow-on). Depth-only stays best: July
native 22.4 against 25.5 for spectral + depth and 36.4 for spectral alone; on the grid
spectral alone is worse than the mean (46.5 against 38.5). Under displacement, spectral
alone beats the mean in 5 of 36 native scenarios, against 20 of 36 before correction. Most
of the apparent spectral skill under displacement was the NIR-correlated term the
correction removes. The NIR controls are constants under this variant and say nothing.

**The July REFLECTANCE conversion is consistent with September's BASIC product.** Over
stable land, PNeo/Sentinel-2 ratios on 8 July (against 2 July) and 3 September (same
minute) differ by 0.007 (green), 0.046 (red) and 0.009 (NIR), all within the registered
0.05. This checks for gross error only; it is not the DIMAP numerical test.

It does show something v1 missed: **on both dates PNeo reads 18–27% below Sentinel-2 over
land in green and red** (ratios 0.73 and 0.77–0.82; 0.85–0.96 in NIR). v1's own
cross-sensor output has it too (green 0.052 against 0.071 on 3 September). v1 reported
only the correlation (about 0.9) and called the scaling sound, but a correlation cannot see
a gain. Over water PNeo reads higher than Sentinel-2, over land lower. Calibration, band
response and view geometry are all possible causes, none is established, and it applies to
both PNeo products alike.

## Q3: is there spectral information beyond depth?

July leave-one-transect-out, v1 families and settings, every family on the same
QC-consistent units. Two negative controls are added: ridge on near-infrared alone, and
on near-infrared plus depth. Water absorbs NIR within about a metre, so at 7–17 m it
carries no seabed signal; any skill it shows comes from something correlated with cover
and depth: adjacency from the cliffs, glint, position.

MAE in percentage points (v1 in brackets where it differs):

| Candidate | PNeo native (29) | PNeo 10 m (18) | S2 2 Jul, 10 m (14) |
|---|---:|---:|---:|
| Training mean | 34.9 (36.0) | 38.5 | 34.2 |
| **Depth-only ridge** | **22.4** (26.3) | **32.3** | **27.5** |
| Blue/green/red ridge | 38.1 (39.3) | 35.8 | 32.9 |
| Blue/green/red + depth ridge | 25.8 (33.4) | 41.7 | 31.2 |
| Spectral + depth forest | 24.1 (27.7) | 37.8 | 34.0 |
| *Control:* NIR ridge | 35.9 | 36.4 | 28.2 |
| *Control:* NIR + depth ridge | 26.0 | 40.8 | 27.7 |

- **Depth-only stays best on every sensor and scale.** Removing one unit with a single
  masked pixel (`july-shallow-T3-Q2`) moves the native spectral + depth ridge by 7.6
  points. That is how fragile the v1 differences between spectral models were, as the
  review suspected, but the ranking does not change.
- **The 10 m forest is mostly an intercept.** With 11–13 training cells per fold and a
  minimum leaf of five, 241–296 of its 300 trees never split. This holds in v1 too, so
  the grid forest says little about nonlinear skill. The native forest (19–20 training
  units) splits normally.
- **The NIR control matches the visible bands.** Natively, NIR + depth (26.0) ties
  blue/green/red + depth (25.8). On Sentinel-2, NIR alone (28.2) is within a point of
  depth-only and beats every visible-band model. The visible bands show no information
  that NIR, which sees no seabed, does not also show.
- Sentinel-2 on 2 July is exploratory twice over: the scene failed the registered cell
  check, and its four lost cells are all shallow.

**Sensitivities** (PNeo, July). Every position shifted jointly by 2, 5, 10 or 20 m in each
of eight directions, quadrat order reversed, and native windows of 1, 5 and 7 pixels: 36
scenarios natively, 33 on the grid.

| Paired MAE difference (negative is better) | Native: range | Native: scenarios ≤ −2 | 10 m: range | 10 m: scenarios ≤ −2 |
|---|---|---:|---|---:|
| Blue/green/red + depth − depth-only | −6.9 to +9.3 | 5 of 36 | −5.9 to +12.0 | 3 of 33 |
| *Control:* NIR + depth − depth-only | −3.5 to +6.9 | 3 of 36 | −5.5 to +8.7 | 3 of 33 |
| Blue/green/red − training mean | −14.8 to +7.5 | 20 of 36 | −19.7 to +13.2 | 26 of 33 |
| *Control:* NIR − training mean | −12.4 to +7.8 | not tallied | −18.6 to +2.0 | not tallied |

- **Position alone moves the comparison by about ±7 points.** Differences smaller than
  that are not interpretable with these 18–30 units.
- **Where the spectral bands beat the mean, NIR does too, by a similar margin.** Quadrat
  order reversed: −8.4 and −8.4; 20 m south-east on the grid: −19.7 and −18.6; 10 m west
  natively: −5.1 and −12.4. The
  apparent spectral skill under displacement is the common-mode, depth-correlated signal
  E2c found in the attenuation regression, not seabed reflectance.
- No scenario gives the visible bands a consistent advantage over depth-only that the NIR
  control does not share.

**September transfer** (refit on all July units, applied once; exploratory, since v1 already
evaluated September). Every PNeo family still overpredicts by 37–58 points with MAE ≈
bias: the level shift is not a QC artefact. On Sentinel-2, the blue/green/red ridge fails
outright (MAE 143, bias −143), because 31 August reflectance sits outside the 2 July range,
while NIR alone scores 7.0. That number is a coincidence of the seasonal level shift, and it
is why the control is there: on 12 cells, a band that cannot see the seabed can "transfer"
best.

## Injection and recovery: could this design have seen a signal?

A synthetic cover target (July cover permuted within each stratum) was painted into the
real 8 July reflectance as 2.5 m discs at each quadrat's true position. The signal came
from the library's forward model, with E5s endmembers and Sesimbra QAA water. The
quadrats were then extracted at their recorded positions under the Q1 QC, and v1's July
comparison was run. 100 draws per condition. Recovered means the blue/green/red + depth
ridge beats depth-only by at least 2 points in 80% of draws.

Signal at the quadrats for 100% kelp, median rhos change at 486 / 563 nm:

| Kelp against | Shallow (≈ 7–9 m) | Deep (≈ 17 m) |
|---|---|---|
| Understory: turf, crust, red foliose (realistic here) | +0.002 / +0.004 | +0.0001 / +0.0003 |
| Bare rock and sand (optimistic bound) | −0.023 / −0.025 | −0.0014 / −0.0015 |

Share of draws in which spectral + depth beats depth-only by ≥ 2 points (median water;
the more attenuating date gives the same picture):

| Condition | Native ridge | Native forest | 10 m ridge |
|---|---:|---:|---:|
| No signal (null) | 5% | 2% | 5% |
| Understory, no position error | 4% | 0% | 2% |
| Understory, 5–20 m error | 5% | 1% | 3–4% |
| Bare, no position error | 46% | 70% | 9% |
| Bare, 5 m error | 2% | 7% | 7% |
| Bare, 10–20 m error | 2% | 7–8% | 4% |

- **No condition met the recovery criterion, not even the optimistic signal with perfect
  positions.** The
  injected signal does reach the units: natively, the spectral ridge beats the training
  mean in 98–100% of bare-substrate draws (80–87% on the grid). What fails is the step
  beyond depth.
- **v1's negative result is therefore inconclusive as a test of the sensor.** With these
  units, these models and the registered patch model, "no skill beyond depth" was the
  likely outcome even if a strong signal had been there. That is a statement about the
  design, not a demonstration that the sensor cannot detect kelp.
- **The 10 m forest column says little.** Most of its trees never split (see Q3). The
  native forest is the only nonlinear model that was really tested.
- **Kelp over understory is indistinguishable from no signal.** This is consistent with
  E5s (not separable at 1 m or deeper).
- **A 5 m position error removes even the bare-substrate signal.** That depends on the
  registered patch model, discs of 2.5 m. That radius is an assumption about spatial
  coherence, not a consequence of the 0.25 m² quadrat, which can sit in a much larger
  bed. Larger patches would lose less.
- **A likely reason the ridge fails while the forest partly succeeds** (an interpretation,
  not tested here): the signal shrinks about fifteen-fold from the shallow to the deep
  stratum, and an additive ridge cannot represent a slope that depends on depth. Scaling
  for depth is what the water-column-corrected predictors were meant to do, and that
  branch is unavailable.

## Q4: what fixed campaign dates can expect

Every Sentinel-2 L2A overpass in equivalent windows of 2017–2026, screened over the
quadrat area (water from the lidar and DTM) with the library's usability gate: viewing
geometry at least 20° from specular, water 1610 nm ≤ 0.02, no haze, tile cloud ≤ 60%.
"Likely to pass" adds 1610 nm ≤ 0.012, where v1's calibration puts ACOLITE under its 0.01
limit. Neither tier checks quadrat cells, so both bound availability from above: 2 July
2026 is "likely to pass" here and still lost four shallow cells.

Years in 2017–2025 with at least one likely-to-pass overpass:

| Window around the survey dates | June (2 Jun) | July (7–8 Jul) | September (1–2 Sep) |
|---|---:|---:|---:|
| ± 3 days | 1 of 9 | 1 of 9 | 5 of 9 |
| ± 7 days | 2 of 9 | 2 of 9 | 7 of 9 |
| ± 14 days | 4 of 9 | 4 of 9 | 8 of 9 |
| 2026, ± 7 days | 0 | 1 (2 July) | 2 |

- **Under this screen, June and July are routinely fragile, not unlucky.** Within ± 14
  days, 4 of 101 June overpasses and 4 of 113 July overpasses in 2017–2025 passed the
  "likely to pass" tier. That counts scenes that pass v1-style glint limits without
  further processing. It does not bound what better glint processing could recover (r3). 2026, with one
  clean July scene, was an ordinary-to-good year.
- **R080 never qualifies in June or July.** All 94 of its non-cloudy overpasses fail on
  geometry (median 14° from specular). R037 (median 29°) is usable 15 times in ten summers;
  measured glint, which depends on sea state, accounts for most of its failures.
- **September is dependable:** a usable overpass within a week in every year.
- Widening the window trades glint for ecological mismatch. Cover fell from 48% to 9%
  between July and September 2026, so a scene two weeks off a June or July survey is not
  interchangeable with one on the day.

## What this changes in the v1 conclusions

| v1 statement | After r1 |
|---|---|
| Depth-only is the best candidate; every spectral candidate is worse (results.md) | Stands under consistent QC, on both PNeo variants and on Sentinel-2 2 July. The exact margins were fragile: one masked pixel moved a native MAE by 7.6 points |
| PNeo "does not predict cover beyond what depth alone implies" | True of this experiment, but inconclusive about the sensor. Under the registered models and 2.5 m patch model, and with float-scale position error, the design would very probably have returned the same answer had a strong signal been present; with perfect positions, about half the time (injection) |
| "Radiometric scaling sound" from land correlation (optics-readiness.md) | Correlation holds; the gain does not. PNeo reads 18–27% below Sentinel-2 over land in green and red on both dates |
| Removing the residual "is new method work, not a processing setting" (optics-readiness.md) | Partly overtaken: a pinned ACOLITE setting gives a plausible QAA fit on 3 September. Not on 8 July |
| Recovering Sentinel-2 "needs a survey timed to a low-glint overpass, … not a processing change" (scene-selection.md) | Not tested by r1, and not established. r1's registered mask cannot use glinted scenes. Fixed June–July dates get an overpass passing the glint screen within a week in 2 of 9 years. r3's glint-tolerant mask was not accepted under its own rules, which have weaknesses of their own. September is dependable |
| The July Sentinel-2 exclusion (14 of 18 cells) | A coverage rule, as the review said. Run anyway (exploratory), the surviving 14 cells give the same ranking |

**What would make a next test informative.** All of these come from r1's results; none is
tested here.

- **Sample design.** More units per stratum and more than two depth strata.
- **Models whose spectral term can scale with depth.** The water-column-corrected
  predictors, once a July-type scene has a plausible reference fit, or a registered
  interaction model.
- **Quadrat positions known to a few metres.** For example, a GPS on the frame rather
  than a dragged float.
- **If a summer Sentinel-2 comparison is wanted without new processing,** a June/July
  survey day chosen from forecast R037 overpasses would help. This is an option the
  window counts support, not a requirement r1 established.
- **A reference that can carry the signal,** such as cover estimated over a pixel-sized
  footprint rather than a 0.25 m² quadrat. Kelp against understory is at the noise floor
  even under ideal positions, so a quadrat-scale reference is unlikely to show it.

## Reproduce

```bash
PYTHONPATH=. python -m benchmark.revision qc          # Q1 support, v1 reproduction, QAA
PYTHONPATH=. python -m benchmark.revision screen      # water-only screen (Planetary Computer)
PYTHONPATH=. python -m benchmark.revision acolite     # V1 glint runs, about 25 min, 12 GB
PYTHONPATH=. python -m benchmark.revision radiometry  # Q2 measures and conversion check
PYTHONPATH=. python -m benchmark.revision models      # Q3 and sensitivities, raw PNeo
PYTHONPATH=. python -m benchmark.revision models --variant V1
PYTHONPATH=. python -m benchmark.revision injection   # 100 draws, about 5 min on 6 cores
PYTHONPATH=. python -m benchmark.revision windows     # Q4, 2017-2026 (Planetary Computer)
```

Code: `benchmark/revision.py`, tests in `tests/test_revision.py`. Every stage writes a run
under `out/r1/<stage>/` with a manifest carrying the registration hash, both runners'
hashes and the library commit. The registration's two errata (a mis-aligned September
grid in the pre-registration notes; two automatic passes in the PNeo improvement rule)
are recorded in the registration, not repaired.

## Narrowed after a second review (2026-10-02)

No result changed. These statements in the first version went beyond the evidence, and
have been reworded above:

- **"The design could not have detected a realistic signal."** What was shown is that the
  registered ridge and forest, on these units and under the 2.5 m patch model, did not
  reliably recover injected signals. The 10 m forest barely splits, and the additive
  ridge cannot represent a depth-dependent response.
- **The Q4 counts.** They describe overpasses passing a fixed glint screen, not what
  processing could recover. They do not establish that survey dates must move.

The next step these results point to is a detectability envelope, not another habitat
fit: idealised positive controls, a model with explicit depth dependence, a range of
patch sizes and position errors, and the real sampling geometries. See revision-r3.md,
consequence 4.
