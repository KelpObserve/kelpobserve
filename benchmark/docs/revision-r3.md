# Sesimbra habitat benchmark: a glint-tolerant Sentinel-2 mask (r3)

Executed 2026-10-02; interpretation corrected the same day after a second review (see the
end). **No glinted June or July scene passed the registered acceptance rules, so no
summer Sentinel-2 comparison was selected.** That is the registered outcome. It does not
show that processing cannot recover these scenes, that the sensor could not detect
cover, or that fieldwork has to be rescheduled.

- **The mask admits the water.** It raises v1's threshold by each scene's glint above
  clean conditions, and every glinted R037 scene keeps 97–99% of its geometric water;
  9 July passes the cell check (17 of 18). This shows admission, not verified clear water:
  the mask is a brightness threshold with no separate cirrus or shadow screening.
- **Each glinted scene fails at least one acceptance check, most narrowly.** The near
  misses lie inside their own sampling uncertainty, and the checks have weaknesses of
  their own (below).
- **Glinted dates show more small-scale variability and overcorrection.** Neighbouring-
  pixel scatter offshore is 1.4–3 times the clean scenes', and 5–28% of offshore pixels
  have negative NIR after correction, against under 2%. The cause is not identified.
- **On 2 July the registered models did not reliably recover injected signals.** That is
  a statement about those models on 14 cells, not about Sentinel-2 in general.

Numbered r3 because r2 is a separate draft protocol
([`registration/revision-r2.json`](../registration/revision-r2.json), not executed) saved
earlier the same day; this design did not draw on it. Exploratory, like r1: every cover
value had been seen. The design was registered in
[`registration/revision-r3.json`](../registration/revision-r3.json) before any r3 value was
computed or any new scene downloaded.

## Design in one paragraph

All June and July R037 candidates whose water-only 1610 nm screen was at most 0.059 (that
is, likely within ACOLITE's 0.05 glint-correction limit) were processed with r1's SWIR
residual glint correction: five June dates and 9 July, plus 2 July and 31 August as clean
controls and 5 July R080 as a severe one. Three were downloaded for r3.

- **Glint background B:** the median top-of-atmosphere 1610 nm over habitat water at
  least 10 m deep and 100 m from land.
- **Mask threshold:** T = 0.0215 + max(0, B − B_clean), where B_clean is the larger clean
  background (0.0065). Clean scenes keep v1's 0.0215 exactly.
- **Glint cap:** B ≤ 0.05.
- **Acceptance:** a scene passes only if v1's valid-water and cell checks pass and all five
  label-free checks do:
  - A1: offshore NIR;
  - A2: glint leaking into green;
  - A3: near-shore NIR/green depth-slope ratio;
  - A4: offshore green spread (spatial variability, not noise);
  - A5: green depth slope.

  Limits are set from the clean controls.

## Results

| Scene | B | T | Valid water | Cells | A1 NIR ≤ 0.0102 | A2 \|slope\| ≤ 0.2 | A3 ratio ≤ 0.212 | A4 MAD ≤ 0.00236 | A5 slope ×0.5–2 | Accepted |
|---|---:|---:|---:|---:|---|---|---|---|---|---|
| 2 Jun R037 | 0.038 | 0.053 | 98% | 17/17 | 0.0036 | **0.51** | 0.18 | **0.002364** | **2.25×** | no |
| 30 May R037 | 0.034 | 0.049 | 98% | 17/17 | 0.0043 | 0.16 | 0.08 | **0.00263** | 1.6× | no |
| 28 May R037 | 0.023 | 0.038 | 97% | 17/17 | 0.0033 | 0.19 | **0.226** | 0.00134 | 1.8× | no |
| 7 Jun R037 | 0.035 | 0.050 | 98% | 17/17 | 0.0023 | 0.18 | 0.18 | **0.00295** | 1.8× | no |
| 9 Jun R037 | 0.031 | 0.046 | 98% | 17/17 | 0.0032 | 0.197 | **0.37** | **0.00256** | **2.1×** | no |
| 9 Jul R037 | 0.040 | 0.055 | 99% | 17/18 | 0.0019 | −0.14 | **0.39** | 0.00136 | 1.35× | no |
| 2 Jul R037 (clean) | 0.0065 | 0.0215 | 97% | **14/18** | 0.0082 | – | 0.17 | 0.00118 | 0.85× | radiometry yes, cells no (known) |
| 31 Aug R037 (clean) | 0.0065 | 0.0215 | 98% | 14/17 | 0.0030 | – | 0.21 | 0.00113 | 1.15× | yes (September) |
| 5 Jul R080 | **0.107** | – | – | 0/18 | fails everything | | | | | no: beyond the cap |

Bold marks failures. A5 is the ratio to the clean controls' mean green slope.

- **Near misses, reported and not acted on.** Four scenes fail a single check: 30 May
  and 7 June on A4 (1.11× and 1.25× the limit), 28 May on A3 (1.06×), 9 July on A3. 2 June
  fails A4 by 0.0000014, a coin flip, but also A2 and A5. Relaxing a limit now would let
  the outcome choose the rule, as v1's scene-selection report put it for 2 July.
- **The near misses are within sampling uncertainty.** A 100 m spatial-block bootstrap
  (300 resamples, descriptive only) gives 28 May's A3 ratio 0.226 (95% 0.072–0.356, limit
  0.212) and 30 May's A4 0.00263 (0.00201–0.00313, limit 0.00236). The limits come from two
  clean dates and carry uncertainty of their own.
- **What the correction does.** Offshore NIR medians fall to 0.002–0.004, and glint
  leaking into green is under 0.2 on every R037 scene except 2 June. ACOLITE corrected
  every valid offshore pixel (checked on 28 May, 2 June and 9 July), so skipped pixels do
  not explain any failure.
- **QAA: some quantities become plausible, the spectra do not.** With the forced 560 nm
  reference every corrected scene meets E2c's broad criterion: bbp(555) 0.006–0.011 m⁻¹,
  Kd(490) 0.12–0.17 m⁻¹, and on 9 July the automatic rule picks 560 nm. But every R037
  scene, clean controls included, infers absorption near 665 nm of 0.14–0.19 m⁻¹, below
  pure water's 0.42 m⁻¹, and stays unavailable under r1's fuller flags. That needs
  diagnosing before these fits feed a water-column correction.
- **What it leaves, diagnosed only in part:**
  - **Small-scale variability.** Offshore, the scatter between neighbouring pixels (robust
    SD of differences / √2) is 0.0014–0.0029 on glinted dates against 0.0009–0.0010 on the
    clean ones. On glinted dates it accounts for nearly all of A4, whereas clean scenes
    also carry larger-scale gradients.
  - **Its cause is open.** Candidates include SWIR noise carried into the visible bands,
    the bands' different resolution and slightly different viewing times over a moving
    sea surface, and real structure.
  - **Overcorrection.** Negative NIR in 5–28% of offshore pixels; on 7 June also 1.3%
    negative green.
  - **On 9 June and 9 July, a near-shore NIR depth slope** about twice the clean scenes'
    ratio, the common-mode signature r1 found at this reef.
- **Cross-date context** (reported, not scored):
  - The corrected June scenes sit within ±0.003 of 31 August in green over habitat water
    5–20 m; 9 July is 0.006 darker.
  - Their spatial agreement with 31 August (Spearman 0.25–0.54 in green) is at most the
    two clean scenes' agreement with each other (0.55).
  - 2 July, a clean control, is itself 0.009 brighter in green than 31 August.

## Injection on Sentinel-2

The registered fallback scene, since no July scene was selected: 2 July under v1's mask,
14 cells. The design is r1's, with E5s Sentinel-2 endmembers and the synthetic field built
on a 1 m grid and averaged to 10 m. At the real quadrat coordinates each 2.5 m disc
rasterises to 19–21 m² (mean 19.7 against 19.6 exact). Allowing for overlap, the patches
total 541 m² on the 1 m grid and 537 m² at 0.25 m, so the injection is close to its
nominal size.

| Condition | Ridge ≥ 2 points better than depth-only | Forest |
|---|---:|---:|
| No signal (null) | 18% | 4% |
| Understory, any position error | 15–16% | 4% |
| Bare substrate, no position error | 20–21% | 4% |
| Bare substrate, 5–20 m error | 13–16% | 4% |

**The registered models did not reliably recover the injected signals on the 14
retained cells.** That is the whole finding, and it is narrow:

- **Neither model can test the physics properly.** The forest barely functions: with 8–11
  training cells and a minimum leaf of five, 296, 300 and 300 of its 300 trees never split
  in the three folds. The additive ridge cannot represent a spectral response that
  changes strongly with depth.
- **One sampling geometry was tested.** 9 July's 17-cell support was not.
- **The 2.5 m patch radius is an assumption about spatial coherence**, not a consequence
  of quadrat size. A 0.25 m² quadrat can sit in a continuous bed much larger than a pixel.
  Larger patches would give each 10 m cell more of the contrast.
- **The spectral ridge beats the training mean in 88% of null draws**, so stratum-
  correlated reflectance alone produces that much apparent skill.

## Consequences

1. **Under these rules, June and July Sentinel-2 remain unselected.** This was one mask,
   one correction and one acceptance rule. It does not establish that summer scenes
   cannot be recovered, nor that CCMAR's survey dates must move.
2. **r1's campaign-window counts stand as counts of overpasses passing the library
   screen.** No glinted scene was accepted, so the registered "correctable" tier is
   undefined.
3. **Diagnose before correcting again.**
   - Compare variability across spatial scales, NIR, visible and SWIR relationships, and
     distance to shore.
   - Then test a small predefined set of correction variants, with matched pixel support
     and explicit checks against overcorrection and signal loss.
   - SWIR smoothing is one candidate. It should follow the diagnosis, because lowering the
     MAD alone could simply remove structure.
4. **Make recovery tests informative before another habitat fit.**
   - An idealised positive control with the injected physics.
   - A small model with explicit depth dependence.
   - Patch radii from sub-pixel to continuous beds, with positional uncertainty.
   - Both the 14- and 17-cell geometries, with reflectance computed at fine resolution
     and passed through the sensor's spatial response.

   The useful product is a detectability envelope: which combinations of depth, patch
   extent, position accuracy and residual radiometric error allow recovery. Run as r4
   ([revision-r4.md](revision-r4.md)).

## Reproduce

```bash
PYTHONPATH=. python -m benchmark.revision_r3 fetch --env-file <file with CDSE_USERNAME/PASSWORD>
PYTHONPATH=. python -m benchmark.revision_r3 acolite     # three scenes, about 9 min, 5 GB
PYTHONPATH=. python -m benchmark.revision_r3 mask
PYTHONPATH=. python -m benchmark.revision_r3 models      # records that nothing was selected
PYTHONPATH=. python -m benchmark.revision_r3 injection   # 100 draws, about 3 min on 8 cores
PYTHONPATH=. python -m benchmark.revision_r3 windows     # records the undefined tier
```

Code: `benchmark/revision_r3.py`, tests in `tests/test_revision_r3.py`. Runs are under
`out/r3/<stage>/`, each manifest carrying both registrations' hashes and the library
commit.

## Corrected after a second review (2026-10-02)

Registered outcomes are unchanged. These statements in the first version of this report
were wrong or went further than the evidence:

- **The 18% patch-area shortfall** came from a unit-test geometry centred on the grid. At
  the real coordinates the injection is close to nominal (above).
- **"Sentinel-2 at 10 m could not have detected the signal anyway"** extrapolated from
  models that cannot test it, from one patch radius, and from the 14-cell 2 July support
  to geometries that were not tested.
- **"About twice the pixel noise".** A4 is the spread of all offshore pixels: spatial
  variability, not noise. The neighbouring-pixel diagnostic above confirms more
  small-scale variability on glinted dates, but not its cause.
- **"The correction works offshore"** was too broad: QAA's red absorption is inconsistent
  on every scene.
- **The claim that v1's call for low-glint survey timing "stands"** read more into one
  failed variant than it supports.

Weaknesses in the registered checks, recorded in the registration's errata rather than
repaired:

- **A1 has no lower bound.** It would accept severe negative overcorrection: a fixture
  with an NIR median of −0.1 passes the whole rule.
- **The cap compares unlike quantities:** a scene-median top-of-atmosphere value with
  ACOLITE's per-pixel surface-reflectance limit.
- **A3 and A5 partly ask new scenes to reproduce a spatial gradient E2c could not
  explain.** A3 also improves whenever its green denominator grows.
- **A4 measures spatial variability, not noise.**
- **The limits rest on two clean dates,** with no account of their uncertainty.
