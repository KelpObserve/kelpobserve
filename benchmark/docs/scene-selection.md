# Sesimbra habitat benchmark: imagery availability and scene selection

Executed 2026-09-20. **Both Pléiades Neo deliveries pass the registered scene checks.
Sentinel-2 is unavailable for June and July, and one September scene qualifies.** The
benchmark's modelling therefore runs on Pléiades Neo alone, and the matched-sensor
comparison in the plan cannot be performed. No model was fitted and no cover value was
read at any point in this selection.

## Result

| Scene | Orbit | Glint, rhos 1614 nm (limit < 0.01) | Valid water (≥ 80%) | Quadrat cells (≥ 80%) | Verdict |
|---|---|---:|---:|---:|---|
| 2026-05-28 S2C | R037 | 0.0222 | 29.9% | 0/17 | rejected |
| 2026-06-02 S2B | R037 | 0.0310 | 0.0% | 0/17 | rejected |
| 2026-06-02 S2A | R080 | 0.1099 | 0.0% | 0/17 | rejected |
| 2026-07-02 S2B | R037 | **0.0040** | **97.1%** | 14/18 = 77.8% | rejected, quadrat cells |
| 2026-07-05 S2B | R080 | 0.1127 | 0.0% | 0/18 | rejected |
| 2026-07-09 S2A | R037 | 0.0295 | 0.0% | 0/18 | rejected |
| 2026-08-31 S2B | R037 | 0.0027 | 97.8% | 14/17 = 82.4% | **selected, September** |
| 2026-09-03 S2B | R080 | 0.0155 | 91.8% | 10/17 = 58.8% | rejected |
| PNeo 2026-07-08 | – | NIR 0.038 over water, 0.013 at cells | 90.9% | 30/30 | **accepted** |
| PNeo 2026-09-03 | – | NIR 0.011 over water, 0.015 at cells | 88.6% | 25/30 | **accepted** |

A further 14 candidate dates were screened out on measured glint without being
processed. Every candidate and its screen value is listed in
[the portable export](scene-selection.json).

## Why Sentinel-2 fails

Only two Sentinel-2 orbits cover Sesimbra, and in summer both view the sea close to the
specular direction. Measured with ACOLITE's own geometry rasters:

| Orbit and season | Specular angle | Median glint |
|---|---:|---:|
| R080, June to July | 11° | 0.11 |
| R037, June to July | 27–28° | 0.03 |
| R080, 3 September (lower sun) | 26° | 0.0155 |

Glint of that size exceeds the substrate contrast the retrieval depends on at 7–17 m, so
it is not a threshold artefact. A calm sea narrows the glint lobe, which is why two R037
dates (2 July and 31 August) are clean while other R037 dates on the same geometry are
not.

The second failure mode is shallow-water contamination. On every date, the rejected
quadrat cells are the shallow ones nearest the cliffs, where the 20 m SWIR band mixes in
land or surf and ACOLITE's water mask removes them. The shallow stratum carries the
highest cover, so losing it is not a neutral loss of coverage.

## Registered amendment

Under the plan's original four-day matchup window, **every** Sentinel-2 candidate failed,
on glint alone. One amendment was registered before any further processing, and is
recorded in `registration/experiment.json`:

- the matchup window widens from four days to seven, which admits 2 July;
- candidates whose screened glint exceeds 0.02 are recorded as failures without a full
  ACOLITE run. The screen is the median L2A reflectance at 1610 nm over the quadrat area,
  calibrated against ACOLITE on four processed scenes, where ACOLITE reads 0.69 to 0.85 of
  the screen value. Candidates with no L2A product to screen are processed in full.

The amendment changes which scenes are examined, never how a scene is judged, and scene
choice never uses cover-prediction error. The four-day result stands in the record: June,
July and September Sentinel-2 all unavailable.

## The July near-miss, and its sensitivity

2 July 2026 passes on glint and on valid water, and fails only the quadrat-cell check, by
one cell: 14 of 18, where 15 are needed. Two of its four rejected cells sit barely above
the mask threshold, at 0.0219 and 0.0246 against 0.0215; a third, at 0.0874, is clearly
contaminated.

The verdict is therefore fragile, and the alternatives disagree:

| Alternative | 2 July | Effect |
|---|---|---|
| Registered: cells, 1614 nm mask, 80% | 77.8% | fails |
| Cell fraction relaxed to 75% | 77.8% | would pass |
| Counted per quadrat instead of per cell | 80.0% | would pass |
| Mask taken from 10 m NIR instead of SWIR | 61.1% | fails, and is stricter |

These are reported, not acted on. The registered check counts cells because quadrats
sharing a pixel are not independent samples, and changing a threshold on a one-cell margin
would let the outcome choose the rule. Note that the NIR variant, the amendment that might
have looked natural, would have rejected 2 July more firmly.

## Consequences

1. **The modelling runs on Pléiades Neo only.** Both deliveries pass, and the plan's
   native 1.2 m experiment and its 10 m aggregation are unaffected.
2. **The matched-sensor comparison cannot run.** The protocol fits on July and evaluates on
   September; with no July Sentinel-2 scene there is nothing to fit, so the September scene
   cannot be evaluated either.
3. **The 31 August scene is kept for labelled diagnostics only**, such as the cross-sensor
   radiometric comparison, never for fitting.
4. **June Sentinel-2 is unavailable on any date**, so the supplementary June experiment
   lapses. The best June candidate is four times the glint limit.
5. Recovering a Sentinel-2 comparison needs a survey timed to a low-glint overpass, which
   is a new field campaign, not a processing change.
6. *Update 2026-10-02, exploratory ([revision-r1.md](revision-r1.md)).* Across 2017–2025,
   a week either side of the June or July survey dates held a Sentinel-2 overpass likely
   to pass in 2 of 9 years, against 7 of 9 for September; R080 never qualified in June or
   July. SWIR residual glint correction improves the clean scenes but cannot rescue a
   glinted one under the registered mask.
   A glint-tolerant mask (r3, [revision-r3.md](revision-r3.md)) admits their water, but no
   glinted June or July scene passed its registered acceptance checks. Those checks have
   known weaknesses, so this does not show that processing cannot recover such scenes.

## Reproduce

ACOLITE is pinned at `64a02ff` (`20260421.0-86-g64a02ff`), with settings in
`registration/acolite/`.

```bash
PYTHONPATH=. python benchmark/run_benchmark.py glint-screen
PYTHONPATH=. python benchmark/run_benchmark.py fetch --env-file .env --top 2
PYTHONPATH=. python benchmark/run_benchmark.py prepare
PYTHONPATH=. python benchmark/run_benchmark.py scene-qc
PYTHONPATH=. python benchmark/run_benchmark.py export
```

Downloads need `CDSE_USERNAME` and `CDSE_PASSWORD`. The screen reads public Sentinel-2
L2A through the Planetary Computer and needs no credentials.
