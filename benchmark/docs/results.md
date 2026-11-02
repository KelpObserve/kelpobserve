# Sesimbra habitat benchmark v1: results

Executed 2026-09-20. **Pléiades Neo does not predict kelp percentage cover at Sesimbra
beyond what depth alone implies, and a July-fitted model does not transfer to September.**
Depth-only is the best candidate on the development campaign at both spatial scales; every
spectral candidate is worse. On the held-out September campaign every candidate
overpredicts by 35 to 50 percentage points, and the error is almost entirely a level shift.

This is a negative result from a pilot-scale dataset, which the plan counts as valid. It is
not evidence that satellite cover mapping is infeasible in general; it measures what these
deliveries, at this site, support.

*Update 2026-10-02.* An exploratory re-analysis after an external review
([revision-r1.md](revision-r1.md)) keeps the ranking below under consistent QC. In an injection test,
the registered models did not reliably recover injected spectral signals on these units.
Read the negative result as inconclusive about the sensor, not as evidence against it.

## What ran, and what could not

| Branch | Status |
|---|---|
| Uncorrected Pléiades Neo, common 10 m grid and native 1.2 m | **Ran** |
| Sentinel-2 on the common grid | Unavailable: glint and near-shore masking, [scene selection](scene-selection.md) |
| June Sentinel-2 supplement | Unavailable: no qualifying scene on any date |
| Water-column-corrected variants | Unavailable: reference QC fails on both dates, [optics readiness](optics-readiness.md) |

Three of the four registered comparisons are therefore reported as unavailable with
evidence, as the plan requires. What follows is the fourth.

## Development: July, leave-one-transect-out

Thirty eligible quadrats, forming 18 cells on the 10 m grid. All three folds met the
registered minimum training support. Errors are in percentage points of cover.

| Candidate | 10 m grid MAE | RMSE | Bias | Native MAE | RMSE | Bias |
|---|---:|---:|---:|---:|---:|---:|
| Training-set mean | 38.5 | 42.8 | +2.0 | 36.0 | 38.8 | +0.0 |
| **Depth-only ridge** | **32.3** | 40.7 | +1.1 | **26.3** | 34.1 | +0.2 |
| Blue/green/red ridge | 35.9 | 42.3 | +7.7 | 39.3 | 43.5 | −0.7 |
| Blue/green/red + depth ridge | 41.9 | 51.2 | +13.9 | 33.4 | 41.9 | −0.3 |
| Spectral + depth forest | 37.8 | 42.1 | +2.9 | 27.7 | 33.6 | +2.2 |

Paired against depth-only on the same observations, every spectral candidate is worse: by
3.7 to 9.6 points on the 10 m grid, and by 1.4 to 13.0 points natively. Adding reflectance
to the depth model does not help; on the 10 m grid it actively hurts.

Depth-only beats the training mean by 6 points on the grid and 10 natively, which is the
expected consequence of the survey design: the shallow stratum averages about 70% cover and
the deep stratum about 26%, so a model that knows depth reproduces most of the variation
available.

## Held-out: September, sealed fits applied once

The fits were sealed before the September targets were opened; the evaluation stage is the
only code that can read them. Thirteen cells on the 10 m grid and 22 native quadrats were
eligible; one native quadrat is kelp-positive with unrecorded cover and is excluded rather
than counted as zero.

| Candidate | 10 m grid MAE | Bias | Native MAE | Bias |
|---|---:|---:|---:|---:|
| Training-set mean | 41.2 | +41.2 | 42.0 | +39.4 |
| Depth-only ridge | 46.2 | +46.2 | 50.1 | +50.1 |
| Blue/green/red ridge | 37.9 | +34.7 | 44.9 | +42.4 |
| Blue/green/red + depth ridge | 36.5 | +36.5 | 45.4 | +45.4 |
| Spectral + depth forest | 49.2 | +49.2 | 38.5 | +37.3 |

Mean absolute error and bias are nearly identical for every candidate, which means the error
is a level shift rather than misplaced spatial detail. Observed cover fell from a July mean
of 48% to a September mean of 8.6%: of 21 scored September quadrats, nine are zero and eight
are 5%. *Saccorhiza polyschides* is a seasonal alga, and the models carry July's level into
September.

Depth-only, the best candidate in July, is the worst in September, because the steep
depth–cover gradient it learned in July no longer holds. Its September error is far larger
in the shallow stratum (61 points natively) than the deep one (38), which is where the
seasonal decline is largest.

No candidate was selected using these results, and no fit was retuned.

## What this does and does not establish

1. **It does not establish that reflectance carries no cover information here.** The
   reflectance itself is compromised: ACOLITE's PNeo output carries an additive term of
   roughly 0.01 to 0.04 over water, larger than the seabed contrast at these depths. A
   spectral model fed that signal is not a fair test of the physics.
2. **The statistical support is thin.** Eighteen cells, three transects 25 to 53 m apart,
   inside a 120 by 140 m area. Folds are neighbouring transects, not independent sites.
3. **The comparison is against a strong baseline for the wrong reason.** With only two depth
   strata, depth-only is close to a stratum-mean model, and the design confounds depth with
   cover.
4. **September tests transfer across a season, not across space.** Its transects were re-laid
   21 to 102 m from July's and sampled after a large seasonal decline, so the holdout mixes a
   temporal shift with a positional one.
5. Eight of 30 September events remain provisional: four deep quadrats sitting on 5–8 m of
   lidar depth fail the depth check, and four shallow quadrats lack lidar entirely. A ninth
   provisional event is in June. Positions came from a surface float the divers dragged
   along, which fits the 10–30 m offsets these flags imply; the CCMAR query that could have
   supplied corrected positions was closed unsent on 2026-09-24.
6. **The September images looked through different water, not only at less kelp.** This
   was added on 2026-09-27, after the results, from the CCMAR sonde casts. It is context,
   not a registered analysis (`kelp_observe/docs/aquaserve/sesimbra-sonde-context.md`).
   - **The casts match the images:** they were taken on 7 July and 2 September, a day
     before each Pléiades Neo image.
   - **July:** the water was stratified, with 15.2 °C at 15–17 m below a 2.5 °C
     thermocline near 5 m.
   - **September:** the water was mixed and 17.2 °C at the same depth. Chlorophyll
     fluorescence there was about five times July's, and phycoerythrin about three times.
   - **The shallow cast stopped at 5.9 m,** so the shallow stratum is less well
     characterised; there the rise was smaller (0.13 to 0.21 RFU for chlorophyll).

   Less clear water changes the spectral predictors independently of cover. The depth-only
   model is immune to that, yet its September error is among the largest (46 points on the
   grid, 50 natively), so the seasonal decline alone explains most of the level shift. For
   the spectral candidates, a change in water clarity is a plausible second contributor, of
   unknown sign. The sonde has no light sensor, so it cannot be quantified.

## What would change the answer

- Surface correction that works without SWIR, or a SWIR-bearing sensor, so the spectral
  predictors carry seabed signal rather than a residual.
- Surveys timed to low-glint overpasses, which would restore the Sentinel-2 comparison.
- A design with more than two depth strata, and transects separated enough that folds are
  spatially independent.
- Within-season repeat sampling, so a transfer test is not confounded by the seasonal cycle.
- A light sensor (PAR or spectral irradiance) on the sonde casts at each image date, so a
  change in water clarity can be separated from a change in cover.

## Reproduce

```bash
PYTHONPATH=. python benchmark/run_benchmark.py freeze
PYTHONPATH=. python benchmark/run_benchmark.py fit
PYTHONPATH=. python benchmark/run_benchmark.py evaluate
PYTHONPATH=. python benchmark/run_benchmark.py export
```

The input lock pins 96 files: sources, processed imagery, the registration and the coastal
source code. Fitting refuses to run if any of them changed, and evaluation refuses to run
twice. Numbers behind these tables are in
[the portable results export](results.json).
