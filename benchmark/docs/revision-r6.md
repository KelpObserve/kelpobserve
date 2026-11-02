# Sesimbra habitat benchmark: glint correction and date combination against the envelope (r6)

Executed 2026-10-02, registered in
[`registration/revision-r6.json`](../registration/revision-r6.json) before anything was
computed. The design is the user's specification of the same day. The recovery criterion
is r4's, unchanged. No model is fitted to observed cover.

> Across the tested Sesimbra campaign geometries, neither glint correction, two-date
> Sentinel-2 averaging, nor removal of fixed background offsets enabled reliable recovery
> of the simulated kelp-cover signals. The experiments identify strong dependence on
> observation conditions, spatial variation and model formulation. They provide a
> quantified limitation of this campaign benchmark, rather than a general limit on
> satellite kelp detection.

On real backgrounds no tested model reaches the recovery criterion, in any condition, on
either sensor, overall or when its predictions are scored separately within either depth
stratum. The models were fitted on both strata together; a model trained for one stratum
alone was not tested.

- **Enabling PNeo's glint correction does not attenuate an injected perturbation**
  relative to v1's processing, for the pattern tested.
- **It does not unlock recovery.**
  - The corrected water is darker than the injection model predicts for the tested
    substrates under the assumed water.
  - The spread that remains is still larger than recovery tolerates.
  - **Removing its median per stratum does not rescue recovery** (offset control,
    section E, registered as an amendment). The between-unit variation was kept. Of 48
    setting, contrast and model combinations, two improve, eight deteriorate
    substantially and none recover.
- **Combining two dates reduces noise about as independent noise predicts, but two dates
  are not enough.** On the real scenes, part of the unexplained variation does not average
  out.

## A1: does enabling the correction attenuate an injected signal? No, for this pattern

One synthetic pattern (bare substrate, 10 m patches, draw 0) was added to the July product
itself, before any processing, and the product run through ACOLITE with V0 and V1 settings.

| Band, stratum | Pixels | V1 response on V0 response: slope, r | V0 response on injected reflectance |
|---|---:|---|---:|
| Blue, shallow | 6,315 | 1.000, 1.000 | 1.018 |
| Blue, deep | 3,135 | 1.000, 1.000 | 1.013 |
| Green, shallow | 6,279 | 1.000, 1.000 | 1.035 |
| Green, deep | 3,465 | 1.000, 1.000 | 1.053 |

- **The aerosol estimate is unchanged** (0.0781 with and without injection).
- **Every pixel outside the pattern is bit-identical to the original runs,** so ACOLITE is
  deterministic here.
- **The registered criterion is met, for this pattern.** The NIR-referenced correction
  subtracts a term computed from NIR alone, which the injection does not touch. So adding
  the signal after V1 is equivalent to adding it after V0.
- **That is narrower than equivalence with injecting before all processing.** The full
  processing passes a product-level perturbation through with a gain of 1.013–1.053 (last
  column). A signal added directly to V0 or V1 reflectance, as A does, corresponds to a
  product-level perturbation 1–5% smaller. The transfer would differ again for a signal
  entering at the top of the atmosphere.

## A: paired PNeo V0/V1 recovery

Identical units (29 native, 18 cells), identical draws and masks.

| Bare substrate, 10 m patches, no position error | V0 | V1 |
|---|---:|---:|
| v1 additive ridge, native: share (shallow, deep) | 0.45 (0.72, 0.21) | 0.48 (0.74, 0.19) |
| v1 additive ridge, 10 m | 0.43 (0.66, 0.14) | 0.29 (0.34, 0.33) |
| Depth-aware ridge, native | 0.01 | 0.07 |
| Oracle, native: share (median gain) | 0.00 (−21 pp) | 0.00 (−28 pp) |
| Null, additive ridge | 0.05 | 0.05 |

- **Understory stays at the null on both versions** (0.04–0.07).
- **Paired over all 100 real-background conditions, correction changes recovery little
  on the median, in both directions.**
  The additive ridge's share moves by a median of +0.03 (range −0.16 to +0.08), the
  depth-aware ridge's by +0.02, and the oracle's by nothing.
- **The ideal controls on the same units reproduce r4.**
- **Why correction does not help is not settled.** V1 water is darker than either
  modelled substrate (C), but removing that median offset per stratum does not rescue
  recovery (E). The variation within strata remains, and depth-dependent or spatially
  varying mismatch was not tested.

## B: Sentinel-2, single dates and their combination on fixed cells

Common cells: 14 (valid on both 2 and 9 July). Union: 17.

| Bare substrate, 10 m patches, no position error | 2 Jul | 9 Jul | Combined (14) | Combined (17) |
|---|---:|---:|---:|---:|
| Real background, additive ridge | 0.35 | 0.37 | 0.37 | 0.41 |
| Real background, depth-aware ridge | 0.10 | 0.03 | 0.05 | 0.07 |
| Ideal, noise 0.0005 per pixel per date, depth-aware ridge | 0.38 | 0.27 | 0.56 | 0.58 |
| Ideal, noise 0.001, oracle | 0.22 | 0.20 | 0.38 | 0.29 |
| Null, additive ridge | 0.18 | 0.13 | 0.04 | 0.11 |

- **With independent noise, combining two dates helps as expected** (0.38 → 0.56, 0.22 →
  0.38) but stays below 0.8. Two dates cut independent noise by at most √2.
- **On the real scenes the variance gain is about zero** (−0.07 to +0.14, within the null
  rates). That fits r5: about half the unexplained variation is stable across dates, and
  averaging cannot remove it.
- **The coverage gain from 14 to 17 cells is small and inconsistent** (−0.17 to +0.13). The
  extra cells are 9 July's alone, near the cliffs.

## C: background bias at the units (the registered interpretation)

How far the real unit reflectance departs from the injection model's background, in green,
and the bottom-reflectance bias after inversion with the assumed water (median by stratum):

| Units | Shallow, vs understory: rhos gap / bottom bias | Deep, vs understory: rhos gap / bottom bias |
|---|---|---|
| PNeo V0, native | −0.007 / −0.07 | −0.003 / −0.36 (SD 1.2) |
| PNeo V1, native | −0.017 / −0.18 | −0.021 / −3.3 (SD 1.8) |
| S2 2 July | **−0.0007 / −0.012** | +0.003 / +0.59 (SD 0.27) |
| S2 9 July (corrected, glinted) | −0.016 / −0.21 | −0.012 / −2.4 (SD 0.7) |

The understory endmember is 0.034 (PNeo) and 0.041 (S2), so a bottom bias of −0.07 or
worse means the inverted seabed is negative.

- **Only Sentinel-2 on 2 July, in the shallow stratum, matches the model's background**
  closely.
- **Glint-corrected scenes are darker than either tested substrate (understory, bare)
  predicts under the assumed water.** Inverting them gives negative bottom reflectance, and
  at 17 m the inversion multiplies the error about a hundredfold.
- **"Darker" is relative to the assumed water, substrate and observation model.** It does
  not by itself identify an atmospheric-correction error: the water could be less clear
  than assumed, the substrates darker, the corrections could overshoot, or several of
  these at once.
- **This mismatch can strongly impair the inversion.**
  - The offset control (E) removed its median per stratum. That narrows the oracle's
    error in most settings but rescues recovery nowhere.
  - The specific conclusion: removing this median offset alone does not rescue recovery.
  - Depth-dependent or spatially varying mismatch remains, and was not tested.

## E: offset control (registered amendment, 2026-10-02)

An intentionally favourable diagnostic, not an operational correction:

1. **Expected background.** The injection's physics with no kelp is computed for each
   scene, scale and background substrate, through the same 1 m grid, area averaging,
   PSF and unit footprints as the signal.
2. **Offset.** The median real-minus-model residual per depth stratum and band is computed
   on the uninjected background, with no cover read.
3. **Removal.** That offset is subtracted from every unit, keeping the differences between
   units.
4. **Rerun.** The same injections, models, nulls and criterion are run, paired with A and
   B.

**Offsets (green, shallow / deep).** They agree with section C within about 0.003, so averaging
the model across a unit's footprint rather than at its mean depth made little difference
here.

| Setting | Understory | Bare |
|---|---|---|
| PNeo V0, native | −0.0072 / −0.0029 | −0.031 / −0.005 |
| PNeo V1, native | −0.0170 / −0.0207 | −0.048 / −0.023 |
| S2 2 July | −0.0007 / +0.0031 | −0.032 / +0.002 |
| S2 combined (common cells) | −0.0088 / −0.0043 | −0.036 / −0.006 |

**Outcome, by the registered rules** (reference set: 5, 10 and 25 m patches, no position
error):

- **Two improve, eight deteriorate substantially, and none recover.** These are median
  paired changes over the reference set. 38 of the 48 combinations change by less than
  10 points either way.
  - **Improve:** the oracle on PNeo V0 10 m cells (0 → 0.04–0.10), and the additive ridge
    on the combined Sentinel-2 dates for understory (0.07–0.09 → 0.20, against an
    offset-removed null of 0.15).
  - **Deteriorate by 17–43 points:** every one of the eight is the pooled additive ridge
    with bare substrate, one in each setting.
- **The registered rule's labels are a naming weakness** (registration erratum). It
  counts 2 as "improves but fails" and 46 as "barely changes", so its "barely changes"
  includes the eight large decreases. The formal classification is kept as registered.
- **No condition reaches the criterion after offset removal,** overall or when scored
  separately within either stratum.

| Bare substrate, 10 m patches, no error | Before: share (median gain) | After: share (median gain) |
|---|---|---|
| PNeo V0 native, oracle | 0.00 (−21.1 pp) | 0.00 (−12.8 pp) |
| PNeo V1 native, oracle | 0.00 (−28.0 pp) | 0.00 (−9.5 pp) |
| S2 combined, oracle | 0.00 (−35.2 pp) | 0.00 (−14.0 pp) |
| PNeo V0 native, depth-aware ridge | 0.01 | 0.02 |
| PNeo V0 native, additive ridge | 0.45 (shallow 0.72) | 0.05 (shallow 0.20) |
| PNeo V1 native, additive ridge | 0.48 (shallow 0.74) | 0.05 (shallow 0.23) |

- **Removing the median offset changes the size of the physics-based error, not the
  recovery outcome.** It narrows the oracle's deficit against depth-only by
  8–21 points in most settings; on 2 July Sentinel-2 for bare substrate it widens
  instead (−4.8 to −8.4). The oracle still loses to depth-only in almost every draw, in
  the shallow stratum too (shares 0.04–0.26). The depth-aware ridge barely moves.
- **The variation between units within each stratum remains.** On PNeo native units in
  the shallow stratum it is about 0.011 in green against the bare model. A fully
  kelp-covered unit changes reflectance by about 0.023 there.
  - The spread is not wholly sensor noise. It includes real substrate variation, model
    mismatch (including any depth-dependent part) and radiometric error, which this
    simulation counts together as noise against a synthetic target.
- **The additive ridge's synthetic recovery is sensitive to the transformation**
  (0.45–0.48 → 0.05).
  - A per-stratum constant leaves every within-stratum difference intact, but changes
    the relation between strata, and a pooled, regularised model responds to that.
  - A ridge fitted within one stratum would standardise the constant away and give
    unchanged predictions.
  - So the earlier result is best described as partial recovery of synthetic signals
    under the original background, which was sensitive to the pooled model's treatment
    of between-stratum differences.
  - It was never evidence of observed kelp detection, since the targets were synthetic.
    The offset experiment does not retroactively invalidate it.
  - A small paired invariance check (pooled against within-stratum ridge, before and
    after the offset) would test the mechanism. It needs no new envelope.
- **In the terms of the specification,** removing this median offset is insufficient.
  The remaining within-stratum variation, the deep stratum and the pooled design still
  defeat the models tested.

## D: matched spatial windows

Square windows of equal size offshore and in the habitat replace r5's unequal comparison.
The offshore region is narrow, so only 100 m windows are informative (23–29 offshore;
at 300 m only 1–2 qualify). Green, median over windows:

| Scene, unit | Offshore local SD (after NIR) | Habitat bins local SD (after NIR) | D1m met |
|---|---|---|---|
| PNeo 8 Jul V0, native | 0.0010 (0.0006) | 0.0052–0.0074 (0.0014–0.0022) | no |
| S2 2 Jul, 10 m | 0.0008 (0.0008) | 0.0007–0.0010 (0.0006–0.0008) | yes at 10–30 m; not at 5–10 m (0.78) |

- **Matched windows confirm r5 on both sensors.** PNeo's native-unit habitat varies 5–7
  times more than equal offshore windows (2–4 times after NIR). Sentinel-2's habitat cells
  vary about as much as offshore; slightly more at 5–10 m.
- **D3m (seabed candidate) is not met anywhere.** Sentinel-2's habitat variation is not
  larger than offshore, even though about 40% of it is stable across dates (r5).

## What r6 establishes

1. **Enabling PNeo's NIR glint correction does not attenuate the tested perturbation**
   (A1), but does not put this design within reach of recovery (A). Its corrected water
   is darker than the injection model predicts (C). Removing that mismatch's median per
   stratum does not rescue recovery (E).
2. **Two July Sentinel-2 dates are not enough,** and on the real scenes averaging gains
   little, because much of the variation is stable (B and r5).
3. **The binding constraints are now specific:**
   - **Unit-level variation within strata, above tolerance even after correction and
     after the median offset is removed (E).** It mixes real substrate variation, model
     mismatch and radiometric error; this simulation cannot separate them.
   - **The deep stratum,** where any inversion error is amplified about a hundredfold.
   - **The mismatch between corrected water and the model backgrounds** (C; also QAA's
     red absorption in r3). Removing its median per stratum narrows the physics-based
     error but does not rescue recovery (E). Depth-dependent or spatially varying parts
     were not tested.
4. **The shallow stratum's 0.72–0.74 for bare substrate** (additive ridge, scored on
   shallow units, fitted on both strata) is sensitive to the per-stratum offset (E). No
   configuration tested reaches the registered criterion.
5. **For presenting the project:**

   > Across the tested Sesimbra campaign geometries, neither glint correction, two-date
   > Sentinel-2 averaging, nor removal of fixed background offsets enabled reliable recovery
   > of the simulated kelp-cover signals. The experiments identify strong dependence on
   > observation conditions, spatial variation and model formulation. They provide a
   > quantified limitation of this campaign benchmark, rather than a general limit on
   > satellite kelp detection.

   This is a negative methodological result about this campaign benchmark. It is useful
   as it stands and needs no detection result to support it.

## Reproduce

```bash
PYTHONPATH=. python -m benchmark.revision_r6 check     # copies the product, injects, two ACOLITE runs (~2 min, 1.4 GB)
PYTHONPATH=. python -m benchmark.revision_r6 pneo      # ~3 min on 5 cores
PYTHONPATH=. python -m benchmark.revision_r6 s2        # ~3 min on 5 cores
PYTHONPATH=. python -m benchmark.revision_r6 bias
PYTHONPATH=. python -m benchmark.revision_r6 windows   # ~1 min
PYTHONPATH=. python -m benchmark.revision_r6 offset    # E, the amendment: ~6 min on 8 cores
```

Code: `benchmark/revision_r6.py`, tests in `tests/test_revision_r6.py`. Outputs under
`out/r6/<stage>/`. The injected product copy is `.cache/pneo-r6/` and its ACOLITE outputs
`.cache/acolite-r6/`. The registration's errata record four implementation points: the
300 m window shortfall, the small-sample SD, the first check run stopping on a
floating-point grid comparison, and the √2 limit of two dates.
