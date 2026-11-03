# PNeo glint correction and water assumptions

Development review, 3 November 2026. No cover targets were read or scored.

The July scene supports changing the glint spectral slopes: the existing correction
subtracts too much blue and green signal relative to NIR within the reference blocks.
An empirical correction reduces negative values and the reversed residual relationship
with NIR. Its absolute reflectance remains unqualified. Retaining a fifth-percentile
NIR baseline reduces negative values further but transfers poorly between reference
areas. No tested branch is adopted for cover retrieval.

## Existing processing and implemented changes

The pinned ACOLITE alternative correction clips negative reference NIR to zero and
subtracts that amplitude from visible bands using LUT spectral ratios. It corrects
pixels with reference NIR below 0.1 and assumes zero non-glint NIR background.
The ratios are almost spectrally flat, particularly in July. The reviewed code is
`acolite/acolite_l2r.py`, the `glint_correction_method == 'alternative'` branch.
Its behaviour matches its settings; this review found an unsupported scene assumption,
without identifying an upstream arithmetic defect.

The new `optical_corrections` module requires an explicit NIR baseline, estimates
slopes from spatial reference blocks, checks calibration support, and retains negative
outputs for diagnosis. The review runner applies every variant at native 1.2 m
resolution before the same masked area average to 10 m. It reproduces the saved
V1 blue/green/red/coastal grids to better than 3 × 10⁻⁹ rhos.

Correction and averaging cannot be interchanged: mean(max(NIR, 0)) differs from
max(mean(NIR), 0). The initial implementation stopped at the reproduction check,
before empirical fitting or scoring, with a maximum discrepancy of 0.000475.
The native-resolution clarification was recorded before the completed comparison.

## Reference support and held-out comparison

The previously fixed deep-reference and habitat polygons are retained. Calibration
uses reference cells only, with per-block centred regression and equal weight for
informative 100 m blocks. Each block needs at least 10 cells and a NIR p95−p05
span of 0.002; at least six blocks are required. The western and eastern reference
sets are separated by an omitted 100 m block column and evaluated in both directions.
The centre column remains in the full-reference development preview.

The regression assumes water/background reflectance is approximately constant within
each informative block and that NIR-correlated variation is dominated by surface glint.
Covarying water properties, adjacency or atmospheric residuals can bias these slopes.
Block centring removes between-block offsets but cannot eliminate such within-block effects.

Full-reference spectral slopes, ordered coastal / blue / green / red:

| Date | Existing LUT | Empirical preview |
| --- | --- | --- |
| 20260708 | 0.972 / 0.980 / 0.988 / 0.994 | 0.710 / 0.806 / 0.883 / 0.954 |
| 20260903 | 0.898 / 0.928 / 0.956 / 0.976 | 0.780 / 0.901 / 0.976 / 1.018 |

The following scores use held-out reference cells. Percent negative is blue / green.
The NIR-support percentage counts cells within the training p1–p99 range; every
held-out cell remains in the headline score. NIR support does not establish geographic
or optical validity. `empirical_zero` isolates the slope change from the baseline change.

| Date | Train to test | Test cells | NIR support | Training NIR p05 | V1 negative | Empirical zero negative |
| --- | --- | --- | --- | --- | --- | --- |
| 20260708 | west to east | 1176 | 58.2% | 0.02414 | 0.0% / 4.7% | 0.0% / 2.1% |
| 20260708 | east to west | 1315 | 4.6% | 0.00373 | 56.8% / 85.2% | 0.0% / 0.5% |
| 20260903 | west to east | 1176 | 21.7% | 0.00562 | 0.0% / 0.0% | 0.0% / 0.0% |
| 20260903 | east to west | 1025 | 5.3% | 0.00317 | 0.0% / 0.0% | 0.0% / 0.0% |

July's west-to-east baseline is 0.02414, against 0.00373 in the reverse direction.
The brighter western reference does not contain a comparable low-glint baseline.
Applying its p05 baseline to the eastern reference leaves substantial NIR-coherent
signal. In the reverse direction, only about 5% of test cells fall within the
calibration NIR range. A good-looking full-reference preview cannot resolve this
transfer problem. September's eastern calibration has only six informative blocks,
at the recorded support minimum, with broad slope variation.

The zero-baseline empirical slopes also leave red-band problems. In the full-reference
preview, about 89% of July and 97% of September deep-reference red cells are negative.
Changing the visible slopes therefore does not repair all bands. With the p05 variant,
the full-reference habitat preview still has negative red values on about 17% of
September cells. These previews share calibration data and are not held-out evidence.

## Water assumptions

Frozen r7 scaled winter non-water absorption and particle backscatter with one factor
to match KD490. A single attenuation value cannot identify both components. The new
`water_constraints` module sets backscatter from BBP443, then solves absorption to
match KD490. It rejects infeasible KD targets, wavelength extrapolation and interpolation
through unphysical reference components. It reports CDM conflicts explicitly.

Nominal coefficients at the PNeo raster wavelengths, in m⁻¹:

| Date | Band | r7 a | KD+BBP a | r7 bb | KD+BBP bb |
| --- | --- | --- | --- | --- | --- |
| 20260708 | blue | 0.04961 | 0.06530 | 0.00822 | 0.00218 |
| 20260708 | green | 0.07919 | 0.08581 | 0.00630 | 0.00131 |
| 20260903 | blue | 0.12350 | 0.17653 | 0.02188 | 0.00454 |
| 20260903 | green | 0.11039 | 0.13279 | 0.01760 | 0.00327 |

These candidates retain the winter spectral shapes. KD490 and BBP443 constrain their
amplitudes; neither supplies a locally measured summer spectrum. The winter reference
has negative non-water absorption at 704 nm under the current pure-water table. The
old band extension clips it silently; the new helper discloses it and refuses requests
that depend on it. Blue/green interpolation does not use that point, so this finding
does not explain the present blue/green failure. Coastal 439 nm also lies outside the
442–704 nm reference support and is excluded from the water-state construction.

The product-qualified Copernicus pixel is 0.891 km from the reef, with daily 1 km
support. Same-day retrievals do not establish simultaneous or spatially uniform reef
water. Provider percentage uncertainties on July KD490/BBP/CDM are approximately
36/20/57%; September values are 58/40/60%. The nine combinations of nominal and
±provider error are sensitivity cases, without probability weights. The products
share retrieval inputs, so independence of their errors is unestablished.

CDM443 is a lower bound on total non-water absorption. July's nominal constructed
non-water absorption (0.0636 m⁻¹) is below nominal CDM (0.0795 m⁻¹), but above its
error-adjusted lower bound (0.0342 m⁻¹). `input_constraints_consistent` therefore means
consistent with that declared tolerance; `cdm_nominal_conflict` preserves the nominal
disagreement. All three low-KD July sensitivity states conflict even with that lower
bound. SPM supplies context only; no local mass-specific optical conversion was assumed.

## Optical closure

Closure checks both blue and green against the finite-depth Lee signal between black
and 0.6-reflectance beds, allowing ±0.001 rhos. The tolerance is a retained diagnostic,
without a calibrated probability interpretation. Water states are fixed from the
products, without fitting image residuals. The reference depth and water level remain
assumptions. These checks do not independently validate either the water products or
the corrected imagery.

Fraction of held-out reference cells compatible with the nominal KD+BBP water state:

| Date | Train to test | Existing V1 | Empirical zero | Empirical p05 |
| --- | --- | --- | --- | --- |
| 20260708 | west to east | 59.6% | 44.3% | 0.0% |
| 20260708 | east to west | 1.1% | 1.0% | 0.1% |
| 20260903 | west to east | 0.2% | 0.3% | 0.0% |
| 20260903 | east to west | 0.4% | 0.1% | 0.0% |

The p05 candidate is close to zero compatibility in both scenes. July empirical-zero
closure differs strongly between the two reference areas; September remains very low.
The p1 and p10 baseline sensitivities do not resolve this. No joint glint/water choice
has been selected by its best closure fraction. The complete sensitivity grid is in
[water-closure.csv](water-closure.csv), including the conflicting water states.

## Next processing step

Retain the slope correction as a development candidate and constrain the additive
baseline independently. Suitable evidence would be a demonstrably low-glint water
reference with comparable water properties, checked across dates and viewing geometry,
or a separately qualified atmospheric/surface-reflection model. A reference percentile
alone is insufficient in these scenes.

Evaluate the aerosol and land-adjacency treatment next, with the same spatial reference
split and unchanged masks. The present experiment changes residual glint only. Both
original branches used adjacency correction disabled. Sky illumination near the cliffs
also remains unqualified. A favourable blue/green regression cannot resolve these terms.

Before another cover run, fix a radiometric branch using held-out optical evidence and
carry water-shape and representativeness uncertainty through it. Retain original r7
results and its failed compatibility outcome. No new habitat map or cover score follows
from this review.

## Reproduction and files

Run `kelpobserve pneo-optical-review --output /absolute/path/to/new-output` using the
existing local cache configuration. The destination must be new. Requires the archived
optical diagnostic and native PNeo V0 rasters; no download is performed.

- [design.json](design.json): development choices, prior information and the implementation clarification.
- [evaluation.csv](evaluation.csv): all bands, variants, held-out directions and preview regions.
- [calibration.json](calibration.json) and [spatial-assignments.csv](spatial-assignments.csv): slopes and support.
- [water-states.json](water-states.json): coefficients, CDM checks and rejected assumptions.
- [candidate-grid-diagnostics.csv.gz](candidate-grid-diagnostics.csv.gz): full-reference p05 preview, explicitly unvalidated.
- [provenance.json](provenance.json), [verification.json](verification.json) and [outputs.json](outputs.json): hashes and checks.

The regression method follows the NIR differential-glint approach of
[Hedley et al. (2005)](https://doi.org/10.1080/01431160500034086), with this experiment's
recorded block estimator and percentile baseline. Its correction equation is also
described in the [ESA Sen2Coral tutorial](https://eo4society.esa.int/wp-content/uploads/2021/04/OCEA05_Sen2Coral_GBR_Tutorial.pdf).
Those sources do not validate the assumptions for these PNeo scenes.
