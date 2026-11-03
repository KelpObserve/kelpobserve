Sesimbra r7: predictions and optical attributes at scored quadrats
================================================================

Source: benchmark/reports/r7-local-cover-20261102-r3
Export date: 2026-11-02

The main inspection table contains 91 quadrat/image associations across 50
distinct field events. It uses V1 radiometry, Sentinel-2 at 10 m and Pleiades
Neo at the registered native extraction support. Every closest fit in this
table failed the nominal reflectance-compatibility rule. Compatible cover
estimates are empty. The closest-fit numbers remain available for diagnosis.

                              Scored footprints   Associated quadrats
Sentinel-2, July                      14                    24
Sentinel-2, September                 12                    20
Pleiades Neo native, July             29                    29
Pleiades Neo native, September        18                    18

Several quadrats share one Sentinel-2 footprint. These repeated predictions
must not be counted as independent observations or used to recompute the
registered footprint-weighted metrics. Use footprints-all.csv for scoring.
Pleiades Neo native extraction takes the mean of a 3 x 3 multispectral
pixel window, approximately 3.6 m across. GeoJSON points mark the quadrat
locations; they do not represent the full prediction footprint.

Files
-----

scored-quadrats-primary.csv
  Main inspection table: one row per scored quadrat and image, at the sensor
  supports above. Includes the diver observation, footprint mean, closest
  fit, optical diagnostics, sensitivity ranges and baseline predictions.

scored-quadrats-primary.geojson
  The same rows as EPSG:4326 point features. A location can appear once for
  each sensor. The properties identify each image and footprint.

unscored-quadrats-primary.csv
  Thirteen candidate quadrat/image associations excluded from these scores,
  with extraction and analysis exclusions. Missing estimates remain empty.
  This covers r7 candidate events; field events excluded before r7 extraction
  remain documented in the source inventory.

all-quadrat-associations.csv
  All 312 associations, including V0 radiometry, Pleiades Neo 10 m comparisons
  and excluded candidates. Includes the original extraction metadata.

footprints-all.csv
  One row per original unit (228 rows), preserving the scoring support.
  Includes all original unit results and additional optical diagnostics.

optical-inputs-by-scene.csv
  Water coefficients, solar angle, water level, wavelengths and fixed kelp
  and background spectra assigned to each image/variant. These inputs are
  shared across footprints in a scene.

copernicus-context-by-scene.csv
  KD490, BBP443, CDM443 and suspended particulate matter from the qualified
  fixed Copernicus neighbour. Only KD490 entered this retrieval. The other
  quantities are candidate context, with used_in_r7=false. Representativeness
  at the reef is unestablished; these are satellite products, without local
  water sampling. Provider percentage errors are not confidence intervals.

export_quadrats.py
  Reproducible export and numerical verification. It reads the frozen run
  and its input inventory, reproduces diagnostic forward calculations, and
  leaves the run and its retrieval method unchanged.

verification.json
  Source checksums, code checksums, counts and numerical consistency checks.

Reading the columns
-------------------

event_id: field quadrat identifier. For example, july-deep-T3-Q2.
unit_id/source_unit_index: original prediction footprint and source row.
scene/image_date/variant: image and radiometric processing variant.
field_date_workbook/field_date_transects: source dates, preserved separately.
field_date_discrepancy: inventory flag; dates have not been silently reconciled.
field_cover_pct: individual diver observation of Saccorhiza polyschides cover.
field_cover_status: measured cover or an explicit absence-derived/recorded zero.
footprint_observed_mean_cover_pct: mean of usable field observations associated
  with the footprint. The quadrats sample the footprint; their mean is not a
  complete-pixel habitat census.
n_quadrats/n_scored_quadrats: associated and usable field quadrat counts.

closest_fit_cover_pct: spectral mixture percentage giving the smallest mean
  squared blue/green reflectance residual on the frozen 0.1 percentage-point
  grid. A boundary estimate of zero does not establish kelp absence.
compatible_fit_cover_pct: closest fit among fractions satisfying the nominal
  tolerance in both bands. Empty when no fraction qualifies.
fit_status/nominal_compatible: whether a fraction met the declared rule.
boundary: closest fit is at zero or 100 percent.
compatibility_tolerance_rhos: fixed absolute tolerance, 0.001 in each band.
  This diagnostic threshold does not establish cover accuracy.

footprint_water_depth_m: depth used by the retrieval, including acquisition
  water level. It can differ from the individual quadrat's lidar depth and
  the diver's nominal depth label, both retained in the full export.
observed_blue/green/red/nir_rhos: processed image surface reflectance,
  dimensionless. Only blue and green enter the physical cover retrieval.
blue_nm/green_nm: actual raster band centres used for water coefficients.
predicted_blue/green_rhos: forward-model reflectance at the closest fit.
model_minus_observed_blue/green_rhos: signed per-band residual at that fit.
residual_rhos: largest absolute blue/green residual at the closest fit.
model_min/max_blue/green_rhos: per-band envelope across the frozen cover grid,
  available in the full export. Matching each envelope separately does not
  establish a common cover fraction satisfying both bands.

effective_bottom_blue/green_unclipped: diagnostic inversion of the observed
  reflectance with the assigned water and depth. Negative values or values
  above one are retained as evidence of incompatibility. They are not valid
  measured bottom reflectances.

a_blue/green_m_inv: assigned total absorption coefficient, per metre.
bb_blue/green_m_inv: assigned total backscattering coefficient, per metre.
kd490_m_inv: same-day Copernicus diffuse attenuation coefficient at 490 nm.
kd_blue/green_m_inv: modelled diffuse attenuation at the respective bands.
  Absorption and backscattering use the winter E5s spectral shape scaled to
  KD490. They have not been retrieved independently at each quadrat.
water_input_status: records this shared scene-level assumption.

bottom_attenuation_factor_blue/green: exp(-kb * depth), the model's attenuation
  factor on the subsurface bottom contribution. This is neither an observed
  fraction of bottom light nor a probability of detecting kelp. kb differs
  from the downwelling diffuse attenuation Kd.
black_bottom_floor_blue/green_rhos: predicted water-column signal with zero
  bottom reflectance under the assigned water and depth.
observed_below_floor_blue/green: observed reflectance is below that floor by
  more than the nominal tolerance. A nonnegative bottom cannot repair that
  mismatch with the assigned water and depth held fixed.

sensitivity_low/high_pct: union of compatible cover fractions over the frozen
  sensitivity scenarios; empty if none qualifies. This is not a confidence
  interval. Check sensitivity_compatible_scenarios alongside the range.
all_scenarios_best_fit_min/max_pct: range of closest fits, including rejected
  scenarios. A wide range describes instability, without probability weights.
depth/raw/raw_depth/mean_baseline_cover_pct: stored spatially held-out baseline
  predictions. Missing values indicate insufficient eligible training support.
  Compare models only where both predictions exist at the same footprints.

Example: July deep transect 3, quadrat 2
--------------------------------------

Diver cover was 10 percent. The Sentinel-2 image of 2 July gives a closest fit
of 100 percent; Pleiades Neo on 8 July gives zero percent. Both are incompatible.

Attribute                              Sentinel-2       Pleiades Neo
Water depth used (m)                     19.2940             19.4600
Observed blue rhos                        0.0344624           0.0065149
Observed green rhos                       0.0233912           0.0024827
Largest absolute residual                 0.0113054           0.0175298
Diagnostic bottom blue                    0.2807144          -0.4770199
Diagnostic bottom green                   0.7788689          -0.7540698

The sensors have different acquisition dates, bands and extraction footprints.
This example shows each pipeline output at an associated field location; it
does not isolate a sensor effect. All 162 sensitivity scenarios also fail
compatibility for this example on both sensors.

Verification and reproduction
-----------------------------

The export reproduced all 206 eligible predictions across V0/V1 and both
Pleiades Neo supports, checking cover, residual, compatibility and diagnostic
bottom reflectance against the saved results. It checked the field membership
and sampled cover mean of all 228 candidate units. Eight frozen code hashes
and every consumed frozen input hash were verified. Source run files remained
unchanged.

From the repository root, with its virtual environment:

  .venv/bin/python benchmark/reports/r7-quadrat-inspection-20261102/export_quadrats.py \
    --repo . --output /tmp/r7-quadrat-inspection-new

The output directory must not exist. The script emits the CSVs, GeoJSON and
verification manifest. This explanatory file accompanies the delivered export.
