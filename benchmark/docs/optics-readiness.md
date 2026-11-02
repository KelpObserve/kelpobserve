# Sesimbra habitat benchmark: water-column correction readiness

Executed 2026-09-20. **The water-column-corrected branch is unavailable on both Pléiades
Neo dates.** The QAA fit over the registered offshore reference region fails its
plausibility checks on 8 July and 3 September, and an independent cross-sensor comparison
shows why: ACOLITE's PNeo surface reflectance carries a large additive term over water.
The uncorrected models are unaffected. No cover value was read for either test.

## QAA fit over the registered reference region

The region was registered from depth and distance alone before its reflectance was ever
examined: 0.287 km², 51–110 m deep, at least 1 km from land and from the habitat AOI.
Coverage is not the problem — each fit used about 199,000 pixels, far above the 1,000
required.

| Quantity | 8 July | 3 September | Expected for clear coastal water |
|---|---:|---:|---|
| K_d at 486 nm (m⁻¹) | 1.32 | 1.26 | 0.05–0.2 |
| Absorption at 439 nm (m⁻¹) | 0.89 | 0.78 | about 0.01 |
| Particle backscatter at 555 nm (m⁻¹) | 0.126 | 0.120 | 0.001–0.01 |
| Reference wavelength chosen | 655 nm | 655 nm | a blue or green band |

Both fits raise the same flags: a red reference wavelength, implausible particle
backscatter, and absorption below the pure-water floor at 723 and 826 nm. Values like these
are the signature of a bright, spectrally flat signal over optically deep water, where the
water-leaving contribution should be near zero. QAA has no way to represent that except as
extreme absorption and backscatter.

Under the registered rule — at least 1,000 reference pixels and no plausibility flags — the
corrected branch is therefore **unavailable for both dates**.

## Cross-sensor comparison

Pléiades Neo aggregated to the common 10 m grid, against Sentinel-2, over water inside the
habitat AOI and over stable land three cells inland. Median surface reflectance:

| Pair | Band | PNeo over water | S2 over water | Difference | Correlation over water | Correlation over land |
|---|---|---:|---:|---:|---:|---:|
| 31 Aug S2B (clean, selected) vs 3 Sep PNeo | green | 0.0245 | 0.0184 | **+0.0057** | 0.08 | 0.88 |
| | red | 0.0177 | 0.0077 | **+0.0094** | 0.15 | 0.88 |
| | NIR | 0.0191 | 0.0059 | **+0.0120** | 0.30 | 0.91 |
| 3 Sep S2B (same minute, failed checks) vs 3 Sep PNeo | green | 0.0245 | 0.0302 | −0.0061 | 0.40 | 0.91 |
| | NIR | 0.0191 | 0.0195 | −0.0016 | 0.32 | 0.93 |

Two things follow. Over land the sensors track each other closely (r about 0.9), so
geolocation, aggregation and radiometric scaling are sound. Over water, against the clean
Sentinel-2 scene, PNeo sits 0.006 to 0.012 higher, and its NIR of 0.019 is far above the
0.006 that clear water gives. The same-minute pair looks closer only because that
Sentinel-2 scene is itself glinted.

The pairs are three days apart in the first case and glint-affected in the second, so
these differences are indicative rather than a calibration. They agree in magnitude and
sign with the QAA failure, which is the point.

## Why this happens

Pléiades Neo has no SWIR band. ACOLITE's dark-spectrum fitting is far less constrained
without it, and the correction cannot separate a spectrally flat surface term from the
water signal. The 8 July acquisition also sits in near-specular geometry. The residual is
about 0.01 to 0.04 in reflectance, which is the same size as, or larger than, the seabed
contrast the inversion is meant to measure at 7–17 m.

## Consequences

1. **The optical-correction comparison in the plan cannot run.** The question of whether
   water-column correction improves prediction is unanswerable with these deliveries, and is
   recorded as unavailable with this evidence.
2. **The uncorrected models still run**: training-set mean, depth-only, blue/green/red,
   blue/green/red with bathymetry, and the small forest, on July with September held back.
3. **Their spectral predictors carry this additive term.** Any apparent spectral skill has
   to be weighed against the depth-only baseline, which the plan already requires.
4. The residual is not a threshold artefact and cannot be tuned away by relaxing a check.
   Removing it needs either a SWIR-bearing sensor or a defensible surface correction fitted
   to these scenes, which is new method work, not a processing setting.
5. *Update 2026-10-02, exploratory ([revision-r1.md](revision-r1.md)).* ACOLITE's own
   NIR-referenced glint correction gives a plausible QAA fit on 3 September (bbp(555)
   0.006 m⁻¹, Kd(490) 0.24 m⁻¹) but over-corrects 8 July, so the corrected branch still
   has no July fit. The land comparison above also hides a gain gap: PNeo reads 18–27%
   below Sentinel-2 over land in green and red on both dates.

## Reproduce

```bash
PYTHONPATH=. python benchmark/run_benchmark.py optics-reference
PYTHONPATH=. python benchmark/run_benchmark.py cross-sensor
```

Outputs land under `out/optics-reference/` and
`out/cross-sensor/`, each with its own manifest and code provenance.
