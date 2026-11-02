# Sesimbra habitat benchmark: what the unexplained variation is (r5)

Executed 2026-10-02, registered in
[`registration/revision-r5.json`](../registration/revision-r5.json) before any measure was
computed. Label-free: no cover value is read. Attribution here is by association, not
proof of cause.

r4 found that reflectance at the model units varies, within each depth stratum, more than
cover recovery can tolerate. **That variation is different on each sensor and date:**

- **PNeo, 8 July (r4's PNeo basis): mostly sea-surface glint texture.** It is coherent
  with NIR, sits at scales under about 5 m, and is unrelated to distance from shore.
  Removing it with NIR cuts the unit spread about five-fold. What remains is still 1.5–2
  times the offshore level, largest at 14–20 m depth.
- **PNeo, 3 September: a smaller near-shore term.** It is coherent with NIR in the
  shallows and rises toward the shore.
- **Sentinel-2, glint-corrected: habitat cells vary about as much as offshore water.**
  About 40% of the depth- and NIR-removed residual in 5–20 m water repeats across a week,
  and across 60 days. It weakens with depth and is absent offshore: a stable,
  seabed-like pattern, alongside date-specific variation of similar size.
- **Sentinel-2's per-pixel SWIR correction adds little fine-scale noise.** The fine-scale
  glint left on glinted dates was already in the visible bands and is poorly matched by
  the 20 m SWIR reference. On one date, a one-pixel along-track offset matches it better.

## Regions and scales

- **Offshore:** the registered PNeo reference region, 51–110 m deep and at least 1 km from
  land, so there is no seabed signal.
- **Habitat bins:** 5–10 m (the shallow stratum), 10–14, 14–20 (the deep stratum) and
  20–30 m.
- **Unit scale,** matching r4: PNeo 3 × 3 means (native units) and 10 m cells, and
  Sentinel-2 10 m pixels.
- **Versions:** V0 is v1's processing; V1 adds r1's glint correction.

## Results

**Unit-scale spread of green** (robust SD; after the slash, with the NIR-explained part
removed):

| Scene, unit | Offshore | 5–10 m | 10–14 m | 14–20 m | 20–30 m |
|---|---|---|---|---|---|
| PNeo 8 Jul V0, native | 0.0028 / 0.0010 | 0.0095 / 0.0016 | 0.0143 / 0.0020 | 0.0156 / 0.0020 | 0.0109 / 0.0021 |
| PNeo 8 Jul V1, native | 0.0013 | 0.0019 | 0.0024 | 0.0029 | 0.0026 |
| PNeo 8 Jul V0, 10 m | 0.0104 / 0.0008 | 0.0106 / 0.0011 | 0.0124 / 0.0012 | 0.0142 / 0.0011 | 0.0085 / 0.0009 |
| PNeo 3 Sep V0, native | 0.0016 / 0.0010 | 0.0046 / 0.0027 | 0.0028 / 0.0021 | 0.0018 / 0.0015 | 0.0017 / 0.0015 |
| S2 2 Jul V1 | 0.0018 / 0.0014 | 0.0018 / 0.0012 | 0.0012 / 0.0012 | 0.0011 / 0.0011 | 0.0011 / 0.0010 |
| S2 9 Jul V1 | 0.0020 / 0.0015 | 0.0024 / 0.0015 | 0.0018 / 0.0016 | 0.0016 / 0.0017 | 0.0017 / 0.0015 |
| S2 31 Aug V1 | 0.0017 / 0.0016 | 0.0020 / 0.0012 | 0.0014 / 0.0013 | 0.0012 / 0.0011 | 0.0010 / 0.0009 |

**PNeo, 8 July.**

- **The variation is NIR-coherent.** NIR explains 86–94% of the native-unit variance of
  green in the habitat bins, and 96–99% at 10 m.
- **It is fine-scale.** In the 14–20 m bin, green's semivariance reaches most of its
  plateau by 5 m (spread 0.0053 at 1.2 m, 0.0080 at 4.8 m), with NIR identical.
- **It does not depend on distance to shore:** an R² gain of 0.002 for NIR.
- **Reading:** wave-facet glint under the near-specular July geometry. It is three to five
  times stronger over the habitat than offshore; why is not resolved here (shelter, slicks
  or local geometry).
- **After NIR correction (V1),** native-unit spread falls about five-fold. The 14–20 m bin
  is still 2.2 times the offshore level, at every lag from 1.2 m. At 17 m the seabed
  cannot produce that much (r4's signal for 100% kelp against bare substrate there is
  0.0015), so imperfectly scaled glint is the more likely remainder.

**PNeo, 3 September.**

- **Near-shore only.** The habitat varies more than offshore only in the shallow bins.
  There it is NIR-coherent (R² 0.58 natively, 0.83 at 10 m) and NIR rises toward the
  shore by 0.0065 per tenfold decrease in distance.
- **Reading:** a near-shore common mode, adjacency or surf. It narrowly misses the
  registered rule D2, with an R² gain of 0.03–0.045 against the 0.05 required.

**Sentinel-2, glint-corrected.**

- **Habitat cells vary about as much as offshore water** (D1 met on the clean dates).
- **Cross-date correlation of the depth- and NIR-removed residuals** (green, 10 m):

| Pair | Offshore | 5–10 m | 10–14 m | 14–20 m | 20–30 m |
|---|---:|---:|---:|---:|---:|
| 2 Jul–9 Jul (7 d) | −0.04 | 0.38 | 0.42 | 0.38 | 0.24 |
| 28 May–2 Jun (5 d) | −0.11 | 0.49 | 0.44 | 0.32 | 0.11 |
| 2 Jun–9 Jun (7 d) | 0.29 | 0.28 | 0.18 | 0.25 | 0.04 |
| 30 May–7 Jun (8 d) | 0.05 | 0.08 | 0.31 | 0.40 | 0.14 |
| 2 Jul–31 Aug (60 d) | 0.45 | 0.35 | 0.52 | 0.47 | 0.21 |

- **A stable spatial pattern exists in 5–20 m water,** fades by 20–30 m, and is mostly
  absent offshore over a week. That is how seabed structure should behave.
  - Its size: with r ≈ 0.4 on residuals of SD 0.0011–0.0015, the stable part is about
    0.0007 in green per 10 m pixel. The date-specific part is of similar size.
  - Persistent near-shore water features would behave the same way, so this is
    seabed-like, not shown to be seabed.
  - The 60-day pair also correlates offshore (0.45), so some large-scale pattern
    persists in open water as well.
- **Visible reflectance rises away from the shore at fixed depth** (R² gain 0.14–0.30),
  with no matching NIR trend: part of the same stable structure.
- **The registered seabed rule D3 is not met anywhere.** Its size criterion compares the
  habitat bins with a 5.7 km offshore region that carries larger-scale gradients
  (registration erratum).
- **In uncorrected Sentinel-2 (V0) the near-shore common mode is clear** (D2 met on 28 May
  and 2 July). The SWIR correction removes it.

**PNeo across dates** (8 July–3 September, 10 m, 57 days): residual correlations are
0.02–0.19 (V0) and −0.07–0.27 (V1). The two PNeo dates share little stable structure at
10 m, unlike Sentinel-2 over the same span. Their residuals are dominated by
date-specific terms: glint in July, the near-shore term in September.

**How well the glint correction matches (M6)**, offshore, glint reference against
uncorrected green:

| Scene | Fine-scale r | Coarse r | Best shift | Correction's share of fine variance |
|---|---:|---:|---|---:|
| S2 28 May (glinted) | 0.63 | 0.99 | none | 0.26 |
| S2 2 Jun (glinted) | 0.31 | 0.91 | one 20 m pixel along track, r 0.51 | 0.13 |
| S2 9 Jul (glinted) | 0.28 | 0.94 | none | 0.08 |
| S2 2 Jul / 31 Aug (clean) | 0.48 / 0.44 | 0.96 / 0.93 | none / small | 0.24 / 0.26 |
| PNeo 8 Jul / 3 Sep | 0.82 / 0.78 | 0.99 / 0.98 | none | 2.45 / 1.59 |

- **SWIR captures glint well at coarse scale** (r 0.91–0.99) and poorly at its own 20 m
  pixel scale (0.28–0.63).
- **On 9 July, D4 is met as registered** (fine r < 0.3, coarse r > 0.7), nominally
  warranting a smoothed glint estimate. But the correction term carries only 8% of the
  corrected image's fine-scale variance (V0 0.00109, V1 0.00110). The fine-scale glint
  was already in the visible bands and the SWIR band does not see it. Smoothing would
  therefore leave it in place (registration erratum).
- **On 2 June a one-pixel along-track shift improves the match by 0.2.** Sentinel-2's
  bands view the surface at slightly different times along track, and the sea surface
  moves in between.
- **PNeo's NIR follows the visible fine structure closely** (0.78–0.82), so an
  NIR-referenced correction is well matched at PNeo's scales.

## What this means for r4 and next steps

1. **r4's PNeo envelope was run on uncorrected 8 July reflectance,** whose unit spread is
   mostly glint texture that NIR removes. The depth-aware model and the oracle were never
   tested on the corrected PNeo reflectance (V1), whose unit spread is about five times
   smaller. Repeating r4 on V1 is the most direct next test. Done as r6 ([revision-r6.md](revision-r6.md)): it did not reach the
   envelope either, because the corrected water is biased dark.
2. **On Sentinel-2, about half the unexplained variation is date-specific** and could
   shrink with multi-date averaging. The other half is stable and seabed-like, so it is
   either habitat signal or confounding structure. Nothing here says which; only a
   reference at pixel scale could.
3. **For glinted Sentinel-2 dates the next correction variant is band matching,** not
   smoothing. Candidates are co-registering the SWIR glint estimate along track, or
   taking glint from the 10 m NIR band at the depths where its seabed signal is
   negligible. That is a hypothesis for a new registration.
4. **For PNeo, July's near-specular geometry is the main problem,** not the sensor. The
   September near-shore term is smaller but real.

## Reproduce

```bash
PYTHONPATH=. python -m benchmark.revision_r5 run   # about 1 min, 3 GB memory
```

Code: `benchmark/revision_r5.py`, tests in `tests/test_revision_r5.py`. Output:
`out/r5/diagnosis/` (`diagnosis.json`: every measure per scene and unit scale, the
cross-date pairs and the four decisions).
