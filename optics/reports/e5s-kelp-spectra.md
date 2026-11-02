# E5s: kelp spectral contrast and detectability

Registration: `registration/e5s-kelp-spectra.json` (4285ff4). Inputs frozen in
`registration/e5s-kelp-spectra-lock.json` (30cdf9e, 648 files). Results:
[`reports/e5s-kelp-spectra/results.json`](e5s-kelp-spectra/results.json). Characterisation:
there is no pass or fail, and the four decision rules were fixed before any contrast was
computed.

**Run history.** The first run (after freeze 035b4f1) stopped before producing any
result. The runner dropped every spectrum with a missing edge band, which removed all
Ria de Vigo spectra, since they end at 700 nm. The fix (fd6739a) keeps a spectrum when
blue, green and red are present. The inputs were re-frozen (30cdf9e) before the run
reported here.

Kelp spectra: Rowan, Gillian, Phinn, Stuart, Smart, Joanna N. and Roelfsema, Chris (2025).
Hyperspectral Reflectance Measurements of Kelps, Macroalgae, and Co-ocurring Materials
Collected in Australia and New Zealand. The University of Queensland. (Dataset) doi:
10.48610/c5f745e. The spectra are used under the dataset's terms (re-use with
acknowledgement) and are not redistributed; only derived statistics appear here.

## Summary

1. **Kelp against bare substrate stays separable down to about 12–14 m** in blue and green
   at Sesimbra's winter water clarity.
   - Under water, single blue and green bands hold AUC ≥ 0.8 to 12 m, after the forward
     model, noise of 0.003 in reflectance and the library's correction.
   - The detection depth from the median contrast is 14 m at 492 and 560 nm, and 12 m at
     443 nm.
   - Across the two Sesimbra dates it ranges from 12 to 18 m at 492 nm. Across the noise
     levels it ranges from 11 to 17 m.
   - The library's "about 20 m at 444 nm" came from a lidar-fitted attenuation that was
     later voided (E2c). With forced-560 QAA attenuation it is 10–15 m at 443 nm.
2. **Kelp is not separable from the understory under water.**
   - The understory classes are red foliose algae, turf and coralline crust. No band or
     index reaches AUC 0.8 at 1 m or deeper in the pooled analysis.
   - The median contrast survives to 2–7 m, but the classes overlap too much for a pixel
     to be assigned.
   - The statement "brown cannot be separated from red" stands (D2). At Sesimbra, where
     kelp grows over turf and coralline crust at 7–17 m, per-pixel kelp mapping is
     outside what these bands can carry.
3. **Kelp is not separable from other brown algae** in the pooled data: AUC 0.61–0.72
   even in air.
4. **The flat-test index for E4t-B is ND(B2,B3)** (D1), (B2 − B3) / (B2 + B3). Kelp's
   green peak puts it lower than bare substrate. AUC is 0.92 in air and 0.73 at 5 m.
   - Under water, single bands (B2: 0.92 at 5 m) do better than any normalised
     difference, because a ratio inherits the faster attenuation of its longer band.
   - A flat test compares an exposed with a submerged scene, so it needs a scale-free
     index. ND(B2,B3) is the best available among the 10 m bands.
5. **Measurement protocol matters as much as class.** In the Hawaii library, the same
   genera measured in situ are 0.4 to 4.7 times as bright as ex situ: turf about 3×,
   *Lobophora* about 3×. Brightness differences between libraries are therefore not
   reliable. The within-library comparisons, and the AUCs that use the full spread, are
   the robust results.

## Data

| Class | Spectra | Taxa | Sources |
|---|---:|---:|---|
| kelp | 111 | 4 | Rowan & Phinn (*Ecklonia* 51, *Macrocystis* 45, *Undaria* 14); Vahtmäe (*Chorda* 1) |
| brown_other | 359 | 19 | Rowan & Phinn 201, Fuller 145, Vahtmäe 9, Glovers 4 |
| red_foliose | 150 | 9 | Fuller 138, Rowan & Phinn 11, Vahtmäe 1 |
| turf | 54 | 7 | Fuller 32, Glovers 9, Rowan & Phinn 8, Vahtmäe 5 |
| coralline_crust | 33 | 4 | Fuller 32, Glovers 1 |
| sand | 51 | 6 | Vigo 40, Rowan & Phinn 6, Vahtmäe 4, Glovers 1 |
| rock | 29 | 5 | Rowan & Phinn 26, Vahtmäe 3 |
| sediment (reported) | 18 | 1 | Vigo lab sediments |
| sediment_film (not scored) | 2 | 1 | Wallis Lake |

- **Bands:** Sentinel-2A (ESA v4.0) B1–B5 and Pléiades Neo 3 B1–B5, from ACOLITE 64a02ff.
- **Weighting:** every taxon has equal weight, so *Durvillaea* (119 spectra) does not set
  the brown median.
- **Coverage:** the Vigo spectra stop at 700 nm, so they have no B5 or Pléiades Neo coastal
  blue value.

## Band values (Sentinel-2, reflectance, taxon-balanced medians)

| Class | B1 443 | B2 493 | B3 560 | B4 665 | B5 704 |
|---|---:|---:|---:|---:|---:|
| kelp | 0.037 | 0.041 | 0.069 | 0.056 | 0.393 |
| brown_other | 0.021 | 0.025 | 0.043 | 0.036 | 0.151 |
| red_foliose | 0.014 | 0.020 | 0.032 | 0.035 | 0.122 |
| turf | 0.024 | 0.030 | 0.045 | 0.051 | 0.143 |
| coralline_crust | 0.014 | 0.021 | 0.044 | 0.037 | 0.136 |
| sand | 0.289 | 0.338 | 0.381 | 0.336 | 0.463 |
| rock | 0.205 | 0.244 | 0.303 | 0.321 | 0.335 |

- **Kelp's high B5** is the red edge of near-surface blades in air. It is attenuated within
  about 1 m under water (two-way k about 1.4 per m).
- **Kelp is brighter in the visible than the other algal classes.**
  - Within Rowan & Phinn alone, measured with one protocol, the same holds against other
    brown algae (about 2–3×, AUC about 0.8 in air).
  - *Undaria* is twice as bright as *Ecklonia* and *Macrocystis*.
  - Rowan & Phinn's filamentous red algae are brighter than kelp.
  - This comparison is descriptive, not registered.

## Separability (AUC, kelp against each class)

Noise 0.003 in reflectance; under water, Sesimbra's median forced-560 QAA IOPs and the
library correction.

| Against | Air, best band | Air, best index | 5 m, best | 12 m, best | 16 m, best |
|---|---|---|---|---|---|
| bare (rock + sand) | B4 0.93 | ND with B5 1.00 | B2 0.92 | B2 0.82 | B2 0.67 |
| sand | B2 0.96 | ND with B5 1.00 | B2 0.95 | B2 0.86 | B2 0.71 |
| rock | B4 0.90 | ND with B5 1.00 | B2 0.88 | B2 0.77 | B2 0.61 |
| red_foliose | B5 0.76 | ND(B3,B4) 0.85 | B3 0.71 | (0.75)* | 0.53 |
| turf | B5 0.68 | ND(B3,B4) 0.90 | ND(B2,B3) 0.59 | (0.69)* | 0.52 |
| coralline_crust | B1 0.75 | ND(B1,B4) 0.90 | B2 0.67 | (0.83)* | 0.51 |
| brown_other | B5 0.72 | ND(B3,B5) 0.66 | B1 0.60 | 0.64 | 0.51 |

\* **An artefact, not signal.** Some indices that use B4 score higher at 12 m than at 5 m.
At 12 m the corrected red band is amplified noise, because its two-way attenuation is
about 0.96 per m. A normalised difference with a noise-dominated denominator then takes a
distribution that depends on the other band. These values do not show separability, and
none of them reach 0.8 in the pooled analysis. The same artefact explains the air-only
result for red foliose algae with ND(B3,B4) at 12 m (see D2).

Pléiades Neo 3 in air gives the same ordering (`separability_pleiades_neo_air`).

## Detection depth (z_max, metres)

From the median contrast, `detectability.z_max_from_k`, two-way k of the forward model.
Sesimbra's median IOPs give k = 0.30, 0.26, 0.26, 0.96 and 1.38 per m at B1–B5.

| Kelp against | B1 | B2 | B3 | B4 |
|---|---:|---:|---:|---:|
| bare | 11.9 | 14.1 | 14.3 | 4.0 |
| sand | 12.5 | 15.0 | 15.2 | 4.1 |
| rock | 11.2 | 13.6 | 14.1 | 4.0 |
| red_foliose | 4.6 | 4.9 | 7.1 | 1.4 |
| coralline_crust | 4.6 | 4.7 | 5.6 | 1.3 |
| turf | 2.6 | 2.4 | 5.4 | 0.0 |
| brown_other | 3.4 | 3.8 | 5.8 | 1.3 |

For kelp against bare substrate at B2:
- **By noise level:** 16.7 m (0.0015), 14.1 m (0.003) and 11.4 m (0.006).
- **By date:** 17.6 m on 24 December 2025 and 11.7 m on 22 February 2026.

For the understory and other brown algae, z_max from the median contrast is an upper
bound: the AUCs show the classes overlap even where their medians differ.

## Decisions

- **D1, the E4t-B flat-test index: ND(B2,B3).**

  | Candidate (10 m bands) | Air | 5 m |
  |---|---:|---:|
  | ND(B2,B3) | 0.92 | 0.73 |
  | ND(B2,B4) | 0.66 | 0.80 |
  | ND(B3,B4) | 0.91 | 0.25 |

  Indices of B5 with B1, B2 or B3 score 0.78–0.82 at 5 m and 1.00 in air, but B5 is
  recorded at 20 m; they are reported, not chosen. E4t-B re-checks the choice with its own water before it
  is registered.
- **D2, brown against red: the statement stands.**
  - Pooled data: no band or index reaches AUC 0.8 at 1 m or deeper, against red foliose
    algae or turf.
  - Air-only data: some indices reach it for red foliose algae, ND(B3,B4) at 12 m and B3
    at 5 m. The 12 m value is the artefact above, and the rule needs both analyses.
- **D3, library defaults, no change in this run.** CC-BY-only medians compared with
  `SUBSTRATE_ALBEDO`:
  - brown: 0.97–1.06× in B1–B3, 1.21× in B4;
  - red: 0.83–1.28×;
  - bare: 1.26–1.53× brighter. The default is conservative, because Iberian seafloor sand
    from Vigo is brighter than the Baltic samples.
- **D4, detectability:** the table above. At Sesimbra, "brown against bare to about 20 m
  at 444 nm" becomes about 12 m at 443 nm and about 14 m at 492 and 560 nm, with
  forced-560 QAA attenuation.

## Measurement geometry (Fuller et al., descriptive)

Ratio of in-situ to ex-situ medians, same genus, Sentinel-2 bands:

| Genus | n ex / in | B1 | B2 | B3 | B4 |
|---|---|---:|---:|---:|---:|
| *Lobophora* (brown) | 3 / 27 | 2.67 | 2.89 | 2.79 | 2.29 |
| Turf | 6 / 23 | 2.88 | 3.02 | 3.44 | 3.15 |
| *Eucheuma* | 3 / 25 | 4.70 | 4.60 | 5.06 | 3.31 |
| *Acanthophora* | 2 / 25 | 3.81 | 3.89 | 4.61 | 2.77 |
| *Dictyosphaeria* | 9 / 49 | 2.01 | 1.78 | 1.89 | 1.27 |
| *Dichotomaria* | 3 / 19 | 0.96 | 1.10 | 1.30 | 1.02 |
| *Symploca* | 9 / 16 | 0.92 | 1.04 | 1.65 | 0.78 |
| *Halimeda* | 5 / 47 | 0.39 | 0.43 | 0.64 | 0.48 |

The data record does not say how the in-situ spectra were measured, and the specimens
differ between the two. Whatever the cause, protocol and specimen differences of 2–5× are
larger than most class contrasts among algae. That is why the pooled AUCs, which include
this spread, are the numbers to rely on.

## What this changes

- **E4t-B:** the flat test uses ND(B2,B3). With bare substrate at up to about 14 m of
  detectability, the fringe at a few metres is well within reach.
- **Plan and README:** the brown-against-bare depth is about 12–14 m at Sesimbra with QAA
  attenuation, not 20 m. The brown-against-red statement stands and now also covers turf
  and coralline crust.
- **KelpObserve at Sesimbra:** kelp over turf and coralline crust at 7–17 m is not
  separable per pixel from these spectra. What the imagery can carry is kelp against
  bare rock or sand, and a change in average reflectance where canopy cover changes a
  lot. The Sesimbra benchmark's negative result is consistent with this.
- **Library:** `sd-data-ingest/oceanstream/coastal/spectra.py` loads all six libraries.
  The Baltic defaults stay: brown and red agree within about 20 %, and bare is
  conservative.

## Limits

- **Kelp spectra** are three Australian and New Zealand genera, one *Chorda*, and no
  Iberian species. *Saccorhiza polyschides* belongs to Tilopteridales, not Laminariales.
  The argument rests on shared brown-algal pigments.
- **The water column is Sesimbra in winter**, two dates of forced-560 QAA, which is
  validated at Key West only.
- **The noise is a single level,** with no spatial structure. Real scenes add adjacency
  and date-to-date error (see E2e).
- **Single-material spectra stand in for mixed canopies.** Upright fronds under water are
  not represented.
- **The understory classes are mostly tropical** (Hawaii, Belize) standing in for
  Atlantic turf and crusts.
