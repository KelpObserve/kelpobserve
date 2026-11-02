# Sesimbra habitat benchmark: detectability envelope (r4)

Executed 2026-10-02, registered in
[`registration/revision-r4.json`](../registration/revision-r4.json) before any simulation
ran. Synthetic cover only: no model was fitted to observed cover.

**On the real scenes, no model recovered a synthetic cover signal in any of the four
sampling geometries, five patch sizes or five position errors tested.** That includes a
physics oracle given perfect positions and patches 25 m across. The pipeline works: on an
ideal background without noise, the depth-aware model and the oracle recover the signal in
80–100% of draws in every geometry. What stops recovery is how much the real reflectance
at the units varies for reasons other than depth. That variation is several times the
noise the recovery tolerates.

## Design

- **Geometries.**
  - Pléiades Neo, 8 July: native units (29) and 10 m cells (18).
  - Sentinel-2 R037: 2 July (14 cells, v1's mask) and 9 July (17 cells, r3's mask).
  - All use glint-corrected Sentinel-2 and the real valid-pixel masks.
- **Signal.**
  - July cover permuted within strata, painted as discs of radius 1.25–25 m.
  - Each transect displaced 0–20 m.
  - Turned into reflectance by the library's forward model on a 1 m grid, with E5s kelp
    against understory or bare substrate and Sesimbra QAA water.
  - Area-averaged to each sensor and blurred by a nominal point-spread function
    (σ = 0.5 pixel).
- **Backgrounds.** The real scene, or an ideal one: forward-model background plus
  independent pixel noise σ = 0–0.004.
- **Models.**
  - Depth-only, the baseline.
  - v1's additive blue/green/red + depth ridge.
  - A depth-aware ridge on Lee bottom reflectance (blue, green) + depth.
  - An oracle: unmixing with the true endmembers, no fitting.
- **Criterion.** "Recovered" means a model's MAE is at least 2 points below depth-only's
  in ≥ 80% of 100 draws.

## Pipeline check (registered): passed

On an ideal background with no noise, bare substrate, 10 m patches and no position error,
the depth-aware ridge and the oracle met the 2-point gain in 100% of draws on PNeo native
units. The envelope below can therefore be read.

## Ideal background: how much noise recovery tolerates

Largest per-pixel noise σ at which a model still recovers (10 m patches, no position
error). "None" means not recovered even without noise.

| Geometry | Contrast | v1 additive ridge | Depth-aware ridge | Oracle |
|---|---|---|---|---|
| PNeo native | bare | 0.004 | 0.002 | 0.004 |
| PNeo native | understory | 0.0005 | 0 only | 0.0005 |
| PNeo 10 m | bare | about 0.004 (0.79–0.87) | 0.002 | 0.004 |
| PNeo 10 m | understory | none (0.79 at 0) | 0.0005 | 0.0005 |
| S2 2 Jul, 14 cells | bare / understory | none (0.11 / 0.27) | 0 only | 0 only |
| S2 9 Jul, 17 cells | bare / understory | none (0.27 / 0.41) | 0 only | 0 only |

- **On Sentinel-2's 14 and 17 cells, the registered additive ridge does not recover the
  signal even without noise.** The depth-aware ridge does (80–98%). The models v1 used
  were a weak test there, as the second review argued.
- **Noise of 0.0005 per pixel already defeats every model on the Sentinel-2 cells.** That
  is about half the pixel-to-pixel scatter of the clean scenes offshore (about 0.001, r3). Near 17 m depth the bottom
  signal is attenuated about a hundredfold: to the inversion, 0.0005 in reflectance is
  about ±50% cover against bare substrate.
- **PNeo tolerates more** because each unit averages 9 pixels natively and about 70 per
  10 m cell.

## Real background: no recovery anywhere

Best share of draws with a 2-point gain over depth-only, any model:

| Geometry | Bare substrate | Understory | Null (no signal) |
|---|---:|---:|---:|
| PNeo native | 0.45 (R ≥ 5 m, no error) | 0.06 | 0.05 |
| PNeo 10 m | 0.44 (R 25 m, no error) | 0.05 | 0.05 |
| S2 2 Jul | 0.37 (R 25 m, no error) | 0.18 | 0.18 |
| S2 9 Jul | 0.41 (R 25 m, no error) | 0.12 | 0.12 |

- **Understory is indistinguishable from the null everywhere.**
- **Bare substrate comes closest, and only through v1's additive ridge** fed the shallow
  stratum's signal. On PNeo native with 10 m patches, scored on shallow units only
  (the model is fitted on both strata) it reaches 0.72, and on deep units 0.21.
- **Patch size and position matter, but cannot rescue it.** Larger patches help up to
  about 5 m (PNeo native: 0.06 at R 1.25 m, 0.45 from R 5 m). Position error removes most
  of it: 0.45, 0.33, 0.08, 0.02 and 0.03 at 0, 2, 5, 10 and 20 m.
- **The oracle fails in every condition,** with a median MAE 5–35 points worse than
  depth-only. Inverting the real reflectance amplifies whatever depth does not explain,
  about a hundredfold at the deep stratum. Because the oracle has perfect knowledge of the
  injected physics, its failure locates the limit in the background, not in support,
  position or model form.
- **The PSF and the IOP mismatch change little.** Box aggregation gives 0.45 against
  0.45 on PNeo native and 0.37 against 0.35 on 2 July. The more attenuating water lowers
  shares by 0–0.14.

**Why the real background defeats it** (descriptive, not registered). The within-stratum
spread of green reflectance across the real units is:

- 0.0048 for PNeo native units and 0.0023 for PNeo 10 m cells;
- 0.0006 on 2 July and 0.0014 on 9 July for Sentinel-2.

At unit level, the ideal ladder lost understory recovery at 0.0002–0.0003 (PNeo native)
and everything at 0.0005 (Sentinel-2). The real spread is at or above those levels. It
combines residual atmosphere, glint, adjacency and real seabed variation. Against a
permuted target, real seabed variation counts as noise too, so this comparison does not
show the variation is all error.

## What this establishes, and what it does not

- **It establishes, for these conditions,** that the unit-level reflectance variation in
  these scenes is larger than any of these models can see past to a cover signal of E5s
  size, at these depths and with these sample sizes. Kelp over understory is out of reach
  everywhere tested. Kelp over bare substrate comes closest, on PNeo, with large patches
  and positions within about 2 m.
- **It does not establish that the sensors cannot detect kelp at Sesimbra.**
  - One kelp spectral class (E5s, Australian and New Zealand kelps).
  - Two depth strata.
  - A disc model of beds.
  - Real seabed structure counted as noise.
  - Leave-one-transect-out on 14–29 units.
- **What would move a design into the envelope** (from these simulations, not tested in
  the field):
  - Lower unit-level residual, for example multi-date averaging or a correction that
    removes the common-mode term.
  - Larger units: PNeo 10 m cells tolerate more noise than native units.
  - More shallow units, where the signal is about 16 times larger than at 17 m.
  - Depth-aware models.
  - Positions known to about 2 m.

The obvious next step is the one the second review placed before any new correction:
diagnose what the unit-level residual is, by spatial scale, band relationship and
distance to shore. Done as r5 ([revision-r5.md](revision-r5.md)): on PNeo 8 July it is mostly
glint texture that NIR removes, which r4 did not test.

## Reproduce

```bash
PYTHONPATH=. python -m benchmark.revision_r4 run --draws 100 --jobs 8   # about 2 min
```

Code: `benchmark/revision_r4.py`, tests in `tests/test_revision_r4.py`. Output:
`out/r4/envelope/` (`draws.csv`, `envelope.json` with every condition, the envelope
tables, the noise ladder and the pipeline check).
