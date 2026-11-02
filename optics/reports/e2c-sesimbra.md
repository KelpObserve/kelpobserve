# E2c at Sesimbra: forced 560 nm QAA against lidar k in the reef region

Registration: `registration/e2c-sesimbra.json` (78923f8). Library sd-data-ingest
db605f7 (attenuation fitted on fine depth only; Sesimbra AOI with the prototype's
pre-registered regions). Inputs frozen in `registration/e2c-sesimbra-lock.json`
(511e05f). Results: `reports/e2c-sesimbra/`.

## Outcome

**Not passed; the lidar comparison (C2) could not be assessed on either date.** QAA
stays experimental (E2 decision D2); the 560 nm reference stays an option.

| Date | Variant | Reference | bbp(555) m⁻¹ | Kd(492) m⁻¹ | C1 | C2 |
|---|---|---|---|---|---|---|
| 2025-12-24 | auto | 665 nm | 0.034 | 0.44 | fail | no trustworthy lidar k |
| 2025-12-24 | forced 560 | 559 nm | 0.010 | 0.13 | pass | no trustworthy lidar k |
| 2026-02-22 | auto | 665 nm | 0.028 | 0.42 | fail | no trustworthy lidar k |
| 2026-02-22 | forced 560 | 559 nm | 0.013 | 0.20 | pass | no trustworthy lidar k |

## Why there was no lidar k: the reef-region "attenuation" is a common-mode gradient

Fits in the dive-reef region on lidar depth (about 3,100–3,300 pixels):

| Band | 2025-12-24 k (r²) | 2026-02-22 k (r²) |
|---|---|---|
| 442 nm | 0.035 (0.46) | 0.145 (0.91) |
| 492 nm | −0.019 (0.63) | 0.110 (0.91) |
| 559 nm | −0.032 (0.92) | 0.071 (0.84) |
| 704 nm | 0.071 (0.90) | 0.178 (0.97) |
| 780 nm | 0.055 (0.95) | 0.165 (0.98) |
| 864 nm | 0.058 (0.89) | 0.165 (0.98) |

Pure water removes any bottom signal in the near-infrared within about a metre, so
the NIR bands cannot decay with depth across 1–20 m. On both dates they fit a depth
trend better than any visible band, with a slope as large as the visible ones. The
library's `depth_correlated_artefact` check therefore voids every band
(`scene_depth_gradient`), as designed: the regression is tracking a gradient that
correlates with depth and is shared by all bands. Over a reef at the foot of the
Arrábida cliffs the likely cause is adjacency light from land, strongest over the
shallowest water nearest the cliffs; the 40 m adjacency buffer does not remove it.

This also bears on the prototype's reef k of 0.108 m⁻¹ at 444 nm in the same
rectangle: on 2026-02-22 the visible slopes are no larger than the NIR ones, so a
visible lidar k at this reef cannot be taken as water attenuation without first
removing the common-mode gradient.

## What E2b and E2c together establish

- At two sites (Key West, Sesimbra) and four dates, automatic QAA switches to the
  red reference with bbp(555) above 0.02; forcing 560 nm gives plausible IOPs every
  time. That is a plausibility result, not an accuracy result.
- Neither Sesimbra region gives an independent attenuation to check QAA against:
  the whole scene mixes coarse depth and water masses (E2b), and the reef region is
  dominated by a depth-correlated common-mode gradient (E2c).

## Options for confirming the 560 nm reference

- A site without steep land next to the lidar: Ria Formosa (low barrier islands,
  lidar over the lagoon and inlets), with regions defined from the lidar before any
  fit.
- Subtracting a NIR-derived common-mode term from the visible bands before the fit,
  as a library change tested first at Key West. E2 fitted only visible bands there,
  so whether Key West passes the NIR artefact check is not yet known.
