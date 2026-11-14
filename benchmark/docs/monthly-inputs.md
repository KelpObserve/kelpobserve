# Sesimbra field inputs and PNEO metadata

Source facts from the September-delivery inspection of 2026-09-19, retained because
they affect target aggregation, spatial alignment and interpretation. The current
evaluation uses the frozen [monthly input bundle](../inputs/monthly/) and the
[separate monthly protocol](monthly/README.md).

## Field observations

The source directory is
`ccmar-diving/sept/` in the data folder.

- `TNA_Aquaserv_KelpData2026.xlsx`, sheet `Folha1`, contains 260 observation
  rows across June, July and September. Multiple rows describe individual kelp
  within a quadrat; they are not independent labelled map samples.
- September occupies rows 213–261: 49 records representing 30 quadrats,
  keyed by date, depth stratum, transect and quadrat number. Twenty quadrats
  contain a kelp species record and ten are explicitly `EMPTY`.
- Nineteen September quadrats have a numeric percentage-cover entry. The ten
  `EMPTY` records have blank percentage cover; one kelp-positive quadrat also
  lacks percentage cover (`K222`, 15 m, transect 2, quadrat 3). Preserve explicit
  absence separately from unknown cover. A generic blank-to-zero conversion
  would turn that positive quadrat into a false absence.
- `Transects kelp AQUASERV Jun-Sep 2026.gpkg`, layer `transects`, contains
  90 point features, 30 per campaign. September has three transects with five
  quadrats at each of two labelled depth strata.
- `Maps kelp AQUASERV Jun-Sep 2026.gpkg`, layer `maps`, contains three
  MULTILINESTRING features labelled June, July and September. These are line
  geometries; they do not by themselves provide labelled habitat polygons.
- The sonde workbook includes 211 September profile rows: 146 under
  `Arrábida_Sep26_17m` and 65 under `Arrábida_Sep26_6m`. This inventory does not
  establish profile quality or optical-model validation.

## Metadata discrepancies requiring resolution

| Source | September date | Depth labels |
|---|---|---|
| Kelp workbook | 2026-09-02 | 6 m and 15 m |
| Transect GeoPackage | `01/set` | 7 m and 17 m |
| September mapping line | Name includes `01/09/26, 15:33:37` | Not established |
| Sonde workbook | Underlying Excel value 46062 = 2026-02-09 | Site names say September, 6 m and 17 m |

The sonde discrepancy is present in the original XLSX XML, not introduced by
extraction. It is consistent with a day/month interpretation problem, but its
correction is not confirmed. June sonde date values also need review. Retain raw
dates, depth labels and source identifiers until an explicit reconciliation is
recorded. Differences in depth labels are not acquisition-time tide measurements.

## Pléiades Neo delivery

Archives are under `pneo/` in the data folder.

| Property | July | September |
|---|---|---|
| Archive prefix | `WO_000598128_1_2` | `WO_000618890_1_2` |
| Acquisition date/time from DIM metadata | 2026-07-08 11:38:15.6 | 2026-09-03 11:33:59.9 |
| Processing level | ORTHO | ORTHO |
| Multispectral product | MS-FS, six bands | MS-FS, six bands |
| Radiometric processing | REFLECTANCE | BASIC |
| Delivered multispectral grid | EPSG:32629, 1.2 m, 4920 × 2817 | Same |

Both contain red, green, blue, NIR, red-edge and deep-blue bands, plus a separate
panchromatic product. Identical grid metadata does not establish subpixel
coregistration, habitat classification resolution or bathymetric visibility.

All 90 supplied quadrat coordinates fall inside both image footprints. All 30
September coordinates have nonzero, unmasked RGB samples in both deliveries.
This establishes raster availability at the points, not usable seabed signal:
cloud, glint, water clarity, positional errors and habitat mixing remain unchecked.

September is one day after the workbook survey date, or two days after the
GeoPackage date. July is one to two days after the corresponding workbook survey
dates. These are promising temporal pairings, subject to confirming field dates.

The two delivery types require their appropriate radiometric conversions before
comparison. The local ACOLITE implementation has separate BASIC and REFLECTANCE
conversion branches in `acolite/pleiades/l1_convert.py`. Reprocess under a pinned,
consistent configuration and verify units and image alignment before treating
July–September differences as ecological change. Supplier REFLECTANCE metadata
alone does not establish water-surface atmospheric correction.

## Use in the current evaluation

The monthly input bundle contains one event per quadrat, preserves positive
missing cover and derives zero only from explicit absence. It retains the source
metadata and quality flags for a September sensitivity analysis. Cover is a
supplied quadrat-level numerical estimate; these files do not establish a formal
observer-error or photographic scoring protocol.

July and September are fitted separately. This avoids pooling the two radiometric
products, but does not correct their water, illumination or spatial-alignment
uncertainties. Repeated transect/quadrat identifiers do not establish that the
same exact seabed footprints were measured in both months.

## Site context

The retained maps provide geographic context, not training labels or model accuracy:
[location](figures/map-1-location.jpg), [bathymetry](figures/map-2-bathymetry.png),
[survey transects](figures/map-3-transects.jpg), and
[habitat context](figures/map-5-habitat.jpg). Background credits are in
[THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md).
