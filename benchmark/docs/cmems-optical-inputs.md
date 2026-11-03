# Sesimbra optical input bundle

Generated: 2026-11-02T16:02:22+00:00

KD490, particulate backscatter (BBP), coloured dissolved and detrital absorption (CDM), and suspended matter (SPM) were downloaded with flags and provider uncertainties. Dataset versions and parts are pinned in the manifest; the source files are immutable local inputs. BBP and CDM refer to 443 nm.

Each row uses the same pixel as the repaired KD490 audit. Values in parentheses are provider percentage error bars. They have no assigned standard deviation or confidence level here.

| Date | KD490 (m⁻¹) | BBP443 (m⁻¹) | CDM443 (m⁻¹) | SPM (g/m³) | All four pass product checks |
|---|---:|---:|---:|---:|---|
| 2025-12-24 | 0.0879767 (35.71%) | 0.0111439 (38.57%) | 0.0512665 (58.62%) | 1.97551 (42.81%) | True |
| 2026-02-22 | 0.241208 (25.65%) | 0.00662852 (17.44%) | 0.178681 (41.72%) | 2.71533 (30.92%) | True |
| 2026-07-02 | 0.0710151 (15.29%) | 0.00408647 (19.79%) | 0.0367391 (36.60%) | 0.801044 (45.82%) | True |
| 2026-07-08 | 0.0784245 (36.48%) | 0.000601314 (19.52%) | 0.0795367 (56.94%) | 0.268153 (60.84%) | True |
| 2026-07-09 | 0.099137 (41.14%) | 0.00343695 (40.98%) | 0.0487983 (62.23%) | 0.361072 (49.76%) | True |
| 2026-08-31 | 0.150157 (10.57%) | 0.00367858 (46.01%) | 0.110645 (40.81%) | 0.647276 (39.92%) | True |
| 2026-09-03 | 0.220088 (58.27%) | 0.00326076 (40.00%) | 0.0297419 (60.00%) | 1.17295 (49.71%) | True |

Passing product checks does not establish that these 1 km pixels describe the reef water. BBP and CDM constrain parts of the optical model, while spectral slopes and phytoplankton absorption remain assumptions. Correlated product errors and mutual consistency must be tested before using these quantities together in inversion. No r7 habitat score was computed.

## Reproduction and source record

The [source manifest](../registration/evidence/cmems-optical-bundle-20261102.json) records the eight version-202603 subsets and their SHA-256 checksums. The [qualification record](../registration/evidence/cmems-optical-inputs-20261102.json) includes source-code hashes, flags, uncertainties and neighbourhood support. The local bundle is stored under `.cache/cmems-optics/20261102-v202603/`.

```bash
kelpobserve optics-inputs audit .cache/cmems-optics/20261102-v202603/manifest.json
```

All seven KD490 values reproduce the previously frozen source values exactly. No previously frozen files or habitat scores were replaced. Product-valid support within 5 km ranges from 48 to 55 pixels; this count does not establish reef representativeness.

The [Copernicus quality document](https://documentation.marine.copernicus.eu/QUID/CMEMS-OC-QUID-009-101to104-111-113-116-118.pdf) describes the algorithms and percentage error bars. Joint use of BBP, CDM and KD490 requires a consistency assessment, including their shared retrieval inputs and correlated errors.

OceanStream's local optical export repairs reproduce these 28 date-variable inputs through its native raster, point and NetCDF paths. Existing production rasters still require regeneration with quality fields and pinned source versions before they can supply equivalent benchmark inputs.
