# Sesimbra Copernicus KD490 matchup audit

Generated: 2026-11-02T15:05:41+00:00

The previous selector accepted finite values on LAND-flagged pixels. The repaired selector requires the reviewed LAND=1 flag schema, flags=0, positive finite KD490 within the encoded valid range, and finite nonnegative uncensored uncertainty in percent. An unknown schema stops the audit. Rejection counts, including overlapping reasons, are in the JSON evidence.

This is a development input audit following inspection of the earlier matchups. The existing 5 km / ±3 day search limits are retained. Selection prefers the nearest day, then the nearest eligible pixel; ties between days favour the earlier day. Uncertainty caps of 20%, 30% and 50% are sensitivity examples evaluated on the fixed candidate. None is adopted as an application acceptance threshold.

## Candidate values

| Date | Role | Previous KD490 / flag | Repaired KD490 | Uncertainty | Distance / day offset | Status |
|---|---|---|---:|---:|---|---|
| 2025-12-24 | winter Sentinel-2 QAA cross-check | 0.0902 / 1 | 0.0880 | 35.71% | 0.891 km / +0 d | candidate_only |
| 2026-02-22 | winter Sentinel-2 QAA cross-check | 0.2714 / 1 | 0.2412 | 25.65% | 0.891 km / +0 d | candidate_only |
| 2026-07-02 | July Sentinel-2 candidate | 0.0708 / 1 | 0.0710 | 15.29% | 0.891 km / +0 d | candidate_only |
| 2026-07-08 | July Pleiades Neo scene | 0.0861 / 1 | 0.0784 | 36.48% | 0.891 km / +0 d | candidate_only |
| 2026-07-09 | July Sentinel-2 candidate | 0.0810 / 1 | 0.0991 | 41.14% | 0.891 km / +0 d | candidate_only |
| 2026-08-31 | late-August Sentinel-2 candidate | 0.1560 / 1 | 0.1502 | 10.57% | 0.891 km / +0 d | candidate_only |
| 2026-09-03 | September Pleiades Neo scene | 0.1527 / 1 | 0.2201 | 58.27% | 0.891 km / +0 d | candidate_only |

KD490 is in m⁻¹. Flag 1 means LAND. Previous values are retained here for error tracing; historical run files and the frozen subset are unchanged.

## Spatial and temporal support

Spatial summaries use product-valid water pixels within 1, 3 and 5 km of the reef and disjoint distance rings. The JSON also reports the eight adjacent grid positions, excluding the selected pixel. These are correlated coastal observations. Their spread has no assigned confidence level.

| Date | Water pixels ≤5 km | KD490 within 1 km (count; median) | KD490 1–3 km / 3–5 km medians | Fixed-pixel valid days / 7 | KD490 ±3-day range |
|---|---:|---|---|---:|---|
| 2025-12-24 | 55 | 1; 0.0880 | 0.0949 / 0.1081 | 6 | 0.0752–0.1040 |
| 2026-02-22 | 55 | 1; 0.2412 | 0.2384 / 0.2192 | 5 | 0.1742–0.3150 |
| 2026-07-02 | 55 | 1; 0.0710 | 0.0747 / 0.0760 | 6 | 0.0710–0.1404 |
| 2026-07-08 | 55 | 1; 0.0784 | 0.0807 / 0.0749 | 6 | 0.0751–0.1157 |
| 2026-07-09 | 53 | 1; 0.0991 | 0.0977 / 0.0921 | 6 | 0.0696–0.1157 |
| 2026-08-31 | 55 | 1; 0.1502 | 0.1590 / 0.1600 | 7 | 0.0931–0.2201 |
| 2026-09-03 | 55 | 1; 0.2201 | 0.2109 / 0.1992 | 6 | 0.0836–0.2201 |

The temporal check keeps the selected coordinate fixed and includes missing or rejected days in the JSON. Same-day daily bins do not establish coincidence with the satellite overpass. Spatial or temporal agreement cannot rule out shared atmospheric, adjacency or shallow-bottom bias.

## Use in the benchmark

7 of 7 dates have product-valid candidates for input sensitivity work. Reef representativeness remains unestablished for every date. The choice of fallback, uncertainty treatment, spectral extension and r7 outcome criteria must be fixed before corrected habitat scoring. KD490 alone cannot identify absorption and backscatter separately. The winter-QAA spectral shape remains an assumption.

The [Copernicus manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-OC-PUM.pdf), II.4.3, describes the uncertainty as a percentage error bar and caps it at 327.67%. This audit treats capped values as unquantified and rejects them. It does not interpret other error bars as standard deviations or 95% intervals.

## Reproduction

```bash
kelpobserve kd490 qualify
```

The command verifies the frozen subset and writes a new JSON evidence file and this report under `out/kd490/`. Each run records the subset hashes, audit code hashes and policy. It performs no network fetch or habitat scoring.

Frozen subset SHA-256: `cb3dfbcb1d97e12f3305a3ab9fb4ef2e0d74a4724c94ca9a42e557824b810462`.
