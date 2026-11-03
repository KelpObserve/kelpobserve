"""Input-only audit of frozen Sesimbra KD490 candidates; no habitat scoring.

Selection stays nearest-day/nearest-pixel after product validity checks. Uncertainty
cutoffs are reported as sensitivity checks on that fixed candidate; they never cause
selection of a more favourable pixel. The reef fallback decision remains open.
"""
from __future__ import annotations

import json
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

import numpy as np

from kelpobserve.attenuation import cmems

POLICY = {
    "version": "cmems-matchups-v1",
    "status": "development input audit; existing matchup values were already inspected",
    "site_lat_lon": list(cmems.SITE),
    "max_km": 5.0,
    "max_days": 3,
    "selection": "nearest product-valid pixel on nearest available day; earlier day wins a date tie",
    "product_checks": "reviewed LAND=1 schema; flags=0; finite positive in-range KD490; "
                      "finite nonnegative in-range uncensored uncertainty in percent",
    "spatial_radii_km": [1.0, 3.0, 5.0],
    "temporal_window_days": 3,
    "uncertainty_sensitivity_caps_pct": [20, 30, 50],
    "uncertainty_cap_status": "diagnostic examples; no cap chosen as a physical acceptance rule",
    "spatial_summary": "eligible water pixels; radii centred on reef; rings disjoint; "
                       "3x3 neighbours exclude selected pixel; no independent replicate assumption",
    "temporal_summary": "fixed selected coordinate within +/-3 days of target; daily bins, "
                        "not exact acquisition-time matchups; no independent replicate assumption",
    "uncertainty_interpretation": "provider error bar in percent; not assigned a standard deviation or CI",
    "application_status": "not qualified as reef IOPs; fallback and r7 scoring remain closed",
}
TARGETS = {
    "2025-12-24": "winter Sentinel-2 QAA cross-check",
    "2026-02-22": "winter Sentinel-2 QAA cross-check",
    "2026-07-02": "July Sentinel-2 candidate",
    "2026-07-08": "July Pleiades Neo scene",
    "2026-07-09": "July Sentinel-2 candidate",
    "2026-08-31": "late-August Sentinel-2 candidate",
    "2026-09-03": "September Pleiades Neo scene",
}
MANUAL = "https://documentation.marine.copernicus.eu/PUM/CMEMS-OC-PUM.pdf"


def stats(values) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return {"count": 0, "min": None, "q25": None, "median": None, "q75": None, "max": None}
    return {"count": int(values.size), **dict(zip(
        ["min", "q25", "median", "q75", "max"],
        map(float, np.quantile(values, [0, .25, .5, .75, 1])), strict=True,
    ))}


def _number(value):
    return float(value) if np.isfinite(value) else None


def _grid(dataset, name, index):
    return dataset[name].isel(time=index).transpose("latitude", "longitude").values


def audit_match(dataset, when: date, max_km: float = 5, max_days: int = 3, site=cmems.SITE) -> dict:
    distance = cmems.search_geometry(dataset, max_km, max_days, site)
    match = cmems.nearest_valid(dataset, when, max_km, max_days, site)
    days = dataset.time.values.astype("datetime64[D]")
    radius = distance <= max_km
    window = []
    previous = None
    for offset in sorted(range(-max_days, max_days + 1), key=lambda value: (abs(value), value)):
        day = np.datetime64(when + timedelta(days=offset), "D")
        positions = np.flatnonzero(days == day)
        if not positions.size:
            window.append({"date": str(day), "day_offset": offset, "status": "missing_day"})
            continue
        index = int(positions[0])
        usable, reasons = cmems.quality_masks(dataset, index)
        kd = _grid(dataset, "KD490", index)
        window.append({
            "date": str(day), "day_offset": offset, "status": "available",
            "grid_pixels_within_radius": int(radius.sum()),
            "eligible_within_radius": int((usable & radius).sum()),
            "rejected_unique_within_radius": int((~usable & radius).sum()),
            "rejection_reasons_overlapping": {name: int((mask & radius).sum()) for name, mask in reasons.items()},
        })
        old_usable = np.isfinite(kd) & radius
        if previous is None and old_usable.any():
            row, column = np.unravel_index(np.where(old_usable, distance, np.inf).argmin(), distance.shape)
            previous = {
                "date": str(day), "day_offset": offset, "distance_km": float(distance[row, column]),
                "latitude": float(dataset.latitude.values[row]), "longitude": float(dataset.longitude.values[column]),
                "kd490": float(kd[row, column]),
                "uncertainty_pct": _number(_grid(dataset, "KD490_uncertainty", index)[row, column]),
                "flags": _number(_grid(dataset, "flags", index)[row, column]),
                "rejection_reasons": [name for name, mask in reasons.items() if mask[row, column]],
            }
    result = {
        "target_date": when.isoformat(), "candidate": match.to_json() if match else None,
        "previous_finite_only_candidate": previous, "search_days": window,
        "qualification": {
            "status": "candidate_only" if match else "unavailable",
            "product_checks": "pass" if match else "no_product_valid_candidate_in_search_window",
            "reef_representativeness": "unestablished",
            "r7_fallback_selected": False,
            "reasons": ["coastal 1 km footprint not independently matched to reef water",
                        "same-day daily bin does not establish an acquisition-time match",
                        "KD490 does not identify absorption and backscatter separately"],
        },
    }
    if match is None:
        result["qualification"]["reasons"].insert(0, "no eligible pixel within fixed space/time limits")
        return result
    index = int(np.flatnonzero(days == np.datetime64(match.date))[0])
    row = int(np.flatnonzero(dataset.latitude.values == match.latitude)[0])
    column = int(np.flatnonzero(dataset.longitude.values == match.longitude)[0])
    usable, _ = cmems.quality_masks(dataset, index)
    kd = _grid(dataset, "KD490", index)
    uncertainty = _grid(dataset, "KD490_uncertainty", index)
    spatial = {}
    inner = 0.0
    for outer in POLICY["spatial_radii_km"]:
        for label, mask in ((f"within_{outer:g}_km", distance <= outer),
                            (f"ring_{inner:g}_{outer:g}_km", (distance > inner) & (distance <= outer))):
            selected = usable & mask
            spatial[label] = {"kd490": stats(kd[selected]), "uncertainty_pct": stats(uncertainty[selected])}
        inner = outer
    neighbours = np.zeros(usable.shape, dtype=bool)
    neighbours[max(row-1, 0):row+2, max(column-1, 0):column+2] = True
    neighbours[row, column] = False
    spatial["adjacent_3x3_excluding_selected"] = {
        "eligible_count": int((usable & neighbours).sum()),
        "kd490": stats(kd[usable & neighbours]),
        "uncertainty_pct": stats(uncertainty[usable & neighbours]),
    }
    result["spatial_support"] = spatial
    result["uncertainty_sensitivity"] = [
        {"cap_pct": cap, "fixed_candidate_meets_cap": bool(match.uncertainty_pct <= cap),
         "eligible_pixels_within_radius_meeting_cap": int((usable & radius & (uncertainty <= cap)).sum())}
        for cap in POLICY["uncertainty_sensitivity_caps_pct"]
    ]
    temporal = []
    for offset in range(-POLICY["temporal_window_days"], POLICY["temporal_window_days"] + 1):
        day = np.datetime64(when + timedelta(days=offset), "D")
        positions = np.flatnonzero(days == day)
        item = {"date": str(day), "day_offset": offset}
        if not positions.size:
            item.update(status="missing_day", kd490=None, uncertainty_pct=None)
        else:
            t = int(positions[0])
            valid, reasons = cmems.quality_masks(dataset, t)
            item.update(
                status="eligible" if valid[row, column] else "rejected",
                kd490=_number(_grid(dataset, "KD490", t)[row, column]),
                uncertainty_pct=_number(_grid(dataset, "KD490_uncertainty", t)[row, column]),
                flags=_number(_grid(dataset, "flags", t)[row, column]),
                rejection_reasons=[name for name, mask in reasons.items() if mask[row, column]],
            )
        temporal.append(item)
    result["temporal_support"] = {
        "latitude": match.latitude, "longitude": match.longitude, "days": temporal,
        "kd490": stats([item["kd490"] for item in temporal if item["status"] == "eligible"]),
        "uncertainty_pct": stats([item["uncertainty_pct"] for item in temporal if item["status"] == "eligible"]),
    }
    if match.day_offset:
        result["qualification"]["reasons"].append(f"candidate is {match.day_offset:+d} days from target")
    return result


def build(evidence_path: Path | None = None, subsets: Path | None = None) -> dict:
    from kelpobserve.attenuation import table

    evidence_path = evidence_path or table.EVIDENCE
    subsets = subsets or table.SUBSETS
    evidence = json.loads(evidence_path.read_text())
    dataset = cmems.open_subset(subsets / evidence["file"], evidence)
    return {
        "schema_version": "1.0", "created_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "policy": POLICY, "subset": evidence,
        "source_sha256": {f"src/kelpobserve/attenuation/{p.name}": cmems.sha256(p)
                          for p in (Path(cmems.__file__), Path(__file__))},
        "references": {"manual": MANUAL, "manual_section": "II.4.3, issue 7.0, March 2026"},
        "targets": [{"role": role, **audit_match(dataset, date.fromisoformat(day))}
                    for day, role in TARGETS.items()],
        "scope": "Frozen subset input audit. No habitat predictions, corrected-branch metrics or label selection.",
    }


def markdown(record: dict) -> str:
    lines = [
        "# Sesimbra Copernicus KD490 matchup audit", "", f"Generated: {record['created_at']}", "",
        ("The previous selector accepted finite values on LAND-flagged pixels. The repaired selector requires "
        "the reviewed LAND=1 flag schema, flags=0, positive finite KD490 within the encoded valid range, and "
        "finite nonnegative uncensored uncertainty in percent. An unknown schema stops the audit. "
        "Rejection counts, including overlapping reasons, are in the JSON evidence."), "",
        ("This is a development input audit following inspection of the earlier matchups. The existing "
        "5 km / ±3 day search limits are retained. Selection prefers the nearest day, then the nearest eligible "
        "pixel; ties between days favour the earlier day. Uncertainty caps of 20%, 30% and 50% are sensitivity "
        "examples evaluated on the fixed candidate. None is adopted as an application acceptance threshold."), "",
        "## Candidate values", "",
        "| Date | Role | Previous KD490 / flag | Repaired KD490 | Uncertainty | Distance / day offset | Status |",
        "|---|---|---|---:|---:|---|---|",
    ]
    for target in record["targets"]:
        match, old = target["candidate"], target["previous_finite_only_candidate"]
        prior = f"{old['kd490']:.4f} / {old['flags']:g}" if old else "none"
        values = (f"{match['kd490']:.4f} | {match['uncertainty_pct']:.2f}% | "
                  f"{match['distance_km']:.3f} km / {match['day_offset']:+d} d") if match else "— | — | —"
        lines.append(f"| {target['target_date']} | {target['role']} | {prior} | {values} | "
                     f"{target['qualification']['status']} |")
    lines += [
        "", ("KD490 is in m⁻¹. Flag 1 means LAND. Previous values are retained here for error tracing; "
        "historical run files and the frozen subset are unchanged."), "",
        "## Spatial and temporal support", "",
        ("Spatial summaries use product-valid water pixels within 1, 3 and 5 km of the reef and disjoint "
        "distance rings. The JSON also reports the eight adjacent grid positions, excluding the selected "
        "pixel. These are correlated coastal observations. Their spread has no assigned confidence level."), "",
        ("| Date | Water pixels ≤5 km | KD490 within 1 km (count; median) | "
        "KD490 1–3 km / 3–5 km medians | Fixed-pixel valid days / 7 | KD490 ±3-day range |"),
        "|---|---:|---|---|---:|---|",
    ]
    def num(value):
        return f"{value:.4f}" if value is not None else "—"
    for target in record["targets"]:
        if not target["candidate"]:
            continue
        spatial, temporal = target["spatial_support"], target["temporal_support"]["kd490"]
        near = spatial["within_1_km"]["kd490"]
        rings = [spatial[key]["kd490"]["median"] for key in ["ring_1_3_km", "ring_3_5_km"]]
        lines.append(f"| {target['target_date']} | {target['candidate']['valid_within_radius']} | "
                     f"{near['count']}; {num(near['median'])} | {num(rings[0])} / {num(rings[1])} | "
                     f"{temporal['count']} | {num(temporal['min'])}–{num(temporal['max'])} |")
    lines += [
        "", ("The temporal check keeps the selected coordinate fixed and includes missing or rejected days "
        "in the JSON. Same-day daily bins do not establish coincidence with the satellite overpass. "
        "Spatial or temporal agreement cannot rule out shared atmospheric, adjacency or shallow-bottom bias."), "",
        "## Use in the benchmark", "",
        (f"{sum(t['candidate'] is not None for t in record['targets'])} of {len(record['targets'])} dates "
        "have product-valid candidates for input sensitivity work. Reef representativeness remains "
        "unestablished for every date. The choice of fallback, uncertainty treatment, spectral extension and "
        "r7 outcome criteria must be fixed before corrected habitat scoring. KD490 alone cannot identify "
        "absorption and backscatter separately. The winter-QAA spectral shape remains an assumption."), "",
        (f"The [Copernicus manual]({MANUAL}), II.4.3, describes the uncertainty as a percentage error bar "
        "and caps it at 327.67%. This audit treats capped values as unquantified and rejects them. "
        "It does not interpret other error bars as standard deviations or 95% intervals."), "",
        "## Reproduction", "", "```bash", "kelpobserve kd490 qualify", "```", "",
        ("The command verifies the frozen subset and writes a new JSON evidence file and this report under "
        "`out/kd490/`. Each run records the subset hashes, audit code hashes and policy. "
        "It performs no network fetch or habitat scoring."), "",
        f"Frozen subset SHA-256: `{record['subset']['sha256']}`.", "",
    ]
    return "\n".join(lines)


def write(record: dict, destination: Path | None = None) -> Path:
    from kelpobserve.config import REPO

    folder = destination or REPO / "out" / "kd490" / datetime.now(UTC).strftime("qualification-%Y%m%dT%H%M%S%fZ")
    folder.mkdir(parents=True, exist_ok=False)
    (folder / "qualification.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    (folder / "qualification.md").write_text(markdown(record))
    return folder
