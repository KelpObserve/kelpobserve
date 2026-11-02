"""Sesimbra habitat benchmark v1: pre-registration stages.

Implements the stages of docs/plan.md that
need no model fitting. Commands are separate, following the component-benchmark
pattern:

- inventory: versioned three-campaign reference dataset with source revisions,
  reconciliation evidence, depth-consistency QC and join classes. Cover targets
  are written outside the event table, and September targets to their own file,
  which fitting never opens.
- reference-region: the geometry-only PNeo deep-water reference region.
- fetch: Sentinel-2 L1C candidates, downloaded in pre-QC rank order.
- prepare: the pinned ACOLITE revision on every downloaded Sentinel-2 candidate and
  both PNeo deliveries, with registered sensor-specific settings. Failed scenes are
  recorded, never dropped.
- scene-qc: the registered scene checks and Sentinel-2 ranking on prepared scenes.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from kelpobserve.coastal.provenance import complete_run, file_digest, new_run_directory

ROOT = Path(__file__).resolve().parents[1]
REGISTRATION = ROOT / "registration"
OUTPUT = ROOT / "out"
E5S_RESULTS = ROOT.parent / "optics" / "reports" / "e5s-kelp-spectra" / "results.json"

PT_MONTHS = {
    "jan": 1,
    "fev": 2,
    "mar": 3,
    "abr": 4,
    "mai": 5,
    "jun": 6,
    "jul": 7,
    "ago": 8,
    "set": 9,
    "out": 10,
    "nov": 11,
    "dez": 12,
}
KEY = ["campaign", "stratum", "transect", "quadrat"]
TARGET_FIELDS = [
    "occurrence",
    "cover_pct",
    "cover_status",
    "plant_rows",
    "juvenile_total",
    "adult_total",
]
JUVENILE, ADULT, COVER = "Juvenile (< 30cm)", "Adult (>= 30cm)", "% cover"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _flag(value) -> bool:
    """True only for real booleans; outer joins leave NaN where a side is missing."""
    return isinstance(value, (bool, np.bool_)) and bool(value)


def _number(value):
    return None if value is None or (isinstance(value, float) and math.isnan(value)) else value


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + "\n")
    temp.replace(path)


def load_spec(path: Path = REGISTRATION / "experiment.json") -> dict:
    spec = json.loads(Path(path).read_text())
    require(spec.get("experiment_id") == "sesimbra-habitat-v1", "Unexpected registration file.")
    return spec


def campaign_of_month(month: int, spec: dict) -> str:
    matches = [name for name, c in spec["campaigns"].items() if c["month"] == month]
    require(len(matches) == 1, f"No registered campaign for month {month}.")
    return matches[0]


def parse_transect_date(text, year: int) -> date:
    """Parse GeoPackage dates such as '01/set' (Portuguese month abbreviations)."""
    day, _, month = str(text).strip().lower().partition("/")
    require(day.isdigit() and month in PT_MONTHS, f"Unrecognised transect date: {text!r}")
    return date(year, PT_MONTHS[month], int(day))


def read_workbook(path, sheet) -> pd.DataFrame:
    rows = pd.read_excel(path, sheet_name=sheet)
    rows.insert(0, "source_row", rows.index + 2)  # Excel row number; header is row 1
    return rows


def aggregate_events(rows: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """One record per surveyed quadrat. Plant rows never become extra samples."""
    target = spec["taxon_mapping"]["target"].strip().lower()
    absence = spec["taxon_mapping"]["absence"].strip().lower()
    taxa = rows["Species"].astype(str).str.strip().str.lower()
    unknown = sorted(set(taxa) - {target, absence})
    require(not unknown, f"Unmapped taxa {unknown}; extend taxon_mapping explicitly.")
    rows = rows.assign(_taxon=taxa, _date=pd.to_datetime(rows["Date"]).dt.date)
    records = []
    for (day, depth, transect, quadrat), group in rows.groupby(
        ["_date", "Depth", "Transect", "Quadrat"], sort=True
    ):
        campaign = campaign_of_month(day.month, spec)
        strata = spec["campaigns"][campaign]["strata"]["workbook"]
        label = str(int(depth))
        require(label in strata, f"Unmapped workbook depth label {label} m in {campaign}.")
        source_rows = [int(r) for r in group["source_row"]]
        empty = group["_taxon"] == absence
        require(
            not (empty.any() and (~empty).any()),
            f"Quadrat has both EMPTY and plant rows: {source_rows}",
        )
        covers = group[COVER].dropna()
        require(len(covers) <= 1, f"Repeated cover measurements in one quadrat: {source_rows}")
        cover = float(covers.iloc[0]) if len(covers) else math.nan
        if empty.all():
            require(
                math.isnan(cover) or cover == 0.0, f"EMPTY record with nonzero cover: {source_rows}"
            )
            status = "absence_derived_zero" if math.isnan(cover) else "absence_recorded_zero"
            cover, occurrence, plants = 0.0, False, 0
        else:
            status = "missing_positive" if math.isnan(cover) else "measured"
            occurrence, plants = True, len(group)
        records.append(
            {
                "campaign": campaign,
                "stratum": strata[label],
                "transect": int(transect),
                "quadrat": int(quadrat),
                "workbook_date": day.isoformat(),
                "workbook_depth_label_m": int(depth),
                "source_rows": json.dumps(source_rows),
                "occurrence": occurrence,
                "cover_pct": cover,
                "cover_status": status,
                "plant_rows": plants,
                "juvenile_total": float(group[JUVENILE].sum()),
                "adult_total": float(group[ADULT].sum()),
            }
        )
    return pd.DataFrame(records)


def _numbered(text, prefix: str) -> int:
    match = re.fullmatch(prefix + r"(\d+)", str(text).strip())
    require(match is not None, f"Unrecognised {prefix} identifier: {text!r}")
    return int(match.group(1))


def read_transects(path, layer: str, spec: dict):
    import geopandas as gpd

    points = gpd.read_file(path, layer=layer)
    records = []
    for fid, row in points.iterrows():
        day = parse_transect_date(row["Date"], spec["year"])
        campaign = campaign_of_month(day.month, spec)
        label = str(row["Depth"]).strip()
        strata = spec["campaigns"][campaign]["strata"]["transects"]
        require(label in strata, f"Unmapped transect depth label {label!r} in {campaign}.")
        records.append(
            {
                "campaign": campaign,
                "stratum": strata[label],
                "transect": _numbered(row["Transect number"], "T"),
                "quadrat": _numbered(row["Quadrat"], "Q"),
                "transects_date": day.isoformat(),
                "transects_depth_label": label,
                "transects_fid": int(fid),
                "geometry": row.geometry,
            }
        )
    return gpd.GeoDataFrame(records, geometry="geometry", crs=points.crs)


def join_events(events: pd.DataFrame, points):
    """Outer join on campaign/stratum/transect/quadrat; expose every non-unique match."""
    import geopandas as gpd

    dup_events = events.duplicated(KEY, keep=False)
    dup_points = points.duplicated(KEY, keep=False)
    left = events.assign(_dup_workbook=dup_events)
    right = pd.DataFrame(points).assign(_dup_transects=dup_points)
    joined = left.merge(right, on=KEY, how="outer", indicator=True)
    issues = []
    for _, row in joined.iterrows():
        found = []
        if row["_merge"] == "right_only":
            found.append("no_workbook_event")
        if row["_merge"] == "left_only":
            found.append("no_transect_point")
        if _flag(row.get("_dup_workbook")):
            found.append("duplicate_workbook_key")
        if _flag(row.get("_dup_transects")):
            found.append("duplicate_transect_key")
        issues.append(found)
    joined["join_issues"] = issues
    joined["event_id"] = [
        f"{c}-{s}-T{int(t)}-Q{int(q)}"
        for c, s, t, q in zip(joined.campaign, joined.stratum, joined.transect, joined.quadrat)
    ]
    joined = joined.drop(columns=["_merge", "_dup_workbook", "_dup_transects"])
    return gpd.GeoDataFrame(joined, geometry="geometry", crs=points.crs)


def check_dates(joined, spec: dict, extra: dict[str, list[str]] | None = None) -> None:
    """Every recorded date must be one of the registered plausible survey dates."""
    for _, row in joined.iterrows():
        plausible = set(spec["campaigns"][row["campaign"]]["plausible_survey_dates"])
        recorded = {row.get("workbook_date"), row.get("transects_date")}
        recorded |= set((extra or {}).get(row["event_id"], []))
        recorded = {d for d in recorded if isinstance(d, str)}
        unregistered = sorted(recorded - plausible)
        require(
            not unregistered,
            f"{row['event_id']}: unregistered survey date(s) {unregistered}; "
            "register them as plausible or reconcile the source.",
        )


def previous_timed_positions(spec: dict) -> pd.DataFrame:
    """Earlier July GeoPackage: tape positions and dive times, matched to Qn by coordinates."""
    import geopandas as gpd

    source = spec["sources"]["previous_transects"]
    timed = pd.DataFrame(gpd.read_file(source["path"], layer=source["timed_layer"])).drop(
        columns="geometry"
    )
    numbered = pd.DataFrame(gpd.read_file(source["path"], layer=source["numbered_layer"])).drop(
        columns="geometry"
    )
    merged = timed.merge(
        numbered, on=["GPS Lat", "GPS Long"], suffixes=("_timed", "_numbered"), how="left"
    )
    records = []
    for _, row in merged.iterrows():
        day = parse_transect_date(row["Date_timed"], spec["year"])
        campaign = campaign_of_month(day.month, spec)
        label = str(row["Depth_timed"]).strip()
        records.append(
            {
                "campaign": campaign,
                "stratum": spec["campaigns"][campaign]["strata"]["transects"][label],
                "transect": _numbered(row["Transect number_timed"], "T"),
                "quadrat": _numbered(row["Quadrat_numbered"], "Q"),
                "tape_position_m": _numbered(row["Quadrat_timed"], "#"),
                "dive_time_local": str(row["Time"]),
                "previous_timed_date": day.isoformat(),
                "previous_numbered_date": parse_transect_date(
                    row["Date_numbered"], spec["year"]
                ).isoformat(),
            }
        )
    return pd.DataFrame(records)


def quadrat_order_evidence(timed: pd.DataFrame) -> list[dict]:
    """Per transect: does Q number follow survey order (time) and descending tape position?"""
    evidence = []
    for (campaign, stratum, transect), group in timed.groupby(["campaign", "stratum", "transect"]):
        group = group.sort_values("quadrat")
        evidence.append(
            {
                "campaign": campaign,
                "stratum": stratum,
                "transect": int(transect),
                "q_to_tape_m": {
                    f"Q{int(q)}": int(t) for q, t in zip(group.quadrat, group.tape_position_m)
                },
                "q_follows_survey_time": bool(group.dive_time_local.is_monotonic_increasing),
                "q1_is_tape_20m": bool(int(group.tape_position_m.iloc[0]) == 20),
            }
        )
    return evidence


def sample_lidar(points, path) -> pd.DataFrame:
    """Containing lidar cell plus the 3x3 neighbourhood range; no interpolation."""
    import rasterio

    with rasterio.open(path) as src:
        projected = points.to_crs(src.crs)
        data = src.read(1).astype("float64")
        if src.nodata is not None and not np.isnan(src.nodata):
            data[data == src.nodata] = np.nan
        values = []
        for geom in projected.geometry:
            if geom is None or geom.is_empty:
                values.append((math.nan, math.nan, math.nan))
                continue
            r, c = src.index(geom.x, geom.y)
            if not (0 <= r < src.height and 0 <= c < src.width):
                values.append((math.nan, math.nan, math.nan))
                continue
            window = data[max(r - 1, 0) : r + 2, max(c - 1, 0) : c + 2]
            finite = window[np.isfinite(window)]
            values.append(
                (
                    float(data[r, c]),
                    float(finite.min()) if finite.size else math.nan,
                    float(finite.max()) if finite.size else math.nan,
                )
            )
    return pd.DataFrame(values, columns=["lidar_depth_m", "lidar_nbr_min_m", "lidar_nbr_max_m"])


def water_level_interval(campaign: str, spec: dict, levels: pd.DataFrame | None):
    """Water-level range above LAT over all plausible survey dates, or None if incomplete."""
    if levels is None:
        return None
    rows = levels[levels.campaign == campaign]
    require((rows.datum.astype(str).str.upper() == "LAT").all(), "Water levels must be above LAT.")
    plausible = set(spec["campaigns"][campaign]["plausible_survey_dates"])
    if not plausible <= set(rows.date.astype(str)):
        return None
    rows = rows[rows.date.astype(str).isin(plausible)]
    return float(rows.min_m.min()), float(rows.max_m.max())


def depth_consistency(label_m, lidar, nbr_min, nbr_max, levels, tolerance_m):
    """Recorded depth minus (lidar + water level): status, residual range, neighbour check."""
    if not np.isfinite(label_m):
        return "no_label", math.nan, math.nan, None
    if not np.isfinite(lidar):
        return "no_lidar", math.nan, math.nan, None
    if levels is None:
        return "pending_water_level", math.nan, math.nan, None
    low, high = label_m - lidar - levels[1], label_m - lidar - levels[0]
    if high >= -tolerance_m and low <= tolerance_m:
        return "pass", low, high, None
    reconciles = bool(
        label_m - nbr_min - levels[0] >= -tolerance_m
        and label_m - nbr_max - levels[1] <= tolerance_m
    )
    return "flag", low, high, reconciles


def classify_join(campaign, issues, depth_status, spec) -> tuple[str, list[str]]:
    reasons = list(issues)
    if reasons:
        return "provisional", reasons
    if depth_status == "flag":
        return "provisional", ["depth_consistency_flag"]
    if depth_status == "no_lidar":
        return "provisional", ["depth_consistency_not_evaluable"]
    if campaign in spec["ccmar_confirmation"].get("confirmed_campaigns", []):
        return (
            ("confirmed", []) if depth_status == "pass" else ("pending", ["water_levels_missing"])
        )
    if depth_status == "pending_water_level":
        return "pending", ["water_levels_missing"]
    return "reconciled_by_rule", []


def workbook_revisions(previous: pd.DataFrame, current: pd.DataFrame, spec: dict) -> dict:
    """Cell- and event-level differences for the campaigns both workbooks contain."""
    months = set(pd.to_datetime(previous["Date"]).dt.month)
    overlap = current[pd.to_datetime(current["Date"]).dt.month.isin(months)].reset_index(drop=True)
    columns = [c for c in previous.columns if c in overlap.columns and c != "source_row"]
    report: dict = {"previous_rows": len(previous), "current_rows_same_campaigns": len(overlap)}
    if len(previous) == len(overlap):
        a, b = previous[columns].reset_index(drop=True), overlap[columns]
        differs = ~((a == b) | (a.isna() & b.isna()))
        report["cells_differing_by_column"] = {c: int(n) for c, n in differs.sum().items() if n}
    else:
        report["cells_differing_by_column"] = "row counts differ; compare events only"
    old, new = aggregate_events(previous, spec), aggregate_events(overlap, spec)
    merged = old.merge(new, on=KEY, how="outer", suffixes=("_previous", "_current"), indicator=True)
    report["events_only_previous"] = int((merged._merge == "left_only").sum())
    report["events_only_current"] = int((merged._merge == "right_only").sum())
    both = merged[merged._merge == "both"]
    for field in ["cover_pct", "plant_rows", "juvenile_total", "adult_total"]:
        x, y = both[f"{field}_previous"], both[f"{field}_current"]
        changed = both[~((x == y) | (x.isna() & y.isna()))]
        report[f"events_with_changed_{field}"] = [
            {
                "campaign": row.campaign,
                "stratum": row.stratum,
                "transect": int(row.transect),
                "quadrat": int(row.quadrat),
                "previous": _number(float(getattr(row, f"{field}_previous"))),
                "current": _number(float(getattr(row, f"{field}_current"))),
            }
            for row in changed.itertuples(index=False)
        ]
    return report


def git_state(paths: list[Path]) -> dict:
    def run(*args):
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True
        ).stdout.strip()

    dirty = run("status", "--porcelain", "--", *[str(p) for p in paths])
    return {"commit": run("rev-parse", "HEAD"), "uncommitted_changes": bool(dirty)}


def code_provenance() -> dict:
    """Runner checksum and git state, captured when a stage starts, not when it ends."""
    return {
        "runner_sha256": file_digest(Path(__file__)),
        "git": git_state([Path(__file__), REGISTRATION]),
    }


def inventory(spec_path: Path = REGISTRATION / "experiment.json", output: Path = OUTPUT) -> Path:
    import geopandas as gpd

    spec = load_spec(spec_path)
    code = code_provenance()
    sources = spec["sources"]
    rows = read_workbook(sources["kelp_workbook"], sources["kelp_sheet"])
    events = aggregate_events(rows, spec)
    points = read_transects(sources["transects"]["path"], sources["transects"]["layer"], spec)
    joined = join_events(events, points)

    timed = previous_timed_positions(spec)
    joined = joined.merge(timed, on=KEY, how="left")
    extra = {
        e: [d for d in (t, n) if isinstance(d, str)]
        for e, t, n in zip(
            joined.event_id, joined.previous_timed_date, joined.previous_numbered_date
        )
    }
    check_dates(joined, spec, extra)

    lidar = sample_lidar(joined, sources["lidar"])
    joined = pd.concat([joined.reset_index(drop=True), lidar], axis=1)
    levels_path = spec["water_levels"].get("path")
    levels = pd.read_csv(REGISTRATION / levels_path) if levels_path else None
    qc = spec["depth_consistency"]
    label_column = {"workbook": "workbook_depth_label_m", "transects": "transects_label_m"}[
        qc["label_source"]
    ]
    joined["transects_label_m"] = joined.transects_depth_label.str.rstrip("m").astype(float)
    statuses, lows, highs, neighbours, classes, reasons = [], [], [], [], [], []
    for _, row in joined.iterrows():
        status, low, high, neighbour = depth_consistency(
            row[label_column],
            row.lidar_depth_m,
            row.lidar_nbr_min_m,
            row.lidar_nbr_max_m,
            water_level_interval(row.campaign, spec, levels),
            qc["tolerance_m"],
        )
        join_class, why = classify_join(row.campaign, row.join_issues, status, spec)
        statuses.append(status)
        lows.append(low)
        highs.append(high)
        neighbours.append(neighbour)
        classes.append(join_class)
        reasons.append(json.dumps(why))
    joined = joined.assign(
        depth_residual_nominal_m=joined[label_column] - joined.lidar_depth_m,
        depth_qc_status=statuses,
        depth_residual_low_m=lows,
        depth_residual_high_m=highs,
        depth_qc_neighbour_reconciles=pd.array(neighbours, dtype="boolean"),
        join_class=classes,
        join_reasons=reasons,
        date_discrepancy=joined.workbook_date != joined.transects_date,
        depth_label_discrepancy=joined.workbook_depth_label_m != joined.transects_label_m,
        join_issues=joined.join_issues.map(json.dumps),
    )
    joined["eligible_registered"] = joined.join_class.isin(["confirmed", "reconciled_by_rule"]) & (
        joined.depth_qc_status == "pass"
    )

    directory = new_run_directory(output / "inventory")
    evaluation = [n for n, c in spec["campaigns"].items() if c["role"] == "evaluation"]
    targets = pd.DataFrame(joined[["event_id", "campaign", *TARGET_FIELDS]])
    held_out = targets.campaign.isin(evaluation)
    targets[~held_out].to_csv(directory / "targets_development.csv", index=False)
    targets[held_out].to_csv(directory / "targets_september.csv", index=False)
    table = gpd.GeoDataFrame(
        joined.drop(columns=TARGET_FIELDS), geometry="geometry", crs=points.crs
    )
    table.to_file(directory / "events.gpkg", layer="events", driver="GPKG")
    table.drop(columns="geometry").assign(
        longitude=table.geometry.x, latitude=table.geometry.y
    ).to_csv(directory / "events.csv", index=False)

    previous = read_workbook(sources["previous_workbook"], sources["kelp_sheet"])
    write_json(directory / "revisions.json", {"workbook": workbook_revisions(previous, rows, spec)})
    reconciliation = {
        "plausible_survey_dates": {
            n: c["plausible_survey_dates"] for n, c in spec["campaigns"].items()
        },
        "recorded_dates": {
            campaign: {
                "workbook": sorted(set(group.workbook_date.dropna())),
                "transects": sorted(set(group.transects_date.dropna())),
                "previous_timed": sorted(set(group.previous_timed_date.dropna())),
                "previous_numbered": sorted(set(group.previous_numbered_date.dropna())),
            }
            for campaign, group in joined.groupby("campaign")
        },
        "depth_labels": {
            f"{campaign}/{stratum}": {
                "workbook_m": sorted({int(v) for v in group.workbook_depth_label_m.dropna()}),
                "transects": sorted(set(group.transects_depth_label.dropna())),
            }
            for (campaign, stratum), group in joined.groupby(["campaign", "stratum"])
        },
        "quadrat_order": {**spec["quadrat_order"], "evidence": quadrat_order_evidence(timed)},
        "ccmar_confirmation": spec["ccmar_confirmation"],
        "water_levels": spec["water_levels"],
        "depth_consistency": qc,
    }
    write_json(directory / "reconciliation.json", reconciliation)

    inputs = [
        Path(sources["kelp_workbook"]),
        Path(sources["transects"]["path"]),
        Path(sources["previous_workbook"]),
        Path(sources["previous_transects"]["path"]),
        Path(sources["lidar"]),
        Path(spec_path),
    ] + ([REGISTRATION / levels_path] if levels_path else [])
    outputs = [
        "events.gpkg",
        "events.csv",
        "targets_development.csv",
        "targets_september.csv",
        "revisions.json",
        "reconciliation.json",
    ]
    summary = {
        "events": int(len(joined)),
        "expected_events": spec["expected_events"],
        "events_by_campaign": joined.campaign.value_counts().sort_index().to_dict(),
        "join_class": joined.join_class.value_counts().sort_index().to_dict(),
        "depth_qc_status": joined.depth_qc_status.value_counts().sort_index().to_dict(),
        "eligible_registered": int(joined.eligible_registered.sum()),
    }
    require(
        len(joined) == spec["expected_events"],
        f"Expected {spec['expected_events']} events; found {len(joined)}. "
        f"Diagnostics kept in {directory} without a completion manifest.",
    )
    complete_run(
        output / "inventory",
        directory,
        {
            "stage": "inventory",
            "created_at": datetime.now(UTC).isoformat(),
            "inputs": {str(p): file_digest(p) for p in inputs},
            "outputs": {name: file_digest(directory / name) for name in outputs},
            **code,
            "evaluation_targets": "targets_september.csv is opened only by the evaluation command.",
            "summary": summary,
        },
    )
    return directory


def survey_window(plausible: list[str], days: int) -> tuple[date, date]:
    """Scene dates within `days` of every plausible survey date."""
    dates = sorted(date.fromisoformat(d) for d in plausible)
    start, end = dates[-1] - timedelta(days=days), dates[0] + timedelta(days=days)
    require(
        start <= end, "No acquisition date lies within the window of every plausible survey date."
    )
    return start, end


def _utc(value: datetime) -> datetime:
    return value if value.tzinfo else value.replace(tzinfo=UTC)


def rank_candidates(candidates, plausible: list[str], pneo_utc: datetime | None) -> list[dict]:
    """Pre-QC ranking: max date offset, then time from PNeo, then identifier.

    Valid-water coverage (the third registered criterion) needs processed imagery
    and is applied after ACOLITE; it never reorders scenes ahead of these keys.
    """
    dates = [date.fromisoformat(d) for d in plausible]
    ranked = []
    for candidate in candidates:
        sensed = _utc(candidate.sensing_datetime)
        ranked.append(
            {
                "candidate": candidate,
                "max_offset_days": max(abs((sensed.date() - d).days) for d in dates),
                "pneo_offset_s": abs((sensed - pneo_utc).total_seconds()) if pneo_utc else None,
            }
        )
    return sorted(
        ranked,
        key=lambda r: (
            r["max_offset_days"],
            math.inf if r["pneo_offset_s"] is None else r["pneo_offset_s"],
            r["candidate"].name,
        ),
    )


def latest_baseline(candidates) -> list:
    """Keep the newest processing baseline of each datatake/tile."""
    best: dict = {}
    for candidate in candidates:
        key = (_utc(candidate.sensing_datetime), candidate.tile)
        baseline = re.search(r"_N(\d{4})_", candidate.name)
        rank = int(baseline.group(1)) if baseline else -1
        if key not in best or rank > best[key][0]:
            best[key] = (rank, candidate)
    return [c for _, c in best.values()]


def read_env_file(path: Path) -> dict[str, str]:
    values = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        values[key.strip().removeprefix("export ").strip()] = value.strip().strip("'\"")
    return values


PLANETARY_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"


def window_days(spec: dict) -> tuple[int, bool]:
    """Matchup window in days, and whether the registered amendment supplies it."""
    s2 = spec["sentinel2"]
    amended = s2.get("window_days_amended")
    return (int(amended), True) if amended else (int(s2["window_days"]), False)


def screen_key(when: datetime, platform: str) -> str:
    """Datatake minute plus platform: the only fields the CDSE and STAC catalogues share."""
    return f"{_utc(when):%Y-%m-%dT%H:%M}|{platform}"


def platform_of(name: str) -> str:
    return f"Sentinel-{name[1]}{name[2]}"


def load_screen(spec: dict) -> tuple[dict, float | None]:
    """Registered glint screen, keyed by datatake minute and platform."""
    path = REGISTRATION / "glint-screen.json"
    threshold = spec["sentinel2"].get("screen_glint_max")
    if not path.is_file() or threshold is None:
        return {}, None
    document = json.loads(path.read_text())
    return {row["key"]: row for row in document["scenes"]}, float(threshold)


def screen_fields(candidate, screen: dict, threshold: float | None) -> dict:
    """Screen value for a CDSE candidate, matched on datatake minute and platform."""
    row = screen.get(screen_key(candidate.sensing_datetime, platform_of(candidate.name)))
    if row is None or threshold is None:
        return {"screen_b11_median": None, "screened_out": False}
    return {"screen_b11_median": row["b11_median"], "screened_out": bool(row["screened_out"])}


def quadrat_window_bounds(buffer_m: float = 250.0) -> tuple[float, float, float, float]:
    """Quadrat bounding box in EPSG:32629, padded. The screen reads only this window."""
    import geopandas as gpd
    from kelpobserve.coastal.provenance import resolve_run

    events = gpd.read_file(resolve_run(OUTPUT / "inventory") / "events.gpkg", layer="events")
    x0, y0, x1, y1 = events.to_crs(32629).total_bounds
    return (x0 - buffer_m, y0 - buffer_m, x1 + buffer_m, y1 + buffer_m)


def glint_screen() -> Path:
    """Median L2A 1610 nm over the quadrat area for every candidate date.

    An image-quality pre-filter registered as an amendment: it never reads cover, and
    it decides only which candidates are worth the full ACOLITE run.
    """
    import geopandas as gpd
    import planetary_computer as pc
    import rasterio
    from pystac_client import Client
    from rasterio.windows import from_bounds

    spec = load_spec()
    code = code_provenance()
    days, amended = window_days(spec)
    bounds = quadrat_window_bounds()
    lon0, lat0, lon1, lat1 = (
        gpd.GeoSeries.from_xy([bounds[0], bounds[2]], [bounds[1], bounds[3]], crs=32629)
        .to_crs(4326)
        .total_bounds
    )
    client = Client.open(PLANETARY_STAC, modifier=pc.sign_inplace)
    seen, rows = set(), []
    for name, campaign in spec["campaigns"].items():
        start, end = survey_window(campaign["plausible_survey_dates"], days)
        found = client.search(
            collections=["sentinel-2-l2a"],
            bbox=[lon0, lat0, lon1, lat1],
            datetime=f"{start.isoformat()}/{end.isoformat()}",
        ).items()
        for item in sorted(found, key=lambda i: i.datetime):
            properties = item.properties
            key = screen_key(item.datetime, properties.get("platform", ""))
            if key in seen:
                continue
            seen.add(key)
            offset = (
                -1000 if str(properties.get("s2:processing_baseline", "00.00")) >= "04.00" else 0
            )
            row = {
                "key": key,
                "campaign": name,
                "sensing_utc": _utc(item.datetime).isoformat(),
                "platform": properties.get("platform"),
                "relative_orbit": properties.get("sat:relative_orbit"),
                "tile_cloud_pct": properties.get("eo:cloud_cover"),
            }
            with rasterio.open(item.assets["B11"].href) as src:
                data = src.read(1, window=from_bounds(*bounds, transform=src.transform)).astype(
                    float
                )
            data[data == 0] = np.nan
            reflectance = (data + offset) / 10000.0
            row["b11_median"] = round(float(np.nanmedian(reflectance)), 4)
            row["b11_p10"] = round(float(np.nanpercentile(reflectance, 10)), 4)
            rows.append(row)
            orbit, median = row["relative_orbit"], row["b11_median"]
            print(f"  {row['sensing_utc'][:16]} R{orbit:03d} B11={median:.4f}", flush=True)
    threshold = float(spec["sentinel2"]["screen_glint_max"])
    for row in rows:
        row["screened_out"] = row["b11_median"] > threshold
    path = REGISTRATION / "glint-screen.json"
    write_json(
        path,
        {
            "experiment_id": spec["experiment_id"],
            "created_at": datetime.now(UTC).isoformat(),
            "amendment": "glint-window-and-screen",
            "window_days": days,
            "window_amended": amended,
            "screen_glint_max": threshold,
            "bounds_epsg32629": [round(v, 1) for v in bounds],
            "collection": "sentinel-2-l2a (Planetary Computer)",
            "definition": spec["amendments"][0]["screen_definition"],
            "scenes": rows,
            **code,
        },
    )
    return path


def fetch(env_file: Path | None, campaigns: list[str] | None, top: int, download: bool) -> None:
    from kelpobserve.coastal.acquire.cdse import (
        CDSECredentials,
        download_scene,
        get_access_token,
        search_scenes,
    )

    spec = load_spec()
    s2 = spec["sentinel2"]
    credentials = CDSECredentials.from_env(read_env_file(env_file) if env_file else None)
    token = get_access_token(credentials)
    cache = ROOT / s2["l1c_cache"]
    days, amended = window_days(spec)
    screen, threshold = load_screen(spec)
    for name, campaign in spec["campaigns"].items():
        if campaigns and name not in campaigns:
            continue
        start, end = survey_window(campaign["plausible_survey_dates"], days)
        pneo = campaign.get("pneo")
        pneo_utc = (
            datetime.fromisoformat(
                spec["sources"]["pneo"][pneo]["acquired_utc"].replace("Z", "+00:00")
            )
            if pneo
            else None
        )
        found = search_scenes(
            tile=s2["tile"], start=start, end=end, product_type=s2["product_type"], token=token
        )
        ranked = rank_candidates(
            latest_baseline(found), campaign["plausible_survey_dates"], pneo_utc
        )
        write_json(
            OUTPUT / "scenes" / f"{name}-candidates.json",
            {
                "campaign": name,
                "window": [start.isoformat(), end.isoformat()],
                "searched_at": datetime.now(UTC).isoformat(),
                "candidates": [
                    {
                        "rank_pre_qc": i + 1,
                        "name": r["candidate"].name,
                        "product_id": r["candidate"].product_id,
                        "sensing_utc": _utc(r["candidate"].sensing_datetime).isoformat(),
                        "tile_cloud_pct": r["candidate"].cloud_cover_pct,
                        "size_mb": r["candidate"].size_mb,
                        "online": r["candidate"].online,
                        "max_offset_days": r["max_offset_days"],
                        "pneo_offset_s": r["pneo_offset_s"],
                        **screen_fields(r["candidate"], screen, threshold),
                    }
                    for i, r in enumerate(ranked)
                ],
                "window_days": days,
                "window_amended": amended,
            },
        )
        print(f"{name}: {len(ranked)} candidates in {start}..{end}", flush=True)
        wanted = [
            r
            for r in ranked
            if not screen_fields(r["candidate"], screen, threshold).get("screened_out")
        ]
        for r in wanted[: top if download else 0]:
            path = download_with_backoff(download_scene, r["candidate"], cache / name, credentials)
            print(f"  downloaded {path.name}", flush=True)


def download_with_backoff(
    download_scene, candidate, directory, credentials, attempts=3, wait_s=120
):
    """CDSE answers 429 when the account's request or concurrency quota is spent."""
    import time

    import requests

    for attempt in range(1, attempts + 1):
        try:
            return download_scene(candidate, directory, credentials=credentials)
        except requests.HTTPError as error:
            if error.response is None or error.response.status_code != 429 or attempt == attempts:
                raise
            delay = int(error.response.headers.get("Retry-After", wait_s))
            print(
                f"  CDSE quota reached (429); retry {attempt}/{attempts - 1} in {delay}s",
                flush=True,
            )
            time.sleep(delay)


def s2_sensing_time(name: str) -> datetime:
    """Datatake start from a Sentinel-2 product name."""
    match = re.search(r"_(\d{8}T\d{6})_", name)
    require(match is not None, f"No sensing time in {name}")
    return datetime.strptime(match.group(1), "%Y%m%dT%H%M%S").replace(tzinfo=UTC)


def water_levels(env_file: Path | None) -> tuple[Path, Path]:
    """Modelled sea level at the site, expressed above the registered chart-datum proxy.

    Writes the survey-day ranges the depth check consumes and the level at each
    acquisition. Both carry the datum proxy and its uncertainty.
    """
    import os

    import copernicusmarine as cm
    import geopandas as gpd
    from kelpobserve.coastal.provenance import resolve_run

    spec = load_spec()
    code = code_provenance()
    config = spec["water_levels"]
    source = config["source"]
    if env_file:
        os.environ.update(read_env_file(env_file))
    x0, y0, x1, y1 = quadrat_window_bounds(0.0)
    centre = gpd.GeoSeries.from_xy([(x0 + x1) / 2], [(y0 + y1) / 2], crs=32629).to_crs(4326)
    lon, lat = float(centre.x.iloc[0]), float(centre.y.iloc[0])
    data = cm.open_dataset(
        dataset_id=source["dataset_id"],
        variables=[source["variable"]],
        minimum_longitude=lon - 0.2,
        maximum_longitude=lon + 0.2,
        minimum_latitude=lat - 0.2,
        maximum_latitude=lat + 0.2,
    )
    grid = data[source["variable"]]
    first = grid.isel(time=0)
    best, best_km = None, math.inf
    for i, cell_lat in enumerate(grid.latitude.values):
        for j, cell_lon in enumerate(grid.longitude.values):
            if not np.isfinite(first.isel(latitude=i, longitude=j).item()):
                continue
            km = math.hypot((cell_lat - lat) * 111.0, (cell_lon - lon) * 87.0)
            if km < best_km:
                best, best_km = (i, j, float(cell_lat), float(cell_lon)), km
    require(best is not None, "No model water cell near the quadrats.")
    series = grid.isel(latitude=best[0], longitude=best[1]).to_series().dropna()
    datum = float(series.min())
    above = series - datum
    station = f"{best[2]:.4f}N {abs(best[3]):.4f}W, {best_km:.1f} km from the quadrats"
    provider = source["provider"]

    rows = []
    for campaign, config_campaign in spec["campaigns"].items():
        for day in config_campaign["plausible_survey_dates"]:
            values = above[day]
            require(len(values) >= 20, f"Incomplete modelled day for {day}.")
            rows.append(
                {
                    "campaign": campaign,
                    "date": day,
                    "min_m": round(float(values.min()), 3),
                    "max_m": round(float(values.max()), 3),
                    "source": provider,
                    "station": station,
                    "datum": "LAT",
                }
            )
    survey_path = REGISTRATION / config["path"]
    pd.DataFrame(rows).to_csv(survey_path, index=False)

    hourly = above.copy()
    hourly.index = pd.to_datetime(hourly.index, utc=True)
    acquisitions = []
    prepared = json.loads((resolve_run(OUTPUT / "prepare") / "manifest.json").read_text())
    for scene in prepared["scenes"]:
        if scene["sensor"] == "pleiades_neo":
            when = datetime.fromisoformat(
                spec["sources"]["pneo"][spec["campaigns"][scene["campaign"]]["pneo"]][
                    "acquired_utc"
                ].replace("Z", "+00:00")
            )
        else:
            when = s2_sensing_time(Path(scene["input"]["archive"]).name)
        window = hourly.reindex(hourly.index.union([when])).interpolate("time")
        acquisitions.append(
            {
                "scene_id": scene["id"],
                "datetime_utc": when.isoformat(),
                "water_level_m": round(float(window.loc[when]), 3),
                "source": provider,
                "station": station,
                "datum": "LAT",
            }
        )
    acquisition_path = REGISTRATION / config["acquisition_path"]
    pd.DataFrame(acquisitions).to_csv(acquisition_path, index=False)

    write_json(
        REGISTRATION / "water-levels-provenance.json",
        {
            "experiment_id": spec["experiment_id"],
            "created_at": datetime.now(UTC).isoformat(),
            "amendment": "modelled-water-levels",
            "dataset_id": source["dataset_id"],
            "variable": source["variable"],
            "cell": {"latitude": best[2], "longitude": best[3], "distance_km": round(best_km, 2)},
            "datum_proxy_model_metres": round(datum, 3),
            "datum_period": [str(series.index.min())[:10], str(series.index.max())[:10]],
            "datum_uncertainty_m": config["datum_uncertainty_m"],
            "hours_used": int(len(series)),
            **code,
        },
    )
    return survey_path, acquisition_path


def reference_mask(valid, depth, habitat, pixel_m: float, rules: dict) -> np.ndarray:
    """Geometry and depth rules only: no reflectance enters this mask."""
    from scipy import ndimage

    land = ~np.isfinite(depth)
    with np.errstate(invalid="ignore"):
        deep = np.isfinite(depth) & (depth >= rules["min_dtm_depth_m"])
    return (
        valid
        & deep
        & (ndimage.distance_transform_edt(valid) * pixel_m >= rules["min_inside_footprint_m"])
        & (ndimage.distance_transform_edt(~land) * pixel_m >= rules["min_distance_from_land_m"])
        & (
            ndimage.distance_transform_edt(~habitat) * pixel_m
            >= rules["min_distance_from_habitat_m"]
        )
    )


def reference_region() -> Path:
    import geopandas as gpd
    import rasterio
    from pyproj import Transformer
    from rasterio import features
    from rasterio.transform import Affine
    from rasterio.warp import Resampling, reproject
    from scipy import ndimage
    from shapely import set_precision
    from shapely.geometry import mapping, shape
    from shapely.ops import transform as shp_transform
    from shapely.ops import unary_union

    spec = load_spec()
    code = code_provenance()
    config = spec["pneo_reference_region"]
    rules = config["rules"]
    masks, grid = [], None
    for acquisition in ("july", "september"):
        source = spec["sources"]["pneo"][acquisition]
        with rasterio.open(f"/vsizip/{source['archive']}/{source['ms_rgb']}") as src:
            factor = int(round(rules["working_resolution_m"] / src.res[0]))
            shape_ = (src.height // factor, src.width // factor)
            band = src.read(1, out_shape=shape_)
            transform = src.transform * src.transform.scale(
                src.width / shape_[1], src.height / shape_[0]
            )
            require(
                grid is None or grid == (transform, shape_, src.crs),
                "PNeo grids differ between dates.",
            )
            grid = (transform, shape_, src.crs)
            masks.append(band > 0)
    transform, _, crs = grid
    pad = int(
        math.ceil(
            max(rules["min_distance_from_land_m"], rules["min_distance_from_habitat_m"])
            / transform.a
        )
    )
    valid = np.pad(masks[0] & masks[1], pad, constant_values=False)
    transform = transform * Affine.translation(-pad, -pad)
    depth = np.full(valid.shape, np.nan, dtype="float32")
    with rasterio.open(spec["sources"]["dtm"]) as src:
        reproject(
            rasterio.band(src, 1),
            depth,
            dst_transform=transform,
            dst_crs=crs,
            resampling=Resampling.nearest,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )
    habitat_aoi = gpd.read_file(ROOT / spec["sources"]["habitat_aoi"]).to_crs(crs)
    habitat = features.rasterize(
        habitat_aoi.geometry, out_shape=valid.shape, transform=transform
    ).astype(bool)
    mask = reference_mask(valid, depth, habitat, transform.a, rules)
    require(
        mask.any(),
        "No area satisfies the PNeo reference rules; the corrected branch is unavailable.",
    )

    to_wgs84 = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform
    labels, count = ndimage.label(mask)
    parts = []
    for part in range(1, count + 1):
        inside = labels == part
        polygon = unary_union(
            [
                shape(g)
                for g, v in features.shapes(
                    inside.astype("uint8"), mask=inside, transform=transform
                )
                if v == 1
            ]
        )
        # Sub-millimetre coordinate precision, as for the six-scene regions.
        geographic = set_precision(shp_transform(to_wgs84, polygon), 1e-9)
        parts.append(
            {
                "type": "Feature",
                "properties": {
                    "role": "deepwater_reference",
                    "sensor": "pleiades_neo",
                    "part_id": f"pneo_ref_{part:02d}",
                    "area_m2": round(polygon.area, 1),
                    "dtm_depth_range_m": [
                        round(float(np.nanmin(depth[inside])), 1),
                        round(float(np.nanmax(depth[inside])), 1),
                    ],
                    "is_provisional": True,
                },
                "geometry": mapping(geographic),
            }
        )
    inputs = [
        spec["sources"]["dtm"],
        ROOT / spec["sources"]["habitat_aoi"],
        spec["sources"]["pneo"]["july"]["archive"],
        spec["sources"]["pneo"]["september"]["archive"],
    ]
    document = {
        "type": "FeatureCollection",
        "name": "pneo_reference",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
        "metadata": {
            "experiment_id": spec["experiment_id"],
            "created_at": datetime.now(UTC).isoformat(),
            "rules": rules,
            "provenance": config["provenance"],
            "land_mask": "EMODnet DTM NaN cells (land and depths at or above datum)",
            "total_area_m2": round(sum(f["properties"]["area_m2"] for f in parts), 1),
            "inputs": {str(p): file_digest(p) for p in inputs},
            **code,
        },
        "features": parts,
    }
    path = REGISTRATION / config["path"]
    write_json(path, document)
    return path


def read_settings_file(path: Path) -> dict[str, str]:
    """KEY=VALUE lines of an ACOLITE settings file; comments and blanks dropped."""
    values = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, sep, value = line.partition("=")
            require(sep == "=", f"Malformed ACOLITE setting in {path}: {line!r}")
            values[key.strip()] = value.strip()
    return values


def acolite_revision(path: Path) -> dict:
    def run(*args):
        return subprocess.run(
            ["git", *args], cwd=path, capture_output=True, text=True
        ).stdout.strip()

    return {
        "commit": run("rev-parse", "HEAD"),
        "describe": run("describe", "--tags"),
        "uncommitted_changes": bool(run("status", "--porcelain")),
    }


def verify_acolite(config: dict) -> dict:
    revision = acolite_revision(Path(config["path"]))
    require(
        revision["commit"] == config["commit"] and not revision["uncommitted_changes"],
        f"ACOLITE checkout {revision} differs from the registered revision {config['commit']}.",
    )
    return revision


def geographic_limit(bounds, crs, pad_deg: float = 0.0) -> tuple[float, float, float, float]:
    """ACOLITE's (south, west, north, east) from projected bounds."""
    from pyproj import Transformer

    west, south, east, north = Transformer.from_crs(
        crs, "EPSG:4326", always_xy=True
    ).transform_bounds(*bounds)
    return (
        round(south - pad_deg, 6),
        round(west - pad_deg, 6),
        round(north + pad_deg, 6),
        round(east + pad_deg, 6),
    )


def s2_scene_id(name: str) -> str:
    """S2A_MSIL1C_20260602T113911_N0512_R080_T29SMC_... -> s2-20260602-S2A-R080."""
    match = re.match(r"(S2[A-D])_MSIL1C_(\d{8})T\d{6}_N\d{4}_(R\d{3})_", name)
    require(match is not None, f"Unrecognised Sentinel-2 product name: {name}")
    platform, day, orbit = match.groups()
    return f"s2-{day}-{platform}-{orbit}"


def extract_pneo(source: dict, root: Path) -> Path:
    """Extract a delivered PNeo archive once; the marker records the archive checksum."""
    import zipfile

    archive = Path(source["archive"])
    product = root / Path(source["ms_rgb"]).parts[0]
    marker = product / ".extracted-from.json"
    digest = file_digest(archive)
    if marker.is_file() and json.loads(marker.read_text()).get("sha256") == digest:
        return product
    root.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as bundle:
        bundle.extractall(root)
    write_json(marker, {"archive": str(archive), "sha256": digest})
    return product


def prepare_scenes(spec: dict) -> list[dict]:
    """Every downloaded Sentinel-2 candidate plus both PNeo deliveries."""
    import geopandas as gpd
    import rasterio

    extent = gpd.read_file(ROOT / spec["sources"]["processing_extent"])
    s2_limit = geographic_limit(extent.total_bounds, extent.crs)
    scenes = []
    cache = ROOT / spec["sentinel2"]["l1c_cache"]
    for campaign in spec["campaigns"]:
        for bundle in sorted((cache / campaign).glob("*.zip")):
            scenes.append(
                {
                    "id": s2_scene_id(bundle.name),
                    "campaign": campaign,
                    "sensor": "sentinel2",
                    "input": bundle,
                    "limit": s2_limit,
                }
            )
    for campaign, config in spec["campaigns"].items():
        acquisition = config.get("pneo")
        if not acquisition:
            continue
        source = spec["sources"]["pneo"][acquisition]
        with rasterio.open(f"/vsizip/{source['archive']}/{source['ms_rgb']}") as src:
            limit = geographic_limit(src.bounds, src.crs, pad_deg=0.005)
        day = source["acquired_utc"][:10].replace("-", "")
        scenes.append(
            {
                "id": f"pneo-{day}",
                "campaign": campaign,
                "sensor": "pleiades_neo",
                "source": source,
                "limit": limit,
            }
        )
    return scenes


def prepare(only: list[str] | None = None) -> Path:
    """Run the pinned ACOLITE on every registered scene; failures are recorded, not hidden."""
    import time

    from kelpobserve.coastal.acolite import correct_scene, output_is_complete

    spec = load_spec()
    code = code_provenance()
    config = spec["acolite"]
    revision = verify_acolite(config)
    settings = {name: REGISTRATION / path for name, path in config["settings"].items()}
    overrides = {
        sensor: read_settings_file(settings[sensor]) for sensor in ("sentinel2", "pleiades_neo")
    }
    directory = new_run_directory(OUTPUT / "prepare")
    records = []
    for scene in prepare_scenes(spec):
        if only and scene["id"] not in only:
            continue
        output = ROOT / config["output_root"] / scene["id"]
        started = time.monotonic()
        record = {k: scene[k] for k in ("id", "campaign", "sensor")} | {"limit": scene["limit"]}
        print(f"[prepare] {scene['id']}", flush=True)
        try:
            if scene["sensor"] == "pleiades_neo":
                source = scene["source"]
                record["input"] = {
                    "archive": source["archive"],
                    "sha256": file_digest(source["archive"]),
                }
                bundle = extract_pneo(source, ROOT / config["pneo_extract_root"])
            else:
                bundle = scene["input"]
                record["input"] = {"archive": str(bundle), "sha256": file_digest(bundle)}
            reused = output_is_complete(output)
            correct_scene(
                bundle,
                output,
                scene["limit"],
                Path(config["path"]),
                template=settings["common"],
                extra=overrides[scene["sensor"]],
                python_executable=config["python"],
            )
            record["status"] = "reused" if reused else "ok"
        except Exception as error:  # recorded per scene; the stage continues
            record["status"] = "failed"
            record["error"] = f"{type(error).__name__}: {error}"[:2000]
        record["elapsed_s"] = round(time.monotonic() - started, 1)
        record["output_dir"] = str(output)
        record["outputs"] = (
            {p.name: file_digest(p) for p in sorted(output.glob("*.tif"))}
            if output.is_dir()
            else {}
        )
        records.append(record)
        write_json(directory / "progress.json", {"scenes": records})
        print(f"[prepare] {scene['id']}: {record['status']} ({record['elapsed_s']} s)", flush=True)
    complete_run(
        OUTPUT / "prepare",
        directory,
        {
            "stage": "prepare",
            "created_at": datetime.now(UTC).isoformat(),
            "acolite": revision,
            "settings": {
                name: {"path": str(p), "sha256": file_digest(p)} for name, p in settings.items()
            },
            **code,
            "scenes": records,
        },
    )
    return directory


def acolite_defaults(acolite_path: Path, sensor_key: str) -> dict[str, str]:
    """Generic ACOLITE defaults overlaid with the sensor's own defaults file."""
    values: dict[str, str] = {}
    for path in (
        acolite_path / "config/defaults.txt",
        acolite_path / f"config/defaults/{sensor_key}.txt",
    ):
        for line in path.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                key, _, value = line.partition("=")
                values[key.strip()] = value.strip()
    return values


def acolite_rasters(directory: Path, kind: str) -> dict[float, Path]:
    """Wavelength -> path for ACOLITE per-band L2R rasters of one kind (rhos or rhot)."""
    return {float(p.stem.rsplit("_", 1)[-1]): p for p in directory.glob(f"*_L2R_{kind}_[0-9]*.tif")}


def nearest_band(bands: dict[float, Path], target_nm: float, tolerance_nm: float | None = 15.0):
    if not bands:
        return None, None
    wavelength = min(bands, key=lambda w: abs(w - target_nm))
    if tolerance_nm is not None and abs(wavelength - target_nm) > tolerance_nm:
        return None, None
    return wavelength, bands[wavelength]


def qc_metrics(rgb, mask_rhot, threshold, water, cells, glint=None, nir=None) -> dict:
    """Registered scene checks on arrays sharing one grid; `cells` are (row, col) pairs."""
    finite = np.all([np.isfinite(band) for band in rgb], axis=0)
    with np.errstate(invalid="ignore"):
        clear = np.isfinite(mask_rhot) & (mask_rhot <= threshold)
    valid = water & finite & clear
    unique = sorted(set(cells))
    inside = [(r, c) for r, c in unique if 0 <= r < valid.shape[0] and 0 <= c < valid.shape[1]]
    metrics = {
        "water_pixels": int(water.sum()),
        "valid_water_pixels": int(valid.sum()),
        "valid_water_fraction": float(valid.sum() / water.sum()) if water.any() else 0.0,
        "quadrat_cells": len(unique),
        "valid_quadrat_cells": int(sum(valid[r, c] for r, c in inside)),
    }
    metrics["quadrat_cell_fraction"] = (
        metrics["valid_quadrat_cells"] / len(unique) if unique else 0.0
    )
    if glint is not None:
        values = glint[water & np.isfinite(glint)]
        metrics["median_rhos_glint_band"] = float(np.median(values)) if values.size else None
    if nir is not None:
        over_water = nir[valid]
        at_cells = [nir[r, c] for r, c in inside if valid[r, c]]
        metrics["median_rhos_nir_valid_water"] = (
            float(np.median(over_water)) if over_water.size else None
        )
        metrics["median_rhos_nir_quadrat_cells"] = float(np.median(at_cells)) if at_cells else None
    return metrics


def qc_sensitivity(valid_swir, valid_nir, cells, relaxed_fraction=0.75, registered_fraction=0.8):
    """Labelled alternatives to the registered quadrat-cell check.

    Reported so the margin is visible; never part of the verdict. `cells` maps a
    (row, col) pixel to the number of quadrats inside it.
    """
    inside = {
        cell: n
        for cell, n in cells.items()
        if 0 <= cell[0] < valid_swir.shape[0] and 0 <= cell[1] < valid_swir.shape[1]
    }
    registered_valid = sum(1 for cell in inside if valid_swir[cell])
    nir_valid = sum(1 for cell in inside if valid_nir[cell]) if valid_nir is not None else None
    quadrats = sum(inside.values())
    quadrats_valid = sum(n for cell, n in inside.items() if valid_swir[cell])
    fraction = registered_valid / len(inside) if inside else 0.0
    result = {
        "note": "Labelled sensitivity; the registered verdict uses cells and the 1614 nm mask.",
        "registered_cell_fraction": round(fraction, 4),
        "would_pass_at_relaxed_cell_fraction": bool(fraction >= relaxed_fraction),
        "relaxed_cell_fraction": relaxed_fraction,
        "per_quadrat_fraction": round(quadrats_valid / quadrats, 4) if quadrats else 0.0,
        "would_pass_per_quadrat": bool(
            quadrats and quadrats_valid / quadrats >= registered_fraction
        ),
    }
    if nir_valid is not None:
        nir_fraction = nir_valid / len(inside) if inside else 0.0
        result |= {
            "nir_mask_cell_fraction": round(nir_fraction, 4),
            "would_pass_with_nir_mask": bool(nir_fraction >= registered_fraction),
        }
    return result


def qc_verdict(metrics: dict, sensor: str, checks: dict) -> tuple[bool, list[str]]:
    reasons = []
    if metrics["valid_water_fraction"] < checks["valid_water_min_fraction"]:
        reasons.append("valid_water_fraction_below_minimum")
    if metrics["quadrat_cell_fraction"] < checks["quadrat_cell_min_fraction"]:
        reasons.append("quadrat_cell_fraction_below_minimum")
    if sensor == "sentinel2":
        glint = metrics.get("median_rhos_glint_band")
        if glint is None:
            reasons.append("glint_band_missing")
        elif glint >= checks["sentinel2_max_median_rhos_1614"]:
            reasons.append("glint_above_maximum")
    return not reasons, reasons


def select_scene(candidates: list[dict]) -> dict:
    """First passing scene in registered rank order; ties wait until every member is processed.

    Each candidate carries name, max_offset_days, pneo_offset_s, processed, passed and
    valid_water_fraction. Valid-water coverage only breaks ties on the first two keys.
    """
    groups: dict = {}
    for candidate in candidates:
        key = (
            candidate["max_offset_days"],
            math.inf if candidate["pneo_offset_s"] is None else candidate["pneo_offset_s"],
        )
        groups.setdefault(key, []).append(candidate)
    for key in sorted(groups):
        group = groups[key]
        pending = [c["name"] for c in group if not c["processed"]]
        if pending:
            return {"status": "undetermined", "reason": "not processed", "waiting_for": pending}
        passing = [c for c in group if c["passed"]]
        if passing:
            best = sorted(passing, key=lambda c: (-c["valid_water_fraction"], c["name"]))[0]
            return {"status": "selected", "name": best["name"]}
    return {"status": "unavailable", "reason": "no candidate passed the registered checks"}


def scene_qc() -> Path:
    import geopandas as gpd
    import rasterio
    from kelpobserve.coastal.provenance import resolve_run
    from rasterio import features, windows
    from rasterio.warp import Resampling, reproject

    spec = load_spec()
    code = code_provenance()
    checks = spec["scene_qc"]
    prepared = json.loads((resolve_run(OUTPUT / "prepare") / "manifest.json").read_text())
    inventory_dir = resolve_run(OUTPUT / "inventory")
    events = gpd.read_file(inventory_dir / "events.gpkg", layer="events")
    habitat = gpd.read_file(ROOT / spec["sources"]["habitat_aoi"])
    acolite_path = Path(spec["acolite"]["path"])
    results = {}
    for record in prepared["scenes"]:
        result = {k: record[k] for k in ("id", "campaign", "sensor", "status")}
        results[record["id"]] = result
        if record["status"] not in ("ok", "reused"):
            result |= {"passed": False, "reasons": ["acolite_failed"]}
            continue
        directory = Path(record["output_dir"])
        rhos, rhot = acolite_rasters(directory, "rhos"), acolite_rasters(directory, "rhot")
        bands = checks["bands_nm"][record["sensor"]]
        chosen = {role: nearest_band(rhos, nm) for role, nm in bands.items()}
        if any(chosen[role][1] is None for role in ("blue", "green", "red")):
            result |= {"passed": False, "reasons": ["incomplete_rgb"]}
            continue
        sample = next(iter(rhos.values()))
        sensor_key = re.match(r"(.+?)_\d{4}_\d{2}_\d{2}_", sample.name).group(1)
        defaults = acolite_defaults(acolite_path, sensor_key)
        mask_wave = float(defaults["l2w_mask_wave"])
        threshold = float(defaults["l2w_mask_threshold"])
        mask_nm, mask_path = nearest_band(rhot, mask_wave, tolerance_nm=None)
        if mask_path is None:
            result |= {"passed": False, "reasons": ["mask_band_missing"]}
            continue
        with rasterio.open(chosen["blue"][1]) as src:
            crs = src.crs
            window = windows.from_bounds(*habitat.to_crs(crs).total_bounds, transform=src.transform)
            window = (
                window.round_offsets()
                .round_lengths()
                .intersection(windows.Window(0, 0, src.width, src.height))
            )
            transform = src.window_transform(window)
            shape = (int(window.height), int(window.width))

        def read(path, window=window):
            with rasterio.open(path) as src:
                data = src.read(1, window=window).astype("float64")
                if src.nodata is not None and not np.isnan(src.nodata):
                    data[data == src.nodata] = np.nan
                return data

        water = np.zeros(shape, bool)
        for source in (spec["sources"]["lidar"], spec["sources"]["dtm"]):
            depth = np.full(shape, np.nan, dtype="float32")
            with rasterio.open(source) as src:
                reproject(
                    rasterio.band(src, 1),
                    depth,
                    dst_transform=transform,
                    dst_crs=crs,
                    resampling=Resampling.nearest,
                    src_nodata=np.nan,
                    dst_nodata=np.nan,
                )
            water |= np.isfinite(depth)
        water &= features.rasterize(
            habitat.to_crs(crs).geometry, out_shape=shape, transform=transform
        ).astype(bool)
        quadrats = events[events.campaign == record["campaign"]].to_crs(crs)
        cells = [
            rasterio.transform.rowcol(transform, g.x, g.y)
            for g in quadrats.geometry
            if g is not None
        ]
        metrics = qc_metrics(
            [read(chosen[r][1]) for r in ("blue", "green", "red")],
            read(mask_path),
            threshold,
            water,
            [(int(r), int(c)) for r, c in cells],
            glint=read(chosen["glint"][1]) if chosen.get("glint", (None, None))[1] else None,
            nir=read(chosen["nir"][1]) if chosen["nir"][1] else None,
        )
        passed, reasons = qc_verdict(metrics, record["sensor"], checks)
        counts: dict = {}
        for r, c in cells:
            counts[(int(r), int(c))] = counts.get((int(r), int(c)), 0) + 1
        nir_rhos = read(chosen["nir"][1]) if chosen["nir"][1] else None
        finite_rgb = np.all(
            [np.isfinite(read(chosen[r][1])) for r in ("blue", "green", "red")], axis=0
        )
        with np.errstate(invalid="ignore"):
            valid_swir = water & finite_rgb & (read(mask_path) <= threshold)
            valid_nir = (
                water & finite_rgb & np.isfinite(nir_rhos) & (nir_rhos <= 0.0215)
                if nir_rhos is not None
                else None
            )
        result |= {
            "sensitivity": qc_sensitivity(valid_swir, valid_nir, counts),
            "passed": passed,
            "reasons": reasons,
            "metrics": metrics,
            "bands_used_nm": {role: w for role, (w, _) in chosen.items()},
            "cloud_mask": {
                "acolite_sensor": sensor_key,
                "rhot_nm": mask_nm,
                "threshold": threshold,
            },
        }

    selections = {}
    for campaign in spec["campaigns"]:
        path = OUTPUT / "scenes" / f"{campaign}-candidates.json"
        if not path.is_file():
            continue
        candidates = []
        for c in json.loads(path.read_text())["candidates"]:
            result = results.get(s2_scene_id(c["name"]))
            processed = result is not None and "passed" in result
            if c.get("screened_out") and not processed:
                results.setdefault(
                    c["name"],
                    {
                        "id": c["name"],
                        "campaign": campaign,
                        "sensor": "sentinel2",
                        "status": "screened_out",
                        "passed": False,
                        "reasons": ["screened_out_glint"],
                        "screen_b11_median": c.get("screen_b11_median"),
                    },
                )
                processed = True
            candidates.append(
                {
                    "name": c["name"],
                    "max_offset_days": c["max_offset_days"],
                    "pneo_offset_s": c["pneo_offset_s"],
                    "processed": processed,
                    "passed": bool(result is not None and result.get("passed")),
                    "valid_water_fraction": (result or {})
                    .get("metrics", {})
                    .get("valid_water_fraction", 0.0),
                }
            )
        selections[campaign] = select_scene(candidates)

    directory = new_run_directory(OUTPUT / "scene-qc")
    write_json(
        directory / "scene_qc.json",
        {"checks": checks, "scenes": results, "sentinel2_selection": selections},
    )
    complete_run(
        OUTPUT / "scene-qc",
        directory,
        {
            "stage": "scene-qc",
            "created_at": datetime.now(UTC).isoformat(),
            "prepare_run": prepared.get("created_at"),
            "inventory_run": str(inventory_dir),
            "outputs": {"scene_qc.json": file_digest(directory / "scene_qc.json")},
            **code,
        },
    )
    return directory


def scene_directories(spec: dict) -> dict[str, dict]:
    """Prepared scene records by id, from the completed prepare run."""
    from kelpobserve.coastal.provenance import resolve_run

    prepared = json.loads((resolve_run(OUTPUT / "prepare") / "manifest.json").read_text())
    return {record["id"]: record for record in prepared["scenes"]}


def common_grid(spec: dict, bounds, reference_scene: str):
    """The registered 10 m grid, cropped to `bounds`. Every S2 scene shares it."""
    import rasterio
    from rasterio import windows

    directory = Path(scene_directories(spec)[reference_scene]["output_dir"])
    sample = next(iter(acolite_rasters(directory, "rhos").values()))
    with rasterio.open(sample) as src:
        require(
            src.res[0] == spec["common_grid"]["pixel_m"], "Reference scene is not on the 10 m grid."
        )
        window = (
            windows.from_bounds(*bounds, transform=src.transform)
            .round_offsets()
            .round_lengths()
            .intersection(windows.Window(0, 0, src.width, src.height))
        )
        return src.window_transform(window), (int(window.height), int(window.width)), src.crs


def aggregate_to_grid(path, transform, shape, crs, min_valid_fraction: float):
    """Area-weighted mean of a finer raster on the common grid, with a valid-fraction gate."""
    import rasterio
    from rasterio.warp import Resampling, reproject

    with rasterio.open(path) as src:
        values = src.read(1).astype("float32")
        if src.nodata is not None and not np.isnan(src.nodata):
            values[values == src.nodata] = np.nan
        mean = np.full(shape, np.nan, dtype="float32")
        reproject(
            values,
            mean,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=crs,
            resampling=Resampling.average,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )
        fraction = np.zeros(shape, dtype="float32")
        reproject(
            np.isfinite(values).astype("float32"),
            fraction,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=crs,
            resampling=Resampling.average,
            src_nodata=None,
            dst_nodata=None,
        )
    mean[fraction < min_valid_fraction] = np.nan
    return mean, fraction


def optics_reference() -> Path:
    """QAA scene IOPs over the registered PNeo reference region, one fit per date."""
    import geopandas as gpd
    import rasterio
    from kelpobserve.coastal.optics.qaa import fit_scene_iops
    from kelpobserve.coastal.optics.water import rhos_to_rrs
    from kelpobserve.coastal.scene import Scene
    from rasterio import features

    spec = load_spec()
    code = code_provenance()
    rules = spec["reference_qc"]
    region = gpd.read_file(REGISTRATION / spec["pneo_reference_region"]["path"])
    results = {}
    for scene_id, record in scene_directories(spec).items():
        if record["sensor"] != "pleiades_neo" or record["status"] not in ("ok", "reused"):
            continue
        scene = Scene.from_acolite_dir(Path(record["output_dir"]), sensor="pleiades_neo")
        with rasterio.open(
            next(iter(acolite_rasters(Path(record["output_dir"]), "rhos").values()))
        ) as src:
            mask = features.rasterize(
                region.to_crs(src.crs).geometry, out_shape=src.shape, transform=src.transform
            ).astype(bool)
        rrs = rhos_to_rrs(scene.rhos) if scene.rhos is not None else scene.rrs_above
        valid = mask & np.all(np.isfinite(rrs), axis=0)
        iops = fit_scene_iops(rrs, scene.wavelengths_nm, valid, scene.solar_zenith_deg)
        reasons = []
        if iops.n_deep_pixels < rules["min_pixels"]:
            reasons.append("insufficient_reference_pixels")
        if rules["require_no_qa_flags"] and iops.qa_flags:
            reasons.append("qaa_plausibility_flags")
        results[scene_id] = {
            "scene": scene_id,
            "campaign": record["campaign"],
            "solar_zenith_deg": round(scene.solar_zenith_deg, 2),
            "reference_pixels": int(iops.n_deep_pixels),
            "reference_wavelength_nm": float(iops.ref_wavelength_nm),
            "wavelengths_nm": [float(w) for w in iops.wavelengths_nm],
            "a": [round(float(v), 5) for v in iops.a],
            "bb": [round(float(v), 5) for v in iops.bb],
            "kd": [round(float(v), 5) for v in iops.kd],
            "a_cdm_443": _number(round(float(iops.a_cdm_443), 5)),
            "bbp_555": _number(round(float(iops.bbp_555), 6)),
            "y": _number(round(float(iops.y), 4)),
            "qa_flags": list(iops.qa_flags),
            "corrected_branch_available": not reasons,
            "reasons": reasons,
        }
        print(
            f"  {scene_id}: pixels={iops.n_deep_pixels} flags={iops.qa_flags or 'none'}", flush=True
        )
    directory = new_run_directory(OUTPUT / "optics-reference")
    write_json(directory / "reference_iops.json", {"rules": rules, "scenes": results})
    complete_run(
        OUTPUT / "optics-reference",
        directory,
        {
            "stage": "optics-reference",
            "created_at": datetime.now(UTC).isoformat(),
            "region": str(REGISTRATION / spec["pneo_reference_region"]["path"]),
            "outputs": {"reference_iops.json": file_digest(directory / "reference_iops.json")},
            **code,
        },
    )
    return directory


def cross_sensor() -> Path:
    """Label-free radiometric comparison of PNeo aggregated to 10 m against Sentinel-2."""
    import geopandas as gpd
    import rasterio
    from rasterio import features
    from rasterio.warp import Resampling, reproject
    from scipy import ndimage

    spec = load_spec()
    code = code_provenance()
    scenes = scene_directories(spec)
    bands = spec["scene_qc"]["bands_nm"]
    pairs = [
        (
            "s2-20260903-S2B-R080",
            "pneo-20260903",
            "same datatake minute; this S2 scene failed the scene checks",
        ),
        (
            "s2-20260831-S2B-R037",
            "pneo-20260903",
            "the selected September S2 scene, three days earlier",
        ),
    ]
    habitat = gpd.read_file(ROOT / spec["sources"]["habitat_aoi"])
    results = []
    for s2_id, pneo_id, note in pairs:
        if s2_id not in scenes or pneo_id not in scenes:
            continue
        pneo_dir, s2_dir = Path(scenes[pneo_id]["output_dir"]), Path(scenes[s2_id]["output_dir"])
        pneo_rhos, s2_rhos = acolite_rasters(pneo_dir, "rhos"), acolite_rasters(s2_dir, "rhos")
        sample = next(iter(pneo_rhos.values()))
        with rasterio.open(sample) as src:
            bounds, pneo_crs = src.bounds, src.crs
        transform, shape, crs = common_grid(spec, bounds, s2_id)
        require(str(crs) == str(pneo_crs), "PNeo and Sentinel-2 are on different projections.")
        depth = np.full(shape, np.nan, dtype="float32")
        with rasterio.open(spec["sources"]["dtm"]) as src:
            reproject(
                rasterio.band(src, 1),
                depth,
                dst_transform=transform,
                dst_crs=crs,
                resampling=Resampling.nearest,
                src_nodata=np.nan,
                dst_nodata=np.nan,
            )
        lidar = np.full(shape, np.nan, dtype="float32")
        with rasterio.open(spec["sources"]["lidar"]) as src:
            reproject(
                rasterio.band(src, 1),
                lidar,
                dst_transform=transform,
                dst_crs=crs,
                resampling=Resampling.nearest,
                src_nodata=np.nan,
                dst_nodata=np.nan,
            )
        water = np.isfinite(depth) | np.isfinite(lidar)
        in_habitat = features.rasterize(
            habitat.to_crs(crs).geometry, out_shape=shape, transform=transform
        ).astype(bool)
        # Stable land: inland of the surf and adjacency zone, three cells from any water.
        land = ndimage.binary_erosion(~water, iterations=3)
        comparison = {"sentinel2": s2_id, "pneo": pneo_id, "note": note, "bands": []}
        for role in ("blue", "green", "red", "nir"):
            s2_nm, s2_path = nearest_band(s2_rhos, bands["sentinel2"][role])
            pneo_nm, pneo_path = nearest_band(pneo_rhos, bands["pleiades_neo"][role])
            if s2_path is None or pneo_path is None:
                continue
            aggregated, fraction = aggregate_to_grid(
                pneo_path, transform, shape, crs, spec["common_grid"]["min_valid_fraction"]
            )
            with rasterio.open(s2_path) as src:
                window = rasterio.windows.from_bounds(
                    *rasterio.transform.array_bounds(*shape, transform), transform=src.transform
                )
                s2_values = src.read(1, window=window, boundless=True, fill_value=np.nan).astype(
                    "float32"
                )
            row = {"role": role, "sentinel2_nm": s2_nm, "pneo_nm": pneo_nm}
            for label, mask in (("water_in_habitat", water & in_habitat), ("stable_land", land)):
                both = mask & np.isfinite(aggregated) & np.isfinite(s2_values)
                if both.sum() < 50:
                    row[label] = {"pixels": int(both.sum()), "note": "too few pixels"}
                    continue
                a, b = aggregated[both], s2_values[both]
                row[label] = {
                    "pixels": int(both.sum()),
                    "pneo_median": round(float(np.median(a)), 5),
                    "sentinel2_median": round(float(np.median(b)), 5),
                    "median_difference": round(float(np.median(a - b)), 5),
                    "robust_scatter": round(
                        float(np.median(np.abs((a - b) - np.median(a - b)))), 5
                    ),
                    "correlation": round(float(np.corrcoef(a, b)[0, 1]), 3),
                }
            comparison["bands"].append(row)
        results.append(comparison)
        print(f"  {s2_id} vs {pneo_id}: {len(comparison['bands'])} bands compared", flush=True)
    directory = new_run_directory(OUTPUT / "cross-sensor")
    write_json(
        directory / "cross_sensor.json",
        {
            "note": "Radiometric comparison only. No field target is read, and neither"
            " correction is tuned.",
            "grid": spec["common_grid"],
            "pairs": results,
        },
    )
    complete_run(
        OUTPUT / "cross-sensor",
        directory,
        {
            "stage": "cross-sensor",
            "created_at": datetime.now(UTC).isoformat(),
            "outputs": {"cross_sensor.json": file_digest(directory / "cross_sensor.json")},
            **code,
        },
    )
    return directory


def eligible_events(spec: dict):
    """Registered-eligible events with their positions, from the completed inventory."""
    import geopandas as gpd
    from kelpobserve.coastal.provenance import resolve_run

    directory = resolve_run(OUTPUT / "inventory")
    events = gpd.read_file(directory / "events.gpkg", layer="events")
    return events[events.eligible_registered.astype(bool)].copy(), directory


def acquisition_levels(spec: dict) -> dict[str, float]:
    path = REGISTRATION / spec["water_levels"]["acquisition_path"]
    table = pd.read_csv(path)
    return dict(zip(table.scene_id, table.water_level_m))


def pneo_scene_for(spec: dict, campaign: str) -> str:
    acquisition = spec["campaigns"][campaign]["pneo"]
    require(acquisition, f"No PNeo delivery registered for {campaign}.")
    return f"pneo-{spec['sources']['pneo'][acquisition]['acquired_utc'][:10].replace('-', '')}"


def build_units(spec: dict, scale: str, campaign: str) -> pd.DataFrame:
    """Model units for one campaign: 10 m cells, or quadrats for the native scale."""
    import rasterio
    from rasterio.warp import Resampling, reproject

    events, _ = eligible_events(spec)
    events = events[events.campaign == campaign]
    require(len(events) > 0, f"No eligible {campaign} events.")
    scene_id = pneo_scene_for(spec, campaign)
    record = scene_directories(spec)[scene_id]
    rasters = acolite_rasters(Path(record["output_dir"]), "rhos")
    bands = spec["models"]["bands_nm"]
    level = acquisition_levels(spec)[scene_id]
    rows = []
    if scale == "common_grid_10m":
        margin = 100.0
        with rasterio.open(next(iter(rasters.values()))) as src:
            crs = src.crs
        projected = events.to_crs(crs)
        bounds = (
            projected.total_bounds[0] - margin,
            projected.total_bounds[1] - margin,
            projected.total_bounds[2] + margin,
            projected.total_bounds[3] + margin,
        )
        transform, shape, crs = common_grid(spec, bounds, spec["selected_grid_scene"])
        stack = {}
        for role, nm in bands.items():
            _, path = nearest_band(rasters, nm)
            require(path is not None, f"{scene_id} has no {role} band near {nm} nm.")
            stack[role], _ = aggregate_to_grid(
                path, transform, shape, crs, spec["common_grid"]["min_valid_fraction"]
            )
        depth = np.full(shape, np.nan, dtype="float32")
        with rasterio.open(spec["sources"]["lidar"]) as src:
            reproject(
                rasterio.band(src, 1),
                depth,
                dst_transform=transform,
                dst_crs=crs,
                resampling=Resampling.average,
                src_nodata=np.nan,
                dst_nodata=np.nan,
            )
        cells: dict = {}
        for _, event in projected.iterrows():
            r, c = rasterio.transform.rowcol(transform, event.geometry.x, event.geometry.y)
            cells.setdefault((int(r), int(c)), []).append(event)
        for (r, c), members in sorted(cells.items()):
            strata = {m.stratum for m in members}
            transects = {int(m.transect) for m in members}
            row = {
                "unit_id": f"{campaign}-cell-{r}-{c}",
                "campaign": campaign,
                "row": r,
                "column": c,
                "transect": sorted(transects)[0],
                "mixed_transects": len(transects) > 1,
                "stratum": sorted(strata)[0],
                "mixed_strata": len(strata) > 1,
                "n_quadrats": len(members),
                "event_ids": json.dumps([m.event_id for m in members]),
                "depth": float(depth[r, c]) + level if np.isfinite(depth[r, c]) else math.nan,
            }
            for role in bands:
                row[role] = float(stack[role][r, c])
            rows.append(row)
    else:
        with rasterio.open(next(iter(rasters.values()))) as src:
            crs = src.crs
        projected = events.to_crs(crs)
        windows = {}
        for role, nm in bands.items():
            _, path = nearest_band(rasters, nm)
            require(path is not None, f"{scene_id} has no {role} band near {nm} nm.")
            windows[role] = rasterio.open(path)
        for _, event in projected.iterrows():
            row = {
                "unit_id": event.event_id,
                "campaign": campaign,
                "transect": int(event.transect),
                "mixed_transects": False,
                "stratum": event.stratum,
                "mixed_strata": False,
                "n_quadrats": 1,
                "event_ids": json.dumps([event.event_id]),
                "depth": float(event.lidar_depth_m) + level
                if np.isfinite(event.lidar_depth_m)
                else math.nan,
            }
            for role, src in windows.items():
                r, c = src.index(event.geometry.x, event.geometry.y)
                patch = src.read(
                    1, window=((r - 1, r + 2), (c - 1, c + 2)), boundless=True, fill_value=np.nan
                ).astype("float64")
                row[role] = float(patch.mean()) if np.isfinite(patch).all() else math.nan
            rows.append(row)
        for src in windows.values():
            src.close()
    units = pd.DataFrame(rows)
    units["scale"] = scale
    units["scene"] = scene_id
    units["water_level_m"] = level
    return units


def attach_targets(units: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    """Unit target: the mean cover of its member events. Missing cover is dropped, not zeroed."""
    usable = targets[
        targets.cover_status.isin(["measured", "absence_derived_zero", "absence_recorded_zero"])
    ]
    lookup = dict(zip(usable.event_id, usable.cover_pct))
    covers, counts = [], []
    for ids in units.event_ids:
        values = [lookup[e] for e in json.loads(ids) if e in lookup]
        covers.append(float(np.mean(values)) if values else math.nan)
        counts.append(len(values))
    return units.assign(cover_pct=covers, n_scored_quadrats=counts)


def usable_units(units: pd.DataFrame, predictors: list[str]) -> pd.DataFrame:
    columns = [*predictors, "cover_pct"]
    return units.dropna(subset=columns)


def fit_family(family: dict, train: pd.DataFrame, spec: dict):
    """Fit one registered family. Standardisation uses training data only."""
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge

    predictors = family["predictors"]
    y = train.cover_pct.to_numpy(dtype=float)
    if family["kind"] == "mean":
        return {"kind": "mean", "value": float(y.mean()), "n_train": int(len(y))}
    x = train[predictors].to_numpy(dtype=float)
    if family["kind"] == "ridge":
        centre, scale = x.mean(axis=0), x.std(axis=0)
        scale = np.where(scale > 1e-12, scale, 1.0)
        model = Ridge(alpha=spec["models"]["ridge"]["alpha"]).fit((x - centre) / scale, y)
        return {
            "kind": "ridge",
            "predictors": predictors,
            "centre": centre.tolist(),
            "scale": scale.tolist(),
            "coef": model.coef_.tolist(),
            "intercept": float(model.intercept_),
            "n_train": int(len(y)),
        }
    forest = RandomForestRegressor(**spec["models"]["forest"]).fit(x, y)
    return {"kind": "forest", "predictors": predictors, "model": forest, "n_train": int(len(y))}


def predict_family(fit: dict, frame: pd.DataFrame) -> np.ndarray:
    if fit["kind"] == "mean":
        return np.full(len(frame), fit["value"], dtype=float)
    x = frame[fit["predictors"]].to_numpy(dtype=float)
    if fit["kind"] == "ridge":
        z = (x - np.array(fit["centre"])) / np.array(fit["scale"])
        return z @ np.array(fit["coef"]) + fit["intercept"]
    return np.asarray(fit["model"].predict(x), dtype=float)


def score(predicted: np.ndarray, observed: np.ndarray) -> dict:
    """Errors in percentage points, raw and clipped to the 0-100 range."""
    clipped = np.clip(predicted, 0.0, 100.0)
    error = predicted - observed
    return {
        "n": int(len(observed)),
        "mae_pp": round(float(np.mean(np.abs(error))), 2),
        "rmse_pp": round(float(np.sqrt(np.mean(error**2))), 2),
        "bias_pp": round(float(np.mean(error)), 2),
        "mae_clipped_pp": round(float(np.mean(np.abs(clipped - observed))), 2),
        "clipping_frequency": round(float(np.mean(predicted != clipped)), 3),
    }


def library_root() -> Path:
    """Where the coastal library lives."""
    import kelpobserve.coastal as coastal

    return Path(coastal.__file__).resolve().parent


def targets_path(inventory_dir: Path, role: str, stage: str) -> Path:
    """Target file for a stage. Fitting cannot reach the held-out campaign's targets."""
    require(role in ("development", "holdout"), f"Unknown target role {role!r}.")
    if stage == "fit" and role == "holdout":
        raise ValueError(
            "The fit stage cannot open the held-out targets. Evaluation opens them, "
            "after the fitted parameters are sealed."
        )
    name = "targets_development.csv" if role == "development" else "targets_september.csv"
    return inventory_dir / name


def verify_lock(path: Path) -> dict:
    """Every frozen input must still hash to its registered value."""
    lock = json.loads(path.read_text())
    for name, expected in lock["files"].items():
        candidate = Path(name)
        require(candidate.is_file(), f"Frozen input is missing: {name}")
        require(file_digest(candidate) == expected, f"Frozen input changed: {name}")
    return {
        "path": str(path),
        "sha256": file_digest(path),
        "frozen_at": lock["frozen_at"],
        "n_files": len(lock["files"]),
    }


def freeze() -> Path:
    """Pin every input, output and source file the fit will depend on."""
    from kelpobserve.coastal.provenance import resolve_run

    spec = load_spec()
    code = code_provenance()
    lock_path = REGISTRATION / "input-lock.json"
    require(not lock_path.exists(), "Input lock exists; preserve the registered experiment.")
    inventory_dir = resolve_run(OUTPUT / "inventory")
    paths = set(REGISTRATION.rglob("*.json")) | set(REGISTRATION.rglob("*.csv"))
    paths |= set(REGISTRATION.rglob("*.txt")) | set(REGISTRATION.rglob("*.md"))
    paths |= set(library_root().rglob("*.py"))
    paths.add(Path(__file__))
    paths |= {
        inventory_dir / name
        for name in (
            "events.gpkg",
            "targets_development.csv",
            "targets_september.csv",
            "reconciliation.json",
            "revisions.json",
            "manifest.json",
        )
    }
    for key in ("kelp_workbook", "previous_workbook", "lidar", "dtm"):
        paths.add(Path(spec["sources"][key]))
    paths.add(Path(spec["sources"]["transects"]["path"]))
    paths.add(Path(spec["sources"]["previous_transects"]["path"]))
    for acquisition in spec["sources"]["pneo"].values():
        paths.add(Path(acquisition["archive"]))
    for campaign in spec["campaigns"]:
        if spec["campaigns"][campaign].get("pneo"):
            directory = Path(scene_directories(spec)[pneo_scene_for(spec, campaign)]["output_dir"])
            paths |= set(directory.glob("*.tif")) | set(directory.glob("acolite_settings.txt"))
    paths = {p for p in paths if p.is_file() and p != lock_path}
    write_json(
        lock_path,
        {
            "schema_version": "1.0",
            "frozen_at": datetime.now(UTC).isoformat(),
            "inventory_run": str(inventory_dir),
            "files": {str(p.resolve()): file_digest(p) for p in sorted(paths)},
            **code,
        },
    )
    return lock_path


def fit() -> Path:
    """Leave-one-transect-out on the development campaign, then seal the fits."""
    import joblib
    from kelpobserve.coastal.provenance import resolve_run

    spec = load_spec()
    code = code_provenance()
    lock = verify_lock(REGISTRATION / "input-lock.json")
    design = spec["evaluation"]
    inventory_dir = resolve_run(OUTPUT / "inventory")
    targets = pd.read_csv(targets_path(inventory_dir, "development", stage="fit"))
    directory = new_run_directory(OUTPUT / "fit")
    summary = {}
    for scale in spec["models"]["scales"]:
        units = attach_targets(build_units(spec, scale, design["development_campaign"]), targets)
        units.to_csv(directory / f"units-{scale}.csv", index=False)
        folds, predictions = [], []
        for transect in sorted(units.transect.unique()):
            held = units[units.transect == transect]
            rest = units[units.transect != transect]
            record = {"held_out_transect": int(transect), "held_out_units": int(len(held))}
            per_stratum = rest.groupby("stratum").size().to_dict()
            enough = (
                len(usable_units(rest, ["blue", "green", "red", "depth"]))
                >= design["min_training_units"]
                and all(v >= design["min_training_units_per_stratum"] for v in per_stratum.values())
                and len(per_stratum) > 1
            )
            record["training_units_by_stratum"] = {k: int(v) for k, v in per_stratum.items()}
            record["sufficient"] = bool(enough)
            if not enough:
                record["reason"] = "insufficient training support"
                folds.append(record)
                continue
            for family in spec["models"]["families"]:
                train = usable_units(rest, family["predictors"])
                test = usable_units(held, family["predictors"])
                if len(test) == 0 or len(train) < design["min_training_units"]:
                    continue
                fitted = fit_family(family, train, spec)
                for unit_id, observed, predicted in zip(
                    test.unit_id, test.cover_pct, predict_family(fitted, test)
                ):
                    predictions.append(
                        {
                            "family": family["id"],
                            "unit_id": unit_id,
                            "transect": int(transect),
                            "stratum": test.set_index("unit_id").stratum.get(unit_id),
                            "observed": float(observed),
                            "predicted": float(predicted),
                        }
                    )
            folds.append(record)
        frame = pd.DataFrame(predictions)
        frame.to_csv(directory / f"loto-predictions-{scale}.csv", index=False)
        scores = {}
        for family in spec["models"]["families"]:
            rows = frame[frame.family == family["id"]]
            if rows.empty:
                scores[family["id"]] = {"n": 0, "note": "no fold produced predictions"}
                continue
            entry = score(rows.predicted.to_numpy(), rows.observed.to_numpy())
            entry["per_transect_mae_pp"] = {
                str(int(t)): round(float(np.mean(np.abs(g.predicted - g.observed))), 2)
                for t, g in rows.groupby("transect")
            }
            scores[family["id"]] = entry
        baseline = frame[frame.family == "depth_ridge"].set_index("unit_id")
        for family in spec["models"]["families"]:
            rows = frame[frame.family == family["id"]].set_index("unit_id")
            shared = rows.index.intersection(baseline.index)
            if family["id"] == "depth_ridge" or shared.empty:
                continue
            difference = np.abs(rows.loc[shared].predicted - rows.loc[shared].observed) - np.abs(
                baseline.loc[shared].predicted - baseline.loc[shared].observed
            )
            scores[family["id"]]["paired_vs_depth_only"] = {
                "n": int(len(shared)),
                "mean_absolute_error_difference_pp": round(float(difference.mean()), 2),
                "median_absolute_error_difference_pp": round(float(difference.median()), 2),
                "note": "negative means better than depth-only",
            }
        sealed = {}
        for family in spec["models"]["families"]:
            train = usable_units(units, family["predictors"])
            if len(train) < design["min_training_units"]:
                sealed[family["id"]] = {"unavailable": "insufficient training units"}
                continue
            sealed[family["id"]] = fit_family(family, train, spec)
        joblib.dump(sealed, directory / f"fits-{scale}.joblib")
        summary[scale] = {
            "units": int(len(units)),
            "units_with_target": int(units.cover_pct.notna().sum()),
            "quadrats": int(units.n_quadrats.sum()),
            "folds": folds,
            "loto_scores": scores,
            "sealed_fits_sha256": file_digest(directory / f"fits-{scale}.joblib"),
        }
    write_json(
        directory / "fit.json", {"development": design["development_campaign"], "scales": summary}
    )
    complete_run(
        OUTPUT / "fit",
        directory,
        {
            "stage": "fit",
            "created_at": datetime.now(UTC).isoformat(),
            "lock": lock,
            "targets_used": str(targets_path(inventory_dir, "development", stage="fit")),
            "outputs": {"fit.json": file_digest(directory / "fit.json")},
            "sealed": {scale: summary[scale]["sealed_fits_sha256"] for scale in summary},
            **code,
        },
    )
    return directory


def evaluate() -> Path:
    """Apply the sealed fits to the held-out campaign. One shot, never retuned."""
    import joblib
    from kelpobserve.coastal.provenance import resolve_run

    spec = load_spec()
    code = code_provenance()
    lock = verify_lock(REGISTRATION / "input-lock.json")
    design = spec["evaluation"]
    require(
        not (OUTPUT / "evaluate" / "latest.json").exists(),
        "An evaluation already exists; preserve the first test of the held-out campaign.",
    )
    fit_dir = resolve_run(OUTPUT / "fit")
    manifest = json.loads((fit_dir / "manifest.json").read_text())
    require(manifest.get("complete"), "The fit run is incomplete.")
    inventory_dir = resolve_run(OUTPUT / "inventory")
    targets = pd.read_csv(targets_path(inventory_dir, "holdout", stage="evaluate"))
    directory = new_run_directory(OUTPUT / "evaluate")
    summary = {}
    for scale in spec["models"]["scales"]:
        sealed_path = fit_dir / f"fits-{scale}.joblib"
        require(
            file_digest(sealed_path) == manifest["sealed"][scale],
            f"Sealed fits changed for {scale}.",
        )
        sealed = joblib.load(sealed_path)
        units = attach_targets(build_units(spec, scale, design["holdout_campaign"]), targets)
        units.to_csv(directory / f"units-{scale}.csv", index=False)
        rows, scores = [], {}
        for family in spec["models"]["families"]:
            fitted = sealed.get(family["id"], {})
            if "unavailable" in fitted:
                scores[family["id"]] = {"unavailable": fitted["unavailable"]}
                continue
            test = usable_units(units, family["predictors"])
            if test.empty:
                scores[family["id"]] = {"n": 0, "note": "no usable held-out units"}
                continue
            predicted = predict_family(fitted, test)
            entry = score(predicted, test.cover_pct.to_numpy(dtype=float))
            entry["per_transect_mae_pp"] = {
                str(int(t)): round(float(np.mean(np.abs(p - o))), 2)
                for t, p, o in [
                    (
                        t,
                        predicted[test.transect.to_numpy() == t],
                        test.cover_pct.to_numpy()[test.transect.to_numpy() == t],
                    )
                    for t in sorted(test.transect.unique())
                ]
            }
            entry["per_stratum_mae_pp"] = {
                str(st): round(
                    float(
                        np.mean(
                            np.abs(
                                predicted[test.stratum.to_numpy() == st]
                                - test.cover_pct.to_numpy()[test.stratum.to_numpy() == st]
                            )
                        )
                    ),
                    2,
                )
                for st in sorted(test.stratum.unique())
            }
            scores[family["id"]] = entry
            for unit_id, observed, value in zip(test.unit_id, test.cover_pct, predicted):
                rows.append(
                    {
                        "family": family["id"],
                        "unit_id": unit_id,
                        "observed": float(observed),
                        "predicted": float(value),
                    }
                )
        frame = pd.DataFrame(rows)
        frame.to_csv(directory / f"predictions-{scale}.csv", index=False)
        if not frame.empty and "depth_ridge" in set(frame.family):
            baseline = frame[frame.family == "depth_ridge"].set_index("unit_id")
            for family in spec["models"]["families"]:
                block = frame[frame.family == family["id"]].set_index("unit_id")
                shared = block.index.intersection(baseline.index)
                if family["id"] == "depth_ridge" or shared.empty or "n" not in scores[family["id"]]:
                    continue
                difference = np.abs(
                    block.loc[shared].predicted - block.loc[shared].observed
                ) - np.abs(baseline.loc[shared].predicted - baseline.loc[shared].observed)
                scores[family["id"]]["paired_vs_depth_only"] = {
                    "n": int(len(shared)),
                    "mean_absolute_error_difference_pp": round(float(difference.mean()), 2),
                    "note": "negative means better than depth-only",
                }
        summary[scale] = {
            "units": int(len(units)),
            "units_with_target": int(units.cover_pct.notna().sum()),
            "units_without_target": int(units.cover_pct.isna().sum()),
            "quadrats": int(units.n_quadrats.sum()),
            "scores": scores,
        }
    write_json(
        directory / "evaluation.json",
        {
            "holdout": design["holdout_campaign"],
            "no_selection_on_holdout": design["no_selection_on_holdout"],
            "scales": summary,
        },
    )
    complete_run(
        OUTPUT / "evaluate",
        directory,
        {
            "stage": "evaluate",
            "created_at": datetime.now(UTC).isoformat(),
            "lock": lock,
            "fit_run": str(fit_dir),
            "outputs": {"evaluation.json": file_digest(directory / "evaluation.json")},
            **code,
        },
    )
    return directory


def export_scene_selection() -> Path:
    """Portable scene-selection evidence for docs/coastal, with local paths dropped."""
    from kelpobserve.coastal.provenance import resolve_run

    spec = load_spec()
    code = code_provenance()
    qc_dir = resolve_run(OUTPUT / "scene-qc")
    report = json.loads((qc_dir / "scene_qc.json").read_text())
    screen = json.loads((REGISTRATION / "glint-screen.json").read_text())
    prepared = json.loads((resolve_run(OUTPUT / "prepare") / "manifest.json").read_text())
    scenes = []
    for scene in report["scenes"].values():
        scenes.append({k: v for k, v in scene.items() if k not in ("output_dir", "input")})
    results_path = ROOT / "docs/results.json"
    fit_pointer, evaluate_pointer = (
        OUTPUT / "fit" / "latest.json",
        OUTPUT / "evaluate" / "latest.json",
    )
    if fit_pointer.is_file() and evaluate_pointer.is_file():
        fit_report = json.loads((resolve_run(OUTPUT / "fit") / "fit.json").read_text())
        holdout = json.loads((resolve_run(OUTPUT / "evaluate") / "evaluation.json").read_text())
        write_json(
            results_path,
            {
                "experiment_id": spec["experiment_id"],
                "exported_at": datetime.now(UTC).isoformat(),
                "models": spec["models"],
                "evaluation_design": spec["evaluation"],
                "unavailable_branches": {
                    "sentinel2": "June and July unavailable; see scene-selection.md",
                    "water_column_correction": "Reference QC fails on both dates; see"
                    " optics-readiness.md",
                },
                "development": fit_report,
                "holdout": holdout,
                **code,
            },
        )
    path = ROOT / "docs/scene-selection.json"
    write_json(
        path,
        {
            "experiment_id": spec["experiment_id"],
            "exported_at": datetime.now(UTC).isoformat(),
            "amendment": spec["amendments"][0],
            "checks": report["checks"],
            "acolite": prepared["acolite"],
            "screen": {
                "definition": screen["definition"],
                "screen_glint_max": screen["screen_glint_max"],
                "scenes": screen["scenes"],
            },
            "scenes": scenes,
            "sentinel2_selection": report["sentinel2_selection"],
            **code,
        },
    )
    return path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("inventory", help="Build the versioned reference dataset.")
    commands.add_parser("reference-region", help="Write the geometry-only PNeo reference region.")
    fetcher = commands.add_parser("fetch", help="Rank and download Sentinel-2 L1C candidates.")
    fetcher.add_argument(
        "--env-file", type=Path, help="File holding CDSE_USERNAME and CDSE_PASSWORD."
    )
    fetcher.add_argument("--campaign", action="append", help="Limit to a campaign; repeatable.")
    fetcher.add_argument("--top", type=int, default=1, help="Download this many top-ranked scenes.")
    fetcher.add_argument(
        "--list-only", action="store_true", help="Rank candidates without downloading."
    )
    preparer = commands.add_parser("prepare", help="Run the pinned ACOLITE on every scene.")
    preparer.add_argument("--scene", action="append", help="Limit to a scene id; repeatable.")
    commands.add_parser("scene-qc", help="Apply the registered scene checks and ranking.")
    commands.add_parser("glint-screen", help="Registered glint pre-filter over candidate dates.")
    commands.add_parser("export", help="Write portable scene-selection evidence to docs/coastal.")
    levels = commands.add_parser("water-levels", help="Modelled sea level above the datum proxy.")
    levels.add_argument("--env-file", type=Path, help="File holding Copernicus Marine credentials.")
    commands.add_parser(
        "optics-reference", help="QAA fit over the registered PNeo reference region."
    )
    commands.add_parser("cross-sensor", help="Compare PNeo aggregated to 10 m against Sentinel-2.")
    commands.add_parser("freeze", help="Pin every registered input before fitting.")
    commands.add_parser("fit", help="Leave-one-transect-out on July, then seal the fits.")
    commands.add_parser("evaluate", help="Apply the sealed fits to September, once.")
    args = parser.parse_args()
    if args.command == "freeze":
        print(f"Frozen {len(json.loads(freeze().read_text())['files'])} files")
    elif args.command == "fit":
        directory = fit()
        print(json.dumps(json.loads((directory / "fit.json").read_text()), indent=2)[:4000])
    elif args.command == "evaluate":
        directory = evaluate()
        print(json.dumps(json.loads((directory / "evaluation.json").read_text()), indent=2)[:4000])
    elif args.command == "optics-reference":
        directory = optics_reference()
        report = json.loads((directory / "reference_iops.json").read_text())
        for scene in report["scenes"].values():
            available, why = scene["corrected_branch_available"], scene["reasons"]
            print(f"{scene['scene']}: corrected branch available = {available} {why}")
    elif args.command == "cross-sensor":
        directory = cross_sensor()
        print(f"Cross-sensor comparison written to {directory}")
    elif args.command == "water-levels":
        survey, acquisition = water_levels(args.env_file)
        print(pd.read_csv(survey).to_string(index=False))
        print(pd.read_csv(acquisition).to_string(index=False))
    elif args.command == "export":
        print(f"Scene-selection evidence written to {export_scene_selection()}")
    elif args.command == "glint-screen":
        path = glint_screen()
        rows = json.loads(path.read_text())["scenes"]
        kept = [r for r in rows if not r["screened_out"]]
        print(f"{len(rows)} candidate dates screened; {len(kept)} within the glint limit")
        for r in kept:
            print(
                f"  keep {r['sensing_utc'][:16]} R{r['relative_orbit']:03d} B11={r['b11_median']}"
            )
    elif args.command == "scene-qc":
        directory = scene_qc()
        report = json.loads((directory / "scene_qc.json").read_text())
        for scene in report["scenes"].values():
            print(f"{scene['id']}: passed={scene['passed']} {scene['reasons']}")
        print(json.dumps(report["sentinel2_selection"], indent=2))
    elif args.command == "prepare":
        directory = prepare(args.scene)
        for scene in json.loads((directory / "manifest.json").read_text())["scenes"]:
            print(f"{scene['id']}: {scene['status']} {scene.get('error', '')}")
    elif args.command == "inventory":
        directory = inventory()
        print(
            json.dumps(json.loads((directory / "manifest.json").read_text())["summary"], indent=2)
        )
        print(f"Inventory written to {directory}")
    elif args.command == "reference-region":
        path = reference_region()
        meta = json.loads(path.read_text())["metadata"]
        print(f"PNeo reference region: {meta['total_area_m2'] / 1e6:.3f} km2 -> {path}")
    else:
        fetch(args.env_file, args.campaign, args.top, not args.list_only)


if __name__ == "__main__":
    main()
