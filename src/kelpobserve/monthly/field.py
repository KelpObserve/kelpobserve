"""CCMAR quadrat aggregation, coordinate joins and lidar sampling."""

from __future__ import annotations

import json
import math
import re
from datetime import date

import numpy as np
import pandas as pd

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

JUVENILE, ADULT, COVER = "Juvenile (< 30cm)", "Adult (>= 30cm)", "% cover"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _flag(value) -> bool:
    """True only for real booleans; outer joins leave NaN where a side is missing."""
    return isinstance(value, (bool, np.bool_)) and bool(value)


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
    for (day, depth, transect, quadrat), group in rows.groupby(["_date", "Depth", "Transect", "Quadrat"], sort=True):
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
            require(math.isnan(cover) or cover == 0.0, f"EMPTY record with nonzero cover: {source_rows}")
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
