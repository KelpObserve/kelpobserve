"""Habitat context for the Sesimbra site, from published maps and this campaign's field data.

Presentation and context only: nothing here feeds the registered experiment or any
re-analysis. Sources:

- EMODnet Seabed Habitats, survey map PT010000: seabed habitats and associated
  communities of the Luiz Saldanha Marine Park (EUNIS 2012, determined 2013).
- EMODnet Seabed Habitats, EUSeaMap 2025: the Europe-wide broad-scale predictive map
  (EUNIS 2019).
- CCMAR's GPS mapping tracks and quadrats from the 2026 campaigns, and the EMODnet lidar.

    PYTHONPATH=. python -m benchmark.habitat_map fetch      # download and record provenance
    PYTHONPATH=. python -m benchmark.habitat_map summary    # the numbers the reports quote
"""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from kelpobserve.coastal.provenance import file_digest

from benchmark import run_benchmark as rb

CACHE = rb.ROOT / ".cache" / "emodnet-habitats"
WFS = "https://ows.emodnet-seabedhabitats.eu/geoserver/{workspace}/wfs"
# Lon/lat box for the Europe-wide map: the Setúbal peninsula and its shelf.
REGION = (-9.35, 38.30, -8.85, 38.60)
LAYERS = {
    "park": {"workspace": "emodnet_open_maplibrary", "name": "pt010000", "bbox": None},
    "euseamap": {"workspace": "emodnet_open", "name": "eusm2025_eunis2019_full", "bbox": REGION},
}
SACCORHIZA = "A3.31_PT16"
MAPPING_TRACKS = "Maps kelp AQUASERV Jun-Sep 2026.gpkg"

# Broad groups of the park map, by EUNIS 2012 code prefix, in drawing and legend order.
GROUPS = {
    "saccorhiza": {"label": "Saccorhiza polyschides on infralittoral rock", "colour": "#d97706"},
    "infralittoral_rock": {"label": "Infralittoral rock, other communities", "colour": "#4d9a3a"},
    "littoral_rock": {"label": "Littoral rock", "colour": "#9ca3af"},
    "circalittoral_rock": {"label": "Circalittoral rock", "colour": "#7c5cbf"},
    "coarse": {"label": "Coarse sediment", "colour": "#e6b95c"},
    "sand": {"label": "Sand and muddy sand", "colour": "#f4e3a3"},
    "mud": {"label": "Mud and sandy mud", "colour": "#b9a08a"},
    "mixed": {"label": "Mixed sediment", "colour": "#d8c3e6"},
}
PREFIXES = {
    "A1": "littoral_rock",
    "A3": "infralittoral_rock",
    "A4": "circalittoral_rock",
    "A5.1": "coarse",
    "A5.2": "sand",
    "A5.3": "mud",
    "A5.4": "mixed",
}
# EUNIS 2012 names of the standard units the park map uses.
EUNIS_2012 = {
    "A3.1": "High energy infralittoral rock",
    "A3.2": "Moderate energy infralittoral rock",
    "A3.3": "Low energy infralittoral rock",
    "A4.1": "High energy circalittoral rock",
    "A4.2": "Moderate energy circalittoral rock",
    "A4.3": "Low energy circalittoral rock",
    "A5.13": "Infralittoral coarse sediment",
    "A5.14": "Circalittoral coarse sediment",
    "A5.15": "Deep circalittoral coarse sediment",
    "A5.23": "Infralittoral fine sand",
    "A5.25": "Circalittoral fine sand",
    "A5.26": "Circalittoral muddy sand",
    "A5.27": "Deep circalittoral sand",
    "A5.35": "Circalittoral sandy mud",
    "A5.36": "Circalittoral fine mud",
    "A5.37": "Deep circalittoral mud",
    "A5.44": "Circalittoral mixed sediments",
    "A5.45": "Deep circalittoral mixed sediments",
}
# Colours for the communities mapped on the reef at the dive site.
COMMUNITY_COLOURS = {
    "A3.31_PT16": "#d97706",
    "A3.15_PT12": "#ec7fa9",
    "A3.15_PT13": "#b5179e",
    "A3.14_PT11": "#f4a3b5",
    "A3.16_PT14_PT15": "#3b6fb6",
}


# ---------------------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------------------


def layer_url(layer: dict) -> str:
    url = (
        f"{WFS.format(workspace=layer['workspace'])}?service=WFS&version=2.0.0"
        f"&request=GetFeature&typeNames={layer['workspace']}:{layer['name']}"
        "&outputFormat=application/json&srsName=EPSG:4326"
    )
    if layer["bbox"]:
        url += "&bbox=" + ",".join(str(value) for value in layer["bbox"]) + ",EPSG:4326"
    return url


def layer_path(key: str) -> Path:
    return CACHE / f"{key}.geojson"


def fetch(force: bool = False) -> Path:
    """Download both EMODnet layers once, with where and when they came from."""
    CACHE.mkdir(parents=True, exist_ok=True)
    record_path = CACHE / "provenance.json"
    record = json.loads(record_path.read_text()) if record_path.exists() else {}
    for key, layer in LAYERS.items():
        path = layer_path(key)
        if path.exists() and key in record and not force:
            continue
        url = layer_url(layer)
        request = urllib.request.Request(url, headers={"User-Agent": "kelpobserve"})
        with urllib.request.urlopen(request, timeout=300) as response:
            path.write_bytes(response.read())
        record[key] = {
            "url": url,
            "fetched_at": datetime.now(UTC).isoformat(),
            "sha256": file_digest(path),
            "features": len(json.loads(path.read_text())["features"]),
        }
    record_path.write_text(json.dumps(record, indent=2))
    return CACHE


def read_layer(key: str):
    import geopandas as gpd

    fetch()
    return gpd.read_file(layer_path(key))


# ---------------------------------------------------------------------------------------
# The park map, classified
# ---------------------------------------------------------------------------------------


def group_of(code: str) -> str:
    if code.startswith(SACCORHIZA):
        return "saccorhiza"
    for prefix, group in PREFIXES.items():
        if code.startswith(prefix):
            return group
    raise ValueError(f"No habitat group for EUNIS code {code}")


def community_name(code: str, note: str | None) -> str:
    """The community a polygon is mapped as: the proposed unit's name, or the EUNIS name."""
    if isinstance(note, str) and note:
        match = re.search(
            r"(?:new EUNIS unit|replacement of EUNIS unit [A-Z0-9.]+)\s*:\s*(.+?)\. PMLS", note
        )
        if match:
            return match.group(1).replace("[", "").replace("]", "")
    return EUNIS_2012.get(code.rstrip("*"), code)


def park_map():
    """PT010000 with a group, a readable name and an area for every polygon."""
    frame = read_layer("park")
    frame = frame.assign(
        code=frame.hab_type.str.rstrip("*"),
        group=[group_of(code) for code in frame.hab_type],
        name=[community_name(code, note) for code, note in zip(frame.hab_type, frame.tran_com)],
        area_ha=frame.to_crs(32629).area / 1e4,
    )
    return frame[["code", "group", "name", "area_ha", "geometry"]]


def euseamap():
    """EUSeaMap 2025 over the region, one row per polygon with its EUNIS 2019 class."""
    frame = read_layer("euseamap")
    frame = frame.rename(columns={"eunis2019d": "name", "eunis2019c": "code"})
    return frame[["code", "name", "substrate", "biozone", "energy", "geometry"]]


def track_campaign(name: str) -> str:
    """A track is named after its month ("June 2026") or its recording time (dd/mm/yy)."""
    months = {6: "june", 7: "july", 9: "september"}
    stamp = re.match(r"\d{2}/(\d{2})/\d{2}", name)
    if stamp:
        return months[int(stamp.group(1))]
    return name.split()[0].lower()


def mapping_tracks(spec: dict):
    """CCMAR's GPS mapping tracks, one line per campaign."""
    import geopandas as gpd

    path = Path(spec["sources"]["transects"]["path"]).parent / MAPPING_TRACKS
    frame = gpd.read_file(path, layer="maps")
    frame = frame.assign(campaign=[track_campaign(name) for name in frame.Name])
    return frame[["campaign", "geometry"]]


# ---------------------------------------------------------------------------------------
# What the maps say at the quadrats
# ---------------------------------------------------------------------------------------


def lidar_depths(spec: dict, geometry) -> np.ndarray:
    """Lidar depths (chart datum) of the cells inside a lon/lat geometry."""
    import rasterio
    from rasterio import features

    with rasterio.open(spec["sources"]["lidar"]) as source:
        depth = source.read(1)
        inside = features.geometry_mask(
            [geometry.__geo_interface__], depth.shape, source.transform, invert=True
        )
    values = depth[inside]
    return values[np.isfinite(values)]


def summary(events: pd.DataFrame) -> dict:
    """The numbers the report and the story map quote, from the two maps and the quadrats."""
    import geopandas as gpd

    spec = rb.load_spec()
    park = park_map()
    broad = euseamap().to_crs(32629)
    quadrats = gpd.GeoDataFrame(
        events, geometry=gpd.points_from_xy(events.longitude, events.latitude), crs=4326
    )
    joined = gpd.sjoin(quadrats, park[["code", "name", "geometry"]], how="left", predicate="within")
    joined["name"] = joined["name"].fillna("outside the mapped polygons")
    counts = joined.groupby(["campaign", "stratum", "name"]).size().rename("quadrats").reset_index()
    kelp = park[park.group == "saccorhiza"]
    kelp_utm = kelp.to_crs(32629)
    points = quadrats.to_crs(32629).geometry
    to_kelp = np.array([kelp_utm.distance(point).min() for point in points])
    to_broad = np.array([broad.distance(point).min() for point in points])
    nearest = broad.iloc[int(broad.distance(points.iloc[0]).argmin())]
    rock = broad[broad.substrate.str.contains("Rock")]
    depths = np.concatenate([lidar_depths(spec, geometry) for geometry in kelp.geometry])
    areas = park.groupby("group").area_ha.sum()
    return rb_clean(
        {
            "park_map": {
                "polygons": int(len(park)),
                "area_ha": round(float(park.area_ha.sum())),
                "area_by_group_ha": {
                    group: round(float(areas.get(group, 0.0)), 1) for group in GROUPS
                },
                "saccorhiza_polygons": int(len(kelp)),
            },
            "saccorhiza_lidar_depth_m": {
                "cells": int(depths.size),
                "p10": round(float(np.percentile(depths, 10)), 1),
                "median": round(float(np.median(depths)), 1),
                "p90": round(float(np.percentile(depths, 90)), 1),
            },
            "quadrats": {
                "by_community": counts.to_dict("records"),
                "on_saccorhiza_class": int((joined.code == SACCORHIZA).sum()),
                "distance_to_saccorhiza_m": {
                    "min": round(float(to_kelp.min())),
                    "median": round(float(np.median(to_kelp))),
                    "max": round(float(to_kelp.max())),
                },
            },
            "euseamap": {
                "quadrats_inside_a_polygon": int((to_broad == 0).sum()),
                "distance_to_nearest_polygon_m": {
                    "min": round(float(to_broad.min())),
                    "median": round(float(np.median(to_broad))),
                    "max": round(float(to_broad.max())),
                },
                "nearest_class": nearest["name"],
                "distance_to_nearest_rock_m": round(
                    float(min(rock.distance(point).min() for point in points))
                ),
            },
        }
    )


def rb_clean(value):
    from benchmark import revision as rv

    return rv.clean(value)


# ---------------------------------------------------------------------------------------
# Layers for the story map
# ---------------------------------------------------------------------------------------


def as_features(frame, tolerance: float, digits: int = 5) -> dict:
    """A GeoDataFrame as a compact FeatureCollection, simplified in degrees."""
    from shapely import set_precision

    features = []
    for row in frame.itertuples(index=False):
        geometry = set_precision(row.geometry.simplify(tolerance), 10**-digits)
        if geometry.is_empty:
            continue
        properties = {
            key: (round(value, 2) if isinstance(value, float) else value)
            for key, value in row._asdict().items()
            if key != "geometry"
        }
        features.append(
            {"type": "Feature", "properties": properties, "geometry": geometry.__geo_interface__}
        )
    return {"type": "FeatureCollection", "features": features}


def site_layers(spec: dict) -> dict[str, dict]:
    """The park map, the broad-scale map and CCMAR's tracks, sized for the web."""
    park = park_map()
    park = park.assign(
        colour=[
            COMMUNITY_COLOURS.get(code, GROUPS[group]["colour"])
            for code, group in zip(park.code, park.group)
        ]
    )
    broad = euseamap()
    broad = broad.assign(geometry=broad.geometry.make_valid())
    broad = broad[broad.code != "Na"].dissolve(by="name", as_index=False)[
        ["name", "substrate", "biozone", "geometry"]
    ]
    return {
        "habitat-park": as_features(park, 0, digits=6),
        "habitat-labels": habitat_labels(park),
        "habitat-euseamap": as_features(broad, 0.0003, digits=4),
        "mapping-tracks": as_features(mapping_tracks(spec), 0.00001, digits=6),
    }


def habitat_labels(park) -> dict:
    """Place local reef labels inside the mapped communities, without refining boundaries."""
    from shapely.geometry import box
    from shapely.ops import unary_union

    window = box(-9.205, 38.407, -9.194, 38.4125)
    features = []
    for code in COMMUNITY_COLOURS:
        clipped = [g.intersection(window) for g in park.loc[park.code == code, "geometry"]]
        geometry = unary_union([g for g in clipped if not g.is_empty])
        parts = list(geometry.geoms) if hasattr(geometry, "geoms") else [geometry]
        polygons = [g for g in parts if g.geom_type == "Polygon" and not g.is_empty]
        if not polygons:
            continue
        point = max(polygons, key=lambda g: g.area).representative_point()
        features.append({
            "type": "Feature", "properties": {"code": code},
            "geometry": point.__geo_interface__,
        })
    return {"type": "FeatureCollection", "features": features}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    fetcher = commands.add_parser("fetch", help="Download the EMODnet layers.")
    fetcher.add_argument("--force", action="store_true", help="Download again.")
    commands.add_parser("summary", help="Print the numbers the reports quote.")
    arguments = parser.parse_args()
    if arguments.command == "fetch":
        print(f"EMODnet habitat layers in {fetch(arguments.force)}")
    else:
        from benchmark import report_maps

        print(json.dumps(summary(report_maps.events()), indent=2))


if __name__ == "__main__":
    main()
