"""Maps for the technical report, and map layers for the story map.

Presentation only: reads registered inputs and completed run outputs, fits nothing and
computes no result. Figures go to docs/figures/; layers go to a directory of your choice
(the KelpObserve site keeps them under site/public/sesimbra/).

    PYTHONPATH=. python -m benchmark.report_maps layers --site <directory>
    PYTHONPATH=. python -m benchmark.report_maps figures
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import shutil
import subprocess
import urllib.request
from html import escape
from pathlib import Path
from string import Template
from urllib.parse import unquote, urlsplit

import numpy as np
import pandas as pd
from kelpobserve.coastal.provenance import resolve_run

from benchmark import habitat_map as hm
from benchmark import revision as rv
from benchmark import revision_r5 as r5
from benchmark import run_benchmark as rb

FIGURES = rb.ROOT / "docs" / "figures"
REPORT = rb.ROOT / "docs" / "technical-report.md"
TILE_CACHE = rb.ROOT / ".cache" / "tiles"
USER_AGENT = "kelpobserve report maps"

BASEMAPS = {
    "oceanstream": {
        "tiles": "https://basemap.oceanstream.io/v5/xyz/{z}/{x}/{y}.webp",
        "max_zoom": 14,
        "credit": "Basemap: OceanStream (GEBCO, EMODnet, Copernicus DEM, OpenStreetMap)",
    },
    "esri": {
        "tiles": "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/"
        "MapServer/tile/{z}/{y}/{x}",
        "max_zoom": 19,
        "credit": "Imagery: Esri, Maxar, Earthstar Geographics and the GIS user community",
    },
    "emodnet": {
        "wms": "https://ows.emodnet-bathymetry.eu/wms",
        "layer": "emodnet:mean_atlas_land",
        "credit": "Bathymetry: EMODnet Bathymetry Consortium",
    },
}

# Sentinel-2 scenes with ACOLITE output, in date order. The three June dates added by r3
# exist only under the glint-corrected variant; top-of-atmosphere reflectance is the same
# in both.
S2_SCENES = [
    "s2-20260528-S2C-R037",
    "s2-20260530-S2A-R037",
    "s2-20260602-S2B-R037",
    "s2-20260602-S2A-R080",
    "s2-20260607-S2C-R037",
    "s2-20260609-S2A-R037",
    "s2-20260702-S2B-R037",
    "s2-20260705-S2B-R080",
    "s2-20260709-S2A-R037",
    "s2-20260831-S2B-R037",
    "s2-20260903-S2B-R080",
]
# The overlay window, EPSG:32629: Cabo Espichel to the Arrábida coast, 12 x 9 km.
OVERLAY_BOUNDS = (477000.0, 4246500.0, 489000.0, 4255500.0)
CONTOURS_LIDAR = (5.0, 10.0, 14.0, 20.0)
CONTOURS_DTM = (30.0, 50.0, 75.0, 100.0)
CAMPAIGNS = ("june", "july", "september")


# ---------------------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------------------


def transformer(source: int, target: int):
    from pyproj import Transformer

    return Transformer.from_crs(source, target, always_xy=True)


def lonlat(x, y) -> tuple:
    return transformer(32629, 4326).transform(x, y)


def mercator(lon, lat) -> tuple:
    return transformer(4326, 3857).transform(lon, lat)


def rounded(value, digits: int = 6):
    if isinstance(value, (list, tuple)):
        return [rounded(v, digits) for v in value]
    return round(float(value), digits)


def feature(geometry: dict, **properties) -> dict:
    return {"type": "Feature", "properties": rv.clean(properties), "geometry": geometry}


def collection(features: list[dict]) -> dict:
    return {"type": "FeatureCollection", "features": features}


def polygon_utm(left, bottom, right, top) -> dict:
    ring = [(left, bottom), (right, bottom), (right, top), (left, top), (left, bottom)]
    return {"type": "Polygon", "coordinates": [[rounded(lonlat(x, y)) for x, y in ring]]}


# ---------------------------------------------------------------------------------------
# Survey events and units
# ---------------------------------------------------------------------------------------


def events() -> pd.DataFrame:
    """All 90 survey events with their cover, from the inventory run."""
    directory = resolve_run(rb.OUTPUT / "inventory")
    frame = pd.read_csv(directory / "events.csv")
    targets = pd.concat(
        [
            pd.read_csv(directory / "targets_development.csv"),
            pd.read_csv(directory / "targets_september.csv"),
        ]
    )
    frame = frame.merge(targets[["event_id", "cover_pct", "cover_status"]], on="event_id")
    x, y = transformer(4326, 32629).transform(frame.longitude.values, frame.latitude.values)
    return frame.assign(x=x, y=y)


def quadrat_features(frame: pd.DataFrame) -> dict:
    features = []
    for row in frame.itertuples():
        features.append(
            feature(
                {"type": "Point", "coordinates": rounded([row.longitude, row.latitude], 7)},
                event_id=row.event_id,
                campaign=row.campaign,
                stratum=row.stratum,
                transect=int(row.transect),
                quadrat=int(row.quadrat),
                cover_pct=row.cover_pct,
                cover_status=row.cover_status,
                lidar_depth_m=None if pd.isna(row.lidar_depth_m) else round(row.lidar_depth_m, 1),
                eligible=bool(row.eligible_registered),
                depth_check=row.depth_qc_status,
            )
        )
    return collection(features)


def transect_features(frame: pd.DataFrame) -> dict:
    features = []
    for (campaign, stratum, transect), group in frame.groupby(["campaign", "stratum", "transect"]):
        group = group.sort_values("quadrat")
        line = [rounded([lon, lat], 7) for lon, lat in zip(group.longitude, group.latitude)]
        features.append(
            feature(
                {"type": "LineString", "coordinates": line},
                campaign=campaign,
                stratum=stratum,
                transect=int(transect),
                mean_cover_pct=round(float(group.cover_pct.mean()), 1),
            )
        )
    return collection(features)


def cell_features() -> dict:
    """July 10 m cells of the common grid, with which scenes kept each cell (r1 QC, r3 mask)."""
    spec = rb.load_spec()
    qc = resolve_run(rv.OUTPUT / "qc")
    pneo = pd.read_csv(qc / "units-pleiades_neo-july-grid.csv").set_index("unit_id")
    sentinel = pd.read_csv(qc / "units-sentinel2-july-grid.csv").set_index("unit_id")
    frame = events()
    july = frame[(frame.campaign == "july") & frame.eligible_registered]
    ninth = tolerant_cells(spec, july)
    features = []
    for unit_id, unit in pneo.iterrows():
        members = july[july.event_id.isin(json.loads(unit.event_ids))]
        left = math.floor(members.x.iloc[0] / 10) * 10
        bottom = math.floor(members.y.iloc[0] / 10) * 10
        if not ((members.x // 10 * 10 == left) & (members.y // 10 * 10 == bottom)).all():
            raise ValueError(f"{unit_id}: its quadrats do not share one 10 m cell")
        features.append(
            feature(
                polygon_utm(left, bottom, left + 10, bottom + 10),
                unit_id=unit_id,
                stratum=unit.stratum,
                quadrats=int(unit.n_quadrats),
                mean_cover_pct=round(float(unit.cover_pct), 1),
                depth_m=round(float(unit.depth), 1),
                pneo_8jul=unit.reason == "ok",
                s2_2jul=sentinel.loc[unit_id, "reason"] == "ok",
                s2_9jul=bool(ninth[(left, bottom)]),
            )
        )
    return collection(features)


def tolerant_cells(spec: dict, july: pd.DataFrame) -> dict:
    """Which July cells r3's glint-tolerant mask keeps on 9 July (17 of 18)."""
    import rasterio

    mask = json.loads((resolve_run(rb.OUTPUT / "r3" / "mask") / "mask.json").read_text())
    scene = mask["scenes"]["s2-20260709-S2A-R037"]
    directory = rv.ACOLITE_R1 / "V1" / "s2-20260709-S2A-R037"
    _, swir = rb.nearest_band(rb.acolite_rasters(directory, "rhot"), 1610.0)
    visible = [
        rb.nearest_band(rb.acolite_rasters(directory, "rhos"), nm)[1] for nm in (490, 560, 665)
    ]
    cells = {(math.floor(x / 10) * 10, math.floor(y / 10) * 10) for x, y in zip(july.x, july.y)}
    points = [(left + 5, bottom + 5) for left, bottom in cells]
    with rasterio.open(swir) as source:
        glint = np.array([v[0] for v in source.sample(points)])
    finite = np.ones(len(points), dtype=bool)
    for path in visible:
        with rasterio.open(path) as source:
            finite &= np.isfinite(np.array([v[0] for v in source.sample(points)]))
    kept = dict(zip(cells, finite & (glint <= scene["threshold"])))
    if sum(kept.values()) != scene["checks"]["valid_quadrat_cells"]:
        raise ValueError("9 July cells do not reproduce r3's count")
    return kept


# ---------------------------------------------------------------------------------------
# Areas, depth contours and depth bins
# ---------------------------------------------------------------------------------------


def pneo_footprint(spec: dict):
    """Where the 8 July Pléiades Neo image holds data, as a polygon in EPSG:32629."""
    import rasterio
    from rasterio import features as rio_features
    from shapely.geometry import Polygon, shape

    directory = rv.scene_dir(spec, "pneo-20260708")
    path = next(iter(rb.acolite_rasters(directory, "rhot").values()))
    step = 8
    with rasterio.open(path) as source:
        shape_out = (source.height // step, source.width // step)
        values = source.read(1, out_shape=shape_out)
        transform = source.transform * source.transform.scale(
            source.width / shape_out[1], source.height / shape_out[0]
        )
    filled = np.isfinite(values) & (values != 0)
    pieces = [
        shape(geometry)
        for geometry, value in rio_features.shapes(
            filled.astype("uint8"), mask=filled, transform=transform
        )
        if value == 1
    ]
    largest = max(pieces, key=lambda piece: piece.area)
    return Polygon(largest.exterior).simplify(15)


def area_features(spec: dict) -> dict:
    import geopandas as gpd

    features = []
    habitat = gpd.read_file(rb.ROOT / spec["sources"]["habitat_aoi"])
    features.append(
        feature(
            habitat.geometry.iloc[0].__geo_interface__,
            role="habitat_aoi",
            label="Habitat area of interest (provisional rectangle)",
        )
    )
    extent = gpd.read_file(rb.ROOT / spec["sources"]["processing_extent"])
    features.append(
        feature(
            extent.geometry.iloc[0].__geo_interface__,
            role="processing_extent",
            label="Sentinel-2 processing extent",
        )
    )
    region = gpd.read_file(rb.REGISTRATION / spec["pneo_reference_region"]["path"])
    for geometry in region.geometry:
        features.append(
            feature(
                geometry.__geo_interface__,
                role="offshore_reference",
                label="Offshore reference region, 51-110 m deep",
            )
        )
    footprint = pneo_footprint(spec)
    features.append(
        feature(
            gpd.GeoSeries([footprint], crs=32629).to_crs(4326).iloc[0].__geo_interface__,
            role="pneo_footprint",
            label="Pléiades Neo image footprint, 8 July",
            area_km2=round(footprint.area / 1e6, 1),
        )
    )
    screen = json.loads((rb.REGISTRATION / "glint-screen.json").read_text())
    features.append(
        feature(
            polygon_utm(*screen["bounds_epsg32629"]),
            role="screen_area",
            label="Glint screen area: the quadrats plus 250 m",
        )
    )
    features.append(
        feature(
            polygon_utm(*OVERLAY_BOUNDS),
            role="overlay_window",
            label="Sentinel-2 image window",
        )
    )
    return collection(features)


def contour_lines(path: Path, levels, source_name: str, clip=None) -> list[dict]:
    """Depth contours of a lon/lat raster as line features (depth positive down)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import rasterio
    from shapely.geometry import LineString, box

    with rasterio.open(path) as source:
        depth = source.read(1)
        transform = source.transform
    rows, columns = depth.shape
    longitude = transform.c + transform.a * (np.arange(columns) + 0.5)
    latitude = transform.f + transform.e * (np.arange(rows) + 0.5)
    figure, axis = plt.subplots()
    drawn = axis.contour(longitude, latitude, depth, levels=list(levels))
    features = []
    for level, segments in zip(drawn.levels, drawn.allsegs):
        for segment in segments:
            if len(segment) < 8:
                continue
            line = LineString(segment)
            if clip is not None:
                line = line.intersection(box(*clip))
            parts = getattr(line, "geoms", [line])
            for part in parts:
                if part.geom_type != "LineString" or len(part.coords) < 8:
                    continue
                features.append(
                    feature(
                        {"type": "LineString", "coordinates": rounded(list(part.coords), 5)},
                        depth_m=float(level),
                        source=source_name,
                    )
                )
    plt.close(figure)
    return features


def contour_features(spec: dict) -> dict:
    lidar = Path(spec["sources"]["lidar"])
    dtm = Path(spec["sources"]["dtm"])
    return collection(
        contour_lines(lidar, CONTOURS_LIDAR, "EMODnet HR lidar, 2011")
        + contour_lines(dtm, CONTOURS_DTM, "EMODnet DTM", clip=(-9.30, 38.34, -9.08, 38.47))
    )


def depth_bin_features(spec: dict) -> dict:
    """r5's habitat depth bins as polygons, from the lidar (chart datum, LAT)."""
    import rasterio
    from rasterio import features as rio_features
    from shapely.geometry import shape
    from shapely.ops import unary_union

    with rasterio.open(spec["sources"]["lidar"]) as source:
        depth = source.read(1)
        transform = source.transform
    features = []
    for index, (name, (low, high)) in enumerate(r5.BINS.items(), start=1):
        with np.errstate(invalid="ignore"):
            inside = (depth >= low) & (depth < high)
        shapes = [
            shape(geometry)
            for geometry, value in rio_features.shapes(
                inside.astype("uint8"), mask=inside, transform=transform
            )
            if value == 1
        ]
        merged = unary_union(shapes).simplify(0.00008)
        features.append(
            feature(merged.__geo_interface__, bin=name, order=index, low_m=low, high_m=high)
        )
    return collection(features)


# ---------------------------------------------------------------------------------------
# Sentinel-2 scenes: facts and image overlays
# ---------------------------------------------------------------------------------------


def s2_directory(spec: dict, scene_id: str) -> Path:
    prepared = rb.scene_directories(spec)
    if scene_id in prepared:
        return Path(prepared[scene_id]["output_dir"])
    return rv.ACOLITE_R1 / "V1" / scene_id


def true_colour(spec: dict, scene_id: str, bounds=OVERLAY_BOUNDS) -> np.ndarray:
    """Top-of-atmosphere red, green, blue over the window, one fixed stretch for all dates."""
    import rasterio
    from rasterio.windows import from_bounds

    bands = rb.acolite_rasters(s2_directory(spec, scene_id), "rhot")
    layers = []
    for nm in (665.0, 560.0, 490.0):
        _, path = rb.nearest_band(bands, nm)
        with rasterio.open(path) as source:
            window = (
                from_bounds(*bounds, transform=source.transform).round_offsets().round_lengths()
            )
            layers.append(source.read(1, window=window))
    rgb = np.dstack(layers)
    # The same linear stretch on every date, so brightness differences are real.
    low = np.array([0.02, 0.04, 0.07])
    high = np.array([0.22, 0.22, 0.25])
    scaled = np.clip((rgb - low) / (high - low), 0, 1) ** 0.8
    return np.where(np.isfinite(scaled), scaled, 0)


def scene_facts() -> list[dict]:
    """Per Sentinel-2 scene: v1's checks and r3's glint-tolerant mask, as recorded."""
    selection = json.loads((rb.ROOT / "docs" / "scene-selection.json").read_text())
    v1 = {scene["id"]: scene for scene in selection["scenes"]}
    mask = json.loads((resolve_run(rb.OUTPUT / "r3" / "mask") / "mask.json").read_text())
    screen = json.loads((resolve_run(rv.OUTPUT / "screen") / "screen.json").read_text())
    screened = {
        (scene["sensing_utc"][:10], scene["relative_orbit"]): scene["water_median"]
        for scene in screen["scenes"]
    }
    facts = []
    for scene_id in S2_SCENES:
        _, stamp, platform, orbit = scene_id.split("-")
        date = f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:]}"
        record = {
            "id": scene_id,
            "date": date,
            "platform": f"Sentinel-2{platform[-1]}",
            "orbit": orbit,
            "clean": scene_id in r5.CLEAN,
            # r1's water-only screen: median L2A reflectance at 1610 nm over the quadrat area.
            "screen_1610": screened.get((date, int(orbit[1:]))),
        }
        if scene_id in v1:
            metrics = v1[scene_id]["metrics"]
            record["v1"] = {
                "glint_rhos_1614": metrics["median_rhos_glint_band"],
                "valid_water_fraction": metrics["valid_water_fraction"],
                "quadrat_cells": metrics["quadrat_cells"],
                "valid_quadrat_cells": metrics["valid_quadrat_cells"],
                "passed": v1[scene_id]["passed"],
            }
        if scene_id in mask["scenes"]:
            scene = mask["scenes"][scene_id]
            checks = scene["acceptability"]["checks"]
            record["r3"] = {
                "glint_background_rhot_1610": scene["background"]["B"],
                "threshold": scene["threshold"],
                "within_cap": scene["within_cap"],
                "valid_water_fraction": scene["checks"]["valid_water_fraction"],
                "quadrat_cells": scene["checks"]["quadrat_cells"],
                "valid_quadrat_cells": scene["checks"]["valid_quadrat_cells"],
                "failed_checks": [name for name, passed in checks.items() if passed is False],
                "accepted": scene["usable"],
            }
        facts.append(record)
    return rv.clean(facts)


def overpass_history() -> list[dict]:
    """r1's screen of every Sentinel-2 overpass in the campaign windows, 2017-2026."""
    frame = pd.read_csv(resolve_run(rv.OUTPUT / "windows") / "overpasses.csv")
    frame = frame[frame.offset_days.abs() <= 14]
    return rv.clean(
        [
            {
                "date": row.sensing_utc[:10],
                "year": int(row.year),
                "campaign": row.campaign,
                "offset_days": int(row.offset_days),
                "orbit": f"R{int(row.relative_orbit):03d}",
                "glint_angle_deg": row.glint_angle_deg,
                "swir_glint": row.swir_glint,
                "verdict": row.verdict,
                "likely_pass": bool(row.likely_pass),
            }
            for row in frame.itertuples()
        ]
    )


def candidates_2026() -> list[dict]:
    """Every 2026 candidate in the matchup windows, with r1's water-only screen value."""
    screen = json.loads((resolve_run(rv.OUTPUT / "screen") / "screen.json").read_text())
    return rv.clean(
        [
            {
                "date": scene["sensing_utc"][:10],
                "campaign": scene["campaign"],
                "orbit": f"R{int(scene['relative_orbit']):03d}",
                "screen_1610": scene["water_median"],
                "screened_out": scene["water_screened_out"],
            }
            for scene in screen["scenes"]
        ]
    )


def overlay_corners(bounds=OVERLAY_BOUNDS) -> list:
    left, bottom, right, top = bounds
    corners = [(left, top), (right, top), (right, bottom), (left, bottom)]
    return [rounded(lonlat(x, y)) for x, y in corners]


# ---------------------------------------------------------------------------------------
# Numbers the story map shows, copied from the run outputs
# ---------------------------------------------------------------------------------------


def story_numbers(frame: pd.DataFrame) -> dict:
    results = json.loads((rb.ROOT / "docs" / "results.json").read_text())
    envelope = json.loads(
        (resolve_run(rb.OUTPUT / "r4" / "envelope") / "envelope.json").read_text()
    )
    windows = json.loads((resolve_run(rv.OUTPUT / "windows") / "windows.json").read_text())
    scored = frame[frame.cover_pct.notna()]
    cover = {
        campaign: {
            "all_events": {
                "n": int(len(group)),
                "mean": round(float(group.cover_pct.mean()), 1),
                "shallow": round(float(group[group.stratum == "shallow"].cover_pct.mean()), 1),
                "deep": round(float(group[group.stratum == "deep"].cover_pct.mean()), 1),
            },
            "eligible": {
                "n": int(group.eligible_registered.sum()),
                "mean": round(float(group[group.eligible_registered].cover_pct.mean()), 1),
            },
        }
        for campaign, group in scored.groupby("campaign")
    }
    scores = {
        "july": {
            scale: {
                family: {key: score[key] for key in ("n", "mae_pp", "bias_pp")}
                for family, score in block["loto_scores"].items()
            }
            for scale, block in results["development"]["scales"].items()
        },
        "september": {
            scale: {
                family: {key: score[key] for key in ("n", "mae_pp", "bias_pp")}
                for family, score in block["scores"].items()
            }
            for scale, block in results["holdout"]["scales"].items()
        },
    }
    recovery = {}
    for key, condition in envelope["conditions"].items():
        if condition["block"] not in ("A", "B", "B_null"):
            continue
        recovery[key] = {
            "geometry": condition["geometry"],
            "background": condition["background"],
            "sigma": condition["sigma"],
            "contrast": condition["contrast"],
            "radius": condition["radius"],
            "error": condition["error"],
            "shares": {
                model: condition[model]["share_ge_2pp"]
                for model in ("spectral_depth_ridge", "corrected_depth_ridge", "oracle_unmixing")
            },
        }
    return rv.clean(
        {
            "cover": cover,
            "scores": scores,
            "recovery": recovery,
            "units_per_geometry": envelope["units_per_geometry"],
            "habitat": hm.summary(frame),
            "windows": windows["summary"],
            "windows_by_year": [
                {
                    key: row[key]
                    for key in ("h_days", "year", "campaign", "overpasses", "likely_pass")
                }
                for row in windows["table"]
            ],
        }
    )


# ---------------------------------------------------------------------------------------
# Layers stage
# ---------------------------------------------------------------------------------------


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, separators=(",", ":"), ensure_ascii=False))


TEMPLATE = Path(__file__).with_name("report_template.html")
BANNER = "overlays/s2-20260705-S2B-R080.webp"
SUPPORT_REPORTS = {
    "repository-readme": rb.ROOT.parent / "README.md",
    **{name: rb.ROOT / "docs" / f"{name}.md" for name in (
        "plan", "september-inputs", "scene-selection", "optics-readiness", "results",
        "revision-r1", "revision-r3", "revision-r4", "revision-r5", "revision-r6",
    )},
    "e5s-kelp-spectra": rb.ROOT.parent / "optics" / "reports" / "e5s-kelp-spectra.md",
}


def document_links(body: str, source: Path, published: dict[Path, str]) -> str:
    """Resolve the published supplement; leave unavailable repository references as text."""
    def replace(match):
        href, label = match.groups()
        url = urlsplit(href)
        if url.scheme or url.netloc or not url.path or href.startswith("/sesimbra/"):
            return match.group(0)
        target = (source.parent / unquote(url.path)).resolve()
        if target not in published:
            return label
        suffix = (f"?{url.query}" if url.query else "") + (f"#{url.fragment}" if url.fragment else "")
        return f'<a href="{escape(published[target] + suffix, quote=True)}">{label}</a>'

    return re.sub(r'<a\s+href="([^"]+)">(.*?)</a>', replace, body, flags=re.DOTALL)


def markdown_html(source: Path) -> str:
    return subprocess.run(
        ["pandoc", str(source), "--from", "gfm", "--to", "html5", "--wrap=none"],
        capture_output=True, text=True, check=True,
    ).stdout


def supporting_reports(site: Path) -> dict[Path, str]:
    """Publish only the selected evidence documents, without changing historical sources."""
    destination = site / "supporting"
    destination.mkdir(parents=True, exist_ok=True)
    published = {source.resolve(): f"supporting/{name}.html" for name, source in SUPPORT_REPORTS.items()}
    published[REPORT.resolve()] = "technical-report.html"
    nested ={source: f"../{url}" for source, url in published.items()}
    styles = re.search(r"<style>(.*?)</style>", TEMPLATE.read_text(), flags=re.DOTALL).group(1)
    records = []
    for name, source in SUPPORT_REPORTS.items():
        body = document_links(markdown_html(source), source, nested)
        source_copy = destination / f"{name}.md"
        shutil.copyfile(source, source_copy)
        html_path = destination / f"{name}.html"
        html_path.write_text(
            '<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>{escape(name)} — Sesimbra supporting report</title><style>{styles}</style></head><body>'
            '<nav class="topbar"><a href="../technical-report.html">Current technical report</a>'
            '<a class="story" href="../">Story map →</a></nav>'
            '<main style="max-width:950px;margin:0 auto;padding:2rem 1.5rem"><article>'
            '<p><strong>Supporting evidence document.</strong> This is a snapshot of a '
            'repository report. Historical interpretations may be superseded by the '
            '<a href="../technical-report.html">current synthesis</a>. Repository code, '
            'registrations and input data are not included in this web supplement; '
            'unavailable repository links are shown as text.</p>'
            f'<p><a href="{name}.md" download>Download source Markdown</a> · '
            '<a href="manifest.json">File hashes</a></p>'
            f'{body}</article></main></body></html>\n'
        )
        records.append({"source": source.relative_to(rb.ROOT.parent).as_posix(),
                        "markdown": source_copy.name, "html": html_path.name,
                        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                        "html_sha256": hashlib.sha256(html_path.read_bytes()).hexdigest()})
    write_json(destination / "manifest.json", {
        "scope": "Supporting report snapshots. Source Markdown is unchanged.",
        "reports": records,
    })
    return published


def captioned(image: str, caption: str) -> str:
    """A figure whose caption leads with its label, "Map 1." or "Chart 2."."""
    label, _, rest = caption.partition(". ")
    return f"<figure>{image}<figcaption><strong>{label}.</strong> {rest}</figcaption></figure>"


def journal_body(html: str) -> dict:
    """Shape pandoc's output for the page: masthead text, contents, abstract and figures."""
    title = re.search(r"<h1[^>]*>(.*?)</h1>", html, flags=re.DOTALL)
    heading = re.sub(r"<[^>]+>", "", title.group(1)) if title else "Technical report"
    html = html[title.end() :] if title else html
    first = re.search(r"<p>(.*?)</p>", html, flags=re.DOTALL)
    dek = re.sub(r"<[^>]+>", "", first.group(1)) if first else ""
    html = html[first.end() :] if first else html
    # Images followed by an italic paragraph become captioned figures.
    html = re.sub(
        r"<p>(<img[^>]+>)</p>\s*<p><em>(.*?)</em></p>",
        lambda m: captioned(m.group(1), m.group(2)),
        html,
        flags=re.DOTALL,
    )
    html = re.sub(
        r"<table>(.*?)</table>",
        r'<div class="table-wrap"><table>\1</table></div>',
        html,
        flags=re.DOTALL,
    )
    # The first section is the summary; it reads as an abstract.
    sections = re.split(r"(?=<h2 )", html)
    if len(sections) > 1 and "summary" in sections[1][:80].lower():
        sections[1] = f'<section class="abstract">{sections[1]}</section>'
    html = "".join(sections)
    toc = "".join(
        f'<a href="#{slug}">{re.sub(r"<[^>]+>", "", text)}</a>'
        for slug, text in re.findall(r'<h2 id="([^"]+)"[^>]*>(.*?)</h2>', html, flags=re.DOTALL)
    )
    return {"heading": heading, "dek": dek, "toc": toc, "body": html}


def report_page(site: Path) -> Path:
    """The technical report as one journal-style HTML page beside the story map. Needs pandoc."""
    published = supporting_reports(site)
    html = document_links(markdown_html(REPORT), REPORT, published)
    parts = journal_body(html.replace('src="figures/', 'src="maps/'))
    path = site / "technical-report.html"
    page = Template(TEMPLATE.read_text()).substitute(title=parts["heading"], banner=BANNER, **parts)
    path.write_text(page)
    return path


def photoquadrats(spec: dict, frame: pd.DataFrame, sources: list[Path], site: Path) -> dict:
    """Copy published photoquadrats beside the story map and index them by survey event.

    Each source is a directory written by kelp_observe's publish_photoquadrats, holding
    photoquadrats.json with full/ and thumb/ images. A campaign with no source is listed
    as pending, so the page can show placeholders.
    """
    index = {campaign: {"status": "pending", "items": {}} for campaign in CAMPAIGNS}
    known = set(frame.event_id)
    for source in sources:
        published = json.loads((source / "photoquadrats.json").read_text())
        campaign = published["campaign"].split("_")[0]
        strata = spec["campaigns"][campaign]["strata"]["workbook"]
        items = {}
        for item in published["items"]:
            stratum = strata[str(item["depth_bucket_m"])]
            event_id = f"{campaign}-{stratum}-T{item['transect_int']}-Q{item['quadrat_int']}"
            if event_id not in known:
                raise ValueError(f"Photoquadrat {item['source_filename']} matches no survey event")
            record = {"taken": item["exif_datetime"]}
            for size in ("full", "thumb"):
                target = site / "photoquadrats" / campaign / size / f"{event_id}.jpg"
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(source / item[f"{size}_relpath"], target)
                record[size] = str(target.relative_to(site))
            items[event_id] = record
        index[campaign] = {"status": "available", "items": items}
    return index


def layers(site: Path, photo_sources: list[Path] | None = None) -> Path:
    """Write what the story map serves: layers, photos, scene images, maps and the report."""
    from PIL import Image

    spec = rb.load_spec()
    frame = events()
    data = site / "data"
    for name, features in hm.site_layers(spec).items():
        write_json(data / f"{name}.geojson", features)
    write_json(data / "photoquadrats.json", photoquadrats(spec, frame, photo_sources or [], site))
    write_json(data / "quadrats.geojson", quadrat_features(frame))
    write_json(data / "transects.geojson", transect_features(frame))
    write_json(data / "cells.geojson", cell_features())
    write_json(data / "areas.geojson", area_features(spec))
    write_json(data / "depth-contours.geojson", contour_features(spec))
    write_json(data / "depth-bins.geojson", depth_bin_features(spec))
    facts = scene_facts()
    for record in facts:
        image = (true_colour(spec, record["id"]) * 255).astype("uint8")
        name = f"overlays/{record['id']}.webp"
        (site / "overlays").mkdir(parents=True, exist_ok=True)
        Image.fromarray(image).save(site / name, quality=82, method=6)
        record["image"] = name
    write_json(
        data / "scenes.json",
        {
            "corners": overlay_corners(),
            "scenes": facts,
            "candidates_2026": candidates_2026(),
            "screen_limit": 0.02,
            "history": overpass_history(),
        },
    )
    write_json(data / "numbers.json", story_numbers(frame))
    (site / "maps").mkdir(parents=True, exist_ok=True)
    for figure in sorted(FIGURES.glob("*.*")):
        shutil.copy(figure, site / "maps" / figure.name)
    report_page(site)
    return site


# ---------------------------------------------------------------------------------------
# Basemaps for the static figures
# ---------------------------------------------------------------------------------------


def fetch(url: str, cache: Path) -> bytes:
    if cache.exists():
        return cache.read_bytes()
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read()
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_bytes(content)
    return content


def tile_mosaic(name: str, west, south, east, north, zoom: int):
    """Tiles covering a lon/lat box as one image, with its Web Mercator extent."""
    from PIL import Image

    def index(lon, lat):
        n = 2**zoom
        x = (lon + 180) / 360 * n
        y = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n
        return int(x), int(y)

    x0, y0 = index(west, north)
    x1, y1 = index(east, south)
    mosaic = Image.new("RGB", ((x1 - x0 + 1) * 256, (y1 - y0 + 1) * 256))
    for x in range(x0, x1 + 1):
        for y in range(y0, y1 + 1):
            url = BASEMAPS[name]["tiles"].format(z=zoom, x=x, y=y)
            content = fetch(url, TILE_CACHE / name / str(zoom) / f"{x}-{y}")
            tile = Image.open(io.BytesIO(content)).convert("RGB")
            mosaic.paste(tile, ((x - x0) * 256, (y - y0) * 256))
    world = 20037508.342789244
    size = 2 * world / 2**zoom
    extent = (
        -world + x0 * size,
        -world + (x1 + 1) * size,
        world - (y1 + 1) * size,
        world - y0 * size,
    )
    return np.asarray(mosaic), extent


def wms_image(west, south, east, north, width: int = 1800):
    """One EMODnet bathymetry image for a lon/lat box, with its Web Mercator extent."""
    from PIL import Image

    left, bottom = mercator(west, south)
    right, top = mercator(east, north)
    height = round(width * (top - bottom) / (right - left))
    url = (
        f"{BASEMAPS['emodnet']['wms']}?service=WMS&version=1.3.0&request=GetMap"
        f"&layers={BASEMAPS['emodnet']['layer']}&styles=&crs=EPSG:3857"
        f"&bbox={left:.1f},{bottom:.1f},{right:.1f},{top:.1f}"
        f"&width={width}&height={height}&format=image/png"
    )
    name = f"{west}_{south}_{east}_{north}_{width}.png"
    content = fetch(url, TILE_CACHE / "emodnet" / name)
    return np.asarray(Image.open(io.BytesIO(content)).convert("RGB")), (left, right, bottom, top)


# ---------------------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------------------

INK = "#12303d"
COVER_MAP = "YlGn"
STRATUM_MARKER = {"shallow": "o", "deep": "s"}


def save_options(path: Path) -> dict:
    """Photographic maps are saved as JPEG to keep the repository small."""
    return {"pil_kwargs": {"quality": 88}} if path.suffix == ".jpg" else {}


def new_map(west, south, east, north, width_in: float, columns: int = 1):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    left, bottom = mercator(west, south)
    right, top = mercator(east, north)
    height_in = width_in / columns * (top - bottom) / (right - left)
    figure, axes = plt.subplots(
        1, columns, figsize=(width_in, height_in + 0.5), squeeze=False, constrained_layout=True
    )
    for axis in axes[0]:
        axis.set_xlim(left, right)
        axis.set_ylim(bottom, top)
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
    return figure, list(axes[0])


def scale_bar(axis, metres: float, latitude: float, label: str | None = None, colour="white"):
    """A bar of `metres` on the ground (Web Mercator stretches distance by 1/cos(lat))."""
    import matplotlib.patheffects as effects

    left, right = axis.get_xlim()
    bottom, top = axis.get_ylim()
    length = metres / math.cos(math.radians(latitude))
    x0 = left + 0.05 * (right - left)
    y0 = bottom + 0.06 * (top - bottom)
    outline = [
        effects.withStroke(linewidth=3, foreground="black" if colour == "white" else "white")
    ]
    axis.plot([x0, x0 + length], [y0, y0], color=colour, linewidth=2.5, path_effects=outline)
    text = label or (f"{metres / 1000:g} km" if metres >= 1000 else f"{metres:g} m")
    axis.text(
        x0 + length / 2,
        y0 + 0.015 * (top - bottom),
        text,
        color=colour,
        ha="center",
        va="bottom",
        fontsize=9,
        path_effects=outline,
    )


def credit(axis, text: str, colour="white"):
    import matplotlib.patheffects as effects

    axis.text(
        0.995,
        0.008,
        text,
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.5,
        color=colour,
        path_effects=[
            effects.withStroke(linewidth=2, foreground="black" if colour == "white" else "white")
        ],
    )


def label(axis, lon, lat, text, colour="white", size=9, **style):
    import matplotlib.patheffects as effects

    x, y = mercator(lon, lat)
    axis.text(
        x,
        y,
        text,
        color=colour,
        fontsize=size,
        ha=style.pop("ha", "center"),
        va=style.pop("va", "center"),
        path_effects=[
            effects.withStroke(linewidth=2.5, foreground="black" if colour == "white" else "white")
        ],
        **style,
    )


def draw_geometry(axis, geometry: dict, **style):
    """Outline a GeoJSON polygon, multipolygon or line in Web Mercator."""
    kind = geometry["type"]
    if kind == "Polygon":
        rings = [geometry["coordinates"][0]]
    elif kind == "MultiPolygon":
        rings = [polygon[0] for polygon in geometry["coordinates"]]
    elif kind == "LineString":
        rings = [geometry["coordinates"]]
    else:
        raise ValueError(kind)
    for ring in rings:
        ring = np.asarray(ring)
        x, y = mercator(ring[:, 0], ring[:, 1])
        axis.plot(x, y, **style)


def areas_by_role(spec: dict) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = {}
    for item in area_features(spec)["features"]:
        grouped.setdefault(item["properties"]["role"], []).append(item["geometry"])
    return grouped


# ---------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------


def figure_location(spec: dict, frame: pd.DataFrame) -> Path:
    """Map 1: where the site is, and what each sensor covers (OceanStream basemap)."""
    west, south, east, north = -9.42, 38.10, -8.78, 38.60
    figure, (axis,) = new_map(west, south, east, north, 8.0)
    image, extent = tile_mosaic("oceanstream", west, south, east, north, 12)
    axis.imshow(image, extent=extent, interpolation="bilinear")
    roles = areas_by_role(spec)
    for geometry in roles["processing_extent"]:
        draw_geometry(axis, geometry, color=INK, linewidth=1.4, linestyle="--")
    for geometry in roles["pneo_footprint"]:
        draw_geometry(axis, geometry, color="#c2410c", linewidth=1.8)
    x, y = mercator(frame.longitude.mean(), frame.latitude.mean())
    axis.plot(x, y, marker="*", markersize=15, color="#facc15", markeredgecolor=INK)
    label(axis, -9.203, 38.385, "Dive site", colour=INK, size=9, weight="bold")
    label(axis, -9.10, 38.447, "Sesimbra", colour=INK, size=9)
    label(axis, -9.243, 38.418, "Cabo\nEspichel", colour=INK, size=8, ha="right")
    label(axis, -8.90, 38.535, "Setúbal", colour=INK, size=9)
    label(axis, -8.99, 38.49, "Serra da Arrábida", colour=INK, size=8, style="italic")
    label(axis, -9.33, 38.30, "Sentinel-2\nprocessing extent\n36 × 45 km", colour=INK, size=8.5)
    label(axis, -9.19, 38.455, "Pléiades Neo image\n11 km²", colour="#c2410c", size=8.5)
    scale_bar(axis, 10000, 38.35, colour=INK)
    credit(axis, BASEMAPS["oceanstream"]["credit"], colour=INK)
    inset = axis.inset_axes([0.72, 0.03, 0.26, 0.30])
    wide, wide_extent = tile_mosaic("oceanstream", -12.5, 35.5, -5.0, 44.5, 6)
    inset.imshow(wide, extent=wide_extent, interpolation="bilinear")
    left, bottom = mercator(-12.5, 35.5)
    right, top = mercator(-5.0, 44.5)
    inset.set_xlim(left, right)
    inset.set_ylim(bottom, top)
    inset.plot(x, y, marker="*", markersize=10, color="#facc15", markeredgecolor=INK)
    inset.set_xticks([])
    inset.set_yticks([])
    path = FIGURES / "map-1-location.jpg"
    figure.savefig(path, dpi=200, **save_options(path))
    return path


def figure_bathymetry(spec: dict, frame: pd.DataFrame) -> Path:
    """Map 2: depth at the site, the strata and the offshore reference (EMODnet)."""
    west, south, east, north = -9.250, 38.382, -9.140, 38.436
    figure, (axis,) = new_map(west, south, east, north, 9.0)
    image, extent = wms_image(west, south, east, north)
    axis.imshow(image, extent=extent, interpolation="bilinear")
    colours = {5.0: "#0e7490", 10.0: "#0e7490", 14.0: "#155e75", 20.0: "#164e63"}
    for item in contour_features(spec)["features"]:
        depth = item["properties"]["depth_m"]
        lidar = depth in colours
        draw_geometry(
            axis,
            item["geometry"],
            color=colours.get(depth, "#1e3a8a"),
            linewidth=0.9 if lidar else 1.1,
            linestyle="-" if lidar else ":",
        )
    roles = areas_by_role(spec)
    for geometry in roles["habitat_aoi"]:
        draw_geometry(axis, geometry, color="#15803d", linewidth=1.6)
    for geometry in roles["offshore_reference"]:
        draw_geometry(axis, geometry, color="#b91c1c", linewidth=1.4)
    for geometry in roles["pneo_footprint"]:
        draw_geometry(axis, geometry, color="#c2410c", linewidth=1.2, linestyle="--")
    x, y = mercator(frame.longitude.values, frame.latitude.values)
    axis.scatter(x, y, s=5, color="#facc15", edgecolor=INK, linewidth=0.3, zorder=5)
    label(axis, -9.2005, 38.4125, "90 quadrats", colour=INK, size=8.5, weight="bold", ha="left")
    label(axis, -9.2085, 38.4215, "Habitat area\n(provisional)", colour="#15803d", size=8)
    label(axis, -9.168, 38.4035, "Offshore reference\n51–110 m deep", colour="#b91c1c", size=8)
    label(axis, -9.160, 38.4335, "Pléiades Neo image", colour="#c2410c", size=8)
    label(
        axis,
        -9.166,
        38.3855,
        "Lidar contours 5, 10, 14 and 20 m (solid); DTM 30, 50, 75 and 100 m (dotted)",
        colour=INK,
        size=7.5,
    )
    scale_bar(axis, 1000, 38.41, colour=INK)
    credit(
        axis, BASEMAPS["emodnet"]["credit"] + "; contours from EMODnet HR lidar and DTM", colour=INK
    )
    path = FIGURES / "map-2-bathymetry.png"
    figure.savefig(path, dpi=200, **save_options(path))
    return path


def figure_transects(spec: dict, frame: pd.DataFrame) -> Path:
    """Map 3: quadrat cover in each campaign, on aerial imagery (Esri World Imagery)."""
    import matplotlib.pyplot as plt

    margin = 0.00045
    west, east = frame.longitude.min() - margin, frame.longitude.max() + margin
    south, north = frame.latitude.min() - margin * 0.8, frame.latitude.max() + margin * 0.8
    figure, axes = new_map(west, south, east, north, 13.0, columns=3)
    image, extent = tile_mosaic("esri", west, south, east, north, 19)
    cells = cell_features()["features"]
    for axis, campaign in zip(axes, CAMPAIGNS):
        axis.imshow(image, extent=extent, interpolation="bilinear")
        if campaign == "july":
            for cell in cells:
                draw_geometry(axis, cell["geometry"], color="white", linewidth=0.6, alpha=0.8)
        group = frame[frame.campaign == campaign]
        for stratum, marker in STRATUM_MARKER.items():
            part = group[group.stratum == stratum]
            scored = part[part.cover_pct.notna()]
            x, y = mercator(scored.longitude.values, scored.latitude.values)
            points = axis.scatter(
                x,
                y,
                c=scored.cover_pct,
                cmap=COVER_MAP,
                vmin=0,
                vmax=100,
                s=62,
                marker=marker,
                edgecolor=np.where(scored.eligible_registered, "black", "#dc2626"),
                linewidth=np.where(scored.eligible_registered, 0.6, 1.4),
                zorder=5,
            )
        scored = group[group.cover_pct.notna()]
        eligible = scored[scored.eligible_registered]
        title = f"{campaign.capitalize()} 2026\nmean cover {scored.cover_pct.mean():.0f}%"
        title += f" ({len(scored)} scored quadrats"
        if len(eligible) < len(scored):
            title += f"; {eligible.cover_pct.mean():.1f}% over the {len(eligible)} eligible"
        axis.set_title(title + ")", fontsize=10, color=INK)
        scale_bar(axis, 50, 38.409)
        credit(axis, BASEMAPS["esri"]["credit"])
    axes[0].scatter(
        [], [], marker="o", s=50, color="lightgrey", edgecolor="black", label="Shallow stratum"
    )
    axes[0].scatter(
        [], [], marker="s", s=50, color="lightgrey", edgecolor="black", label="Deep stratum"
    )
    axes[0].scatter(
        [],
        [],
        marker="o",
        s=50,
        color="lightgrey",
        edgecolor="#dc2626",
        linewidth=1.4,
        label="Provisional position",
    )
    axes[0].legend(loc="upper left", fontsize=8, framealpha=0.9)
    axes[1].plot([], [], color="white", linewidth=1, label="10 m grid cells with July quadrats")
    axes[1].legend(
        loc="upper left", fontsize=8, framealpha=0.5, facecolor="black", labelcolor="white"
    )
    bar = figure.colorbar(points, ax=axes, orientation="horizontal", fraction=0.04, aspect=50)
    bar.set_label("Kelp (Saccorhiza polyschides) cover in the quadrat, %")
    path = FIGURES / "map-3-transects.jpg"
    figure.savefig(path, dpi=170, bbox_inches="tight", **save_options(path))
    plt.close(figure)
    return path


def figure_glint(spec: dict, frame: pd.DataFrame) -> Path:
    """Map 4: every processed Sentinel-2 overpass in true colour, one stretch for all."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    facts = {record["id"]: record for record in scene_facts()}
    left, bottom, right, top = OVERLAY_BOUNDS
    figure, axes = plt.subplots(3, 4, figsize=(13, 8.4), constrained_layout=True)
    site = (frame.x.mean(), frame.y.mean())
    for axis, scene_id in zip(axes.flat, S2_SCENES):
        record = facts[scene_id]
        axis.imshow(true_colour(spec, scene_id), extent=(left, right, bottom, top))
        axis.plot(
            *site,
            marker="o",
            markersize=9,
            markerfacecolor="none",
            markeredgecolor="#facc15",
            markeredgewidth=1.6,
        )
        screen = record["screen_1610"]
        first = f"{record['date']} · {record['orbit']}"
        first += " · not screened" if screen is None else f" · screen {screen:.3f}"
        kept = record.get("v1", {}).get("valid_quadrat_cells")
        second = (
            "screened out before processing"
            if kept is None
            else f"registered mask keeps {kept} of {record['v1']['quadrat_cells']} cells"
        )
        axis.set_title(
            f"{first}\n{second}",
            fontsize=8.5,
            color="#166534" if record["clean"] else "#9a3412",
        )
        axis.set_xticks([])
        axis.set_yticks([])
    legend = axes.flat[-1]
    legend.axis("off")
    legend.text(
        0.0,
        0.95,
        "Top-of-atmosphere true colour,\none fixed stretch on every date.\n\n"
        "Yellow ring: the dive site.\n\n"
        "Screen: median reflectance at\n1610 nm over water at the site;\nthe limit is 0.02.\n\n"
        "Green titles: the two clean dates,\n2 July and 31 August, both orbit\n"
        "R037 on a calm sea.\n\n"
        "Window 12 × 9 km.\n"
        "Contains modified Copernicus\nSentinel data 2026.",
        va="top",
        fontsize=9,
        color=INK,
    )
    path = FIGURES / "map-4-glint.jpg"
    figure.savefig(path, dpi=170, **save_options(path))
    plt.close(figure)
    return path


FAMILY_LABELS = {
    "training_mean": "Training mean",
    "depth_ridge": "Depth only",
    "spectral_ridge": "Blue, green, red",
    "spectral_depth_ridge": "Blue, green, red + depth",
    "spectral_depth_forest": "Spectral + depth forest",
}
MODEL_LABELS = {
    "spectral_depth_ridge": "Additive ridge (v1)",
    "corrected_depth_ridge": "Depth-aware ridge",
    "oracle_unmixing": "Oracle unmixing",
}
GEOMETRY_LABELS = {
    "G1_pneo_native": "PNeo native\n29 units",
    "G2_pneo_10m": "PNeo 10 m\n18 cells",
    "G3_s2_2jul": "S2 2 July\n14 cells",
    "G4_s2_9jul": "S2 9 July\n17 cells",
}


def figure_models(numbers: dict) -> Path:
    """Chart 1: v1's errors, July leave-one-transect-out and the September holdout."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(11, 3.8), constrained_layout=True, sharey=True)
    scales = {"native_pneo": "Native 1.2 m (3 × 3 pixels)", "common_grid_10m": "10 m grid"}
    positions = np.arange(len(FAMILY_LABELS))
    for axis, campaign in zip(axes, ("july", "september")):
        for offset, (scale, name) in zip((-0.2, 0.2), scales.items()):
            scores = numbers["scores"][campaign][scale]
            values = [scores[family]["mae_pp"] for family in FAMILY_LABELS]
            bars = axis.barh(
                positions + offset,
                values,
                height=0.38,
                color="#0f766e" if scale == "native_pneo" else "#d4a03c",
                label=name,
            )
            axis.bar_label(bars, fmt="%.1f", fontsize=8, padding=2)
        axis.set_yticks(positions, list(FAMILY_LABELS.values()))
        axis.set_xlim(0, 72)
        axis.set_xlabel("Mean absolute error, percentage points of cover")
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("July: leave-one-transect-out", fontsize=10, color=INK)
    axes[1].set_title("September: sealed July fits, applied once", fontsize=10, color=INK)
    axes[0].invert_yaxis()
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    path = FIGURES / "chart-1-models.png"
    figure.savefig(path, dpi=170, **save_options(path))
    plt.close(figure)
    return path


def figure_recovery(numbers: dict) -> Path:
    """Chart 2: r4's recovery shares, ideal against real background, and position error."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    recovery = numbers["recovery"]
    figure, axes = plt.subplots(1, 3, figsize=(13, 3.9), constrained_layout=True)
    colours = {
        "spectral_depth_ridge": "#64748b",
        "corrected_depth_ridge": "#0f766e",
        "oracle_unmixing": "#d4a03c",
    }
    positions = np.arange(len(GEOMETRY_LABELS))
    panels = {
        "Ideal background, no noise": "A|{}|ideal|bare|R10.0|d0.0|median|psf|s0.0",
        "Real scene": "B|{}|real|bare|R10.0|d0.0|median|psf",
    }
    for axis, (title, pattern) in zip(axes, panels.items()):
        for offset, (model, name) in zip((-0.27, 0.0, 0.27), MODEL_LABELS.items()):
            values = [
                recovery[pattern.format(geometry)]["shares"][model] for geometry in GEOMETRY_LABELS
            ]
            axis.bar(positions + offset, values, width=0.26, color=colours[model], label=name)
        axis.axhline(0.8, color="#b91c1c", linewidth=1, linestyle="--")
        axis.text(3.45, 0.82, "criterion 0.8", color="#b91c1c", fontsize=8, ha="right")
        axis.set_xticks(positions, list(GEOMETRY_LABELS.values()), fontsize=8)
        axis.set_ylim(0, 1.05)
        axis.set_title(
            f"{title}\nkelp on bare substrate, 10 m patches, exact positions",
            fontsize=9.5,
            color=INK,
        )
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Share of 100 draws beating depth-only by ≥ 2 points")
    axes[1].legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.74))
    errors = (0.0, 2.0, 5.0, 10.0, 20.0)
    for radius, shade in zip((1.25, 2.5, 5.0, 10.0, 25.0), np.linspace(0.35, 1.0, 5)):
        values = [
            recovery[f"B|G1_pneo_native|real|bare|R{radius}|d{error}|median|psf"]["shares"][
                "spectral_depth_ridge"
            ]
            for error in errors
        ]
        axes[2].plot(
            range(len(errors)), values, marker="o", color=plt.cm.YlGn(shade), label=f"{radius:g} m"
        )
    axes[2].axhline(0.8, color="#b91c1c", linewidth=1, linestyle="--")
    axes[2].set_xticks(range(len(errors)), [f"{error:g}" for error in errors])
    axes[2].set_xlabel("Position error of each transect, m")
    axes[2].set_ylim(0, 1.05)
    axes[2].set_title(
        "Real scene, PNeo native, additive ridge:\npatch radius and position error",
        fontsize=9.5,
        color=INK,
    )
    axes[2].legend(title="Patch radius", fontsize=8, title_fontsize=8, frameon=False)
    axes[2].spines[["top", "right"]].set_visible(False)
    path = FIGURES / "chart-2-recovery.png"
    figure.savefig(path, dpi=170, **save_options(path))
    plt.close(figure)
    return path


TRACK_COLOURS = {"june": "#4ade80", "july": "#38bdf8", "september": "#fb7185"}


def fill_polygons(axis, frame, colours, alpha: float, edge="white", width=0.3, zorder=2):
    """Fill lon/lat polygons on a Web Mercator axis."""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as Outline

    patches = []
    faces = []
    for geometry, colour in zip(frame.to_crs(3857).geometry, colours):
        for polygon in getattr(geometry, "geoms", [geometry]):
            vertices = []
            codes = []
            for ring in [polygon.exterior, *polygon.interiors]:
                points = np.asarray(ring.coords)
                vertices.append(points)
                codes += (
                    [Outline.MOVETO] + [Outline.LINETO] * (len(points) - 2) + [Outline.CLOSEPOLY]
                )
            patches.append(PathPatch(Outline(np.concatenate(vertices), codes)))
            faces.append(colour)
    axis.add_collection(
        PatchCollection(
            patches, facecolor=faces, edgecolor=edge, linewidth=width, alpha=alpha, zorder=zorder
        )
    )


def figure_habitat(spec: dict, frame: pd.DataFrame) -> Path:
    """Map 5: published habitat maps and this campaign's field data, park-wide and at the site."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle

    park = hm.park_map()
    tracks = hm.mapping_tracks(spec)
    wide = (-9.245, 38.388, -8.935, 38.492)
    near = (-9.2095, 38.4050, -9.1945, 38.4125)
    figure = plt.figure(figsize=(13, 10.6), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=[1.0, 1.25], width_ratios=[2.6, 1])
    top = figure.add_subplot(grid[0, :])
    bottom = figure.add_subplot(grid[1, 0])
    key = figure.add_subplot(grid[1, 1])
    key.axis("off")

    for axis, box in ((top, wide), (bottom, near)):
        left, low = mercator(box[0], box[1])
        right, high = mercator(box[2], box[3])
        axis.set_xlim(left, right)
        axis.set_ylim(low, high)
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])

    image, extent = tile_mosaic("oceanstream", *wide, 13)
    top.imshow(image, extent=extent, interpolation="bilinear")
    order = list(hm.GROUPS)[::-1]
    ordered = park.assign(rank=park.group.map(order.index)).sort_values("rank")
    fill_polygons(top, ordered, [hm.GROUPS[group]["colour"] for group in ordered.group], 0.88)
    kelp = park[park.group == "saccorhiza"]
    fill_polygons(
        top,
        kelp,
        [hm.GROUPS["saccorhiza"]["colour"]] * len(kelp),
        1.0,
        edge=INK,
        width=0.9,
        zorder=3,
    )
    x, y = mercator(frame.longitude.mean(), frame.latitude.mean())
    top.plot(x, y, marker="*", markersize=14, color="#facc15", markeredgecolor=INK, zorder=5)
    left, low = mercator(near[0], near[1])
    right, high = mercator(near[2], near[3])
    top.add_patch(
        Rectangle(
            (left, low),
            right - left,
            high - low,
            fill=False,
            edgecolor=INK,
            linewidth=1.3,
            zorder=5,
        )
    )
    label(top, -9.10, 38.4525, "Sesimbra", colour=INK, size=9)
    label(top, -9.227, 38.4205, "Cabo Espichel", colour=INK, size=8)
    label(top, -9.202, 38.4015, "Dive site", colour=INK, size=8.5, weight="bold")
    top.set_title(
        "a. Luiz Saldanha Marine Park: mapped seabed habitats (EMODnet Seabed Habitats, map"
        " PT010000, 2013)",
        fontsize=10.5,
        color=INK,
        loc="left",
    )
    scale_bar(top, 2000, 38.44, colour=INK)
    credit(top, BASEMAPS["oceanstream"]["credit"], colour=INK)

    image, extent = tile_mosaic("esri", *near, 17)
    bottom.imshow(image, extent=extent, interpolation="bilinear")
    colours = [
        hm.COMMUNITY_COLOURS.get(code, hm.GROUPS[group]["colour"])
        for code, group in zip(park.code, park.group)
    ]
    fill_polygons(bottom, park, colours, 0.62, edge="white", width=0.5)
    for item in contour_features(spec)["features"]:
        if item["properties"]["depth_m"] in (10.0, 20.0):
            draw_geometry(
                bottom, item["geometry"], color="white", linewidth=0.8, linestyle="--", zorder=4
            )
    for track in tracks.to_crs(3857).itertuples():
        for line in getattr(track.geometry, "geoms", [track.geometry]):
            xs, ys = line.xy
            bottom.plot(xs, ys, color=TRACK_COLOURS[track.campaign], linewidth=1.8, zorder=5)
    x, y = mercator(frame.longitude.values, frame.latitude.values)
    bottom.scatter(x, y, s=9, color="white", edgecolor=INK, linewidth=0.4, zorder=6)
    bottom.set_title(
        "b. The dive site: mapped communities, CCMAR's 2026 mapping tracks and the 90 quadrats",
        fontsize=10.5,
        color=INK,
        loc="left",
    )
    scale_bar(bottom, 200, 38.409)
    credit(bottom, BASEMAPS["esri"]["credit"])

    handles = [Patch(color=group["colour"], label=group["label"]) for group in hm.GROUPS.values()]
    first = key.legend(
        handles=handles,
        title="a. Habitat groups",
        loc="upper left",
        fontsize=8.5,
        title_fontsize=9.5,
        frameon=False,
        alignment="left",
    )
    key.add_artist(first)
    names = park.drop_duplicates("code").set_index("code").name
    handles = [
        Patch(color=colour, label=names[code].replace(" on infralittoral rock", ""))
        for code, colour in hm.COMMUNITY_COLOURS.items()
    ]
    handles += [
        Line2D([], [], color=colour, linewidth=2, label=f"Mapping track, {campaign.capitalize()}")
        for campaign, colour in TRACK_COLOURS.items()
    ]
    handles += [
        Line2D([], [], color="grey", linewidth=1, linestyle="--", label="Lidar depth 10 and 20 m"),
        Line2D(
            [],
            [],
            marker="o",
            color="none",
            markerfacecolor="white",
            markeredgecolor=INK,
            markersize=5,
            label="Quadrat",
        ),
    ]
    key.legend(
        handles=handles,
        title="b. At the dive site (other colours as in a)",
        loc="lower left",
        fontsize=8.5,
        title_fontsize=9.5,
        frameon=False,
        alignment="left",
    )
    path = FIGURES / "map-5-habitat.jpg"
    figure.savefig(path, dpi=170, **save_options(path))
    plt.close(figure)
    return path


def figures() -> list[Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    spec = rb.load_spec()
    frame = events()
    numbers = story_numbers(frame)
    return [
        figure_location(spec, frame),
        figure_bathymetry(spec, frame),
        figure_transects(spec, frame),
        figure_glint(spec, frame),
        figure_habitat(spec, frame),
        figure_models(numbers),
        figure_recovery(numbers),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    layer_parser = commands.add_parser("layers", help="Write the story map's layers.")
    layer_parser.add_argument("--site", type=Path, required=True)
    layer_parser.add_argument(
        "--photoquadrats",
        type=Path,
        action="append",
        default=[],
        help="A directory written by kelp_observe's publish_photoquadrats; repeat per campaign.",
    )
    commands.add_parser("figures", help="Render the report's static maps.")
    arguments = parser.parse_args()
    if arguments.command == "layers":
        print(f"layers written to {layers(arguments.site, arguments.photoquadrats)}")
    else:
        for path in figures():
            print(path)


if __name__ == "__main__":
    main()
