"""Exploratory re-analysis r3: a glint-tolerant Sentinel-2 mask.

Registered in registration/revision-r3.json before any of its results were computed. v1's
mask (top-of-atmosphere 1610 nm <= 0.0215) reads glint as cloud and removes every water
pixel on glinted summer dates. r3 raises that threshold by each scene's glint above clean
conditions, keeps v1's threshold exactly on clean scenes, and accepts a scene only if its
glint-corrected water looks like the clean scenes' on label-free measures. Reuses r1's
readers, units and models (benchmark/revision.py); writes under out/r3/.

    PYTHONPATH=. python -m benchmark.revision_r3 fetch --env-file PATH  # three new L1C bundles
    PYTHONPATH=. python -m benchmark.revision_r3 acolite                # V1 on the new scenes
    PYTHONPATH=. python -m benchmark.revision_r3 mask                   # B, T, A1-A5, selection
    PYTHONPATH=. python -m benchmark.revision_r3 models                 # conditional on selection
    PYTHONPATH=. python -m benchmark.revision_r3 injection
    PYTHONPATH=. python -m benchmark.revision_r3 windows
"""

from __future__ import annotations

import argparse
import json
import math
import re
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from kelpobserve.coastal.provenance import file_digest, new_run_directory, resolve_run

from benchmark import revision as rv
from benchmark import run_benchmark as rb

REGISTRATION_R3 = rb.REGISTRATION / "revision-r3.json"
OUTPUT = rb.OUTPUT / "r3"
L1C_R3 = rb.ROOT / ".cache" / "l1c-r3"
NEW_SCENES = {
    "s2-20260530-S2A-R037": ("june", "S2A_MSIL1C_20260530T112131"),
    "s2-20260607-S2C-R037": ("june", "S2C_MSIL1C_20260607T112111"),
    "s2-20260609-S2A-R037": ("june", "S2A_MSIL1C_20260609T112131"),
}
CLEAN = ("s2-20260702-S2B-R037", "s2-20260831-S2B-R037")
CAP = 0.05
A2_MAX_SLOPE = 0.2
MIN_REGION_PIXELS = 1000
MIN_HABITAT_PIXELS = 500


def load_registration() -> dict:
    return json.loads(REGISTRATION_R3.read_text())


def provenance() -> dict:
    return rv.provenance() | {
        "registration_r3": {"path": str(REGISTRATION_R3), "sha256": file_digest(REGISTRATION_R3)},
        "r3_runner_sha256": file_digest(Path(__file__)),
    }


def finish(stage: str, directory: Path, outputs: list[str], started: dict) -> Path:
    from kelpobserve.coastal.provenance import complete_run

    complete_run(
        OUTPUT / stage,
        directory,
        {
            "stage": f"r3-{stage}",
            "created_at": datetime.now(UTC).isoformat(),
            "exploratory": True,
            "outputs": {name: file_digest(directory / name) for name in outputs},
            **started,
        },
    )
    return directory


def scene_list(registration: dict) -> list[str]:
    scenes = registration["scenes"]
    ordered = [*scenes["june"], *scenes["july"], *scenes["controls"]["clean"]]
    ordered.append(scenes["controls"]["severe"].split()[0])
    return list(dict.fromkeys(ordered))


def scene_date(scene_id: str) -> date:
    return datetime.strptime(scene_id.split("-")[1], "%Y%m%d").date()


def campaign_of(scene_id: str, registration: dict) -> str:
    for campaign in ("june", "july"):
        if scene_id in registration["scenes"][campaign]:
            return campaign
    return {6: "june", 7: "july"}.get(scene_date(scene_id).month, "september")


# ---------------------------------------------------------------------------------------
# Stage: fetch and ACOLITE for the three new scenes
# ---------------------------------------------------------------------------------------


def fetch(env_file: Path) -> Path:
    from kelpobserve.coastal.acquire.cdse import (
        CDSECredentials,
        download_scene,
        get_access_token,
        search_scenes,
    )

    started = provenance()
    spec = rb.load_spec()
    credentials = CDSECredentials.from_env(rb.read_env_file(env_file))
    token = get_access_token(credentials)
    records = []
    for scene_id, (campaign, prefix) in NEW_SCENES.items():
        day = scene_date(scene_id)
        found = search_scenes(
            tile=spec["sentinel2"]["tile"],
            start=day,
            end=day + timedelta(days=1),
            product_type=spec["sentinel2"]["product_type"],
            token=token,
        )
        matching = rb.latest_baseline([c for c in found if c.name.startswith(prefix)])
        rb.require(len(matching) == 1, f"{scene_id}: {len(matching)} products match {prefix}")
        path = rb.download_with_backoff(download_scene, matching[0], L1C_R3 / campaign, credentials)
        records.append(
            {
                "id": scene_id,
                "name": matching[0].name,
                "path": str(path),
                "sha256": file_digest(path),
            }
        )
        print(f"  {scene_id}: {path.name}", flush=True)
    directory = new_run_directory(OUTPUT / "fetch")
    rv.write(directory / "fetch.json", {"scenes": records})
    return finish("fetch", directory, ["fetch.json"], started)


def acolite() -> Path:
    import time

    import geopandas as gpd
    from kelpobserve.coastal.acolite import correct_scene

    started = provenance()
    spec = rb.load_spec()
    config = spec["acolite"]
    revision = rb.verify_acolite(config)
    settings = {name: rb.REGISTRATION / path for name, path in config["settings"].items()}
    extra = rb.read_settings_file(settings["sentinel2"]) | rv.GLINT_VARIANT["sentinel2"]
    extent = gpd.read_file(rb.ROOT / spec["sources"]["processing_extent"])
    limit = rb.geographic_limit(extent.total_bounds, extent.crs)
    fetched = json.loads((resolve_run(OUTPUT / "fetch") / "fetch.json").read_text())["scenes"]
    directory = new_run_directory(OUTPUT / "acolite")
    records = []
    for scene in fetched:
        output = rv.ACOLITE_R1 / "V1" / scene["id"]
        began = time.monotonic()
        record = {"id": scene["id"], "input": scene["path"], "output_dir": str(output)}
        print(f"[acolite V1] {scene['id']}", flush=True)
        try:
            correct_scene(
                Path(scene["path"]),
                output,
                limit,
                Path(config["path"]),
                template=settings["common"],
                extra=extra,
                python_executable=config["python"],
            )
            record["status"] = "ok"
            record["effective_settings"] = rv.effective_settings(output)
        except Exception as error:  # recorded per scene
            record["status"] = "failed"
            record["error"] = f"{type(error).__name__}: {error}"[:2000]
        record["elapsed_s"] = round(time.monotonic() - began, 1)
        records.append(record)
        rv.write(directory / "acolite.json", {"acolite": revision, "scenes": records})
        print(
            f"[acolite V1] {scene['id']}: {record['status']} ({record['elapsed_s']} s)", flush=True
        )
    return finish("acolite", directory, ["acolite.json"], started)


# ---------------------------------------------------------------------------------------
# The glint-tolerant mask and the acceptability measures
# ---------------------------------------------------------------------------------------


def study_bounds(spec: dict) -> tuple[float, float, float, float]:
    """Habitat AOI and PNeo reference region, padded by 500 m (r1's radiometry window)."""
    import geopandas as gpd

    region = gpd.read_file(rb.REGISTRATION / spec["pneo_reference_region"]["path"]).to_crs(32629)
    habitat = gpd.read_file(rb.ROOT / spec["sources"]["habitat_aoi"]).to_crs(32629)
    a, b = region.total_bounds, habitat.total_bounds
    return (
        min(a[0], b[0]) - 500,
        min(a[1], b[1]) - 500,
        max(a[2], b[2]) + 500,
        max(a[3], b[3]) + 500,
    )


def rasterize(frame, window: rv.Window) -> np.ndarray:
    from rasterio import features

    return features.rasterize(
        frame.to_crs(window.crs).geometry, out_shape=window.shape, transform=window.transform
    ).astype(bool)


def background_water(spec: dict, window: rv.Window, habitat: np.ndarray) -> np.ndarray:
    """Habitat water at least 10 m deep and at least 100 m from any non-water pixel."""
    from rasterio.warp import Resampling
    from scipy import ndimage

    dtm = rv.raster_on(
        spec["sources"]["dtm"], window.transform, window.shape, window.crs, Resampling.nearest
    )
    depth = np.where(np.isfinite(window.lidar), window.lidar, dtm)
    pixel = abs(window.transform.a)
    far = ndimage.distance_transform_edt(window.water) * pixel >= 100.0
    with np.errstate(invalid="ignore"):
        return window.water & habitat & far & (depth >= 10.0)


def glint_background(window: rv.Window, background: np.ndarray) -> dict:
    values = window.mask_rhot[background & np.isfinite(window.mask_rhot)]
    if values.size == 0:
        return {"B": None, "iqr": None, "pixels": 0}
    q25, q50, q75 = np.percentile(values, [25, 50, 75])
    return {"B": float(q50), "iqr": float(q75 - q25), "pixels": int(values.size)}


def tolerant_threshold(base: float, background: float, clean: float) -> float:
    """v1's threshold raised by the scene's glint above clean conditions."""
    return base + max(0.0, background - clean)


def slope(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) == 0:
        return math.nan
    return float(np.polyfit(x[ok], y[ok], 1)[0])


def mad(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(np.median(np.abs(values - np.median(values)))) if values.size else math.nan


def measures(window: rv.Window, region: np.ndarray, habitat: np.ndarray) -> dict:
    """A1-A5 inputs over the valid pixels of one scene."""
    valid = window.valid()
    in_region, in_habitat = region & valid, habitat & valid
    green, nir, z = window.bands["green"], window.bands["nir"], window.lidar
    with np.errstate(invalid="ignore"):
        common = in_habitat & np.isfinite(z) & (z >= 5) & (z <= 25)
    green_slope = (
        slope(z[common], green[common]) if common.sum() >= MIN_HABITAT_PIXELS else math.nan
    )
    nir_slope = slope(z[common], nir[common]) if common.sum() >= MIN_HABITAT_PIXELS else math.nan
    enough = in_region.sum() >= MIN_REGION_PIXELS
    return {
        "region_pixels": int(in_region.sum()),
        "habitat_pixels_5_25m": int(common.sum()),
        "A1_nir_median": float(np.median(nir[in_region])) if enough else math.nan,
        "A2_glint_slope": slope(window.mask_rhot[in_region], green[in_region])
        if enough
        else math.nan,
        "A3_ratio": abs(nir_slope) / abs(green_slope) if green_slope else math.nan,
        "A4_green_mad": mad(green[in_region]) if enough else math.nan,
        "A5_green_slope": green_slope,
        "nir_slope": nir_slope,
    }


def acceptability(m: dict, clean: list[dict], glinted: bool) -> dict:
    """A1-A5 against the clean controls. A missing measure is a failure."""

    def finite(value):
        return value is not None and math.isfinite(value)

    limit_nir = max(c["A1_nir_median"] for c in clean) + 0.002
    limit_ratio = max(c["A3_ratio"] for c in clean)
    limit_mad = 2 * max(c["A4_green_mad"] for c in clean)
    clean_slopes = [c["A5_green_slope"] for c in clean]
    mean_slope = float(np.mean(clean_slopes))
    same_sign = all(np.sign(s) == np.sign(mean_slope) for s in clean_slopes)
    checks = {
        "A1": bool(finite(m["A1_nir_median"]) and m["A1_nir_median"] <= limit_nir),
        "A2": bool(finite(m["A2_glint_slope"]) and abs(m["A2_glint_slope"]) <= A2_MAX_SLOPE)
        if glinted
        else "not applicable",
        "A3": bool(finite(m["A3_ratio"]) and m["A3_ratio"] <= limit_ratio),
        "A4": bool(finite(m["A4_green_mad"]) and m["A4_green_mad"] <= limit_mad),
        "A5": bool(
            same_sign
            and finite(m["A5_green_slope"])
            and np.sign(m["A5_green_slope"]) == np.sign(mean_slope)
            and 0.5 <= m["A5_green_slope"] / mean_slope <= 2.0
        ),
    }
    passed = all(v == "not applicable" or v is True for v in checks.values())
    limits = {
        "A1_max": limit_nir,
        "A2_max_abs": A2_MAX_SLOPE,
        "A3_max": limit_ratio,
        "A4_max": limit_mad,
        "A5_band": [0.5 * mean_slope, 2.0 * mean_slope],
    }
    return {"checks": checks, "limits": limits, "accepted": bool(passed)}


def scene_checks(spec, window: rv.Window, habitat: np.ndarray, events) -> dict:
    """v1's valid-water and quadrat-cell checks under the window's mask."""
    import rasterio

    valid = window.valid()
    water = window.water & habitat
    cells = {
        rasterio.transform.rowcol(window.transform, g.x, g.y)
        for g in events.to_crs(window.crs).geometry
        if g is not None
    }
    inside = [(r, c) for r, c in cells if 0 <= r < window.shape[0] and 0 <= c < window.shape[1]]
    fraction = float((valid & water).sum() / water.sum()) if water.any() else 0.0
    cell_fraction = sum(bool(valid[r, c]) for r, c in inside) / len(cells) if cells else 0.0
    checks = spec["scene_qc"]
    return {
        "valid_water_fraction": round(fraction, 4),
        "quadrat_cells": len(cells),
        "valid_quadrat_cells": int(sum(bool(valid[r, c]) for r, c in inside)),
        "quadrat_cell_fraction": round(cell_fraction, 4),
        "passed": bool(
            fraction >= checks["valid_water_min_fraction"]
            and cell_fraction >= checks["quadrat_cell_min_fraction"]
        ),
    }


def cross_date(window: rv.Window, reference: rv.Window, habitat: np.ndarray) -> dict:
    """Median difference and Spearman correlation against another date, same grid."""
    from scipy.stats import spearmanr

    if window.shape != reference.shape:
        return {"note": "grids differ"}
    with np.errstate(invalid="ignore"):
        both = (
            habitat
            & window.valid()
            & reference.valid()
            & np.isfinite(window.lidar)
            & (window.lidar >= 5)
            & (window.lidar <= 20)
        )
    out = {"pixels": int(both.sum())}
    if both.sum() < 50:
        return out
    for role in ("blue", "green"):
        a, b = window.bands[role][both], reference.bands[role][both]
        out[role] = {
            "median_difference": float(np.median(a - b)),
            "spearman": float(spearmanr(a, b).statistic),
        }
    return out


def c1(fit: dict) -> bool | None:
    """E2c's criterion C1: reference 555-565 nm, 0 < bbp(555) < 0.02, 0.01 <= Kd <= 2."""
    if "bbp_555" not in fit:
        return None
    return bool(
        555 <= fit["reference_nm"] <= 565
        and 0 < fit["bbp_555"] < 0.02
        and 0.01 <= fit["kd_490"] <= 2.0
    )


def rank_key(row: dict, campaign: str, spec: dict) -> tuple:
    """v1's ranking: maximum offset from the plausible dates, PNeo time, valid water, name."""
    plausible = [
        date.fromisoformat(d) for d in spec["campaigns"][campaign]["plausible_survey_dates"]
    ]
    day = scene_date(row["scene"])
    offset = max(abs((day - d).days) for d in plausible)
    pneo = spec["campaigns"][campaign].get("pneo")
    if pneo:
        acquired = datetime.fromisoformat(
            spec["sources"]["pneo"][pneo]["acquired_utc"].replace("Z", "+00:00")
        )
        pneo_offset = abs((row["sensing_utc"] - acquired).total_seconds())
    else:
        pneo_offset = math.inf
    return (offset, pneo_offset, -row["valid_water_fraction"], row["scene"])


def sensing_time(directory: Path) -> datetime:
    """Sensing time from ACOLITE's output names, e.g. S2B_MSI_2026_07_02_11_30_44_T29SMC_..."""
    sample = next(iter(rb.acolite_rasters(directory, "rhos").values())).name
    match = re.search(r"_(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_", sample)
    return datetime(*map(int, match.groups()), tzinfo=UTC)


def tolerant_windows(spec: dict, registration: dict) -> tuple[dict, dict, float, float]:
    """Every r3 scene read once, with B per scene and the clean background."""
    bounds = study_bounds(spec)
    import geopandas as gpd

    habitat_frame = gpd.read_file(rb.ROOT / spec["sources"]["habitat_aoi"])
    windows, backgrounds = {}, {}
    for scene_id in scene_list(registration):
        window = rv.read_window(spec, scene_id, bounds, "V1")
        habitat = rasterize(habitat_frame, window)
        windows[scene_id] = window
        backgrounds[scene_id] = glint_background(window, background_water(spec, window, habitat))
    _, base, _ = rv.mask_band(spec, rv.scene_dir(spec, CLEAN[0], "V1"))
    clean = max(backgrounds[s]["B"] for s in CLEAN)
    return windows, backgrounds, clean, base


def mask_stage() -> Path:
    import geopandas as gpd

    started = provenance()
    spec = rb.load_spec()
    registration = load_registration()
    windows, backgrounds, clean_b, base = tolerant_windows(spec, registration)
    region_frame = gpd.read_file(rb.REGISTRATION / spec["pneo_reference_region"]["path"])
    habitat_frame = gpd.read_file(rb.ROOT / spec["sources"]["habitat_aoi"])
    events = gpd.read_file(resolve_run(rb.OUTPUT / "inventory") / "events.gpkg", layer="events")
    thresholds = {s: tolerant_threshold(base, backgrounds[s]["B"], clean_b) for s in windows}
    masked = {s: w.with_threshold(thresholds[s]) for s, w in windows.items()}
    rows = {}
    for scene_id, tolerant in masked.items():
        b, threshold = backgrounds[scene_id]["B"], thresholds[scene_id]
        region, habitat = rasterize(region_frame, tolerant), rasterize(habitat_frame, tolerant)
        campaign = campaign_of(scene_id, registration)
        rows[scene_id] = {
            "scene": scene_id,
            "campaign": campaign,
            "sensing_utc": sensing_time(rv.scene_dir(spec, scene_id, "V1")),
            "background": backgrounds[scene_id],
            "threshold": threshold,
            "within_cap": bool(b <= CAP),
            "glinted": bool(b > clean_b),
            "checks": scene_checks(spec, tolerant, habitat, events[events.campaign == campaign]),
            "measures": measures(tolerant, region, habitat),
            "cross_date": {ref: cross_date(tolerant, masked[ref], habitat) for ref in CLEAN},
        }
        print(f"  {scene_id}: B={b:.4f} T={threshold:.4f} {rows[scene_id]['checks']}", flush=True)
    clean = [rows[s]["measures"] for s in CLEAN]
    for scene_id, row in rows.items():
        row["acceptability"] = acceptability(row["measures"], clean, row["glinted"])
        row["usable"] = bool(
            row["within_cap"] and row["checks"]["passed"] and row["acceptability"]["accepted"]
        )
        row["qaa"] = rv.qaa_fits(spec, scene_id, "V1", threshold=row["threshold"])
        for fit in row["qaa"].values():
            fit["e2c_c1"] = c1(fit)
        row["valid_water_fraction"] = row["checks"]["valid_water_fraction"]
    selection = {}
    for campaign in ("june", "july"):
        candidates = [r for r in rows.values() if r["scene"] in registration["scenes"][campaign]]
        usable = sorted(
            (r for r in candidates if r["usable"]), key=lambda r: rank_key(r, campaign, spec)
        )
        accepted = [
            r["scene"] for r in candidates if r["within_cap"] and r["acceptability"]["accepted"]
        ]
        selection[campaign] = {
            "selected": usable[0]["scene"] if usable else None,
            "usable_in_rank_order": [r["scene"] for r in usable],
            "radiometrically_accepted": accepted,
        }
    directory = new_run_directory(OUTPUT / "mask")
    rv.write(
        directory / "mask.json",
        {
            "clean_background": clean_b,
            "base_threshold": base,
            "scenes": rows,
            "selection": selection,
        },
    )
    print(json.dumps(selection, indent=1), flush=True)
    return finish("mask", directory, ["mask.json"], started)


def latest_mask() -> dict:
    return json.loads((resolve_run(OUTPUT / "mask") / "mask.json").read_text())


# ---------------------------------------------------------------------------------------
# Stage: models (conditional) and injection
# ---------------------------------------------------------------------------------------


def scene_units(spec, scene_id, campaign, threshold, *, targets=None, bands=None, window=None):
    """r1 grid units on one Sentinel-2 scene under a given mask threshold."""
    events = rv.campaign_events(spec, campaign)
    level = rb.acquisition_levels(spec).get(scene_id, 0.0)
    if window is None:
        window = rv.read_window(spec, scene_id, rv.event_bounds(events, 250.0), "V1")
        window = window.with_threshold(threshold)
    units = rv.grid_units(spec, events, window, level, campaign, bands=bands)
    target = (
        rv.targets_for("september" if campaign == "september" else "july")
        if targets is None
        else targets
    )
    units = rb.attach_targets(units, target)
    return units.assign(scene=scene_id, water_level_m=level), window


def models_stage() -> Path:
    started = provenance()
    spec = rb.load_spec()
    mask = latest_mask()
    family_list = rv.families(spec)
    selection = mask["selection"]
    result: dict = {"selection": selection, "july": {}, "june": {}}
    directory = new_run_directory(OUTPUT / "models")
    september, _ = scene_units(spec, CLEAN[1], "september", mask["scenes"][CLEAN[1]]["threshold"])
    july_scenes = selection["july"]["usable_in_rank_order"]
    for rank, scene_id in enumerate(july_scenes):
        train, _ = scene_units(spec, scene_id, "july", mask["scenes"][scene_id]["threshold"])
        block, predictions = rv.model_block(train, september, spec, family_list)
        block["role"] = "selected" if rank == 0 else "secondary"
        block["support"] = rv.support(train)
        result["july"][scene_id] = block
        predictions.to_csv(directory / f"loto-{scene_id}.csv", index=False)
        print(f"  july {scene_id}: {rv._headline(block['july_loto'])}", flush=True)
    june = selection["june"]["selected"]
    if june:
        units, _ = scene_units(spec, june, "june", mask["scenes"][june]["threshold"])
        block, predictions = rv.model_block(units, None, spec, family_list)
        block["support"] = rv.support(units)
        result["june"][june] = block
        predictions.to_csv(directory / f"loto-{june}.csv", index=False)
        print(f"  june {june}: {rv._headline(block['july_loto'])}", flush=True)
    rv.write(directory / "models.json", result)
    names = ["models.json", *sorted(p.name for p in directory.glob("loto-*.csv"))]
    return finish("models", directory, names, started)


def e5s_sentinel2() -> dict:
    """E5s pooled Sentinel-2 endmembers at B2-B4 and Sesimbra IOPs at the same bands."""
    path = rb.E5S_RESULTS
    e5s = json.loads(path.read_text())
    pooled = e5s["class_bands_sentinel2"]["pooled"]
    bands = ("B2", "B3", "B4")

    def class_bands(name):
        return np.array([pooled[name]["bands"][b]["median"] for b in bands])

    understory = np.mean(
        [class_bands(n) for n in ("turf", "coralline_crust", "red_foliose")], axis=0
    )
    iops = {
        name: {
            "a": np.array(e5s["iops"][name]["a"][1:4]),
            "bb": np.array(e5s["iops"][name]["bb"][1:4]),
        }
        for name in ("median", "2026-02-22")
    }
    return {
        "path": str(path),
        "sha256": file_digest(path),
        "kelp": class_bands("kelp"),
        "background": {"understory": understory, "bare": class_bands("bare")},
        "iops": iops,
    }


def fine_field(window: rv.Window, xy_true, cover, strata, factor: int = 10) -> np.ndarray:
    """r1's synthetic field on a grid `factor` times finer, averaged back to the window."""
    from types import SimpleNamespace

    from rasterio import Affine

    fine = SimpleNamespace(
        transform=window.transform * Affine.scale(1 / factor),
        shape=(window.shape[0] * factor, window.shape[1] * factor),
    )
    field = rv.synthetic_field(fine, xy_true, cover, strata)
    return field.reshape(window.shape[0], factor, window.shape[1], factor).mean(axis=(1, 3))


def injection_draw(seed, spec, events, window, level, conditions, e5s, family_list, targets):
    rng = np.random.default_rng(seed)
    observed = targets.set_index("event_id").loc[events.event_id, "cover_pct"].to_numpy(dtype=float)
    strata = events.stratum.to_numpy()
    permuted = observed.copy()
    for s in np.unique(strata):
        idx = np.where((strata == s) & np.isfinite(observed))[0]
        permuted[idx] = observed[rng.permutation(idx)]
    synthetic = pd.DataFrame(
        {"event_id": events.event_id, "cover_pct": permuted, "cover_status": "measured"}
    )
    bearings = {t: rng.uniform(0, 2 * math.pi) for t in sorted(events.transect.unique())}
    depth = window.lidar + level
    rows = []
    for condition in conditions:
        d = condition["displacement_m"]
        shift = np.array(
            [[d * math.sin(bearings[t]), d * math.cos(bearings[t])] for t in events.transect]
        )
        bands = window.bands
        if condition["amplitude"]:
            field = fine_field(window, rv.positions(events, shift), permuted / 100.0, strata)
            delta = rv.signal_rhos(
                field,
                depth,
                e5s["kelp"],
                e5s["background"][condition["background"]],
                e5s["iops"][condition["water"]],
                window.solar_zenith_deg,
            )
            bands = dict(window.bands)
            for role, change in zip(rv.RGB, delta):
                bands[role] = window.bands[role] + change
        units = rb.attach_targets(
            rv.grid_units(spec, events, window, level, "july", bands=bands), synthetic
        )
        frame, _ = rv.loto(units, spec, family_list)
        if frame.empty:
            continue
        mae = frame.assign(e=(frame.predicted - frame.observed).abs()).groupby("family").e.mean()
        rows.append(
            {
                "draw": seed,
                "condition": condition["id"],
                "scale": "grid",
                "units": int(len(rv.complete_units(units))),
            }
            | {f"mae_{k}": round(float(v), 3) for k, v in mae.items()}
        )
    return rows


def injection_stage(draws: int = 100, jobs: int = 8) -> Path:
    from joblib import Parallel, delayed

    started = provenance()
    spec = rb.load_spec()
    mask = latest_mask()
    selected = mask["selection"]["july"]["selected"]
    scene_id = selected or CLEAN[0]
    threshold = mask["scenes"][scene_id]["threshold"] if selected else mask["base_threshold"]
    events = rv.campaign_events(spec, "july")
    level = rb.acquisition_levels(spec).get(scene_id, 0.0)
    window = rv.read_window(spec, scene_id, rv.event_bounds(events, 250.0), "V1").with_threshold(
        threshold
    )
    e5s = e5s_sentinel2()
    family_list = rv.families(spec, controls=False)
    targets = rv.targets_for("july")
    conditions = [{"id": "null", "amplitude": 0, "displacement_m": 0}]
    for background in ("understory", "bare"):
        for water in ("median", "2026-02-22"):
            for d in (0, 5, 10, 20):
                conditions.append(
                    {
                        "id": f"{background}/{water}/d{d}",
                        "amplitude": 1,
                        "background": background,
                        "water": water,
                        "displacement_m": d,
                    }
                )
    rows = Parallel(n_jobs=jobs, verbose=5)(
        delayed(injection_draw)(
            42 + i, spec, events, window, level, conditions, e5s, family_list, targets
        )
        for i in range(draws)
    )
    frame = pd.DataFrame([r for block in rows for r in block])
    directory = new_run_directory(OUTPUT / "injection")
    frame.to_csv(directory / "draws.csv", index=False)
    summary = {}
    for condition, group in frame.groupby("condition"):
        gain = group.mae_depth_ridge - group.mae_spectral_depth_ridge
        summary[condition] = {
            "draws": int(len(group)),
            "median_units": float(group.units.median()),
            "median_gain_pp": round(float(gain.median()), 2),
            "share_gain_ge_2pp": round(float((gain >= 2.0).mean()), 3),
            "recovered": bool((gain >= 2.0).mean() >= 0.8),
            "forest_share_gain_ge_2pp": round(
                float((group.mae_depth_ridge - group.mae_spectral_depth_forest >= 2.0).mean()), 3
            ),
            "spectral_vs_mean_share_ge_2pp": round(
                float((group.mae_training_mean - group.mae_spectral_ridge >= 2.0).mean()), 3
            ),
        }
    rv.write(
        directory / "injection.json",
        {
            "scene": scene_id,
            "threshold": threshold,
            "e5s": {"path": e5s["path"], "sha256": e5s["sha256"]},
            "conditions": conditions,
            "summary": summary,
        },
    )
    return finish("injection", directory, ["draws.csv", "injection.json"], started)


# ---------------------------------------------------------------------------------------
# Stage: windows re-tally (no new reads)
# ---------------------------------------------------------------------------------------


def windows_stage() -> Path:
    started = provenance()
    mask = latest_mask()
    r1_screen = json.loads((resolve_run(rv.OUTPUT / "screen") / "screen.json").read_text())
    by_day_orbit = {
        (row["sensing_utc"][:10], int(row["relative_orbit"])): row["water_median"]
        for row in r1_screen["scenes"]
    }
    accepted = [
        s
        for s, row in mask["scenes"].items()
        if row["glinted"] and row["within_cap"] and row["acceptability"]["accepted"]
    ]
    values = [
        by_day_orbit.get((scene_date(s).isoformat(), int(s.split("-R")[1]))) for s in accepted
    ]
    values = [v for v in values if v is not None]
    overpasses = pd.read_csv(resolve_run(rv.OUTPUT / "windows") / "overpasses.csv")
    result: dict = {"accepted_glinted_scenes": accepted, "G": max(values) if values else None}
    if values:
        g = max(values)
        only_glint = overpasses.reasons.fillna("[]").map(
            lambda r: set(json.loads(r.replace("'", '"'))) <= {"glint_measured"}
        )
        correctable = only_glint & overpasses.swir_glint.le(g) & overpasses.glint_angle_deg.ge(20)
        frame = overpasses.assign(correctable=correctable)
        table = {}
        for h in (3, 7, 14):
            inside = frame[(frame.offset_days <= h) & (frame.year < 2026)]
            for campaign, group in inside.groupby("campaign"):
                per_year = group.groupby("year")
                table[f"{campaign} +/-{h} d"] = {
                    "years_with_correctable": int((per_year.correctable.sum() > 0).sum()),
                    "years_with_likely_pass": int((per_year.likely_pass.sum() > 0).sum()),
                    "years": int(per_year.ngroups),
                }
        result["table"] = table
    directory = new_run_directory(OUTPUT / "windows")
    rv.write(directory / "windows.json", result)
    return finish("windows", directory, ["windows.json"], started)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    fetcher = commands.add_parser("fetch", help="Download the three new L1C bundles from CDSE.")
    fetcher.add_argument("--env-file", type=Path, required=True)
    commands.add_parser("acolite", help="ACOLITE V1 on the new scenes.")
    commands.add_parser("mask", help="Glint-tolerant mask, scene checks, A1-A5, selection.")
    commands.add_parser("models", help="Exploratory models on selected and accepted scenes.")
    injector = commands.add_parser("injection", help="Injection and recovery on Sentinel-2.")
    injector.add_argument("--draws", type=int, default=100)
    injector.add_argument("--jobs", type=int, default=8)
    commands.add_parser("windows", help="Re-tally r1's campaign windows with the correctable tier.")
    args = parser.parse_args()
    if args.command == "fetch":
        directory = fetch(args.env_file)
    elif args.command == "acolite":
        directory = acolite()
    elif args.command == "mask":
        directory = mask_stage()
    elif args.command == "models":
        directory = models_stage()
    elif args.command == "injection":
        directory = injection_stage(args.draws, args.jobs)
    else:
        directory = windows_stage()
    print(f"r3 {args.command} written to {directory}")


if __name__ == "__main__":
    main()
