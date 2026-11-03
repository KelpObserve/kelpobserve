"""Label-free radiometry, terrain and water-state diagnosis for the four r7 scenes."""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from kelpobserve import config
from kelpobserve.attenuation import bands
from kelpobserve.coastal.inversion import lee
from kelpobserve.coastal.optics import water
from kelpobserve.local_cover import clean, digest, legacy_readers
from pyproj import Geod, Transformer
from rasterio import features, windows
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROLE_WAVES = {"coastal": 443, "blue": 490, "green": 560, "red": 665, "nir": 835}


def write(path, obj):
    path.write_text(json.dumps(clean(obj), indent=2, allow_nan=False) + "\n")


def quantiles(a):
    a = np.asarray(a)
    a = a[np.isfinite(a)]
    return dict(zip(("p10", "median", "p90"), np.quantile(a, [0.1, 0.5, 0.9]))) if len(a) else {}


def warp_array(a, transform, crs, target_transform, shape, method=Resampling.nearest):
    result = np.full(shape, np.nan, np.float32)
    reproject(
        np.asarray(a, np.float32),
        result,
        src_transform=transform,
        src_crs=crs,
        src_nodata=np.nan,
        dst_transform=target_transform,
        dst_crs="EPSG:32629",
        dst_nodata=np.nan,
        resampling=method,
    )
    return result


def raster_grid(path, transform, shape, method=Resampling.nearest):
    with rasterio.open(path) as src:
        a = src.read(1, masked=True).filled(np.nan)
        return warp_array(a, src.transform, src.crs, transform, shape, method)


def optical_signal(a, bb, sza, depths, rho=0.0):
    """Vectorized Lee water + bounded-bottom signal, dimensionless above-water rhos."""
    a, bb, depths = np.asarray(a), np.asarray(bb), np.asarray(depths)
    u = bb / (a + bb)
    ic = 1 / np.cos(np.deg2rad(lee._subsurface_solar_zenith_deg(sza)))
    kc, kb = (ic + lee._du_column(u)) * (a + bb), (ic + lee._du_bottom(u)) * (a + bb)
    sub = (0.089 * u + 0.125 * u**2) * (1 - np.exp(-depths[..., None] * kc))
    sub += rho / np.pi * np.exp(-depths[..., None] * kb)
    return np.pi * water.rrs_subsurface_to_above(sub)


def candidate_states(e5s, match, wavelengths, sza):
    """Externally constrained scenarios; image reflectance is never fitted."""
    shape = bands.site_shape(e5s)
    kd0, bb0, cdm = [match["variables"][k] for k in ("KD490", "BBP", "CDM")]
    nominal = bands.at_bands(kd0["value"], wavelengths, sza, e5s)
    rows = [
        {
            "state": "r7_nominal",
            "a": nominal["a"],
            "bb": nominal["bb"],
            "kd490_target": kd0["value"],
            "valid_solution": True,
            "candidate": False,
        }
    ]
    bbshape443 = float(np.interp(443, bands.E5S_IOP_NM, shape[1]))
    an443 = float(np.interp(443, bands.E5S_IOP_NM, shape[0]))
    cdm_lower = cdm["value"] * max(0, 1 - cdm["uncertainty_pct"] / 100)
    for ks in (-1, 0, 1):
        for bs in (-1, 0, 1):
            kd_target = kd0["value"] * (1 + ks * kd0["uncertainty_pct"] / 100)
            bb_target = bb0["value"] * (1 + bs * bb0["uncertainty_pct"] / 100)
            bscale = bb_target / bbshape443
            bb490 = water.bb_water(np.array([490])) + bscale * np.interp(490, bands.E5S_IOP_NM, shape[1])

            def gap(ascale, bb490=bb490, kd_target=kd_target):
                a490 = water.a_water(np.array([490])) + ascale * np.interp(490, bands.E5S_IOP_NM, shape[0])
                return float(bands.kd(a490, bb490, sza)[0]) - kd_target

            row = {
                "state": f"KD_BBP_k{ks:+d}_b{bs:+d}",
                "candidate": True,
                "kd_step": ks,
                "bbp_step": bs,
                "kd490_target": kd_target,
                "bbp443_target": bb_target,
                "cdm443_lower": cdm_lower,
            }
            if kd_target <= 0 or bb_target <= 0 or gap(0) > 0 or gap(10000) < 0:
                rows.append({**row, "valid_solution": False, "reason": "KD490_infeasible_for_BBP"})
                continue
            ascale = brentq(gap, 0, 10000, xtol=1e-12)
            a = water.a_water(np.array(wavelengths)) + ascale * np.interp(wavelengths, bands.E5S_IOP_NM, shape[0])
            bb = water.bb_water(np.array(wavelengths)) + bscale * np.interp(wavelengths, bands.E5S_IOP_NM, shape[1])
            rows.append(
                {
                    **row,
                    "a": a,
                    "bb": bb,
                    "valid_solution": True,
                    "a_nonwater443": ascale * an443,
                    "cdm_conflict": ascale * an443 < cdm_lower,
                    "absorption_shape_scale": ascale,
                    "backscatter_shape_scale": bscale,
                }
            )
    return rows


def terrain_model(paths, copdem, centre):
    # 8 km square covers every 3 km ray from the field points.
    x, y = centre
    transform = from_origin(np.floor((x - 4000) / 10) * 10, np.ceil((y + 4000) / 10) * 10, 10, 10)
    shape = (800, 800)
    context = raster_grid(copdem, transform, shape, Resampling.bilinear)
    raw_path = paths.require("data") / "external/590_HR_Lidar_Sul.nc"
    with xr.open_dataset(raw_path) as ds:
        subset = ds.elevation.sel(lon=slice(-9.26, -9.15), lat=slice(38.37, 38.46)).load()
        dx, dy = float(np.diff(subset.lon)[0]), float(np.diff(subset.lat)[0])
        raw_transform = from_origin(float(subset.lon[0]) - dx / 2, float(subset.lat[-1]) + dy / 2, dx, dy)
        raw = warp_array(subset.values[::-1], raw_transform, "EPSG:4326", transform, shape)
    terrain = np.where(np.isfinite(raw), raw, context)
    terrain = np.maximum(terrain, 0)
    return terrain, np.isfinite(raw), transform, raw_path


def ray_horizons(terrain, transform, xy, azimuths, water_level, raised=0.0):
    """Maximum elevation angle along true-azimuth rays at 10 m spacing to 3 km."""
    xy = np.asarray(xy)
    inverse = ~transform
    to_ll = Transformer.from_crs(32629, 4326, always_xy=True)
    to_xy = Transformer.from_crs(4326, 32629, always_xy=True)
    lon, lat = to_ll.transform(*np.mean(xy, axis=0))
    geod = Geod(ellps="WGS84")
    angles, coverage = [], []
    distances = np.arange(10.0, 3000.1, 10)
    for azimuth in np.asarray(azimuths):
        lo, la, _ = geod.fwd(lon, lat, float(azimuth), 1000)
        xx, yy = to_xy.transform(lo, la)
        x0, y0 = to_xy.transform(lon, lat)
        direction = np.array([xx - x0, yy - y0])
        direction /= np.linalg.norm(direction)
        xs = xy[:, 0, None] + direction[0] * distances
        ys = xy[:, 1, None] + direction[1] * distances
        cols, rows = inverse * (xs, ys)
        heights = map_coordinates(terrain, [rows - 0.5, cols - 0.5], order=0, mode="constant", cval=np.nan)
        # Stress raises land, leaving open-water context at zero.
        heights = heights + np.where(heights > 0, raised, 0)
        h = np.rad2deg(np.arctan2(heights - water_level, distances))
        angles.append(np.maximum(0, np.nanmax(h, axis=1)))
        coverage.append(np.isfinite(heights).mean(axis=1))
    return np.array(angles).T, np.array(coverage).T


def geometry_audit(events, scene_info, terrain, transform):
    rows = []
    angles = np.arange(0, 360, 5)
    for key, scene in scene_info.items():
        points = events[events.campaign == scene["campaign"]]
        xy = np.column_stack([points.geometry.x, points.geometry.y])
        horizon, coverage = ray_horizons(terrain, transform, xy, angles, scene["water_level_m"])
        sun, _ = ray_horizons(terrain, transform, xy, [scene["saa"]], scene["water_level_m"])
        stress, water_only_stress, shifted_over_land = [], [], []
        for dx, dy in [(0, 0), (20, 0), (-20, 0), (0, 20), (0, -20)]:
            h, _ = ray_horizons(terrain, transform, xy + [dx, dy], [scene["saa"]], scene["water_level_m"], raised=5)
            stress.append(h[:, 0])
            cols, raster_rows = (~transform) * (xy[:, 0] + dx, xy[:, 1] + dy)
            at_target = map_coordinates(terrain, [raster_rows - 0.5, cols - 0.5], order=0, mode="constant", cval=np.nan)
            on_land = at_target > scene["water_level_m"]
            shifted_over_land.append(on_land)
            water_only_stress.append(np.where(on_land, np.nan, h[:, 0]))
        stress = np.max(stress, axis=0)
        water_only_stress = np.nanmax(water_only_stress, axis=0)
        shifted_over_land = np.sum(shifted_over_land, axis=0)
        elevation = 90 - scene["sza"]
        for i, (_, event) in enumerate(points.iterrows()):
            rows.append(
                {
                    "scene": key,
                    "event_id": event.event_id,
                    "campaign": event.campaign,
                    "x": xy[i, 0],
                    "y": xy[i, 1],
                    "solar_azimuth_deg": scene["saa"],
                    "solar_elevation_deg": elevation,
                    "sun_horizon_deg": sun[i, 0],
                    "sun_clearance_deg": elevation - sun[i, 0],
                    "direct_shadow": bool(sun[i, 0] >= elevation),
                    "stress_sun_horizon_deg": stress[i],
                    "stress_clearance_deg": elevation - 1 - stress[i],
                    "stress_shadow": bool(stress[i] >= elevation - 1),
                    "stress_land_positions": int(shifted_over_land[i]),
                    "water_only_stress_horizon_deg": water_only_stress[i],
                    "water_only_stress_shadow": bool(water_only_stress[i] >= elevation - 1)
                    if np.isfinite(water_only_stress[i])
                    else None,
                    "sky_view_factor": np.mean(np.cos(np.deg2rad(horizon[i])) ** 2),
                    "maximum_horizon_deg": np.max(horizon[i]),
                    "ray_coverage_min": np.min(coverage[i]),
                }
            )
    return pd.DataFrame(rows)


def extract_scene(rb, rv, spec, scene_id, variant, bounds, transform, shape, inputs):
    win = rv.read_window(spec, scene_id, bounds, variant)
    valid = win.valid()
    directory = rv.scene_dir(spec, scene_id, variant)
    rasters = {kind: rb.acolite_rasters(directory, kind) for kind in ("rhos", "rhot")}
    mask_path = rv.mask_band(spec, directory)[0]
    inputs.add(mask_path)
    sample = next(iter(rasters["rhos"].values()))
    with rasterio.open(sample) as src:
        tags = src.tags()
    meta = {
        k: tags.get("NC_GLOBAL#" + k)
        for k in [
            "sza",
            "saa",
            "vza",
            "vaa",
            "isodate",
            "ac_model",
            "ac_aot_550",
            "ac_aot_550_mean",
            "adjacency_correction",
            "dsf_residual_glint_correction",
            "dsf_residual_glint_correction_method",
        ]
    }
    for k in ("sza", "saa", "vza", "vaa"):
        meta[k] = float(meta[k])
    fraction = warp_array(valid.astype(float), win.transform, win.crs, transform, shape, Resampling.average)
    result = {"valid": fraction >= 0.8, "valid_fraction": fraction}
    wavelengths = {}
    roles = {
        **ROLE_WAVES,
        "red": 655 if scene_id.startswith("pneo") else 665,
        "nir": 825 if scene_id.startswith("pneo") else 835,
    }
    if not scene_id.startswith("pneo"):
        roles["swir"] = 1610
    for kind, table in rasters.items():
        for role, target in roles.items():
            wavelength, path = rb.nearest_band(table, target)
            # PNeo has an actual coastal-blue band; fail rather than duplicate blue.
            assert abs(wavelength - target) < 35
            inputs.add(path)
            with rasterio.open(path) as src:
                bbox = rasterio.transform.array_bounds(*win.shape, win.transform)
                w = windows.from_bounds(*bbox, transform=src.transform).round_offsets().round_lengths()
                a = rv.read_clipped(path, w)
            assert a.shape == win.shape
            result[f"{kind}_{role}"] = warp_array(
                np.where(valid, a, np.nan), win.transform, win.crs, transform, shape, Resampling.average
            )
            wavelengths[role] = wavelength
    for path in directory.glob("acolite_run_*_l2r_settings.txt"):
        inputs.add(path)
    return result, wavelengths, meta


def block_table(frame):
    frame = frame.copy()
    frame["block_x"] = np.floor(frame.x / 100).astype(int)
    frame["block_y"] = np.floor(frame.y / 100).astype(int)
    keys = ["scene", "region", "block_x", "block_y"]
    columns = [c for c in frame.select_dtypes(include="number").columns if c not in keys]
    blocks = frame.groupby(keys)[columns].median().reset_index()
    blocks["n_cells"] = frame.groupby(keys).size().values
    return blocks


def shore_associations(blocks):
    rows = []
    for (scene, region), g in blocks.groupby(["scene", "region"]):
        if region != "habitat":
            continue
        for metric in ["rhos_blue_V0", "rhos_blue_V1", "rhos_nir_V0", "removed_blue"]:
            part = g[[metric, "depth_m", "shore_m"]].dropna()
            if len(part) < 8:
                continue
            z = part.depth_m.to_numpy()
            base = np.column_stack([np.ones(len(z)), z, z * z])
            x = np.log10(np.maximum(part.shore_m.to_numpy(), 10))
            full = np.column_stack([base, x])
            y = part[metric].to_numpy()
            r0 = y - base @ np.linalg.lstsq(base, y, rcond=None)[0]
            coef = np.linalg.lstsq(full, y, rcond=None)[0]
            r1 = y - full @ coef
            total = np.sum((y - y.mean()) ** 2)
            rows.append(
                {
                    "scene": scene,
                    "metric": metric,
                    "blocks": len(part),
                    "delta_r2_distance_after_depth": (np.sum(r0 * r0) - np.sum(r1 * r1)) / total if total else None,
                    "slope_per_log10_distance": coef[-1],
                }
            )
    return pd.DataFrame(rows)


def summarize_radiometry(frame):
    rows = []
    subsets = [("all", frame)]
    for lo, hi in zip([0, 50, 100, 200, 500, 1000], [50, 100, 200, 500, 1000, 100000]):
        subsets.append((f"shore_{lo}_{hi}", frame[(frame.shore_m >= lo) & (frame.shore_m < hi)]))
    for label, part in subsets:
        for (scene, region), g in part.groupby(["scene", "region"]):
            for name in [c for c in g if c.startswith(("rhos_", "rhot_", "removed_", "atmospheric_removal_"))]:
                rows.append(
                    {
                        "scene": scene,
                        "region": region,
                        "subset": label,
                        "metric": name,
                        "cells": len(g),
                        **quantiles(g[name]),
                    }
                )
    return pd.DataFrame(rows)


def closure_for(frame, scenes, states, tolerance):
    rows = []
    for scene, g in frame[frame.region == "deep_reference"].groupby("scene"):
        meta = scenes[scene]
        z = g.depth_m.to_numpy() + meta["water_level_m"]
        for state in states[scene]:
            if not state["valid_solution"]:
                continue
            low = optical_signal(state["a"], state["bb"], meta["sza"], z)
            high = optical_signal(state["a"], state["bb"], meta["sza"], z, rho=0.6)
            for variant in ("V0", "V1"):
                obs = g[[f"rhos_{r}_{variant}" for r in ("blue", "green")]].to_numpy()
                inside = (obs >= low - tolerance) & (obs <= high + tolerance)
                all_ok = inside.all(axis=1)
                for i, band in enumerate(("blue", "green")):
                    rrs_sub = water.rrs_above_to_subsurface(obs[:, i] / np.pi)
                    u_required = np.where(
                        rrs_sub >= 0, (-0.089 + np.sqrt(np.maximum(0, 0.089**2 + 4 * 0.125 * rrs_sub))) / 0.25, np.nan
                    )
                    rows.append(
                        {
                            "scene": scene,
                            "variant": variant,
                            "state": state["state"],
                            "band": band,
                            "cells": len(g),
                            "cdm_conflict": state.get("cdm_conflict"),
                            "observed_median": np.median(obs[:, i]),
                            "water_floor_median": np.median(low[:, i]),
                            "observed_minus_floor_median": np.median(obs[:, i] - low[:, i]),
                            "below_floor_fraction": np.mean(obs[:, i] < low[:, i] - tolerance),
                            "within_bounded_bottom_fraction": np.mean(inside[:, i]),
                            "both_bands_within_fraction": np.mean(all_ok),
                            "max_bottom_contribution": np.max(high[:, i] - low[:, i]),
                            "nominal_u": state["bb"][i] / (state["a"][i] + state["bb"][i]),
                            "deep_implied_u_median": np.nanmedian(u_required),
                        }
                    )
    return pd.DataFrame(rows)


def figures(output, terrain, terrain_transform, geometry, radiometry, closure):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8, 7))
    bounds = rasterio.transform.array_bounds(*terrain.shape, terrain_transform)
    ax.imshow(
        terrain, extent=[bounds[0], bounds[2], bounds[1], bounds[3]], origin="upper", cmap="terrain", vmin=0, vmax=220
    )
    g = geometry.drop_duplicates("event_id")
    sc = ax.scatter(g.x, g.y, c=g.sky_view_factor, cmap="viridis", vmin=0.5, vmax=1, edgecolors="black", s=38)
    ax.set_xlim(g.x.min() - 700, g.x.max() + 700)
    ax.set_ylim(g.y.min() - 500, g.y.max() + 900)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.set(
        xlabel="Easting (EPSG:32629, m)",
        ylabel="Northing (m)",
        title="Sesimbra: field locations and surrounding terrain",
    )
    fig.colorbar(sc, ax=ax, label="Horizontal diffuse sky-view factor (modelled)")
    fig.tight_layout()
    fig.savefig(output / "terrain-and-quadrats.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharey=True)
    for ax, scene in zip(axes.flat, sorted(radiometry.scene.unique())):
        for region, linestyle in [("habitat", "-"), ("deep_reference", "--")]:
            g = radiometry[
                (radiometry.scene == scene) & (radiometry.region == region) & (radiometry.subset == "all")
            ].set_index("metric")
            for variant, colour in [("V0", "#b35b3a"), ("V1", "#247b89")]:
                vals = [g.loc[f"rhos_{b}_{variant}", "median"] for b in ("coastal", "blue", "green", "red", "nir")]
                ax.plot(range(5), vals, marker="o", color=colour, linestyle=linestyle, label=f"{region}, {variant}")
        c = closure[(closure.scene == scene) & (closure.variant == "V1") & (closure.state == "r7_nominal")].set_index(
            "band"
        )
        ax.plot(
            [1, 2],
            [c.loc[b, "water_floor_median"] for b in ("blue", "green")],
            color="black",
            marker="x",
            ls=":",
            label="r7 deep-water floor",
        )
        ax.set_xticks(range(5), ["Coastal", "Blue", "Green", "Red", "NIR"])
        ax.set_title(scene)
        ax.set_ylabel("Dimensionless rhos")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "paired-radiometry.png", dpi=160)
    plt.close(fig)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    for name in ("design.json", "diagnose.py", "copdem30-context.tif"):
        shutil.copy2(HERE / name, output / name)
    design = json.loads((HERE / "design.json").read_text())
    repo, paths = config.REPO, config.load()
    source = repo / "benchmark/reports/r7-local-cover-20261102-r3"
    preserved = {str(p): digest(p) for p in source.iterdir() if p.is_file()}
    inputs = {HERE / "design.json", HERE / "copdem30-context.tif", source / "scene-info.json", source / "units.csv"}
    info = json.loads((source / "scene-info.json").read_text())
    lock = json.loads((source / "input-lock.json").read_text())["files"]
    events_path = next(Path(p) for p in lock if p.endswith("/events.gpkg"))
    events = gpd.read_file(events_path).to_crs(32629)
    events = events[events.eligible_registered & events.campaign.isin(["july", "september"])].copy()
    # No field-cover targets, cover values or prediction errors are read.
    usecols = [
        "index",
        "unit_id",
        "campaign",
        "event_ids",
        "scene",
        "variant",
        "sensor",
        "scale",
        "eligible",
        "depth",
        "blue",
        "green",
        "red",
        "nir",
        "analysis_exclusion",
    ]
    units = pd.read_csv(source / "units.csv", usecols=usecols)
    spectra_path = repo / "optics/reports/e5s-kelp-spectra/results.json"
    audit_path = repo / "benchmark/registration/evidence/cmems-optical-inputs-20261102.json"
    spectra = json.loads(spectra_path.read_text())
    water_match = {t["date"]: t for t in json.loads(audit_path.read_text())["targets"]}
    inputs.update([events_path, spectra_path, audit_path])
    scene_metadata, states, tables = {}, {}, []
    with legacy_readers(paths) as (rb, rv, spec):
        reference_path = rb.REGISTRATION / spec["pneo_reference_region"]["path"]
        habitat_path = rb.ROOT / spec["sources"]["habitat_aoi"]
        reference = gpd.read_file(reference_path).to_crs(32629)
        habitat = gpd.read_file(habitat_path).to_crs(32629)
        inputs.update([reference_path, habitat_path, Path(spec["sources"]["lidar"]), Path(spec["sources"]["dtm"])])
        b = np.array([reference.total_bounds, habitat.total_bounds])
        bounds = (
            np.floor(b[:, 0].min() / 10) * 10 - 100,
            np.floor(b[:, 1].min() / 10) * 10 - 100,
            np.ceil(b[:, 2].max() / 10) * 10 + 100,
            np.ceil(b[:, 3].max() / 10) * 10 + 100,
        )
        transform = from_origin(bounds[0], bounds[3], 10, 10)
        shape = (round((bounds[3] - bounds[1]) / 10), round((bounds[2] - bounds[0]) / 10))
        refmask = features.rasterize(reference.geometry, out_shape=shape, transform=transform).astype(bool)
        habmask = features.rasterize(habitat.geometry, out_shape=shape, transform=transform).astype(bool)
        fine = raster_grid(spec["sources"]["lidar"], transform, shape, Resampling.average)
        coarse = raster_grid(spec["sources"]["dtm"], transform, shape)
        depth = np.where(np.isfinite(fine), fine, coarse)
        wet = rv.water_on(spec, transform, shape, "EPSG:32629")
        shore = distance_transform_edt(wet) * 10
        rr, cc = np.indices(shape)
        xx = transform.c + (cc + 0.5) * 10
        yy = transform.f - (rr + 0.5) * 10
        for scene_id in design["scenes"]:
            print("Extracting", scene_id, flush=True)
            pair = {}
            meta_pair = {}
            for variant in ("V0", "V1"):
                pair[variant], wavelengths, meta_pair[variant] = extract_scene(
                    rb, rv, spec, scene_id, variant, bounds, transform, shape, inputs
                )
            assert meta_pair["V0"]["ac_model"] == meta_pair["V1"]["ac_model"]
            for key in ("ac_aot_550", "ac_aot_550_mean"):
                assert meta_pair["V0"][key] == meta_pair["V1"][key]
            common = pair["V0"]["valid"] & pair["V1"]["valid"]
            original = info[scene_id + ":V1"]
            meta = {
                **meta_pair["V1"],
                "campaign": original["campaign"],
                "date": original["date"],
                "water_level_m": original["acquisition_level_m"],
                "wavelengths": wavelengths,
                "variants": meta_pair,
            }
            scene_metadata[scene_id] = meta
            states[scene_id] = candidate_states(
                spectra, water_match[meta["date"]], [wavelengths[b] for b in ("blue", "green")], meta["sza"]
            )
            for region, region_mask in [("habitat", habmask), ("deep_reference", refmask)]:
                mask = common & region_mask & np.isfinite(depth) & (depth > 0)
                data = {
                    "scene": scene_id,
                    "region": region,
                    "x": xx[mask],
                    "y": yy[mask],
                    "depth_m": depth[mask],
                    "shore_m": shore[mask],
                }
                for variant in ("V0", "V1"):
                    for key, a in pair[variant].items():
                        if key.startswith(("rhos_", "rhot_")):
                            data[key + "_" + variant] = a[mask]
                for role in wavelengths:
                    data["removed_" + role] = data["rhos_" + role + "_V0"] - data["rhos_" + role + "_V1"]
                    data["atmospheric_removal_" + role] = data["rhot_" + role + "_V0"] - data["rhos_" + role + "_V0"]
                    assert np.allclose(
                        data["rhot_" + role + "_V0"], data["rhot_" + role + "_V1"], equal_nan=True, atol=1e-7
                    )
                tables.append(pd.DataFrame(data))
            del pair
    frame = pd.concat(tables, ignore_index=True)
    frame.to_csv(output / "paired-grid-cells.csv.gz", index=False)
    blocks = block_table(frame)
    blocks.to_csv(output / "spatial-blocks.csv", index=False)
    shore_associations(blocks).to_csv(output / "shore-associations.csv", index=False)
    summary = summarize_radiometry(frame)
    summary.to_csv(output / "radiometry-summary.csv", index=False)
    closure = closure_for(frame, scene_metadata, states, design["water"]["tolerance_rhos"])
    closure.to_csv(output / "deep-water-closure.csv", index=False)
    write(output / "water-states.json", states)
    write(output / "scene-metadata.json", scene_metadata)
    # Also quantify the minimum possible signal at the original r7 support.
    footprint_rows = []
    for u in units[units.eligible].itertuples():
        state = states[u.scene][0]
        meta = scene_metadata[u.scene]
        floor = optical_signal(state["a"], state["bb"], meta["sza"], u.depth)
        footprint_rows.append(
            {
                "scene": u.scene,
                "variant": u.variant,
                "scale": u.scale,
                "unit_id": u.unit_id,
                "event_ids": u.event_ids,
                "depth_m": u.depth,
                "blue": u.blue,
                "green": u.green,
                "blue_floor": floor[0],
                "green_floor": floor[1],
                "blue_below_floor": u.blue < floor[0] - 0.001,
                "green_below_floor": u.green < floor[1] - 0.001,
            }
        )
    pd.DataFrame(footprint_rows).to_csv(output / "original-footprint-water-floor.csv", index=False)
    print("Computing terrain horizons", flush=True)
    centre = np.array([events.geometry.x.mean(), events.geometry.y.mean()])
    terrain, raw_coverage, terrain_transform, raw_path = terrain_model(paths, HERE / "copdem30-context.tif", centre)
    inputs.add(raw_path)
    geometry = geometry_audit(events, scene_metadata, terrain, terrain_transform)
    scored_pairs = set()
    for u in units[units.eligible & (units.variant == "V1")].itertuples():
        if (u.sensor == "sentinel2" and u.scale == "10m") or (u.sensor == "pleiades_neo" and u.scale == "native"):
            scored_pairs.update((u.scene, e) for e in json.loads(u.event_ids))
    geometry["on_primary_scored_support"] = [(r.scene, r.event_id) in scored_pairs for r in geometry.itertuples()]
    geometry.to_csv(output / "quadrat-geometry.csv", index=False)
    geom_index = geometry.set_index(["scene", "event_id"])
    linked = []
    for row in footprint_rows:
        g = geom_index.loc[[(row["scene"], event) for event in json.loads(row["event_ids"])]]
        linked.append(
            {
                **row,
                "sky_view_mean_over_associated_quadrats": g.sky_view_factor.mean(),
                "minimum_sun_clearance_deg": g.sun_clearance_deg.min(),
                "any_associated_quadrat_shadow": bool(g.direct_shadow.any()),
                "blue_minus_water_floor": row["blue"] - row["blue_floor"],
                "green_minus_water_floor": row["green"] - row["green_floor"],
            }
        )
    pd.DataFrame(linked).to_csv(output / "footprint-radiometry-and-geometry.csv", index=False)
    with rasterio.open(
        output / "terrain-context.tif",
        "w",
        driver="GTiff",
        width=terrain.shape[1],
        height=terrain.shape[0],
        count=2,
        dtype="float32",
        crs="EPSG:32629",
        transform=terrain_transform,
        nodata=np.nan,
        compress="deflate",
    ) as dst:
        dst.write(terrain.astype("float32"), 1)
        dst.write(raw_coverage.astype("float32"), 2)
        dst.set_band_description(1, "Elevation: coastal lidar LAT, CopDEM EGM2008 context")
        dst.set_band_description(2, "Coastal lidar coverage flag")
    figures(output, terrain, terrain_transform, geometry, summary, closure)
    assert {str(p): digest(p) for p in source.iterdir() if p.is_file()} == preserved
    write(
        output / "provenance.json",
        {
            "created_at": datetime.now(UTC).isoformat(),
            "code_sha256": digest(Path(__file__)),
            "design_sha256": digest(HERE / "design.json"),
            "inputs_sha256": {str(p): digest(p) for p in inputs},
            "source_r7_unchanged": True,
            "cover_targets_read": False,
            "cover_retrieval_or_scores_changed": False,
            "paired_cells": len(frame),
            "spatial_blocks": len(blocks),
            "terrain_context_url": "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N38_00_W010_00_DEM/Copernicus_DSM_COG_10_N38_00_W010_00_DEM.tif",
        },
    )
    print("Completed", output, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    run(p.parse_args().output.resolve())
