"""Exploratory r4: a detectability envelope for this design.

Registered in registration/revision-r4.json before any simulation ran. Synthetic cover
(July cover permuted within strata) is painted as patches of several radii at displaced
positions, turned into reflectance by the forward model on a 1 m grid, aggregated to each
sensor through a nominal point-spread function, and added to the real scene or to an
ideal background. Four sampling geometries, a depth-aware model and a physics oracle are
compared against depth-only. Reuses r1's readers and units (benchmark/revision.py);
writes under out/r4/.

    PYTHONPATH=. python -m benchmark.revision_r4 run [--draws 100] [--jobs 8]
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from kelpobserve.coastal.provenance import file_digest, new_run_directory, resolve_run

from benchmark import revision as rv
from benchmark import run_benchmark as rb

REGISTRATION_R4 = rb.REGISTRATION / "revision-r4.json"
OUTPUT = rb.OUTPUT / "r4"
RADII = (1.25, 2.5, 5.0, 10.0, 25.0)
ERRORS = (0.0, 2.0, 5.0, 10.0, 20.0)
SIGMAS = (0.0, 0.0005, 0.001, 0.002, 0.004)
PSF_SIGMA_PX = 0.5
FINE_M = 1.0
FIELD_MARGIN_M = 60.0
WINDOW_MARGIN_M = 150.0
MODELS = ("spectral_depth_ridge", "corrected_depth_ridge", "oracle_unmixing")
FAMILIES = [
    {"id": "training_mean", "kind": "mean", "predictors": []},
    {"id": "depth_ridge", "kind": "ridge", "predictors": ["depth"]},
    {
        "id": "spectral_depth_ridge",
        "kind": "ridge",
        "predictors": ["blue", "green", "red", "depth"],
    },
    {
        "id": "corrected_depth_ridge",
        "kind": "ridge",
        "predictors": ["rb_blue", "rb_green", "depth"],
    },
]


def provenance() -> dict:
    return rv.provenance() | {
        "registration_r4": {"path": str(REGISTRATION_R4), "sha256": file_digest(REGISTRATION_R4)},
        "r4_runner_sha256": file_digest(Path(__file__)),
    }


# ---------------------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------------------


def coefficients(a: float, bb: float, sza: float) -> tuple[float, float, float]:
    """Deep-water subsurface rrs and the column and bottom attenuation of Lee's model."""
    from kelpobserve.coastal.inversion import lee

    u = bb / (a + bb)
    inv_cos = 1.0 / math.cos(math.radians(lee._subsurface_solar_zenith_deg(sza)))
    beam = a + bb
    k_c = (inv_cos + lee._du_column(u)) * beam
    k_b = (inv_cos + lee._du_bottom(u)) * beam
    return 0.089 * u + 0.125 * u * u, float(k_c), float(k_b)


def surface_rhos(rho_b, depth, a: float, bb: float, sza: float):
    """rhos for a bottom reflectance at a depth: pi x lee.forward_model, one band."""
    from kelpobserve.coastal.optics import water

    deep, k_c, k_b = coefficients(a, bb, sza)
    with np.errstate(over="ignore", invalid="ignore"):
        sub = deep * (1.0 - np.exp(-k_c * depth)) + (rho_b / np.pi) * np.exp(-k_b * depth)
    return np.pi * water.rrs_subsurface_to_above(sub)


def bottom_reflectance(rhos, depth, a: float, bb: float, sza: float):
    """The inverse of surface_rhos: lee.invert_scene's closed form, without its clipping."""
    from kelpobserve.coastal.optics import water

    deep, k_c, k_b = coefficients(a, bb, sza)
    sub = water.rrs_above_to_subsurface(np.asarray(rhos) / np.pi)
    with np.errstate(over="ignore", invalid="ignore"):
        return np.pi * (sub - deep * (1.0 - np.exp(-k_c * depth))) * np.exp(k_b * depth)


def unmix(rho_hat: list, kelp, background) -> np.ndarray:
    """Oracle cover (%) from blue and green bottom reflectance and the true endmembers."""
    fractions = [(rho_hat[b] - background[b]) / (kelp[b] - background[b]) for b in range(2)]
    return 100.0 * np.clip(np.mean(fractions, axis=0), 0.0, 1.0)


# ---------------------------------------------------------------------------------------
# Synthetic cover on a fine grid, and its passage to each sensor
# ---------------------------------------------------------------------------------------


@dataclass
class Fine:
    transform: object
    shape: tuple[int, int]
    centres: np.ndarray  # (n, 2) cell centres, EPSG:32629


def fine_grid(events) -> Fine:
    from rasterio.transform import from_origin

    x0, y0, x1, y1 = events.total_bounds
    left, top = math.floor(x0 - FIELD_MARGIN_M), math.ceil(y1 + FIELD_MARGIN_M)
    width = int(math.ceil(x1 + FIELD_MARGIN_M - left) / FINE_M)
    height = int(math.ceil(top - (y0 - FIELD_MARGIN_M)) / FINE_M)
    transform = from_origin(left, top, FINE_M, FINE_M)
    rows, cols = np.indices((height, width))
    xs, ys = transform * (cols + 0.5, rows + 0.5)
    return Fine(transform, (height, width), np.column_stack([xs.ravel(), ys.ravel()]))


def cover_field(fine: Fine, xy_true, cover, strata, radius: float) -> np.ndarray:
    """Fractional cover: discs of `radius` at true positions, else the nearest stratum mean."""
    from scipy.spatial import cKDTree

    distance, nearest = cKDTree(xy_true).query(fine.centres)
    means = {s: float(np.nanmean(cover[strata == s])) for s in np.unique(strata)}
    background = np.array([means[s] for s in strata])
    own = np.where(np.isfinite(cover), cover, background)
    field = np.where(distance <= radius, own[nearest], background[nearest])
    return field.reshape(fine.shape)


def to_sensor(values: np.ndarray, fine: Fine, window: rv.Window, psf: bool) -> np.ndarray:
    """Area-average a fine-grid signal onto the sensor grid, then the nominal PSF."""
    from rasterio.warp import Resampling, reproject
    from scipy.ndimage import gaussian_filter

    out = np.zeros(window.shape, dtype="float32")
    reproject(
        values.astype("float32"),
        out,
        src_transform=fine.transform,
        src_crs="EPSG:32629",
        dst_transform=window.transform,
        dst_crs=window.crs,
        resampling=Resampling.average,
    )
    return gaussian_filter(out, PSF_SIGMA_PX, mode="nearest") if psf else out


# ---------------------------------------------------------------------------------------
# Setup: geometries, endmembers, water
# ---------------------------------------------------------------------------------------


@dataclass
class Scene:
    name: str
    window: rv.Window
    level: float
    fine_depth: np.ndarray
    pixel_depth: np.ndarray


def setup(spec: dict) -> dict:
    from rasterio.warp import Resampling

    events = rv.campaign_events(spec, "july")
    bounds = rv.event_bounds(events, WINDOW_MARGIN_M)
    r3_mask = json.loads((resolve_run(rb.OUTPUT / "r3" / "mask") / "mask.json").read_text())
    levels = rb.acquisition_levels(spec)
    windows = {
        "pneo": ("pneo-20260708", rv.read_window(spec, "pneo-20260708", bounds, "V0")),
        "s2_2jul": (
            "s2-20260702-S2B-R037",
            rv.read_window(spec, "s2-20260702-S2B-R037", bounds, "V1").with_threshold(0.0215),
        ),
        "s2_9jul": (
            "s2-20260709-S2A-R037",
            rv.read_window(spec, "s2-20260709-S2A-R037", bounds, "V1").with_threshold(
                r3_mask["scenes"]["s2-20260709-S2A-R037"]["threshold"]
            ),
        ),
    }
    fine = fine_grid(events)
    lidar_fine = rv.raster_on(
        spec["sources"]["lidar"], fine.transform, fine.shape, "EPSG:32629", Resampling.nearest
    )
    scenes = {}
    for key, (scene_id, window) in windows.items():
        level = levels[scene_id]
        dtm = rv.raster_on(
            spec["sources"]["dtm"], window.transform, window.shape, window.crs, Resampling.nearest
        )
        scenes[key] = Scene(
            scene_id,
            window,
            level,
            lidar_fine + level,
            np.where(np.isfinite(window.lidar), window.lidar, dtm) + level,
        )
    e5s = rv.e5s_inputs()  # PNeo endmembers; E5s IOPs on S2 B2-B4
    from benchmark import revision_r3

    s2 = revision_r3.e5s_sentinel2()
    targets = rv.targets_for("july").set_index("event_id")
    return {
        "events": events,
        "observed": targets.loc[events.event_id, "cover_pct"].to_numpy(dtype=float),
        "fine": fine,
        "scenes": scenes,
        "endmembers": {
            "pneo": {"kelp": e5s["kelp"], **e5s["background"]},
            "s2_2jul": {"kelp": s2["kelp"], **s2["background"]},
            "s2_9jul": {"kelp": s2["kelp"], **s2["background"]},
        },
        "iops": {k: {q: np.asarray(v) for q, v in w.items()} for k, w in s2["iops"].items()},
        "sources": {"e5s_pneo": e5s["sha256"], "e5s_s2": s2["sha256"]},
    }


GEOMETRIES = {
    "G1_pneo_native": ("pneo", "native"),
    "G2_pneo_10m": ("pneo", "grid"),
    "G3_s2_2jul": ("s2_2jul", "grid"),
    "G4_s2_9jul": ("s2_9jul", "grid"),
}


def conditions() -> list[dict]:
    """Every registered condition."""
    out = []
    for g in GEOMETRIES:
        for contrast in ("understory", "bare"):
            for sigma in SIGMAS:
                out.append(
                    dict(
                        block="A",
                        geometry=g,
                        background="ideal",
                        sigma=sigma,
                        contrast=contrast,
                        radius=10.0,
                        error=0.0,
                        water="median",
                        psf=True,
                    )
                )
            for radius in RADII:
                for error in ERRORS:
                    out.append(
                        dict(
                            block="B",
                            geometry=g,
                            background="real",
                            sigma=None,
                            contrast=contrast,
                            radius=radius,
                            error=error,
                            water="median",
                            psf=True,
                        )
                    )
            out.append(
                dict(
                    block="C_iop",
                    geometry=g,
                    background="real",
                    sigma=None,
                    contrast=contrast,
                    radius=10.0,
                    error=0.0,
                    water="2026-02-22",
                    psf=True,
                )
            )
        out.append(
            dict(
                block="B_null",
                geometry=g,
                background="real",
                sigma=None,
                contrast=None,
                radius=None,
                error=0.0,
                water="median",
                psf=True,
            )
        )
        for radius in (2.5, 10.0):
            out.append(
                dict(
                    block="C_psf",
                    geometry=g,
                    background="real",
                    sigma=None,
                    contrast="bare",
                    radius=radius,
                    error=0.0,
                    water="median",
                    psf=False,
                )
            )
    for c in out:
        parts = [
            c["block"],
            c["geometry"],
            c["background"],
            c["contrast"] or "none",
            f"R{c['radius']}",
            f"d{c['error']}",
            c["water"],
            "psf" if c["psf"] else "box",
        ]
        if c["sigma"] is not None:
            parts.append(f"s{c['sigma']}")
        c["id"] = "|".join(str(p) for p in parts)
    return out


# ---------------------------------------------------------------------------------------
# One draw
# ---------------------------------------------------------------------------------------


def evaluate(units: pd.DataFrame, spec: dict, kelp, background, iops, sza) -> dict:
    """LOTO MAE of the fitted families and the oracle's MAE, overall and by stratum."""
    assumed = iops["median"]
    rho = [
        bottom_reflectance(
            units[role].to_numpy(), units.depth.to_numpy(), assumed["a"][b], assumed["bb"][b], sza
        )
        for b, role in enumerate(("blue", "green"))
    ]
    units = units.assign(rb_blue=rho[0], rb_green=rho[1])
    if background is not None:
        units = units.assign(oracle=unmix(rho, kelp, background))
    frame, _ = rv.loto(units, spec, FAMILIES)
    if frame.empty:
        return {}
    data = rv.complete_units(units)
    if background is not None:
        oracle = data.set_index("unit_id").loc[frame.unit_id.unique()]
        frame = pd.concat(
            [
                frame,
                pd.DataFrame(
                    {
                        "family": "oracle_unmixing",
                        "unit_id": oracle.index,
                        "stratum": oracle.stratum.to_numpy(),
                        "observed": oracle.cover_pct.to_numpy(),
                        "predicted": oracle.oracle.to_numpy(),
                    }
                ),
            ]
        )
    frame = frame.assign(error=(frame.predicted - frame.observed).abs())
    out = {"units": int(frame[frame.family == "depth_ridge"].unit_id.nunique())}
    for family, group in frame.groupby("family"):
        out[f"mae_{family}"] = round(float(group.error.mean()), 3)
        for stratum, part in group.groupby("stratum"):
            out[f"mae_{family}_{stratum}"] = round(float(part.error.mean()), 3)
    return out


def draw(seed: int, spec: dict, setup_: dict, condition_list: list[dict]) -> list[dict]:
    rng = np.random.default_rng(seed)
    events, observed = setup_["events"], setup_["observed"]
    strata = events.stratum.to_numpy()
    permuted = observed.copy()
    for s in np.unique(strata):
        idx = np.where((strata == s) & np.isfinite(observed))[0]
        permuted[idx] = observed[rng.permutation(idx)]
    synthetic = pd.DataFrame(
        {"event_id": events.event_id, "cover_pct": permuted, "cover_status": "measured"}
    )
    bearings = {t: rng.uniform(0, 2 * math.pi) for t in sorted(events.transect.unique())}
    noise_rng = np.random.default_rng(10_000 + seed)
    fine, scenes = setup_["fine"], setup_["scenes"]
    fields: dict = {}
    signals: dict = {}
    rows = []
    for c in condition_list:
        scene_key, scale = GEOMETRIES[c["geometry"]]
        scene = scenes[scene_key]
        window, sza = scene.window, scene.window.solar_zenith_deg
        members = setup_["endmembers"][scene_key]
        truth = setup_["iops"][c["water"]]
        # Signal on the sensor grid, cached per field, scene, contrast, water and PSF.
        delta = [np.zeros(window.shape, "float32")] * 3
        if c["contrast"] is not None:
            key = (c["radius"], c["error"], scene_key, c["contrast"], c["water"], c["psf"])
            if key not in signals:
                fkey = (c["radius"], c["error"])
                if fkey not in fields:
                    shift = np.array(
                        [
                            [c["error"] * math.sin(bearings[t]), c["error"] * math.cos(bearings[t])]
                            for t in events.transect
                        ]
                    )
                    fields[fkey] = cover_field(
                        fine, rv.positions(events, shift), permuted / 100.0, strata, c["radius"]
                    )
                field = fields[fkey]
                kelp, bg = members["kelp"], members[c["contrast"]]
                signals[key] = [
                    to_sensor(
                        np.nan_to_num(
                            surface_rhos(
                                field * kelp[b] + (1 - field) * bg[b],
                                scene.fine_depth,
                                truth["a"][b],
                                truth["bb"][b],
                                sza,
                            )
                            - surface_rhos(
                                np.full(fine.shape, bg[b]),
                                scene.fine_depth,
                                truth["a"][b],
                                truth["bb"][b],
                                sza,
                            )
                        ),
                        fine,
                        window,
                        c["psf"],
                    )
                    for b in range(3)
                ]
            delta = signals[key]
        if c["background"] == "real":
            base = [window.bands[r] for r in rv.RGB]
            nir = window.bands["nir"]
        else:
            bg = members[c["contrast"]]
            base = [
                surface_rhos(
                    np.full(window.shape, bg[b]),
                    scene.pixel_depth,
                    truth["a"][b],
                    truth["bb"][b],
                    sza,
                )
                + noise_rng.normal(0.0, c["sigma"], window.shape)
                for b in range(3)
            ]
            nir = np.zeros(window.shape)
        bands = {role: base[b] + delta[b] for b, role in enumerate(rv.RGB)} | {"nir": nir}
        maker = rv.native_units if scale == "native" else rv.grid_units
        units = rb.attach_targets(
            maker(spec, events, window, scene.level, "july", bands=bands), synthetic
        )
        background = members[c["contrast"]] if c["contrast"] else None
        result = evaluate(units, spec, members["kelp"], background, setup_["iops"], sza)
        rows.append({"draw": seed, "condition": c["id"]} | result)
    return rows


# ---------------------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------------------


def summarise(frame: pd.DataFrame, condition_list: list[dict]) -> dict:
    by_id = {c["id"]: c for c in condition_list}
    out = {}
    for cid, group in frame.groupby("condition"):
        entry = {k: v for k, v in by_id[cid].items() if k != "id"}
        entry["draws"] = int(len(group))
        entry["median_units"] = float(group.units.median())
        for model in MODELS:
            column = f"mae_{model}"
            if column not in group:
                continue
            gain = group.mae_depth_ridge - group[column]
            share = float((gain >= 2.0).mean())
            entry[model] = {
                "share_ge_2pp": round(share, 3),
                "median_gain_pp": round(float(gain.median()), 2),
                "recovered": bool(share >= 0.8),
            }
            for stratum in ("shallow", "deep"):
                a, b = f"mae_depth_ridge_{stratum}", f"{column}_{stratum}"
                if a in group and b in group:
                    g = group[a] - group[b]
                    entry[model][f"{stratum}_share_ge_2pp"] = round(float((g >= 2.0).mean()), 3)
        out[cid] = entry
    return out


def envelope_tables(summary: dict) -> dict:
    """Recovered (R, error) cells per geometry, contrast and model on the real background."""
    tables: dict = {}
    for entry in summary.values():
        if entry["block"] != "B":
            continue
        for model in MODELS:
            key = f"{entry['geometry']}|{entry['contrast']}|{model}"
            table = tables.setdefault(key, {})
            table[f"R{entry['radius']:g}"] = table.get(f"R{entry['radius']:g}", {}) | {
                f"d{entry['error']:g}": entry[model]["share_ge_2pp"]
            }
    return tables


def noise_tolerance(summary: dict) -> dict:
    out: dict = {}
    for entry in summary.values():
        if entry["block"] != "A":
            continue
        for model in MODELS:
            key = f"{entry['geometry']}|{entry['contrast']}|{model}"
            out.setdefault(key, {})[f"sigma{entry['sigma']:g}"] = entry[model]["share_ge_2pp"]
    return out


def pipeline_check(summary: dict) -> dict:
    entry = next(
        e
        for e in summary.values()
        if e["block"] == "A"
        and e["geometry"] == "G1_pneo_native"
        and e["contrast"] == "bare"
        and e["sigma"] == 0.0
    )
    shares = {m: entry[m]["share_ge_2pp"] for m in ("corrected_depth_ridge", "oracle_unmixing")}
    return {"shares": shares, "passed": bool(all(v >= 0.95 for v in shares.values()))}


def run(draws: int = 100, jobs: int = 8) -> Path:
    from joblib import Parallel, delayed
    from kelpobserve.coastal.provenance import complete_run

    started = provenance()
    spec = rb.load_spec()
    setup_ = setup(spec)
    condition_list = conditions()
    rows = Parallel(n_jobs=jobs, verbose=5)(
        delayed(draw)(42 + i, spec, setup_, condition_list) for i in range(draws)
    )
    frame = pd.DataFrame([r for block in rows for r in block])
    summary = summarise(frame, condition_list)
    directory = new_run_directory(OUTPUT / "envelope")
    frame.to_csv(directory / "draws.csv", index=False)
    rv.write(
        directory / "envelope.json",
        {
            "pipeline_check": pipeline_check(summary),
            "units_per_geometry": {
                g: float(frame[frame.condition.str.contains(f"|{g}|", regex=False)].units.median())
                for g in GEOMETRIES
            },
            "envelope_real_background": envelope_tables(summary),
            "noise_tolerance_ideal_background": noise_tolerance(summary),
            "conditions": summary,
            "sources": setup_["sources"],
        },
    )
    complete_run(
        OUTPUT / "envelope",
        directory,
        {
            "stage": "r4-envelope",
            "created_at": datetime.now(UTC).isoformat(),
            "exploratory": True,
            "draws": draws,
            "outputs": {n: file_digest(directory / n) for n in ("draws.csv", "envelope.json")},
            **started,
        },
    )
    return directory


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    runner = commands.add_parser("run", help="Simulate every registered condition.")
    runner.add_argument("--draws", type=int, default=100)
    runner.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    print(f"r4 written to {run(args.draws, args.jobs)}")


if __name__ == "__main__":
    main()
