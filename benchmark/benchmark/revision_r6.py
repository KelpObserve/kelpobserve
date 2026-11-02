"""Exploratory r6: glint correction and multi-date combination against the recovery envelope.

Registered in registration/revision-r6.json before anything was computed. Stages:

    PYTHONPATH=. python -m benchmark.revision_r6 check     # A1: inject, then ACOLITE V0 and V1
    PYTHONPATH=. python -m benchmark.revision_r6 pneo      # A: paired PNeo V0/V1 recovery
    PYTHONPATH=. python -m benchmark.revision_r6 s2        # B: Sentinel-2 dates, combined
    PYTHONPATH=. python -m benchmark.revision_r6 bias      # C: background bias at the units
    PYTHONPATH=. python -m benchmark.revision_r6 windows   # D: matched spatial windows
    PYTHONPATH=. python -m benchmark.revision_r6 offset    # E: offset control (amendment)

Reuses r4's physics, fields and models (benchmark/revision_r4.py) and r5's unit-scale fields
(benchmark/revision_r5.py). Writes under out/r6/.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from kelpobserve.coastal.provenance import complete_run, file_digest, new_run_directory, resolve_run

from benchmark import revision as rv
from benchmark import revision_r3 as r3
from benchmark import revision_r4 as r4
from benchmark import revision_r5 as r5
from benchmark import run_benchmark as rb

REGISTRATION_R6 = rb.REGISTRATION / "revision-r6.json"
OUTPUT = rb.OUTPUT / "r6"
ACOLITE_R6 = rb.ROOT / ".cache" / "acolite-r6"
PNEO_R6 = rb.ROOT / ".cache" / "pneo-r6"
S2_SIGMAS = (0.0, 0.0005, 0.001)
PNEO_ROLE_BAND = {"red": 1, "green": 2, "blue": 3}  # RGB band file order (DIMAP BAND_INDEX)


def provenance() -> dict:
    return rv.provenance() | {
        "registration_r6": {"path": str(REGISTRATION_R6), "sha256": file_digest(REGISTRATION_R6)},
        "r6_runner_sha256": file_digest(Path(__file__)),
    }


def finish(stage: str, directory: Path, outputs: list[str], started: dict) -> Path:
    complete_run(
        OUTPUT / stage,
        directory,
        {
            "stage": f"r6-{stage}",
            "created_at": datetime.now(UTC).isoformat(),
            "exploratory": True,
            "outputs": {n: file_digest(directory / n) for n in outputs},
            **started,
        },
    )
    return directory


# ---------------------------------------------------------------------------------------
# Shared setup and signal
# ---------------------------------------------------------------------------------------


def restrict(window: rv.Window, mask: np.ndarray) -> rv.Window:
    """The same window, valid only where `mask` (and its own validity) holds."""
    return replace(window, mask_ok=window.mask_ok & mask)


def setup(spec: dict) -> dict:
    s = r4.setup(spec)
    pneo_v0 = s["scenes"]["pneo"]
    bounds = rv.event_bounds(s["events"], r4.WINDOW_MARGIN_M)
    v1_window = rv.read_window(spec, "pneo-20260708", bounds, "V1")
    common = pneo_v0.window.valid() & v1_window.valid()
    s["pneo"] = {
        "V0": replace(pneo_v0, window=restrict(pneo_v0.window, common)),
        "V1": replace(pneo_v0, name="pneo-20260708-V1", window=restrict(v1_window, common)),
    }
    w2, w9 = s["scenes"]["s2_2jul"].window, s["scenes"]["s2_9jul"].window
    rb.require(w2.transform == w9.transform and w2.shape == w9.shape, "S2 dates differ in grid.")
    return s


def permutation(seed: int, s: dict):
    rng = np.random.default_rng(seed)
    events, observed = s["events"], s["observed"]
    strata = events.stratum.to_numpy()
    permuted = observed.copy()
    for stratum in np.unique(strata):
        idx = np.where((strata == stratum) & np.isfinite(observed))[0]
        permuted[idx] = observed[rng.permutation(idx)]
    bearings = {t: rng.uniform(0, 2 * math.pi) for t in sorted(events.transect.unique())}
    synthetic = pd.DataFrame(
        {"event_id": events.event_id, "cover_pct": permuted, "cover_status": "measured"}
    )
    return permuted, strata, bearings, synthetic


class Signals:
    """Fine fields and sensor-grid signals for one draw, cached."""

    def __init__(self, s: dict, permuted, strata, bearings):
        self.s, self.permuted, self.strata, self.bearings = s, permuted, strata, bearings
        self.fields, self.cache = {}, {}

    def field(self, radius: float, error: float) -> np.ndarray:
        key = (radius, error)
        if key not in self.fields:
            events = self.s["events"]
            shift = np.array(
                [
                    [error * math.sin(self.bearings[t]), error * math.cos(self.bearings[t])]
                    for t in events.transect
                ]
            )
            self.fields[key] = r4.cover_field(
                self.s["fine"],
                rv.positions(events, shift),
                self.permuted / 100.0,
                self.strata,
                radius,
            )
        return self.fields[key]

    def delta(
        self, scene, members, contrast, radius, error, water="median", psf=True, grid=None
    ) -> list:
        key = (scene.name, contrast, radius, error, water, psf, id(grid))
        if key not in self.cache:
            field = self.field(radius, error)
            truth = self.s["iops"][water]
            kelp, bg = members["kelp"], members[contrast]
            sza = scene.window.solar_zenith_deg
            target = grid or scene.window
            self.cache[key] = [
                r4.to_sensor(
                    np.nan_to_num(
                        r4.surface_rhos(
                            field * kelp[b] + (1 - field) * bg[b],
                            scene.fine_depth,
                            truth["a"][b],
                            truth["bb"][b],
                            sza,
                        )
                        - r4.surface_rhos(
                            np.full(field.shape, bg[b]),
                            scene.fine_depth,
                            truth["a"][b],
                            truth["bb"][b],
                            sza,
                        )
                    ),
                    self.s["fine"],
                    target,
                    psf,
                )
                for b in range(3)
            ]
        return self.cache[key]


def ideal_background(scene, members, contrast, sigma, rng, iops) -> list:
    """Forward-model background substrate at each pixel's depth, plus independent noise."""
    return [
        r4.surface_rhos(
            np.full(scene.window.shape, members[contrast][b]),
            scene.pixel_depth,
            iops["a"][b],
            iops["bb"][b],
            scene.window.solar_zenith_deg,
        )
        + rng.normal(0.0, sigma, scene.window.shape)
        for b in range(3)
    ]


def units_for(spec, events, window: rv.Window, level: float, scale: str, bands: dict, synthetic):
    maker = rv.native_units if scale == "native" else rv.grid_units
    return rb.attach_targets(maker(spec, events, window, level, "july", bands=bands), synthetic)


def with_nir(bands: list, window: rv.Window) -> dict:
    return {role: bands[b] for b, role in enumerate(rv.RGB)} | {
        "nir": np.nan_to_num(window.bands["nir"])
    }


# ---------------------------------------------------------------------------------------
# A: paired PNeo V0/V1
# ---------------------------------------------------------------------------------------


def pneo_conditions() -> list[dict]:
    out = []
    for version in ("V0", "V1"):
        for scale in ("native", "grid"):
            for contrast in ("understory", "bare"):
                for radius in r4.RADII:
                    for error in r4.ERRORS:
                        out.append(
                            dict(
                                version=version,
                                scale=scale,
                                background="real",
                                contrast=contrast,
                                radius=radius,
                                error=error,
                                sigma=None,
                            )
                        )
            out.append(
                dict(
                    version=version,
                    scale=scale,
                    background="real",
                    contrast=None,
                    radius=None,
                    error=0.0,
                    sigma=None,
                )
            )
    for scale in ("native", "grid"):
        for contrast in ("understory", "bare"):
            for sigma in r4.SIGMAS:
                out.append(
                    dict(
                        version="ideal",
                        scale=scale,
                        background="ideal",
                        contrast=contrast,
                        radius=10.0,
                        error=0.0,
                        sigma=sigma,
                    )
                )
    for c in out:
        c["id"] = "|".join(
            str(c[k])
            for k in ("version", "scale", "background", "contrast", "radius", "error", "sigma")
        )
    return out


def pneo_draw(seed: int, spec: dict, s: dict, conditions: list[dict]) -> list[dict]:
    permuted, strata, bearings, synthetic = permutation(seed, s)
    signals = Signals(s, permuted, strata, bearings)
    noise = np.random.default_rng(10_000 + seed)
    members, events, rows = s["endmembers"]["pneo"], s["events"], []
    for c in conditions:
        scene = s["pneo"]["V0" if c["version"] == "ideal" else c["version"]]
        window, sza = scene.window, scene.window.solar_zenith_deg
        if c["background"] == "ideal":
            base = ideal_background(
                scene, members, c["contrast"], c["sigma"], noise, s["iops"]["median"]
            )
        else:
            base = [window.bands[r] for r in rv.RGB]
        if c["contrast"] is None:
            bands = base
        else:
            delta = signals.delta(
                s["scenes"]["pneo"], members, c["contrast"], c["radius"], c["error"]
            )
            bands = [base[b] + delta[b] for b in range(3)]
        units = units_for(
            spec, events, window, scene.level, c["scale"], with_nir(bands, window), synthetic
        )
        background = members[c["contrast"]] if c["contrast"] else None
        rows.append(
            {"draw": seed, "condition": c["id"]}
            | r4.evaluate(units, spec, members["kelp"], background, s["iops"], sza)
        )
    return rows


# ---------------------------------------------------------------------------------------
# B: Sentinel-2 dates and their combination
# ---------------------------------------------------------------------------------------

CONFIGS = ("S1_2jul_common", "S2_9jul_common", "S3_combined_common", "S4_combined_union")


def s2_conditions() -> list[dict]:
    out = []
    for config in CONFIGS:
        for contrast in ("understory", "bare"):
            for radius in r4.RADII:
                for error in r4.ERRORS:
                    out.append(
                        dict(
                            config=config,
                            background="real",
                            contrast=contrast,
                            radius=radius,
                            error=error,
                            sigma=None,
                        )
                    )
            for sigma in S2_SIGMAS:
                out.append(
                    dict(
                        config=config,
                        background="ideal",
                        contrast=contrast,
                        radius=10.0,
                        error=0.0,
                        sigma=sigma,
                    )
                )
        out.append(
            dict(
                config=config, background="real", contrast=None, radius=None, error=0.0, sigma=None
            )
        )
    for c in out:
        c["id"] = "|".join(
            str(c[k]) for k in ("config", "background", "contrast", "radius", "error", "sigma")
        )
    return out


def combine(windows: list[rv.Window], bands_by_date: list[list], require_all: bool):
    """Per-pixel mean over the dates on which a pixel is valid (all dates if require_all)."""
    valid = [w.valid() for w in windows]
    mask = np.logical_and.reduce(valid) if require_all else np.logical_or.reduce(valid)
    combined = []
    for b in range(3):
        stack = np.stack([np.where(v, bands[b], np.nan) for v, bands in zip(valid, bands_by_date)])
        with np.errstate(invalid="ignore"):
            mean = np.nanmean(stack, axis=0) if not require_all else stack.mean(axis=0)
        combined.append(np.where(mask, mean, np.nan))
    return combined, mask


def cell_levels(spec, s: dict) -> pd.Series:
    """Per union cell, the mean water level of the dates on which the cell is valid."""
    events = s["events"]
    out = {}
    for key in ("s2_2jul", "s2_9jul"):
        scene = s["scenes"][key]
        units = rv.grid_units(spec, events, scene.window, 0.0, "july")
        for uid, ok in zip(units.unit_id, units.reason == "ok"):
            if ok:
                out.setdefault(uid, []).append(scene.level)
    return pd.Series({uid: float(np.mean(v)) for uid, v in out.items()})


def s2_draw(
    seed: int, spec: dict, s: dict, conditions: list[dict], levels: pd.Series
) -> list[dict]:
    permuted, strata, bearings, synthetic = permutation(seed, s)
    signals = Signals(s, permuted, strata, bearings)
    noise = np.random.default_rng(20_000 + seed)
    events, members = s["events"], s["endmembers"]["s2_2jul"]
    d2, d9 = s["scenes"]["s2_2jul"], s["scenes"]["s2_9jul"]
    common = d2.window.valid() & d9.window.valid()
    rows = []
    for c in conditions:
        dated = []
        for scene in (d2, d9):
            if c["background"] == "ideal":
                base = ideal_background(
                    scene, members, c["contrast"], c["sigma"], noise, s["iops"]["median"]
                )
            else:
                base = [scene.window.bands[r] for r in rv.RGB]
            if c["contrast"] is not None:
                delta = signals.delta(scene, members, c["contrast"], c["radius"], c["error"])
                base = [base[b] + delta[b] for b in range(3)]
            dated.append(base)
        config = c["config"]
        if config == "S1_2jul_common":
            window, bands, level, sza = (
                restrict(d2.window, common),
                dated[0],
                d2.level,
                d2.window.solar_zenith_deg,
            )
        elif config == "S2_9jul_common":
            window, bands, level, sza = (
                restrict(d9.window, common),
                dated[1],
                d9.level,
                d9.window.solar_zenith_deg,
            )
        else:
            union = config == "S4_combined_union"
            bands, mask = combine([d2.window, d9.window], dated, require_all=not union)
            window = replace(
                d2.window,
                mask_ok=mask,
                water=d2.window.water | d9.window.water,
                bands=dict(d2.window.bands) | {r: bands[b] for b, r in enumerate(rv.RGB)},
            )
            level = (d2.level + d9.level) / 2
            sza = (d2.window.solar_zenith_deg + d9.window.solar_zenith_deg) / 2
        units = units_for(
            spec,
            events,
            window,
            0.0 if config == "S4_combined_union" else level,
            "grid",
            with_nir(bands, window),
            synthetic,
        )
        if config == "S4_combined_union":
            units = units.assign(
                depth=units.depth + units.unit_id.map(levels).fillna((d2.level + d9.level) / 2)
            )
        background = members[c["contrast"]] if c["contrast"] else None
        rows.append(
            {"draw": seed, "condition": c["id"]}
            | r4.evaluate(units, spec, members["kelp"], background, s["iops"], sza)
        )
    return rows


# ---------------------------------------------------------------------------------------
# Summaries shared by A and B
# ---------------------------------------------------------------------------------------


def summarise(frame: pd.DataFrame, conditions: list[dict]) -> dict:
    by_id = {c["id"]: c for c in conditions}
    out = {}
    for cid, group in frame.groupby("condition"):
        entry = {k: v for k, v in by_id[cid].items() if k != "id"}
        entry["draws"] = int(len(group))
        entry["median_units"] = float(group.units.median())
        for model in r4.MODELS:
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
                    entry[model][f"{stratum}_recovered"] = bool((g >= 2.0).mean() >= 0.8)
        out[cid] = entry
    return out


def run_simulation(stage: str, draw_fn, conditions, extra_args=(), draws: int = 100, jobs: int = 8):
    from joblib import Parallel, delayed

    started = provenance()
    spec = rb.load_spec()
    s = setup(spec)
    args = [a(spec, s) if callable(a) else a for a in extra_args]
    rows = Parallel(n_jobs=jobs, verbose=5)(
        delayed(draw_fn)(42 + i, spec, s, conditions, *args) for i in range(draws)
    )
    frame = pd.DataFrame([r for block in rows for r in block])
    summary = summarise(frame, conditions)
    directory = new_run_directory(OUTPUT / stage)
    frame.to_csv(directory / "draws.csv", index=False)
    return directory, frame, summary, started, s


def pneo_stage(draws=100, jobs=8) -> Path:
    conditions = pneo_conditions()
    directory, frame, summary, started, s = run_simulation(
        "pneo", pneo_draw, conditions, (), draws, jobs
    )
    paired = {}
    for cid, entry in summary.items():
        if entry["version"] != "V1":
            continue
        twin = cid.replace("V1|", "V0|", 1)
        if twin in summary:
            paired[cid.replace("V1|", "", 1)] = {
                m: round(entry[m]["share_ge_2pp"] - summary[twin][m]["share_ge_2pp"], 3)
                for m in r4.MODELS
                if m in entry and m in summary[twin]
            }
    units = {
        f"{v}|{scale}": float(frame[frame.condition.str.startswith(f"{v}|{scale}|")].units.median())
        for v in ("V0", "V1")
        for scale in ("native", "grid")
    }
    rv.write(
        directory / "pneo.json",
        {"units": units, "paired_v1_minus_v0_share": paired, "conditions": summary},
    )
    return finish("pneo", directory, ["draws.csv", "pneo.json"], started)


def s2_stage(draws=100, jobs=8) -> Path:
    conditions = s2_conditions()
    directory, frame, summary, started, s = run_simulation(
        "s2", s2_draw, conditions, (cell_levels,), draws, jobs
    )
    gains = {}
    for cid, entry in summary.items():
        if entry["config"] != "S3_combined_common":
            continue
        rest = cid.split("|", 1)[1]
        s1, s2_, s4 = (
            summary.get(f"{k}|{rest}")
            for k in ("S1_2jul_common", "S2_9jul_common", "S4_combined_union")
        )
        gains[rest] = {
            m: {
                "variation_gain": round(
                    entry[m]["share_ge_2pp"] - max(s1[m]["share_ge_2pp"], s2_[m]["share_ge_2pp"]), 3
                ),
                "coverage_gain": round(s4[m]["share_ge_2pp"] - entry[m]["share_ge_2pp"], 3),
            }
            for m in r4.MODELS
            if m in entry and s1 and s2_ and s4 and m in s1 and m in s4
        }
    units = {
        k: float(frame[frame.condition.str.startswith(f"{k}|")].units.median()) for k in CONFIGS
    }
    rv.write(directory / "s2.json", {"units": units, "gains": gains, "conditions": summary})
    return finish("s2", directory, ["draws.csv", "s2.json"], started)


# ---------------------------------------------------------------------------------------
# A1: inject before correction and run ACOLITE
# ---------------------------------------------------------------------------------------


def pattern(spec: dict, s: dict, product_grid) -> list:
    """Draw 0, bare, R = 10 m, no position error, on the product grid, no PSF."""
    permuted, strata, bearings, _ = permutation(42, s)
    signals = Signals(s, permuted, strata, bearings)
    return signals.delta(
        s["scenes"]["pneo"],
        s["endmembers"]["pneo"],
        "bare",
        10.0,
        0.0,
        psf=False,
        grid=product_grid,
    )


def check_stage() -> Path:
    import time
    from types import SimpleNamespace

    import rasterio
    from kelpobserve.coastal.acolite import correct_scene

    started = provenance()
    spec = rb.load_spec()
    s = setup(spec)
    config = spec["acolite"]
    rb.verify_acolite(config)
    source = spec["sources"]["pneo"]["july"]
    original = rb.extract_pneo(source, rb.ROOT / config["pneo_extract_root"])
    copy = PNEO_R6 / original.name
    relative = source["ms_rgb"].split("/", 1)[1]
    rgb = copy / relative
    with rasterio.open(original / relative) as src:
        dn = src.read().astype("float64")
        profile = {
            "driver": "GTiff",
            "height": src.height,
            "width": src.width,
            "count": src.count,
            "dtype": "uint16",
            "crs": src.crs,
            "transform": src.transform,
            "compress": "deflate",
            "tiled": True,
        }
        grid = SimpleNamespace(transform=src.transform, shape=(src.height, src.width), crs=src.crs)
    delta = pattern(spec, s, grid)  # blue, green, red
    injected = dn.copy()
    for role, band in PNEO_ROLE_BAND.items():
        d = delta[rv.RGB.index(role)]
        injected[band - 1] = np.where(
            dn[band - 1] > 0, np.clip(np.round(dn[band - 1] + d * 10000), 1, 65535), 0
        )
    injected = injected.astype("uint16")
    if copy.exists():  # an earlier run: reuse only if it holds exactly this injection
        with rasterio.open(rgb) as src:
            rb.require(np.array_equal(src.read(), injected), f"{rgb} holds a different injection.")
        reused = True
    else:
        shutil.copytree(original, copy)
        rgb.unlink()
        with rasterio.open(rgb, "w", **profile) as dst:  # GeoTIFF content under the original name
            dst.write(injected)
        reused = False
    settings = {name: rb.REGISTRATION / path for name, path in config["settings"].items()}
    overrides = rb.read_settings_file(settings["pleiades_neo"])
    limit = next(x for x in rb.prepare_scenes(spec) if x["id"] == "pneo-20260708")["limit"]
    runs = {}
    for version, extra in (("V0", overrides), ("V1", overrides | rv.GLINT_VARIANT["pleiades_neo"])):
        output = ACOLITE_R6 / f"{version}-injected" / "pneo-20260708"
        began = time.monotonic()
        correct_scene(
            copy,
            output,
            limit,
            Path(config["path"]),
            template=settings["common"],
            extra=extra,
            python_executable=config["python"],
        )
        runs[version] = {
            "output": str(output),
            "elapsed_s": round(time.monotonic() - began, 1),
            "effective_glint": rv.effective_settings(output),
        }
    result = compare_check(spec, s, delta, grid, runs)
    directory = new_run_directory(OUTPUT / "check")
    rv.write(directory / "check.json", {"runs": runs, "injection_reused": reused, **result})
    return finish("check", directory, ["check.json"], started)


def compare_check(spec, s, delta, grid, runs) -> dict:
    import rasterio
    from rasterio.warp import Resampling

    base_dirs = {
        "V0": rv.scene_dir(spec, "pneo-20260708", "V0"),
        "V1": rv.scene_dir(spec, "pneo-20260708", "V1"),
    }
    lidar = rv.raster_on(
        spec["sources"]["lidar"], grid.transform, grid.shape, grid.crs, Resampling.nearest
    )
    response, aot = {}, {}
    for version in ("V0", "V1"):
        injected_dir = Path(runs[version]["output"])
        response[version] = {}
        for role in ("blue", "green", "red"):
            nm = spec["scene_qc"]["bands_nm"]["pleiades_neo"][role]
            _, a = rb.nearest_band(rb.acolite_rasters(injected_dir, "rhos"), nm)
            _, b = rb.nearest_band(rb.acolite_rasters(base_dirs[version], "rhos"), nm)
            with rasterio.open(a) as x, rasterio.open(b) as y:
                rb.require(
                    x.transform.almost_equals(grid.transform, 1e-6)
                    and y.transform.almost_equals(grid.transform, 1e-6)
                    and x.shape == grid.shape == y.shape,
                    "Grids differ.",
                )
                response[version][role] = x.read(1).astype("float64") - y.read(1).astype("float64")
                aot[version] = {
                    "injected": float(x.tags().get("NC_GLOBAL#ac_aot_550", "nan")),
                    "original": float(y.tags().get("NC_GLOBAL#ac_aot_550", "nan")),
                }
    out = {"aot": aot, "bands": {}}
    strata = {"shallow": (lidar >= 4) & (lidar < 11), "deep": (lidar >= 13) & (lidar < 21)}
    for role in ("blue", "green", "red"):
        injected = delta[rv.RGB.index(role)]
        signal = np.abs(injected) >= 0.0005
        outside = (injected == 0) & np.isfinite(response["V0"][role])
        entry = {
            "determinism_max_abs_outside_V0": float(
                np.nanmax(np.abs(response["V0"][role][outside]))
            ),
            "determinism_max_abs_outside_V1": float(
                np.nanmax(np.abs(response["V1"][role][outside]))
            ),
        }
        for name, stratum in strata.items():
            sel = (
                signal
                & stratum
                & np.isfinite(response["V0"][role])
                & np.isfinite(response["V1"][role])
            )
            if sel.sum() < 20:
                entry[name] = {"pixels": int(sel.sum())}
                continue
            r0, r1, inj = response["V0"][role][sel], response["V1"][role][sel], injected[sel]
            slope = float(np.polyfit(r0, r1, 1)[0])
            entry[name] = {
                "pixels": int(sel.sum()),
                "slope_v1_on_v0": round(slope, 4),
                "r_v1_v0": round(float(np.corrcoef(r0, r1)[0, 1]), 5),
                "slope_v0_on_injected": round(float(np.polyfit(inj, r0, 1)[0]), 4),
                "median_injected": round(float(np.median(inj)), 5),
            }
        out["bands"][role] = entry
    ok = all(
        0.95 <= out["bands"][r][st].get("slope_v1_on_v0", 0) <= 1.05
        and out["bands"][r][st].get("r_v1_v0", 0) >= 0.99
        for r in ("blue", "green")
        for st in ("shallow", "deep")
    ) and all(abs(aot[v]["injected"] - aot[v]["original"]) < 0.001 for v in aot)
    out["criterion_met"] = bool(ok)
    return out


# ---------------------------------------------------------------------------------------
# C: background bias
# ---------------------------------------------------------------------------------------


def unit_sd(values: np.ndarray) -> float:
    """Robust SD for a stratum's few units (7-15 here); r5's helper needs 20."""
    values = values[np.isfinite(values)]
    if values.size < 5:
        return math.nan
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def bias_stage() -> Path:
    started = provenance()
    spec = rb.load_spec()
    s = setup(spec)
    events, iops = s["events"], s["iops"]["median"]
    cases = {
        "pneo_V0_native": (s["pneo"]["V0"], "native", "pneo"),
        "pneo_V1_native": (s["pneo"]["V1"], "native", "pneo"),
        "pneo_V0_10m": (s["pneo"]["V0"], "grid", "pneo"),
        "pneo_V1_10m": (s["pneo"]["V1"], "grid", "pneo"),
        "s2_2jul": (s["scenes"]["s2_2jul"], "grid", "s2_2jul"),
        "s2_9jul": (s["scenes"]["s2_9jul"], "grid", "s2_9jul"),
    }
    out = {}
    for name, (scene, scale, key) in cases.items():
        members = s["endmembers"][key]
        maker = rv.native_units if scale == "native" else rv.grid_units
        units = maker(spec, events, scene.window, scene.level, "july").dropna(
            subset=["blue", "green", "depth"]
        )
        sza = scene.window.solar_zenith_deg
        entry = {"units": int(len(units))}
        for background in ("understory", "bare"):
            for b, role in enumerate(("blue", "green")):
                model = r4.surface_rhos(
                    members[background][b], units.depth.to_numpy(), iops["a"][b], iops["bb"][b], sza
                )
                gap = units[role].to_numpy() - model
                rho = r4.bottom_reflectance(
                    units[role].to_numpy(), units.depth.to_numpy(), iops["a"][b], iops["bb"][b], sza
                )
                for stratum in ("shallow", "deep"):
                    sel = (units.stratum == stratum).to_numpy()
                    entry.setdefault(stratum, {}).setdefault(background, {})[role] = {
                        "rhos_gap_median": round(float(np.median(gap[sel])), 5),
                        "rhos_gap_robust_sd": round(unit_sd(gap[sel]), 5),
                        "bottom_bias_median": round(
                            float(np.median(rho[sel] - members[background][b])), 4
                        ),
                        "bottom_bias_robust_sd": round(
                            unit_sd(rho[sel] - members[background][b]), 4
                        ),
                        "endmember": round(float(members[background][b]), 4),
                    }
        out[name] = entry
    directory = new_run_directory(OUTPUT / "bias")
    rv.write(directory / "bias.json", out)
    return finish("bias", directory, ["bias.json"], started)


# ---------------------------------------------------------------------------------------
# D: matched spatial windows
# ---------------------------------------------------------------------------------------


def window_stats(f: r5.Field, size_m: float) -> dict:
    """Broad (x-y plane) and local (residual) spread per class over square windows."""
    size = max(3, int(round(size_m / f.pixel_m)))
    rows, cols = np.indices(f.valid.shape)
    classes = {"offshore": (f.regions["offshore"], 0.5)} | {
        b: (f.regions[b], 0.25) for b in r5.BINS
    }
    out = {}
    for name, (mask, minimum) in classes.items():
        broad, local, local_nir, means = [], [], [], []
        for r0 in range(0, mask.shape[0] - size + 1, size):
            for c0 in range(0, mask.shape[1] - size + 1, size):
                sel = mask[r0 : r0 + size, c0 : c0 + size]
                if sel.sum() < minimum * size * size or sel.sum() < 10:
                    continue
                y = f.bands["green"][r0 : r0 + size, c0 : c0 + size][sel]
                x, yy = (
                    cols[r0 : r0 + size, c0 : c0 + size][sel] * f.pixel_m,
                    rows[r0 : r0 + size, c0 : c0 + size][sel] * f.pixel_m,
                )
                columns = [x - x.mean(), yy - yy.mean()]
                if name != "offshore":
                    d = f.depth[r0 : r0 + size, c0 : c0 + size][sel]
                    if not np.isfinite(d).all():
                        continue
                    columns.append(d)
                coef, residual, _ = r5.least_squares(y, columns)
                broad.append(float(np.std(coef[1] * columns[0] + coef[2] * columns[1])))
                local.append(r5.robust_sd(residual))
                nir = f.bands["nir"][r0 : r0 + size, c0 : c0 + size][sel]
                if np.nanstd(nir) > 0:
                    _, residual_nir, _ = r5.least_squares(y, [*columns, nir])
                    local_nir.append(r5.robust_sd(residual_nir))
                means.append(float(np.mean(y)))
        if not local:
            continue

        def q(v):
            v = np.array([x for x in v if x == x])
            return (
                {
                    "median": float(np.median(v)),
                    "q25": float(np.percentile(v, 25)),
                    "q75": float(np.percentile(v, 75)),
                }
                if v.size
                else None
            )

        out[name] = {
            "windows": len(local),
            "broad": q(broad),
            "local": q(local),
            "local_after_nir": q(local_nir),
            "between_window_sd": float(np.std(means)),
        }
    return out


def windows_stage() -> Path:
    started = provenance()
    spec = rb.load_spec()
    registration = json.loads((rb.REGISTRATION / "revision-r5.json").read_text())
    mask = json.loads((resolve_run(rb.OUTPUT / "r3" / "mask") / "mask.json").read_text())
    diagnosis = json.loads(
        (resolve_run(rb.OUTPUT / "r5" / "diagnosis") / "diagnosis.json").read_text()
    )
    bounds = r3.study_bounds(spec)
    plan = [(sid, "V1") for sid in registration["scenes"]["sentinel2"]["V1"]]
    plan += [
        (sid, v) for sid in registration["scenes"]["pleiades_neo"]["dates"] for v in ("V0", "V1")
    ]
    results = {}
    for scene_id, version in plan:
        print(f"  {scene_id} {version}", flush=True)
        window = rv.read_window(spec, scene_id, bounds, version)
        if scene_id.startswith("s2"):
            window = window.with_threshold(mask["scenes"][scene_id]["threshold"])
            fields = {"s2_10m": r5.pixel_field(spec, window, scene_id)}
        else:
            fields = {
                "pneo_native": r5.pneo_native_field(spec, window, scene_id),
                "pneo_10m": r5.grid_field(spec, window, scene_id, bounds),
            }
        for unit, f in fields.items():
            for size in (100.0, 300.0):
                results[f"{scene_id}|{version}|{unit}|{size:g}m"] = window_stats(f, size)
        del window, fields
    decisions = {"D1m": {}, "D3m": {}}
    for key, classes in results.items():
        off = classes.get("offshore")
        if not off:
            continue
        for b in r5.BINS:
            if b in classes:
                ratio = off["local"]["median"] / classes[b]["local"]["median"]
                decisions["D1m"][f"{key}|{b}"] = {
                    "offshore_over_bin": round(ratio, 3),
                    "met": bool(ratio >= 0.8),
                }
    for pair_key, pair in diagnosis["pairs"].items():
        if pair["days"] > 8:
            continue
        for unit, regions in pair["m5"].items():
            for size in (100.0, 300.0):
                for b in r5.BINS:
                    if b not in regions:
                        continue
                    sizes = []
                    for scene in pair["scenes"]:
                        cls = results.get(f"{scene}|{unit}|{size:g}m", {})
                        if (
                            b in cls
                            and "offshore" in cls
                            and cls[b]["local_after_nir"]
                            and cls["offshore"]["local_after_nir"]
                        ):
                            sizes.append(
                                cls[b]["local_after_nir"]["median"]
                                / cls["offshore"]["local_after_nir"]["median"]
                            )
                    r_bin = regions[b]["green"]["r"]
                    r_off = regions.get("offshore", {}).get("green", {}).get("r", math.nan)
                    decisions["D3m"][f"{pair_key}|{unit}|{size:g}m|{b}"] = {
                        "size_ratios": sizes,
                        "r_bin": r_bin,
                        "r_offshore": r_off,
                        "met": bool(
                            len(sizes) == 2
                            and all(x >= 1.5 for x in sizes)
                            and r_bin >= 0.3
                            and (r_off != r_off or r_off < 0.2)
                        ),
                    }
    directory = new_run_directory(OUTPUT / "windows")
    rv.write(directory / "windows.json", {"windows": results, "decisions": decisions})
    return finish("windows", directory, ["windows.json"], started)


# ---------------------------------------------------------------------------------------
# E: offset control (amendment E_offset_control)
# ---------------------------------------------------------------------------------------


def to_sensor_valid(values: np.ndarray, fine, window) -> np.ndarray:
    """Valid-area mean of a fine field on the sensor grid, then the PSF, normalised over
    defined cells."""
    from rasterio.warp import Resampling, reproject
    from scipy.ndimage import gaussian_filter

    out = np.full(window.shape, np.nan, dtype="float32")
    reproject(
        values.astype("float32"),
        out,
        src_transform=fine.transform,
        src_crs="EPSG:32629",
        dst_transform=window.transform,
        dst_crs=window.crs,
        resampling=Resampling.average,
        src_nodata=np.nan,
        dst_nodata=np.nan,
    )
    defined = np.isfinite(out)
    weight = gaussian_filter(defined.astype(float), r4.PSF_SIGMA_PX, mode="nearest")
    blurred = gaussian_filter(np.where(defined, out, 0.0), r4.PSF_SIGMA_PX, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(defined & (weight > 0), blurred / weight, np.nan)


def model_background(s: dict, scene, substrate: str) -> list:
    """The injection's physics for a bare background of one substrate, on the scene's grid."""
    members, iops = s["endmembers"][scene_key(s, scene)], s["iops"]["median"]
    sza = scene.window.solar_zenith_deg
    return [
        to_sensor_valid(
            r4.surface_rhos(
                np.full(s["fine"].shape, members[substrate][b]),
                scene.fine_depth,
                iops["a"][b],
                iops["bb"][b],
                sza,
            ),
            s["fine"],
            scene.window,
        )
        for b in range(3)
    ]


def scene_key(s: dict, scene) -> str:
    for key, candidate in s["scenes"].items():
        if candidate.name == scene.name.replace("-V1", ""):
            return key
    raise KeyError(scene.name)


def offsets_from(real: pd.DataFrame, model: pd.DataFrame) -> dict:
    """Median real-minus-model residual per stratum and band (the b_s of the registration)."""
    merged = real.merge(model[["unit_id", *rv.RGB]], on="unit_id", suffixes=("", "_model"))
    out = {}
    for stratum, group in merged.groupby("stratum"):
        out[str(stratum)] = {}
        for role in rv.RGB:
            e = (group[role] - group[f"{role}_model"]).to_numpy()
            e = e[np.isfinite(e)]
            out[str(stratum)][role] = float(np.median(e)) if e.size else math.nan
    return out


def apply_offsets(units: pd.DataFrame, offsets: dict) -> pd.DataFrame:
    """y* = y - b_s per unit, by the unit's stratum."""
    units = units.copy()
    for role in rv.RGB:
        shift = units.stratum.map(lambda s_: offsets.get(str(s_), {}).get(role, 0.0))
        units[role] = units[role] - shift
    return units


def s2_setting(s: dict, config: str, dated: list[list], levels=None):
    """Window, bands, level and sun zenith for one Sentinel-2 configuration (as in s2_draw)."""
    d2, d9 = s["scenes"]["s2_2jul"], s["scenes"]["s2_9jul"]
    common = d2.window.valid() & d9.window.valid()
    if config == "S1_2jul_common":
        return restrict(d2.window, common), dated[0], d2.level, d2.window.solar_zenith_deg
    if config == "S2_9jul_common":
        return restrict(d9.window, common), dated[1], d9.level, d9.window.solar_zenith_deg
    union = config == "S4_combined_union"
    bands, mask = combine([d2.window, d9.window], dated, require_all=not union)
    window = replace(
        d2.window,
        mask_ok=mask,
        water=d2.window.water | d9.window.water,
        bands=dict(d2.window.bands) | {r: bands[b] for b, r in enumerate(rv.RGB)},
    )
    level = 0.0 if union else (d2.level + d9.level) / 2
    return window, bands, level, (d2.window.solar_zenith_deg + d9.window.solar_zenith_deg) / 2


def setting_units(spec, s, setting: str, bands_override, levels, synthetic=None) -> tuple:
    """Units of one r6 setting from given per-date (or single) bands; also returns sza."""
    events = s["events"]
    # Offsets never read cover: without synthetic targets, units get an empty placeholder.
    target = (
        synthetic
        if synthetic is not None
        else pd.DataFrame(
            {"event_id": events.event_id, "cover_pct": np.nan, "cover_status": "none"}
        )
    )
    if setting.startswith("pneo"):
        _, version, scale = setting.split("|")
        scene = s["pneo"][version]
        units = units_for(
            spec,
            events,
            scene.window,
            scene.level,
            scale,
            with_nir(bands_override, scene.window),
            target,
        )
        return units, scene.window.solar_zenith_deg
    config = setting.split("|")[1]
    window, bands, level, sza = s2_setting(s, config, bands_override)
    units = units_for(spec, events, window, level, "grid", with_nir(bands, window), target)
    if config == "S4_combined_union":
        mean_level = (s["scenes"]["s2_2jul"].level + s["scenes"]["s2_9jul"].level) / 2
        units = units.assign(depth=units.depth + units.unit_id.map(levels).fillna(mean_level))
    return units, sza


SETTINGS = [f"pneo|{v}|{sc}" for v in ("V0", "V1") for sc in ("native", "grid")] + [
    f"s2|{c}" for c in CONFIGS
]


def all_offsets(spec: dict, s: dict, levels) -> dict:
    """b_s per setting and substrate, from the uninjected background. No cover is read."""
    out = {}
    for setting in SETTINGS:
        for substrate in ("understory", "bare"):
            if setting.startswith("pneo"):
                scene = s["pneo"][setting.split("|")[1]]
                real_bands = [scene.window.bands[r] for r in rv.RGB]
                model_bands = model_background(s, scene, substrate)
            else:
                real_bands = [
                    [s["scenes"][k].window.bands[r] for r in rv.RGB] for k in ("s2_2jul", "s2_9jul")
                ]
                model_bands = [
                    model_background(s, s["scenes"][k], substrate) for k in ("s2_2jul", "s2_9jul")
                ]
            real, _ = setting_units(spec, s, setting, real_bands, levels)
            model, _ = setting_units(spec, s, setting, model_bands, levels)
            out[f"{setting}|{substrate}"] = offsets_from(real, model)
    return out


def offset_conditions() -> list[dict]:
    out = []
    for setting in SETTINGS:
        for contrast in ("understory", "bare"):
            for radius in r4.RADII:
                for error in r4.ERRORS:
                    out.append(
                        dict(
                            setting=setting,
                            contrast=contrast,
                            substrate=contrast,
                            radius=radius,
                            error=error,
                        )
                    )
            out.append(
                dict(setting=setting, contrast=None, substrate=contrast, radius=None, error=0.0)
            )
    for c in out:
        c["id"] = "|".join(
            str(c[k]) for k in ("setting", "contrast", "substrate", "radius", "error")
        )
    return out


def offset_draw(
    seed: int, spec: dict, s: dict, conditions: list[dict], offsets: dict, levels
) -> list[dict]:
    permuted, strata, bearings, synthetic = permutation(seed, s)
    signals = Signals(s, permuted, strata, bearings)
    rows = []
    for c in conditions:
        setting = c["setting"]
        if setting.startswith("pneo"):
            scene = s["pneo"][setting.split("|")[1]]
            members = s["endmembers"]["pneo"]
            bands = [scene.window.bands[r] for r in rv.RGB]
            if c["contrast"] is not None:
                delta = signals.delta(
                    s["scenes"]["pneo"], members, c["contrast"], c["radius"], c["error"]
                )
                bands = [bands[b] + delta[b] for b in range(3)]
        else:
            members = s["endmembers"]["s2_2jul"]
            bands = []
            for key in ("s2_2jul", "s2_9jul"):
                scene = s["scenes"][key]
                dated = [scene.window.bands[r] for r in rv.RGB]
                if c["contrast"] is not None:
                    delta = signals.delta(scene, members, c["contrast"], c["radius"], c["error"])
                    dated = [dated[b] + delta[b] for b in range(3)]
                bands.append(dated)
        units, sza = setting_units(spec, s, setting, bands, levels, synthetic)
        units = apply_offsets(units, offsets[f"{setting}|{c['substrate']}"])
        background = members[c["contrast"]] if c["contrast"] else None
        rows.append(
            {"draw": seed, "condition": c["id"]}
            | r4.evaluate(units, spec, members["kelp"], background, s["iops"], sza)
        )
    return rows


def original_shares() -> dict:
    """r6's A and B shares keyed like the offset conditions, for pairing."""
    out = {}
    pneo = json.loads((resolve_run(OUTPUT / "pneo") / "pneo.json").read_text())["conditions"]
    for e in pneo.values():
        if e["background"] == "real" and e["contrast"] is not None:
            out[
                f"pneo|{e['version']}|{e['scale']}|{e['contrast']}|{e['contrast']}|{e['radius']}|{e['error']}"
            ] = e
    s2 = json.loads((resolve_run(OUTPUT / "s2") / "s2.json").read_text())["conditions"]
    for e in s2.values():
        if e["background"] == "real" and e["contrast"] is not None:
            out[f"s2|{e['config']}|{e['contrast']}|{e['contrast']}|{e['radius']}|{e['error']}"] = e
    return out


def classify(after: dict, before: dict, nulls: dict) -> dict:
    """The registered outcome per setting, contrast and model, on the reference set."""
    out = {}
    for setting in SETTINGS:
        for contrast in ("understory", "bare"):
            for model in r4.MODELS:
                keys = [f"{setting}|{contrast}|{contrast}|{r}|0.0" for r in (5.0, 10.0, 25.0)]
                pairs = [(after.get(k, {}).get(model), before.get(k, {}).get(model)) for k in keys]
                pairs = [(a, b) for a, b in pairs if a and b]
                if not pairs:
                    continue
                null = (
                    nulls.get(f"{setting}|None|{contrast}|None|0.0", {})
                    .get(model, {})
                    .get("share_ge_2pp", math.nan)
                )
                appears = (
                    any(a["share_ge_2pp"] >= 0.8 and not b["recovered"] for a, b in pairs)
                    and null < 0.2
                )
                increase = float(
                    np.median([a["share_ge_2pp"] - b["share_ge_2pp"] for a, b in pairs])
                )
                outcome = (
                    "recovery_appears"
                    if appears
                    else ("improves_but_fails" if increase >= 0.10 else "barely_changes")
                )
                out[f"{setting}|{contrast}|{model}"] = {
                    "outcome": outcome,
                    "median_increase": round(increase, 3),
                    "after": [round(a["share_ge_2pp"], 3) for a, _ in pairs],
                    "before": [round(b["share_ge_2pp"], 3) for _, b in pairs],
                    "null_after": null,
                }
    return out


def offset_stage(draws=100, jobs=8) -> Path:
    from joblib import Parallel, delayed

    started = provenance()
    spec = rb.load_spec()
    s = setup(spec)
    levels = cell_levels(spec, s)
    offsets = all_offsets(spec, s, levels)
    conditions = offset_conditions()
    rows = Parallel(n_jobs=jobs, verbose=5)(
        delayed(offset_draw)(42 + i, spec, s, conditions, offsets, levels) for i in range(draws)
    )
    frame = pd.DataFrame([r for block in rows for r in block])
    by_id = {c["id"]: c for c in conditions}
    summary = {}
    for cid, group in frame.groupby("condition"):
        entry = dict(by_id[cid]) | {
            "draws": int(len(group)),
            "median_units": float(group.units.median()),
        }
        for model in r4.MODELS:
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
                    entry[model][f"{stratum}_share_ge_2pp"] = round(
                        float((group[a] - group[b] >= 2.0).mean()), 3
                    )
        summary[cid] = entry
    before = original_shares()
    nulls = {k: v for k, v in summary.items() if v["contrast"] is None}
    paired = {
        k: {
            m: round(v[m]["share_ge_2pp"] - before[k][m]["share_ge_2pp"], 3)
            for m in r4.MODELS
            if m in v and m in before[k]
        }
        for k, v in summary.items()
        if k in before
    }
    directory = new_run_directory(OUTPUT / "offset")
    frame.to_csv(directory / "draws.csv", index=False)
    rv.write(
        directory / "offset.json",
        {
            "offsets": offsets,
            "classification": classify(summary, before, nulls),
            "paired_after_minus_before": paired,
            "conditions": summary,
        },
    )
    return finish("offset", directory, ["draws.csv", "offset.json"], started)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("check", "bias", "windows"):
        commands.add_parser(name)
    for name in ("pneo", "s2", "offset"):
        p = commands.add_parser(name)
        p.add_argument("--draws", type=int, default=100)
        p.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    if args.command == "check":
        directory = check_stage()
    elif args.command == "pneo":
        directory = pneo_stage(args.draws, args.jobs)
    elif args.command == "s2":
        directory = s2_stage(args.draws, args.jobs)
    elif args.command == "offset":
        directory = offset_stage(args.draws, args.jobs)
    elif args.command == "bias":
        directory = bias_stage()
    else:
        directory = windows_stage()
    print(f"r6 {args.command} written to {directory}")


if __name__ == "__main__":
    main()
