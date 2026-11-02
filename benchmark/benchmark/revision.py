"""Exploratory re-analysis r1 of the Sesimbra habitat benchmark.

Registered in registration/revision-r1.json before any result below was computed. It
re-examines v1 after the 2026-10-02 review: QC applied consistently, label-free radiometry
variants, the spectral comparison with its registered sensitivities, an injection-and-
recovery test and multi-year campaign windows. Every cover value was already seen in v1,
so all of it is exploratory. v1 outputs, locks and reports are read, never written;
stages write under out/r1/ and .cache/acolite-r1/.

    PYTHONPATH=. python -m benchmark.revision screen      # Q1: water-only glint screen
    PYTHONPATH=. python -m benchmark.revision qc          # Q1: support, QAA flag classes
    PYTHONPATH=. python -m benchmark.revision acolite     # Q2: glint variants (slow)
    PYTHONPATH=. python -m benchmark.revision radiometry  # Q2: label-free measures
    PYTHONPATH=. python -m benchmark.revision models      # Q3 and its sensitivities
    PYTHONPATH=. python -m benchmark.revision injection   # injection and recovery
    PYTHONPATH=. python -m benchmark.revision windows     # Q4: campaign windows by year
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from dataclasses import dataclass
from datetime import UTC, date, datetime, timedelta
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from kelpobserve.coastal.provenance import complete_run, file_digest, new_run_directory, resolve_run

from benchmark import run_benchmark as rb

REVISION = rb.REGISTRATION / "revision-r1.json"
OUTPUT = rb.OUTPUT / "r1"
ACOLITE_R1 = rb.ROOT / ".cache" / "acolite-r1"
ROLES = ("blue", "green", "red", "nir")
RGB = ("blue", "green", "red")
CONTROL_FAMILIES = [
    {"id": "nir_ridge", "kind": "ridge", "predictors": ["nir"]},
    {"id": "nir_depth_ridge", "kind": "ridge", "predictors": ["nir", "depth"]},
]
PAIRS = [
    ("spectral_depth_ridge", "depth_ridge"),
    ("spectral_depth_forest", "depth_ridge"),
    ("spectral_ridge", "depth_ridge"),
    ("spectral_ridge", "training_mean"),
    ("depth_ridge", "training_mean"),
    ("nir_depth_ridge", "depth_ridge"),
    ("nir_ridge", "training_mean"),
]
COMPASS = {"N": 0, "NE": 45, "E": 90, "SE": 135, "S": 180, "SW": 225, "W": 270, "NW": 315}
SCENES = {
    ("pleiades_neo", "july"): "pneo-20260708",
    ("pleiades_neo", "september"): "pneo-20260903",
    ("sentinel2", "july"): "s2-20260702-S2B-R037",
    ("sentinel2", "september"): "s2-20260831-S2B-R037",
}
GLINT_VARIANT = {
    "sentinel2": {"dsf_residual_glint_correction": "True"},
    "pleiades_neo": {
        "dsf_residual_glint_correction": "True",
        "dsf_residual_glint_correction_method": "alternative",
        "dsf_residual_glint_wave_range": "800,900",
        "glint_mask_rhos_wave": "800",
        "glint_mask_rhos_threshold": "0.1",
    },
}


def load_revision() -> dict:
    return json.loads(REVISION.read_text())


def provenance() -> dict:
    """Registration, runner and library state, captured when a stage starts."""
    library = rb.library_root().parents[1]

    def git(*args):
        return subprocess.run(
            ["git", *args], cwd=library, capture_output=True, text=True
        ).stdout.strip()

    return {
        "registration": {"path": str(REVISION), "sha256": file_digest(REVISION)},
        "revision_runner_sha256": file_digest(Path(__file__)),
        "v1_runner_sha256": file_digest(Path(rb.__file__)),
        "benchmark_git": rb.git_state([Path(__file__), rb.REGISTRATION]),
        "library_git": {
            "commit": git("rev-parse", "HEAD"),
            "uncommitted_changes": bool(git("status", "--porcelain", "--", "src/coastobserve")),
        },
    }


def clean(value):
    """NaN and infinities become null: the v1 JSON writer refuses them."""
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write(path: Path, value) -> None:
    rb.write_json(path, clean(value))


def finish(stage: str, directory: Path, outputs: list[str], started: dict) -> Path:
    complete_run(
        OUTPUT / stage,
        directory,
        {
            "stage": f"r1-{stage}",
            "created_at": datetime.now(UTC).isoformat(),
            "exploratory": True,
            "outputs": {name: file_digest(directory / name) for name in outputs},
            **started,
        },
    )
    return directory


# ---------------------------------------------------------------------------------------
# Reading scenes under the registered valid-water definition
# ---------------------------------------------------------------------------------------


def sensor_of(scene_id: str) -> str:
    return "pleiades_neo" if scene_id.startswith("pneo") else "sentinel2"


def scene_dir(spec: dict, scene_id: str, variant: str = "V0") -> Path:
    if variant == "V0":
        return Path(rb.scene_directories(spec)[scene_id]["output_dir"])
    return ACOLITE_R1 / variant / scene_id


def mask_band(spec: dict, directory: Path) -> tuple[Path, float, float]:
    """rhot at ACOLITE's l2w_mask_wave and its threshold, from the pinned defaults (as v1)."""
    sample = next(iter(rb.acolite_rasters(directory, "rhos").values()))
    key = re.match(r"(.+?)_\d{4}_\d{2}_\d{2}_", sample.name).group(1)
    defaults = rb.acolite_defaults(Path(spec["acolite"]["path"]), key)
    nm, path = rb.nearest_band(
        rb.acolite_rasters(directory, "rhot"), float(defaults["l2w_mask_wave"]), tolerance_nm=None
    )
    return path, float(defaults["l2w_mask_threshold"]), nm


def raster_on(path, transform, shape, crs, resampling) -> np.ndarray:
    import rasterio
    from rasterio.warp import reproject

    out = np.full(shape, np.nan, dtype="float32")
    with rasterio.open(path) as src:
        reproject(
            rasterio.band(src, 1),
            out,
            dst_transform=transform,
            dst_crs=crs,
            resampling=resampling,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )
    return out


def water_on(spec: dict, transform, shape, crs) -> np.ndarray:
    """Lidar or DTM depth present, nearest-neighbour on the target grid (v1 scene-qc)."""
    from rasterio.warp import Resampling

    water = np.zeros(shape, dtype=bool)
    for key in ("lidar", "dtm"):
        water |= np.isfinite(
            raster_on(spec["sources"][key], transform, shape, crs, Resampling.nearest)
        )
    return water


@dataclass
class Window:
    scene_id: str
    sensor: str
    transform: object
    crs: object
    bands: dict
    mask_ok: np.ndarray
    water: np.ndarray
    lidar: np.ndarray
    solar_zenith_deg: float
    #: Top-of-atmosphere reflectance at the mask band, kept so r3 can set its own threshold.
    mask_rhot: np.ndarray | None = None

    @property
    def shape(self) -> tuple[int, int]:
        return self.mask_ok.shape

    def with_threshold(self, threshold: float) -> Window:
        """The same window with the cloud and land mask at another threshold."""
        from dataclasses import replace

        with np.errstate(invalid="ignore"):
            return replace(self, mask_ok=self.mask_rhot <= threshold)

    def valid(self, rule: str = "r1") -> np.ndarray:
        """r1: the registered valid-water definition. v1: finite RGB only, as v1's units."""
        finite = np.all([np.isfinite(self.bands[r]) for r in RGB], axis=0)
        return finite if rule == "v1" else self.water & finite & self.mask_ok


def read_clipped(path, window) -> np.ndarray:
    """A window read that may overhang the raster, NaN outside it.

    rasterio's boundless read goes through a VRT and segfaults in this environment, so
    the overlap is read and placed by hand.
    """
    import rasterio
    from rasterio.windows import Window as RasterWindow

    height, width = int(window.height), int(window.width)
    out = np.full((height, width), np.nan, dtype="float64")
    with rasterio.open(path) as src:
        row0, col0 = int(window.row_off), int(window.col_off)
        r0, c0 = max(row0, 0), max(col0, 0)
        r1, c1 = min(row0 + height, src.height), min(col0 + width, src.width)
        if r1 > r0 and c1 > c0:
            data = src.read(1, window=RasterWindow(c0, r0, c1 - c0, r1 - r0)).astype("float64")
            if src.nodata is not None and not np.isnan(src.nodata):
                data[data == src.nodata] = np.nan
            out[r0 - row0 : r1 - row0, c0 - col0 : c1 - col0] = data
    return out


def read_window(spec: dict, scene_id: str, bounds_utm, variant: str = "V0") -> Window:
    """Rhos at the four roles, the cloud/glint mask, water and lidar over EPSG:32629 bounds."""
    import rasterio
    from rasterio import windows
    from rasterio.warp import Resampling, transform_bounds

    directory = scene_dir(spec, scene_id, variant)
    sensor = sensor_of(scene_id)
    rhos = rb.acolite_rasters(directory, "rhos")
    roles = spec["scene_qc"]["bands_nm"][sensor]
    mask_path, threshold, _ = mask_band(spec, directory)
    sample = next(iter(rhos.values()))
    with rasterio.open(sample) as src:
        bounds = transform_bounds("EPSG:32629", src.crs, *bounds_utm)
        window = windows.from_bounds(*bounds, transform=src.transform).round_offsets()
        window = window.round_lengths()
        transform, crs = src.window_transform(window), src.crs
        sza = float(src.tags().get("NC_GLOBAL#sza", "nan"))

    def read(path):
        return read_clipped(path, window)

    bands = {role: read(rb.nearest_band(rhos, roles[role])[1]) for role in ROLES}
    mask_rhot = read(mask_path)
    with np.errstate(invalid="ignore"):
        mask_ok = mask_rhot <= threshold
    shape = mask_ok.shape
    return Window(
        scene_id,
        sensor,
        transform,
        crs,
        bands,
        mask_ok,
        water_on(spec, transform, shape, crs),
        raster_on(spec["sources"]["lidar"], transform, shape, crs, Resampling.nearest),
        sza,
        mask_rhot,
    )


# ---------------------------------------------------------------------------------------
# Model units
# ---------------------------------------------------------------------------------------


def campaign_events(spec: dict, campaign: str):
    """v1-eligible events of one campaign, in EPSG:32629."""
    events, _ = rb.eligible_events(spec)
    return events[events.campaign == campaign].to_crs(32629).reset_index(drop=True)


def event_bounds(events, margin: float) -> tuple[float, float, float, float]:
    x0, y0, x1, y1 = events.total_bounds
    return (x0 - margin, y0 - margin, x1 + margin, y1 + margin)


def positions(events, shift: np.ndarray | None = None) -> np.ndarray:
    xy = np.column_stack([events.geometry.x.to_numpy(), events.geometry.y.to_numpy()])
    return xy if shift is None else xy + shift


def lidar_at(spec: dict, xy: np.ndarray) -> np.ndarray:
    """Lidar depth of the containing lidar cell at EPSG:32629 points (v1 native rule)."""
    import rasterio
    from rasterio.warp import transform

    with rasterio.open(spec["sources"]["lidar"]) as src:
        lon, lat = transform("EPSG:32629", src.crs, xy[:, 0].tolist(), xy[:, 1].tolist())
        values = np.array([v[0] for v in src.sample(zip(lon, lat))], dtype=float)
        if src.nodata is not None and not np.isnan(src.nodata):
            values[values == src.nodata] = np.nan
    return values


def _reason(window: Window, region, rule: str) -> str:
    if rule == "v1":
        return "nonfinite"
    if not window.water[region].all():
        return "nonwater"
    if not window.mask_ok[region].all():
        return "masked"
    return "nonfinite"


def native_units(
    spec: dict,
    events,
    window: Window,
    level: float,
    campaign: str,
    *,
    size: int = 3,
    shift: np.ndarray | None = None,
    rule: str = "r1",
    bands: dict | None = None,
) -> pd.DataFrame:
    """One unit per quadrat: mean of a size x size window, kept only if every pixel is valid."""
    bands = bands or window.bands
    valid = window.valid(rule)
    half = size // 2
    inverse = ~window.transform
    xy = positions(events, shift)
    depth = lidar_at(spec, xy)
    rows = []
    for (x, y), z, (_, event) in zip(xy, depth, events.iterrows()):
        col, row_f = inverse * (x, y)
        r, c = int(math.floor(row_f)), int(math.floor(col))
        row = {
            "unit_id": event.event_id,
            "campaign": campaign,
            "transect": int(event.transect),
            "stratum": event.stratum,
            "n_quadrats": 1,
            "event_ids": json.dumps([event.event_id]),
            "depth": float(z) + level if np.isfinite(z) else math.nan,
        }
        r0, r1, c0, c1 = r - half, r + half + 1, c - half, c + half + 1
        if r0 < 0 or c0 < 0 or r1 > window.shape[0] or c1 > window.shape[1]:
            row |= {role: math.nan for role in ROLES} | {"reason": "outside", "valid_pixels": 0}
        else:
            region = (slice(r0, r1), slice(c0, c1))
            ok = bool(valid[region].all())
            row |= {role: float(bands[role][region].mean()) if ok else math.nan for role in ROLES}
            row["reason"] = "ok" if ok else _reason(window, region, rule)
            row["valid_pixels"] = int(valid[region].sum())
        rows.append(row)
    return pd.DataFrame(rows)


def _average(values: np.ndarray, window: Window, transform, shape, crs, nodata) -> np.ndarray:
    from rasterio.warp import Resampling, reproject

    out = np.full(shape, np.nan if nodata is not None else 0.0, dtype="float32")
    reproject(
        values.astype("float32"),
        out,
        src_transform=window.transform,
        src_crs=window.crs,
        dst_transform=transform,
        dst_crs=crs,
        resampling=Resampling.average,
        src_nodata=nodata,
        dst_nodata=nodata,
    )
    return out


def grid_units(
    spec: dict,
    events,
    window: Window,
    level: float,
    campaign: str,
    *,
    shift: np.ndarray | None = None,
    rule: str = "r1",
    bands: dict | None = None,
) -> pd.DataFrame:
    """One unit per 10 m cell holding quadrats: valid-area mean, kept at >= 0.8 valid."""
    import rasterio
    from rasterio.warp import Resampling

    bands = bands or window.bands
    minimum = spec["common_grid"]["min_valid_fraction"]
    transform, shape, crs = rb.common_grid(
        spec, event_bounds(events, 100.0), spec["selected_grid_scene"]
    )
    valid = window.valid(rule)
    fraction = _average(valid.astype("float32"), window, transform, shape, crs, None)
    water = _average(window.water.astype("float32"), window, transform, shape, crs, None)
    means = {
        role: _average(np.where(valid, bands[role], np.nan), window, transform, shape, crs, np.nan)
        for role in ROLES
    }
    depth = raster_on(spec["sources"]["lidar"], transform, shape, crs, Resampling.average)
    cells: dict = {}
    for i, (x, y) in enumerate(positions(events, shift)):
        r, c = rasterio.transform.rowcol(transform, x, y)
        cells.setdefault((int(r), int(c)), []).append(i)
    rows = []
    for (r, c), members in sorted(cells.items()):
        group = events.iloc[members]
        inside = 0 <= r < shape[0] and 0 <= c < shape[1]
        ok = inside and fraction[r, c] >= minimum
        strata, transects = set(group.stratum), {int(t) for t in group.transect}
        if not inside:
            reason = "outside"
        elif ok:
            reason = "ok"
        elif rule == "r1" and water[r, c] < minimum:
            reason = "nonwater"
        else:
            reason = "masked" if rule == "r1" else "nonfinite"
        rows.append(
            {
                "unit_id": f"{campaign}-cell-{r}-{c}",
                "campaign": campaign,
                "row": r,
                "column": c,
                "transect": sorted(transects)[0],
                "mixed_transects": len(transects) > 1,
                "stratum": sorted(strata)[0],
                "mixed_strata": len(strata) > 1,
                "n_quadrats": len(members),
                "event_ids": json.dumps(list(group.event_id)),
                "valid_fraction": round(float(fraction[r, c]), 4) if inside else 0.0,
                "water_fraction": round(float(water[r, c]), 4) if inside else 0.0,
                "depth": float(depth[r, c]) + level
                if inside and np.isfinite(depth[r, c])
                else math.nan,
                **{role: float(means[role][r, c]) if ok else math.nan for role in ROLES},
                "reason": reason,
            }
        )
    return pd.DataFrame(rows)


def targets_for(campaign: str) -> pd.DataFrame:
    """Cover targets. Opening September here is exploratory: v1 already evaluated it."""
    directory = resolve_run(rb.OUTPUT / "inventory")
    name = "targets_development.csv" if campaign == "july" else "targets_september.csv"
    return pd.read_csv(directory / name)


def reversed_targets(targets: pd.DataFrame) -> pd.DataFrame:
    """Quadrat k takes the cover of quadrat 6 - k in the same campaign, stratum and transect."""
    import geopandas as gpd

    events = gpd.read_file(resolve_run(rb.OUTPUT / "inventory") / "events.gpkg", layer="events")
    key = {
        e.event_id: (e.campaign, e.stratum, int(e.transect), int(e.quadrat))
        for e in events.itertuples()
    }
    by_key = {v: k for k, v in key.items()}
    source = targets.set_index("event_id")
    out = targets.copy()
    for i, event_id in enumerate(out.event_id):
        campaign, stratum, transect, quadrat = key[event_id]
        partner = by_key.get((campaign, stratum, transect, 6 - quadrat))
        for column in ("cover_pct", "cover_status"):
            out.loc[out.index[i], column] = (
                source.loc[partner, column] if partner in source.index else np.nan
            )
    return out


# ---------------------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------------------


def families(spec: dict, controls: bool = True) -> list[dict]:
    return [*spec["models"]["families"], *(CONTROL_FAMILIES if controls else [])]


def complete_units(units: pd.DataFrame) -> pd.DataFrame:
    """Every family is scored on the same units: all predictors and the target present."""
    return units.dropna(subset=[*ROLES, "depth", "cover_pct"])


def loto(units: pd.DataFrame, spec: dict, family_list: list[dict]):
    """v1 leave-one-transect-out, with v1's fold-sufficiency rule, on complete units."""
    design = spec["evaluation"]
    data = complete_units(units)
    predictions, folds = [], []
    for transect in sorted(data.transect.unique()):
        held, rest = data[data.transect == transect], data[data.transect != transect]
        per_stratum = rest.groupby("stratum").size()
        sufficient = bool(
            len(rest) >= design["min_training_units"]
            and len(per_stratum) > 1
            and (per_stratum >= design["min_training_units_per_stratum"]).all()
        )
        folds.append(
            {
                "held_out_transect": int(transect),
                "held_out_units": int(len(held)),
                "training_units_by_stratum": {k: int(v) for k, v in per_stratum.items()},
                "sufficient": sufficient,
            }
        )
        if not sufficient:
            continue
        for family in family_list:
            fitted = rb.fit_family(family, rest, spec)
            for unit_id, stratum, observed, predicted in zip(
                held.unit_id, held.stratum, held.cover_pct, rb.predict_family(fitted, held)
            ):
                predictions.append(
                    {
                        "family": family["id"],
                        "unit_id": unit_id,
                        "transect": int(transect),
                        "stratum": stratum,
                        "observed": float(observed),
                        "predicted": float(predicted),
                    }
                )
    return pd.DataFrame(predictions), folds


def transfer(train: pd.DataFrame, test: pd.DataFrame, spec: dict, family_list: list[dict]):
    """Fit on every complete training unit, apply once to the test units."""
    train, test = complete_units(train), complete_units(test)
    rows = []
    if len(train) < spec["evaluation"]["min_training_units"] or test.empty:
        return pd.DataFrame(rows)
    for family in family_list:
        fitted = rb.fit_family(family, train, spec)
        for unit_id, stratum, observed, predicted in zip(
            test.unit_id, test.stratum, test.cover_pct, rb.predict_family(fitted, test)
        ):
            rows.append(
                {
                    "family": family["id"],
                    "unit_id": unit_id,
                    "stratum": stratum,
                    "observed": float(observed),
                    "predicted": float(predicted),
                }
            )
    return pd.DataFrame(rows)


def summarise(frame: pd.DataFrame, family_list: list[dict]) -> dict:
    if frame.empty:
        return {"n": 0}
    scores = {}
    for family in family_list:
        rows = frame[frame.family == family["id"]]
        if rows.empty:
            continue
        entry = rb.score(rows.predicted.to_numpy(), rows.observed.to_numpy())
        entry["per_stratum_mae_pp"] = {
            str(s): round(float(np.mean(np.abs(g.predicted - g.observed))), 2)
            for s, g in rows.groupby("stratum")
        }
        if "transect" in rows:
            entry["per_transect_mae_pp"] = {
                str(int(t)): round(float(np.mean(np.abs(g.predicted - g.observed))), 2)
                for t, g in rows.groupby("transect")
            }
        scores[family["id"]] = entry
    paired = {}
    for a, b in PAIRS:
        first = frame[frame.family == a].set_index("unit_id")
        second = frame[frame.family == b].set_index("unit_id")
        shared = first.index.intersection(second.index)
        if shared.empty:
            continue
        difference = (first.loc[shared].predicted - first.loc[shared].observed).abs() - (
            second.loc[shared].predicted - second.loc[shared].observed
        ).abs()
        paired[f"{a} - {b}"] = {
            "n": int(len(shared)),
            "mean_pp": round(float(difference.mean()), 2),
            "median_pp": round(float(difference.median()), 2),
        }
    return {"scores": scores, "paired_mae_difference_pp": paired, "note": "negative is better"}


def support(units: pd.DataFrame) -> dict:
    """Units, kept units and quadrats by stratum, with exclusion reasons."""
    out = {}
    for stratum, group in units.groupby("stratum"):
        kept = group[group.reason == "ok"]
        out[str(stratum)] = {
            "units": int(len(group)),
            "kept": int(len(kept)),
            "quadrats": int(group.n_quadrats.sum()),
            "quadrats_kept": int(kept.n_quadrats.sum()),
            "excluded": {
                k: int(v) for k, v in group[group.reason != "ok"].reason.value_counts().items()
            },
            "excluded_units": list(group[group.reason != "ok"].unit_id),
        }
    return out


# ---------------------------------------------------------------------------------------
# QAA flag classes
# ---------------------------------------------------------------------------------------


def classify_flags(flags: list[str]) -> dict:
    """Registered classes: informational, band-limited (above 700 nm) and fatal."""
    classes = {"informational": [], "band_limited": [], "fatal": []}
    for flag in flags:
        if flag.startswith("a_cdm_443_not_computable"):
            classes["informational"].append(flag)
        elif flag.startswith("a_below_pure_water_floor"):
            nm = float(re.search(r"\(([\d.]+)nm", flag).group(1))
            classes["band_limited" if nm > 700 else "fatal"].append(flag)
        else:
            classes["fatal"].append(flag)
    return classes


def qaa_fits(
    spec: dict, scene_id: str, variant: str = "V0", threshold: float | None = None
) -> dict:
    """The four registered QAA variants over the PNeo reference region.

    ``threshold`` replaces ACOLITE's mask threshold for the masked variants (r3).
    """
    import geopandas as gpd
    import rasterio
    from kelpobserve.coastal.optics.qaa import fit_scene_iops
    from kelpobserve.coastal.optics.water import rhos_to_rrs
    from kelpobserve.coastal.scene import Scene
    from rasterio import features

    directory = scene_dir(spec, scene_id, variant)
    scene = Scene.from_acolite_dir(directory, sensor=sensor_of(scene_id))
    region = gpd.read_file(rb.REGISTRATION / spec["pneo_reference_region"]["path"])
    mask_path, default, _ = mask_band(spec, directory)
    threshold = default if threshold is None else threshold
    with rasterio.open(mask_path) as src:
        inside = features.rasterize(
            region.to_crs(src.crs).geometry, out_shape=src.shape, transform=src.transform
        ).astype(bool)
        with np.errstate(invalid="ignore"):
            clear = src.read(1) <= threshold
    rrs = rhos_to_rrs(scene.rhos) if scene.rhos is not None else scene.rrs_above
    finite = inside & np.all(np.isfinite(rrs), axis=0)
    waves = np.asarray(scene.wavelengths_nm, dtype=float)
    out = {}
    for name, pixels, reference in (
        ("auto_unmasked", finite, None),
        ("auto_masked", finite & clear, None),
        ("forced560_unmasked", finite, 560.0),
        ("forced560_masked", finite & clear, 560.0),
    ):
        try:
            iops = fit_scene_iops(
                rrs, waves, pixels, scene.solar_zenith_deg, force_ref_wavelength_nm=reference
            )
        except ValueError as error:  # too few pixels: the registered "unavailable"
            out[name] = {"pixels": int(pixels.sum()), "available": False, "error": str(error)}
            continue
        classes = classify_flags(list(iops.qa_flags))
        enough = iops.n_deep_pixels >= spec["reference_qc"]["min_pixels"]
        out[name] = {
            "pixels": int(iops.n_deep_pixels),
            "reference_nm": float(iops.ref_wavelength_nm),
            "bbp_555": round(float(iops.bbp_555), 5),
            "kd_490": round(float(iops.kd[int(np.argmin(np.abs(waves - 490)))]), 4),
            "mean_rhos_nir": round(
                float(scene.rhos[int(np.argmin(np.abs(waves - 830)))][pixels].mean()), 5
            )
            if scene.rhos is not None
            else None,
            "flags": classes,
            "available": bool(enough and not classes["fatal"]),
            "v1_rule_available": bool(enough and not iops.qa_flags),
        }
    return out


# ---------------------------------------------------------------------------------------
# Stage: water-only screen (Q1, D2)
# ---------------------------------------------------------------------------------------


def screen() -> Path:
    import planetary_computer as pc
    import rasterio
    from pystac_client import Client
    from rasterio.windows import from_bounds
    from scipy import ndimage

    started = provenance()
    spec = rb.load_spec()
    registered = json.loads((rb.REGISTRATION / "glint-screen.json").read_text())
    bounds = tuple(registered["bounds_epsg32629"])
    threshold = registered["screen_glint_max"]
    client = Client.open(rb.PLANETARY_STAC, modifier=pc.sign_inplace)
    import geopandas as gpd

    lon0, lat0, lon1, lat1 = (
        gpd.GeoSeries.from_xy([bounds[0], bounds[2]], [bounds[1], bounds[3]], crs=32629)
        .to_crs(4326)
        .total_bounds
    )
    rows = []
    for entry in sorted(registered["scenes"], key=lambda r: r["sensing_utc"]):
        day = entry["sensing_utc"][:10]
        items = client.search(
            collections=["sentinel-2-l2a"], bbox=[lon0, lat0, lon1, lat1], datetime=day
        ).items()
        item = next(
            (
                i
                for i in sorted(items, key=lambda i: i.datetime)
                if rb.screen_key(i.datetime, i.properties.get("platform", "")) == entry["key"]
            ),
            None,
        )
        row = {k: entry[k] for k in ("key", "campaign", "sensing_utc", "relative_orbit")}
        row["v1_b11_median"] = entry["b11_median"]
        if item is None:
            rows.append(row | {"error": "item not found"})
            continue
        offset = -1000 if str(item.properties.get("s2:processing_baseline", "0")) >= "04.00" else 0
        with rasterio.open(item.assets["B11"].href) as src:
            window = from_bounds(*bounds, transform=src.transform)
            data = src.read(1, window=window).astype(float)
            transform = src.window_transform(window)
            crs = src.crs
        data[data == 0] = np.nan
        reflectance = (data + offset) / 10000.0
        water = water_on(spec, transform, reflectance.shape, crs)
        interior = ndimage.binary_erosion(water, iterations=1, border_value=1)
        values = {
            "rectangle": reflectance[np.isfinite(reflectance)],
            "water": reflectance[water & np.isfinite(reflectance)],
            "interior_water": reflectance[interior & np.isfinite(reflectance)],
        }
        for name, v in values.items():
            row[f"{name}_median"] = round(float(np.median(v)), 4) if v.size else None
            row[f"{name}_pixels"] = int(v.size)
        row["v1_screened_out"] = bool(entry["screened_out"])
        row["water_screened_out"] = bool(row["water_median"] > threshold)
        row["verdict_changes"] = row["v1_screened_out"] != row["water_screened_out"]
        rows.append(row)
        print(
            f"  {row['sensing_utc'][:16]} R{row['relative_orbit']:03d} "
            f"rect={row['rectangle_median']} water={row['water_median']} "
            f"interior={row['interior_water_median']}",
            flush=True,
        )
    directory = new_run_directory(OUTPUT / "screen")
    write(
        directory / "screen.json",
        {
            "question": "Q1/D2: does a water-only screen change any v1 screen verdict?",
            "threshold": threshold,
            "bounds_epsg32629": list(bounds),
            "scenes": rows,
            "verdict_changes": [r["key"] for r in rows if r.get("verdict_changes")],
        },
    )
    return finish("screen", directory, ["screen.json"], started)


# ---------------------------------------------------------------------------------------
# Stage: QC-consistent support and QAA flag classes (Q1, D1, D3)
# ---------------------------------------------------------------------------------------


def campaign_units(
    spec: dict, sensor: str, campaign: str, scale: str, *, variant: str = "V0", **options
):
    """Units for one sensor, campaign and scale, with targets attached."""
    events = campaign_events(spec, campaign)
    scene_id = SCENES[(sensor, campaign)]
    level = rb.acquisition_levels(spec)[scene_id]
    window = options.pop("window", None) or read_window(
        spec, scene_id, event_bounds(events, 250.0), variant
    )
    target = options.pop("targets", None)
    target = targets_for(campaign) if target is None else target
    if scale == "native":
        units = native_units(spec, events, window, level, campaign, **options)
    else:
        units = grid_units(spec, events, window, level, campaign, **options)
    units = rb.attach_targets(units, target)
    units["scene"], units["scale"], units["water_level_m"] = scene_id, scale, level
    return units


def qc() -> Path:
    started = provenance()
    spec = rb.load_spec()
    fit_dir, eval_dir = resolve_run(rb.OUTPUT / "fit"), resolve_run(rb.OUTPUT / "evaluate")
    v1_units = {
        ("july", "native"): pd.read_csv(fit_dir / "units-native_pneo.csv"),
        ("july", "grid"): pd.read_csv(fit_dir / "units-common_grid_10m.csv"),
        ("september", "native"): pd.read_csv(eval_dir / "units-native_pneo.csv"),
        ("september", "grid"): pd.read_csv(eval_dir / "units-common_grid_10m.csv"),
    }
    directory = new_run_directory(OUTPUT / "qc")
    result: dict = {"support": {}, "v1_reproduction": {}, "qaa": {}}
    for (sensor, campaign), scene_id in SCENES.items():
        scales = ("native", "grid") if sensor == "pleiades_neo" else ("grid",)
        for scale in scales:
            units = campaign_units(spec, sensor, campaign, scale)
            units.to_csv(directory / f"units-{sensor}-{campaign}-{scale}.csv", index=False)
            label = f"{sensor}/{campaign}/{scale}"
            result["support"][label] = {"scene": scene_id, "by_stratum": support(units)}
            if sensor != "pleiades_neo":
                continue
            # The same code under v1's finite-only rule must reproduce v1's units.
            as_v1 = campaign_units(spec, sensor, campaign, scale, rule="v1")
            saved = v1_units[(campaign, scale)].set_index("unit_id")
            mine = as_v1.set_index("unit_id")
            shared = saved.index.intersection(mine.index)
            columns = ["blue", "green", "red", "depth"]
            difference = (mine.loc[shared, columns] - saved.loc[shared, columns]).abs()
            result["v1_reproduction"][label] = {
                "v1_units": int(len(saved)),
                "reproduced_units": int(len(mine)),
                "shared": int(len(shared)),
                "max_abs_difference": {k: float(v) for k, v in difference.max().items()},
                "v1_units_kept_in_r1": int(
                    units.set_index("unit_id").reindex(saved.index).reason.eq("ok").sum()
                ),
            }
    for scene_id in ("pneo-20260708", "pneo-20260903"):
        result["qaa"][scene_id] = qaa_fits(spec, scene_id)
    write(directory / "qc.json", result)
    names = ["qc.json", *sorted(p.name for p in directory.glob("units-*.csv"))]
    return finish("qc", directory, names, started)


# ---------------------------------------------------------------------------------------
# Stage: models (Q3) and the registered sensitivities
# ---------------------------------------------------------------------------------------


def model_block(train_units, test_units, spec, family_list) -> dict:
    july, folds = loto(train_units, spec, family_list)
    result = {
        "july_loto": summarise(july, family_list) | {"folds": folds},
        "july_units": int(len(complete_units(train_units))),
    }
    if test_units is not None:
        september = transfer(train_units, test_units, spec, family_list)
        result["september_transfer"] = summarise(september, family_list)
        result["september_units"] = int(len(complete_units(test_units)))
    return result, july


def models(variant: str = "V0") -> Path:
    started = provenance()
    spec = rb.load_spec()
    family_list = families(spec)
    directory = new_run_directory(OUTPUT / ("models" if variant == "V0" else f"models-{variant}"))
    result: dict = {"variant": variant, "primary": {}, "sensitivity": {}}
    july_events = campaign_events(spec, "july")
    windows = {}
    for (sensor, campaign), scene_id in SCENES.items():
        if variant != "V0" and sensor != "pleiades_neo":
            continue
        events = campaign_events(spec, campaign)
        windows[(sensor, campaign)] = read_window(
            spec, scene_id, event_bounds(events, 250.0), variant
        )
    for sensor in ("pleiades_neo", "sentinel2"):
        scales = ("native", "grid") if sensor == "pleiades_neo" else ("grid",)
        if (sensor, "july") not in windows:
            continue
        for scale in scales:
            train = campaign_units(spec, sensor, "july", scale, window=windows[(sensor, "july")])
            test = campaign_units(
                spec, sensor, "september", scale, window=windows[(sensor, "september")]
            )
            block, predictions = model_block(train, test, spec, family_list)
            label = f"{sensor}/{scale}"
            result["primary"][label] = block
            predictions.to_csv(directory / f"loto-{sensor}-{scale}.csv", index=False)
            print(f"  {label}: {_headline(block['july_loto'])}", flush=True)
    # Sensitivities: PNeo, July leave-one-transect-out, both scales.
    window = windows[("pleiades_neo", "july")]
    july_targets = targets_for("july")
    scenarios: list[tuple[str, dict]] = []
    for distance in (2, 5, 10, 20):
        for name, bearing in COMPASS.items():
            angle = math.radians(bearing)
            offset = np.tile(
                [distance * math.sin(angle), distance * math.cos(angle)], (len(july_events), 1)
            )
            scenarios.append((f"shift_{distance}m_{name}", {"shift": offset}))
    scenarios.append(("reversed_order", {"targets": reversed_targets(july_targets)}))
    for size in (1, 5, 7):
        scenarios.append((f"window_{size}px", {"size": size, "scales": ("native",)}))
    for name, options in scenarios:
        scales = options.pop("scales", ("native", "grid"))
        for scale in scales:
            opts = dict(options)
            if scale == "grid":
                opts.pop("size", None)
            units = campaign_units(spec, "pleiades_neo", "july", scale, window=window, **opts)
            frame, folds = loto(units, spec, family_list)
            entry = summarise(frame, family_list)
            entry["units"] = int(len(complete_units(units)))
            entry["sufficient_folds"] = sum(f["sufficient"] for f in folds)
            result["sensitivity"].setdefault(name, {})[scale] = entry
        print(f"  {name}: done", flush=True)
    result["sensitivity_summary"] = sensitivity_summary(result["sensitivity"])
    write(directory / "models.json", result)
    names = ["models.json", *sorted(p.name for p in directory.glob("loto-*.csv"))]
    return finish("models" if variant == "V0" else f"models-{variant}", directory, names, started)


def _headline(block: dict) -> str:
    scores = block.get("scores", {})
    return ", ".join(f"{k}={v['mae_pp']}" for k, v in scores.items())


def sensitivity_summary(sensitivity: dict) -> dict:
    """Range across scenarios of the key paired differences, per scale."""
    out = {}
    for scale in ("native", "grid"):
        for pair in (
            "spectral_depth_ridge - depth_ridge",
            "spectral_ridge - training_mean",
            "nir_depth_ridge - depth_ridge",
        ):
            values = {
                name: block[scale]["paired_mae_difference_pp"].get(pair, {}).get("mean_pp")
                for name, block in sensitivity.items()
                if scale in block and "paired_mae_difference_pp" in block[scale]
            }
            values = {k: v for k, v in values.items() if v is not None}
            if not values:
                continue
            out[f"{scale}: {pair}"] = {
                "min": min(values.values()),
                "max": max(values.values()),
                "scenarios_better_by_2pp": sorted(k for k, v in values.items() if v <= -2.0),
                "n_scenarios": len(values),
            }
    return out


# ---------------------------------------------------------------------------------------
# Stage: injection and recovery
# ---------------------------------------------------------------------------------------


def e5s_inputs() -> dict:
    """E5s Pleiades Neo endmembers (taxon-balanced medians) and Sesimbra IOPs."""
    path = rb.E5S_RESULTS
    e5s = json.loads(path.read_text())
    bands = ("B2", "B3", "B4")  # 486, 563, 655 nm: the predictor bands

    def class_bands(name):
        return np.array(
            [e5s["class_bands_pleiades_neo"][name]["bands"][b]["median"] for b in bands]
        )

    understory = np.mean(
        [class_bands(n) for n in ("turf", "coralline_crust", "red_foliose")], axis=0
    )
    # E5s IOPs are on Sentinel-2 B1-B5; B2-B4 sit within 10 nm of the PNeo bands.
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


def signal_rhos(cover: np.ndarray, depth: np.ndarray, kelp, background, iops, sza) -> list:
    """rhos change per band for fractional kelp cover over a background, from the forward model."""
    from kelpobserve.coastal.inversion.lee import forward_model

    out = []
    for b in range(3):
        water = SimpleNamespace(a=np.array([iops["a"][b]]), bb=np.array([iops["bb"][b]]))
        mixed = cover * kelp[b] + (1.0 - cover) * background[b]
        z = np.where(np.isfinite(depth), depth, 0.0)
        delta = forward_model(mixed, z, water, sza) - forward_model(
            np.full_like(mixed, background[b]), z, water, sza
        )
        out.append(np.where(np.isfinite(depth), np.pi * delta, 0.0))
    return out


def synthetic_field(window: Window, xy_true, cover, strata) -> np.ndarray:
    """Cover on the window grid: 2.5 m discs at true positions, else nearest stratum mean."""
    from scipy.spatial import cKDTree

    rows, cols = np.indices(window.shape)
    xs, ys = window.transform * (cols + 0.5, rows + 0.5)
    tree = cKDTree(xy_true)
    distance, nearest = tree.query(np.column_stack([xs.ravel(), ys.ravel()]))
    means = {s: float(np.nanmean(cover[strata == s])) for s in np.unique(strata)}
    background = np.array([means[s] for s in strata])
    own = np.where(np.isfinite(cover), cover, background)
    field = np.where(distance <= 2.5, own[nearest], background[nearest])
    return field.reshape(window.shape)


def injection_draw(seed, spec, events, window, level, conditions, e5s, family_list):
    """One draw: a within-stratum permutation of July cover, then every condition."""
    rng = np.random.default_rng(seed)
    targets = targets_for("july").set_index("event_id")
    observed = targets.loc[events.event_id, "cover_pct"].to_numpy(dtype=float)
    strata = events.stratum.to_numpy()
    permuted = observed.copy()
    for s in np.unique(strata):
        idx = np.where((strata == s) & np.isfinite(observed))[0]
        permuted[idx] = observed[rng.permutation(idx)]
    synthetic = pd.DataFrame(
        {"event_id": events.event_id, "cover_pct": permuted, "cover_status": "measured"}
    )
    transects = sorted(events.transect.unique())
    bearings = {t: rng.uniform(0, 2 * math.pi) for t in transects}
    depth_pixels = window.lidar + level
    rows = []
    for condition in conditions:
        d = condition["displacement_m"]
        shift = np.array(
            [[d * math.sin(bearings[t]), d * math.cos(bearings[t])] for t in events.transect]
        )
        xy_true = positions(events, shift)
        field = synthetic_field(window, xy_true, permuted / 100.0, strata)
        if condition["amplitude"] == 0:
            bands = window.bands
        else:
            delta = signal_rhos(
                field,
                depth_pixels,
                e5s["kelp"],
                e5s["background"][condition["background"]],
                e5s["iops"][condition["water"]],
                window.solar_zenith_deg,
            )
            bands = dict(window.bands)
            for role, change in zip(RGB, delta):
                bands[role] = window.bands[role] + change
        for scale in ("native", "grid"):
            maker = native_units if scale == "native" else grid_units
            units = rb.attach_targets(
                maker(spec, events, window, level, "july", bands=bands), synthetic
            )
            frame, _ = loto(units, spec, family_list)
            if frame.empty:
                continue
            mae = (
                frame.assign(e=(frame.predicted - frame.observed).abs()).groupby("family").e.mean()
            )
            rows.append(
                {
                    "draw": seed,
                    "condition": condition["id"],
                    "scale": scale,
                    "units": int(len(complete_units(units))),
                }
                | {f"mae_{k}": round(float(v), 3) for k, v in mae.items()}
            )
    return rows


def injection(draws: int | None = None, jobs: int = 8) -> Path:
    from joblib import Parallel, delayed

    started = provenance()
    spec = rb.load_spec()
    registration = load_revision()["injection_recovery"]
    draws = draws or registration["draws"]
    events = campaign_events(spec, "july")
    scene_id = SCENES[("pleiades_neo", "july")]
    level = rb.acquisition_levels(spec)[scene_id]
    window = read_window(spec, scene_id, event_bounds(events, 250.0))
    e5s = e5s_inputs()
    family_list = families(spec, controls=False)
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
    # Signal size at the quadrats for 100% kelp, for the report.
    depths = lidar_at(spec, positions(events)) + level
    sizes = {}
    for background in ("understory", "bare"):
        for water in ("median", "2026-02-22"):
            delta = signal_rhos(
                np.ones_like(depths),
                depths,
                e5s["kelp"],
                e5s["background"][background],
                e5s["iops"][water],
                window.solar_zenith_deg,
            )
            sizes[f"{background}/{water}"] = {
                str(s): [
                    round(float(np.median(band[events.stratum.to_numpy() == s])), 5)
                    for band in delta
                ]
                for s in ("shallow", "deep")
            }
    rows = Parallel(n_jobs=jobs, verbose=5)(
        delayed(injection_draw)(42 + i, spec, events, window, level, conditions, e5s, family_list)
        for i in range(draws)
    )
    frame = pd.DataFrame([r for block in rows for r in block])
    directory = new_run_directory(OUTPUT / "injection")
    frame.to_csv(directory / "draws.csv", index=False)
    summary = {}
    criterion = 2.0
    for (condition, scale), group in frame.groupby(["condition", "scale"]):
        gain = group.mae_depth_ridge - group.mae_spectral_depth_ridge
        forest = group.mae_depth_ridge - group.mae_spectral_depth_forest
        summary[f"{condition}|{scale}"] = {
            "draws": int(len(group)),
            "median_units": float(group.units.median()),
            "median_gain_pp": round(float(gain.median()), 2),
            "share_gain_ge_2pp": round(float((gain >= criterion).mean()), 3),
            "recovered": bool((gain >= criterion).mean() >= 0.8),
            "forest_share_gain_ge_2pp": round(float((forest >= criterion).mean()), 3),
            "spectral_vs_mean_share_ge_2pp": round(
                float((group.mae_training_mean - group.mae_spectral_ridge >= criterion).mean()), 3
            ),
        }
    write(
        directory / "injection.json",
        {
            "registration": registration,
            "e5s": {
                "path": e5s["path"],
                "sha256": e5s["sha256"],
                "kelp": e5s["kelp"].round(5).tolist(),
                "background": {k: v.round(5).tolist() for k, v in e5s["background"].items()},
                "iops": {
                    k: {q: v.round(5).tolist() for q, v in w.items()}
                    for k, w in e5s["iops"].items()
                },
            },
            "signal_rhos_at_quadrats_100pct_kelp": {
                "bands_nm": [486, 563, 655],
                "median_by_stratum": sizes,
            },
            "conditions": conditions,
            "summary": summary,
        },
    )
    return finish("injection", directory, ["draws.csv", "injection.json"], started)


# ---------------------------------------------------------------------------------------
# Stage: campaign windows across years (Q4)
# ---------------------------------------------------------------------------------------


def water_geometry_lonlat(spec: dict, bounds_utm) -> dict:
    """Water (lidar or DTM) inside the screen rectangle, as GeoJSON in longitude/latitude."""
    import rasterio
    from rasterio.features import shapes
    from rasterio.warp import transform_geom
    from shapely.geometry import mapping, shape
    from shapely.ops import unary_union

    pixel = 5.0
    width = int(math.ceil((bounds_utm[2] - bounds_utm[0]) / pixel))
    height = int(math.ceil((bounds_utm[3] - bounds_utm[1]) / pixel))
    transform = rasterio.transform.from_origin(bounds_utm[0], bounds_utm[3], pixel, pixel)
    water = water_on(spec, transform, (height, width), "EPSG:32629")
    polygons = [
        shape(g)
        for g, v in shapes(water.astype(np.uint8), mask=water, transform=transform)
        if v == 1
    ]
    merged = unary_union(polygons)
    return transform_geom("EPSG:32629", "EPSG:4326", mapping(merged))


def assess_full_resolution(item, bbox, thresholds, water_geometry):
    """usability.assess with decimate=1 (registered pre-run clarification)."""
    from kelpobserve.coastal.acquire.geometry import item_tile_angles
    from kelpobserve.coastal.acquire.planetary import measure_over_water, mgrs_tile
    from kelpobserve.coastal.acquire.usability import SWIR_BAND_ID, SceneUsability, classify

    props = item.properties
    tile_cloud = props.get("eo:cloud_cover")
    base = {
        "sensing_utc": item.datetime,
        "platform": props.get("platform"),
        "tile": mgrs_tile(item),
        "relative_orbit": props.get("sat:relative_orbit"),
        "tile_cloud_pct": None if tile_cloud is None else round(float(tile_cloud), 1),
    }
    if tile_cloud is not None and float(tile_cloud) > thresholds.max_tile_cloud_pct:
        return SceneUsability(
            **base,
            glint_angle_deg=None,
            swir_glint=None,
            over_water_red=None,
            n_water_pixels=0,
            verdict="tile_cloudy",
            reasons=["tile_cloudy"],
        )
    try:
        angle = item_tile_angles(item).glint_angle(SWIR_BAND_ID)
        swir, n_water = measure_over_water(
            item,
            bbox,
            "B11",
            decimate=1,
            water_nir_max=thresholds.glint_water_nir_max,
            water_geometry=water_geometry,
        )
        red, _ = measure_over_water(item, bbox, "B04", decimate=1, water_geometry=water_geometry)
    except Exception as error:  # noqa: BLE001 - one unreadable scene is not fatal
        print(f"  unreadable {item.id}: {error}", flush=True)
        return SceneUsability(
            **base,
            glint_angle_deg=None,
            swir_glint=None,
            over_water_red=None,
            n_water_pixels=0,
            verdict="unreadable",
            reasons=["unreadable"],
        )
    return classify(base, angle, swir, red, n_water, thresholds)


def windows_stage(years=range(2017, 2027)) -> Path:
    import geopandas as gpd
    from kelpobserve.coastal.acquire.planetary import open_catalog
    from kelpobserve.coastal.acquire.usability import UsabilityThresholds, overpasses

    started = provenance()
    spec = rb.load_spec()
    registration = load_revision()["Q4_windows"]
    bounds = tuple(
        json.loads((rb.REGISTRATION / "glint-screen.json").read_text())["bounds_epsg32629"]
    )
    bbox = tuple(
        gpd.GeoSeries.from_xy([bounds[0], bounds[2]], [bounds[1], bounds[3]], crs=32629)
        .to_crs(4326)
        .total_bounds
    )
    geometry = water_geometry_lonlat(spec, bounds)
    thresholds = UsabilityThresholds()
    catalog = open_catalog()
    centres = {"june": ((6, 2), (6, 2)), "july": ((7, 7), (7, 8)), "september": ((9, 1), (9, 2))}
    widest = max(registration["windows"]["h_days"])
    likely_max = 0.012
    records = []
    for year in years:
        for campaign, (first, last) in centres.items():
            start = date(year, *first) - timedelta(days=widest)
            end = date(year, *last) + timedelta(days=widest)
            items = list(
                catalog.search(
                    collections=["sentinel-2-l2a"],
                    bbox=list(bbox),
                    datetime=f"{start.isoformat()}/{end.isoformat()}",
                ).items()
            )
            for item in overpasses(items, bbox):
                usability = assess_full_resolution(item, bbox, thresholds, geometry)
                day = usability.sensing_utc.date()
                offset = max(0, (date(year, *first) - day).days, (day - date(year, *last)).days)
                record = usability.to_dict() | {
                    "year": year,
                    "campaign": campaign,
                    "offset_days": offset,
                    "likely_pass": bool(
                        usability.usable
                        and usability.swir_glint is not None
                        and usability.swir_glint <= likely_max
                    ),
                }
                records.append(record)
            screened = sum(r["year"] == year and r["campaign"] == campaign for r in records)
            print(f"  {year} {campaign}: {screened} overpasses", flush=True)
    frame = pd.DataFrame(records)
    directory = new_run_directory(OUTPUT / "windows")
    frame.to_csv(directory / "overpasses.csv", index=False)
    table = []
    for h in registration["windows"]["h_days"]:
        inside = frame[frame.offset_days <= h]
        for (year, campaign), group in inside.groupby(["year", "campaign"]):
            table.append(
                {
                    "h_days": h,
                    "year": int(year),
                    "campaign": campaign,
                    "overpasses": int(len(group)),
                    "usable": int(group.verdict.eq("usable").sum()),
                    "likely_pass": int(group.likely_pass.sum()),
                    "likely_pass_by_orbit": {
                        f"R{int(o):03d}": int(g.likely_pass.sum())
                        for o, g in group.groupby("relative_orbit")
                    },
                    "verdicts": {k: int(v) for k, v in group.verdict.value_counts().items()},
                }
            )
    summary = {}
    table_frame = pd.DataFrame(table)
    for (h, campaign), group in table_frame.groupby(["h_days", "campaign"]):
        past = group[group.year < 2026]
        summary[f"{campaign} +/-{h} d"] = {
            "years": int(len(past)),
            "years_with_usable": int((past.usable > 0).sum()),
            "years_with_likely_pass": int((past.likely_pass > 0).sum()),
            "median_likely_pass": float(past.likely_pass.median()) if len(past) else None,
            "2026_likely_pass": int(group[group.year == 2026].likely_pass.sum()),
            "2026_usable": int(group[group.year == 2026].usable.sum()),
        }
    write(
        directory / "windows.json",
        {
            "registration": registration,
            "bbox_lonlat": list(bbox),
            "likely_pass_max_b11": likely_max,
            "summary": summary,
            "table": table,
        },
    )
    return finish("windows", directory, ["overpasses.csv", "windows.json"], started)


# ---------------------------------------------------------------------------------------
# Stage: ACOLITE glint variants and label-free radiometry (Q2)
# ---------------------------------------------------------------------------------------


def acolite(scene_ids: list[str] | None = None) -> Path:
    """Run the pinned ACOLITE with the registered V1 glint settings."""
    import time

    from kelpobserve.coastal.acolite import correct_scene

    started = provenance()
    spec = rb.load_spec()
    config = spec["acolite"]
    revision = rb.verify_acolite(config)
    registered = load_revision()["Q2_radiometry"]["scenes"]
    wanted = scene_ids or [*registered["sentinel2"], *registered["pleiades_neo"]]
    settings = {name: rb.REGISTRATION / path for name, path in config["settings"].items()}
    overrides = {s: rb.read_settings_file(settings[s]) for s in ("sentinel2", "pleiades_neo")}
    scenes = {s["id"]: s for s in rb.prepare_scenes(spec)}
    directory = new_run_directory(OUTPUT / "acolite")
    records = []
    for scene_id in wanted:
        scene = scenes[scene_id]
        output = ACOLITE_R1 / "V1" / scene_id
        record = {"id": scene_id, "sensor": scene["sensor"], "output_dir": str(output)}
        began = time.monotonic()
        print(f"[acolite V1] {scene_id}", flush=True)
        try:
            bundle = (
                rb.extract_pneo(scene["source"], rb.ROOT / config["pneo_extract_root"])
                if scene["sensor"] == "pleiades_neo"
                else scene["input"]
            )
            correct_scene(
                bundle,
                output,
                scene["limit"],
                Path(config["path"]),
                template=settings["common"],
                extra=overrides[scene["sensor"]] | GLINT_VARIANT[scene["sensor"]],
                python_executable=config["python"],
            )
            record["status"] = "ok"
            record["effective_settings"] = effective_settings(output)
        except Exception as error:  # recorded per scene
            record["status"] = "failed"
            record["error"] = f"{type(error).__name__}: {error}"[:2000]
        record["elapsed_s"] = round(time.monotonic() - began, 1)
        records.append(record)
        write(directory / "acolite.json", {"acolite": revision, "scenes": records})
        print(f"[acolite V1] {scene_id}: {record['status']} ({record['elapsed_s']} s)", flush=True)
    return finish("acolite", directory, ["acolite.json"], started)


def effective_settings(output: Path) -> dict:
    """Glint keys as ACOLITE itself recorded them for the L2R step."""
    files = sorted(output.glob("acolite_run_*_l2r_settings.txt"))
    if not files:
        return {}
    values = rb.read_settings_file(files[-1]) if files else {}
    return {k: v for k, v in values.items() if "glint" in k}


def radiometry_measures(spec: dict, scene_id: str, variant: str) -> tuple[dict, Window]:
    import geopandas as gpd
    from rasterio import features

    region = gpd.read_file(rb.REGISTRATION / spec["pneo_reference_region"]["path"]).to_crs(32629)
    habitat = gpd.read_file(rb.ROOT / spec["sources"]["habitat_aoi"]).to_crs(32629)
    x0 = min(region.total_bounds[0], habitat.total_bounds[0]) - 500
    y0 = min(region.total_bounds[1], habitat.total_bounds[1]) - 500
    x1 = max(region.total_bounds[2], habitat.total_bounds[2]) + 500
    y1 = max(region.total_bounds[3], habitat.total_bounds[3]) + 500
    window = read_window(spec, scene_id, (x0, y0, x1, y1), variant)
    valid = window.valid()

    def rasterize(frame):
        return features.rasterize(
            frame.to_crs(window.crs).geometry, out_shape=window.shape, transform=window.transform
        ).astype(bool)

    in_region, in_habitat = rasterize(region) & valid, rasterize(habitat) & valid
    green, nir, z = window.bands["green"], window.bands["nir"], window.lidar

    def mad(values):
        return float(np.median(np.abs(values - np.median(values)))) if values.size else math.nan

    common = in_habitat & np.isfinite(z) & (z >= 5) & (z <= 25)
    slopes = {}
    for role, band in (("green", green), ("nir", nir)):
        ok = common & np.isfinite(band)
        slopes[role] = float(np.polyfit(z[ok], band[ok], 1)[0]) if ok.sum() > 50 else math.nan
    shallow = in_habitat & np.isfinite(z) & (z >= 5) & (z <= 20)
    measures = {
        "M1_offshore": {
            "pixels": int(in_region.sum()),
            "nir_median": round(float(np.median(nir[in_region])), 5) if in_region.any() else None,
            "green_median": round(float(np.median(green[in_region])), 5)
            if in_region.any()
            else None,
        },
        "M3_common_mode": {
            "pixels": int(common.sum()),
            "green_slope_per_m": round(slopes["green"], 6),
            "nir_slope_per_m": round(slopes["nir"], 6),
            "ratio": round(abs(slopes["nir"]) / abs(slopes["green"]), 3)
            if slopes["green"]
            else None,
        },
        "M4_contrast": {
            "habitat_pixels": int(shallow.sum()),
            "ratio": round(mad(green[shallow]) / mad(green[in_region]), 3)
            if in_region.any() and shallow.any()
            else None,
        },
        "valid_habitat_water_fraction": round(
            float(
                (rasterize(habitat) & window.water & valid).sum()
                / max(1, (rasterize(habitat) & window.water).sum())
            ),
            4,
        ),
    }
    return measures, window


def conversion_check(spec: dict) -> dict:
    """PNeo/S2 rhos over stable land on the 10 m grid, July against September."""
    import rasterio
    from scipy import ndimage

    pairs = [
        ("july", "pneo-20260708", "s2-20260702-S2B-R037"),
        ("september_same_minute", "pneo-20260903", "s2-20260903-S2B-R080"),
        ("september_clean", "pneo-20260903", "s2-20260831-S2B-R037"),
    ]
    roles = spec["scene_qc"]["bands_nm"]
    out = {}
    for label, pneo_id, s2_id in pairs:
        pneo_dir, s2_dir = scene_dir(spec, pneo_id), scene_dir(spec, s2_id)
        pneo_rhos, s2_rhos = (
            rb.acolite_rasters(pneo_dir, "rhos"),
            rb.acolite_rasters(s2_dir, "rhos"),
        )
        with rasterio.open(next(iter(pneo_rhos.values()))) as src:
            bounds = src.bounds
        transform, shape, crs = rb.common_grid(spec, bounds, spec["selected_grid_scene"])
        land = ndimage.binary_erosion(~water_on(spec, transform, shape, crs), iterations=3)
        row = {}
        for role in ("green", "red", "nir"):
            _, pneo_path = rb.nearest_band(pneo_rhos, roles["pleiades_neo"][role])
            _, s2_path = rb.nearest_band(s2_rhos, roles["sentinel2"][role])
            pneo, _ = rb.aggregate_to_grid(pneo_path, transform, shape, crs, 0.8)
            with rasterio.open(s2_path) as src:
                window = rasterio.windows.from_bounds(
                    *rasterio.transform.array_bounds(*shape, transform), transform=src.transform
                )
            s2 = read_clipped(s2_path, window.round_offsets().round_lengths()).astype("float32")
            both = land & np.isfinite(pneo) & np.isfinite(s2) & (s2 > 0.02)
            a, b = pneo[both], s2[both]
            row[role] = {
                "pixels": int(both.sum()),
                "median_ratio": round(float(np.median(a / b)), 4),
                "median_difference": round(float(np.median(a - b)), 5),
                "correlation": round(float(np.corrcoef(a, b)[0, 1]), 3),
            }
        out[label] = {"pneo": pneo_id, "sentinel2": s2_id, "bands": row}
    verdict = {}
    for role in ("green", "red", "nir"):
        gap = abs(
            out["july"]["bands"][role]["median_ratio"]
            - out["september_same_minute"]["bands"][role]["median_ratio"]
        )
        verdict[role] = {"ratio_gap": round(gap, 4), "within_0.05": bool(gap <= 0.05)}
    out["verdict"] = verdict | {
        "consistent": bool(verdict["red"]["within_0.05"] and verdict["nir"]["within_0.05"]),
        "rule": "July ratio within 0.05 of the September same-minute ratio at red and NIR.",
    }
    return out


def improvement(v0: dict, v1: dict, land_change: float) -> dict:
    """The registered rule. A check with a missing measure (no valid water) is None, and
    the scene is then not assessable rather than improved or not."""

    def pair(measure, key):
        first, second = v0[measure].get(key), v1[measure].get(key)
        return (first, second) if first is not None and second is not None else None

    nir, ratio, contrast = (
        pair("M1_offshore", "nir_median"),
        pair("M3_common_mode", "ratio"),
        pair("M4_contrast", "ratio"),
    )
    checks = {
        "nir_offshore_falls_by_0.002": None if nir is None else bool(nir[0] - nir[1] >= 0.002),
        "common_mode_ratio_falls": None if ratio is None else bool(ratio[1] < ratio[0]),
        "contrast_kept_0.9": None if contrast is None else bool(contrast[1] >= 0.9 * contrast[0]),
        "land_change_le_0.002": bool(land_change <= 0.002),
    }
    if any(v is None for v in checks.values()):
        return checks | {"improves": None, "note": "not assessable: a measure has no valid water"}
    return checks | {"improves": all(checks.values())}


def radiometry() -> Path:
    started = provenance()
    spec = rb.load_spec()
    registered = load_revision()["Q2_radiometry"]["scenes"]
    result: dict = {"scenes": {}}
    for scene_id in [*registered["sentinel2"], *registered["pleiades_neo"]]:
        entry: dict = {}
        windows = {}
        for variant in ("V0", "V1"):
            if not (scene_dir(spec, scene_id, variant)).is_dir():
                entry[variant] = {"status": "missing"}
                continue
            measures, windows[variant] = radiometry_measures(spec, scene_id, variant)
            measures["M2_qaa"] = qaa_fits(spec, scene_id, variant)
            entry[variant] = measures
            print(
                f"  {scene_id} {variant}: {measures['M1_offshore']} {measures['M3_common_mode']}",
                flush=True,
            )
        if len(windows) == 2:
            from scipy import ndimage

            v0, v1 = windows["V0"], windows["V1"]
            pixel = abs(v0.transform.a)
            land = ndimage.binary_erosion(~v0.water, iterations=max(1, round(30 / pixel)))
            both = land & np.isfinite(v0.bands["green"]) & np.isfinite(v1.bands["green"])
            m5 = float(np.median(np.abs(v1.bands["green"][both] - v0.bands["green"][both])))
            entry["M5_land_green_change"] = round(m5, 5)
            a, b = entry["V0"], entry["V1"]
            entry["improvement"] = improvement(a, b, m5)
        result["scenes"][scene_id] = entry
    result["conversion_check"] = conversion_check(spec)
    directory = new_run_directory(OUTPUT / "radiometry")
    write(directory / "radiometry.json", result)
    return finish("radiometry", directory, ["radiometry.json"], started)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("screen", help="Q1: water-only glint screen of the v1 candidates.")
    commands.add_parser("qc", help="Q1: QC-consistent support and QAA flag classes.")
    runner = commands.add_parser("acolite", help="Q2: ACOLITE with the registered glint variant.")
    runner.add_argument("--scene", action="append")
    commands.add_parser("radiometry", help="Q2: label-free measures and the conversion check.")
    modeller = commands.add_parser("models", help="Q3: July comparison, transfer, sensitivities.")
    modeller.add_argument("--variant", default="V0", choices=["V0", "V1"])
    injector = commands.add_parser("injection", help="Injection and recovery.")
    injector.add_argument("--draws", type=int)
    injector.add_argument("--jobs", type=int, default=8)
    commands.add_parser("windows", help="Q4: campaign windows across years.")
    args = parser.parse_args()
    if args.command == "screen":
        directory = screen()
    elif args.command == "qc":
        directory = qc()
    elif args.command == "acolite":
        directory = acolite(args.scene)
    elif args.command == "radiometry":
        directory = radiometry()
    elif args.command == "models":
        directory = models(args.variant)
    elif args.command == "injection":
        directory = injection(args.draws, args.jobs)
    else:
        directory = windows_stage()
    print(f"r1 {args.command} written to {directory}")


if __name__ == "__main__":
    main()
