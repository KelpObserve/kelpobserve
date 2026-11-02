"""Exploratory r5: what is the reflectance variation that depth does not explain?

Registered in registration/revision-r5.json before any measure was computed. Label-free:
no cover value is read. Measures, per scene: spatial scales of the variation (M1), its
coherence with NIR and SWIR (M2), its dependence on distance to shore (M3), its size by
depth bin against offshore water (M4), its stability across dates (M5), and the spatial
scale at which glint correction works (M6). Writes under out/r5/.

    PYTHONPATH=. python -m benchmark.revision_r5 run
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from kelpobserve.coastal.provenance import file_digest, new_run_directory, resolve_run

from benchmark import revision as rv
from benchmark import revision_r3 as r3
from benchmark import run_benchmark as rb

REGISTRATION_R5 = rb.REGISTRATION / "revision-r5.json"
OUTPUT = rb.OUTPUT / "r5"
BINS = {"5-10": (5.0, 10.0), "10-14": (10.0, 14.0), "14-20": (14.0, 20.0), "20-30": (20.0, 30.0)}
LAGS = {"pleiades_neo": (1, 2, 4, 8, 17, 42, 83, 167), "sentinel2": (1, 2, 5, 10, 20)}
VISIBLE = ("blue", "green")
PAIRS = [
    ("sentinel2", "s2-20260702-S2B-R037", "s2-20260709-S2A-R037"),
    ("sentinel2", "s2-20260528-S2C-R037", "s2-20260602-S2B-R037"),
    ("sentinel2", "s2-20260530-S2A-R037", "s2-20260607-S2C-R037"),
    ("sentinel2", "s2-20260602-S2B-R037", "s2-20260609-S2A-R037"),
    ("sentinel2", "s2-20260702-S2B-R037", "s2-20260831-S2B-R037"),
    ("pneo_10m", "pneo-20260708", "pneo-20260903"),
]
CLEAN = ("s2-20260702-S2B-R037", "s2-20260831-S2B-R037")


def provenance() -> dict:
    return rv.provenance() | {
        "registration_r5": {"path": str(REGISTRATION_R5), "sha256": file_digest(REGISTRATION_R5)},
        "r5_runner_sha256": file_digest(Path(__file__)),
    }


# ---------------------------------------------------------------------------------------
# Small statistics
# ---------------------------------------------------------------------------------------


def robust_sd(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size < 20:
        return math.nan
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def semivariance(z: np.ndarray, mask: np.ndarray, lag: int) -> float:
    """Robust semivariance at a lag along rows and columns (Gaussian-consistent)."""
    pieces = []
    for axis in (0, 1):
        a = z[:-lag] if axis == 0 else z[:, :-lag]
        b = z[lag:] if axis == 0 else z[:, lag:]
        ma = mask[:-lag] if axis == 0 else mask[:, :-lag]
        mb = mask[lag:] if axis == 0 else mask[:, lag:]
        ok = ma & mb & np.isfinite(a) & np.isfinite(b)
        pieces.append((a - b)[ok])
    d = np.concatenate(pieces)
    if d.size < 50:
        return math.nan
    return float(0.5 * np.median(d**2) / 0.4549)


def least_squares(y: np.ndarray, columns: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, float]:
    """Coefficients, residuals and R-squared of y on an intercept and the given columns."""
    x = np.column_stack([np.ones_like(y), *columns])
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ coef
    total = np.sum((y - y.mean()) ** 2)
    return coef, residual, float(1 - np.sum(residual**2) / total) if total > 0 else math.nan


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 20 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return math.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


# ---------------------------------------------------------------------------------------
# Fields: one scene at one scale, with its regions
# ---------------------------------------------------------------------------------------


@dataclass
class Field:
    name: str
    sensor: str
    bands: dict  # blue, green, nir and, for Sentinel-2, swir (top of atmosphere)
    valid: np.ndarray
    depth: np.ndarray
    pixel_m: float
    regions: dict = field(default_factory=dict)
    shore_m: np.ndarray | None = None


def regions_for(spec: dict, transform, shape, crs, valid, depth, water) -> tuple[dict, np.ndarray]:
    import geopandas as gpd
    from rasterio import features
    from scipy import ndimage

    def rasterize(path):
        frame = gpd.read_file(path).to_crs(crs)
        return features.rasterize(frame.geometry, out_shape=shape, transform=transform).astype(bool)

    region = rasterize(rb.REGISTRATION / spec["pneo_reference_region"]["path"])
    habitat = rasterize(rb.ROOT / spec["sources"]["habitat_aoi"])
    out = {"offshore": region & valid}
    with np.errstate(invalid="ignore"):
        for name, (lo, hi) in BINS.items():
            out[name] = habitat & valid & (depth >= lo) & (depth < hi)
        out["habitat_5_30"] = habitat & valid & (depth >= 5) & (depth < 30)
    shore = ndimage.distance_transform_edt(water) * abs(transform.a)
    return out, shore


def pixel_field(spec, window: rv.Window, name: str) -> Field:
    bands = {r: window.bands[r] for r in ("blue", "green", "nir")}
    if window.sensor == "sentinel2":
        bands["swir"] = window.mask_rhot
    valid = window.valid() & np.isfinite(bands["nir"])
    f = Field(name, window.sensor, bands, valid, window.lidar, abs(window.transform.a))
    f.regions, f.shore_m = regions_for(
        spec, window.transform, window.shape, window.crs, valid, window.lidar, window.water
    )
    return f


def pneo_native_field(spec, window: rv.Window, name: str) -> Field:
    """3 x 3 moving means, kept where all nine pixels are valid (r1's native unit)."""
    from scipy.ndimage import uniform_filter

    base = pixel_field(spec, window, name)
    share = uniform_filter(base.valid.astype(float), 3, mode="constant")
    full = share > 1 - 1e-9
    bands = {
        r: np.where(full, uniform_filter(np.where(base.valid, v, 0.0), 3, mode="constant"), np.nan)
        for r, v in base.bands.items()
    }
    f = Field(name, base.sensor, bands, full, window.lidar, base.pixel_m)
    f.regions, f.shore_m = regions_for(
        spec, window.transform, window.shape, window.crs, full, window.lidar, window.water
    )
    return f


def grid_field(spec, window: rv.Window, name: str, bounds) -> Field:
    """Valid-area means on v1's 10 m grid, kept at >= 80% valid (r1's 10 m unit)."""
    from rasterio.warp import Resampling

    transform, shape, crs = rb.common_grid(spec, bounds, spec["selected_grid_scene"])
    valid = window.valid() & np.isfinite(window.bands["nir"])
    fraction = rv._average(valid.astype("float32"), window, transform, shape, crs, None)
    keep = fraction >= spec["common_grid"]["min_valid_fraction"]
    bands = {
        r: np.where(
            keep,
            rv._average(
                np.where(valid, window.bands[r], np.nan), window, transform, shape, crs, np.nan
            ),
            np.nan,
        )
        for r in ("blue", "green", "nir")
    }
    water = rv._average(window.water.astype("float32"), window, transform, shape, crs, None) >= 0.5
    depth = rv.raster_on(spec["sources"]["lidar"], transform, shape, crs, Resampling.average)
    f = Field(name, window.sensor, bands, keep, depth, 10.0)
    f.regions, f.shore_m = regions_for(spec, transform, shape, crs, keep, depth, water)
    return f


# ---------------------------------------------------------------------------------------
# Measures
# ---------------------------------------------------------------------------------------


def m1_scales(f: Field) -> dict:
    out = {}
    for region in ("offshore", "5-10", "14-20"):
        mask = f.regions[region]
        out[region] = {
            role: {
                f"{lag * f.pixel_m:g}m": _sqrt(semivariance(f.bands[role], mask, lag))
                for lag in LAGS[f.sensor]
            }
            for role in ("blue", "green", "nir")
        }
    return out


def _sqrt(value: float) -> float:
    return math.sqrt(value) if value == value and value >= 0 else math.nan


def detrended(f: Field, role: str, mask: np.ndarray, use_depth: bool) -> np.ndarray:
    y = f.bands[role][mask]
    if not use_depth or y.size < 20:
        return y - np.nanmedian(y) if y.size else y
    _, residual, _ = least_squares(y, [f.depth[mask]])
    return residual


def m2_coherence(f: Field) -> dict:
    out = {}
    for region in ("offshore", *BINS):
        mask = (
            f.regions[region] & np.isfinite(f.depth) if region != "offshore" else f.regions[region]
        )
        if mask.sum() < 30:
            continue
        use_depth = region != "offshore"
        entry = {"pixels": int(mask.sum())}
        references = ("nir", "swir") if "swir" in f.bands else ("nir",)
        for reference in references:
            if reference == "swir" and not np.isfinite(f.bands["swir"][mask]).all():
                continue
            x = detrended(f, reference, mask, use_depth)
            for role in VISIBLE:
                y = detrended(f, role, mask, use_depth)
                coef, _, r2 = least_squares(y, [x])
                entry[f"{role}_on_{reference}"] = {
                    "r2": round(r2, 3),
                    "slope": round(float(coef[1]), 3),
                }
        out[region] = entry
    return out


def m3_shore(f: Field) -> dict:
    mask = f.regions["habitat_5_30"] & np.isfinite(f.depth) & (f.shore_m > 0)
    if mask.sum() < 50:
        return {"pixels": int(mask.sum())}
    depth, distance = f.depth[mask], np.log10(f.shore_m[mask])
    out = {"pixels": int(mask.sum())}
    for role in ("blue", "green", "nir"):
        y = f.bands[role][mask]
        _, _, r2_depth = least_squares(y, [depth])
        coef, _, r2_both = least_squares(y, [depth, distance])
        out[role] = {
            "distance_coef_per_decade": round(float(coef[2]), 5),
            "r2_gain": round(r2_both - r2_depth, 3),
            "r2_depth_only": round(r2_depth, 3),
        }
    return out


def residual_field(f: Field, role: str, mask: np.ndarray, with_nir: bool, use_depth: bool):
    """Residual of a band after depth (and NIR), as a full array NaN outside the mask."""
    out = np.full(f.valid.shape, np.nan)
    columns = []
    if use_depth:
        columns.append(f.depth[mask])
    if with_nir:
        columns.append(f.bands["nir"][mask])
    y = f.bands[role][mask]
    if y.size < 20:
        return out
    if columns:
        _, residual, _ = least_squares(y, columns)
    else:
        residual = y - np.median(y)
    out[mask] = residual
    return out


def m4_bins(f: Field) -> tuple[dict, dict]:
    """Unit-scale spread by bin and offshore; also returns the residual fields for M5."""
    out, residuals = {}, {}
    for region in ("offshore", *BINS):
        use_depth = region != "offshore"
        mask = f.regions[region] & (np.isfinite(f.depth) if use_depth else True)
        if mask.sum() < 20:
            continue
        entry = {"units": int(mask.sum())}
        for role in VISIBLE:
            total = residual_field(f, role, mask, False, use_depth)
            after = residual_field(f, role, mask, True, use_depth)
            entry[role] = {"sd": robust_sd(total[mask]), "sd_after_nir": robust_sd(after[mask])}
            residuals[(region, role)] = after
        out[region] = entry
    for region, entry in out.items():
        if region == "offshore" or "offshore" not in out:
            continue
        for role in VISIBLE:
            off = out["offshore"][role]
            for key in ("sd", "sd_after_nir"):
                a, b = entry[role][key], off[key]
                entry[role][f"{key}_ratio_to_offshore"] = a / b if b and b == b else math.nan
                entry[role][f"{key}_excess"] = (
                    math.sqrt(max(0.0, a**2 - b**2)) if a == a and b == b else math.nan
                )
    return out, residuals


def m5_pair(res_a: dict, res_b: dict) -> dict:
    out = {}
    for (region, role), a in res_a.items():
        b = res_b.get((region, role))
        if b is None or a.shape != b.shape:
            continue
        both = np.isfinite(a) & np.isfinite(b)
        out.setdefault(region, {})[role] = {
            "r": correlation(a[both], b[both]),
            "units": int(both.sum()),
        }
    return out


def blocks(
    z: np.ndarray, size: int, mask: np.ndarray, offset=(0, 0)
) -> tuple[np.ndarray, np.ndarray]:
    """Block means of `size` x `size`, with block validity requiring every pixel valid."""
    z = z[offset[0] :, offset[1] :]
    mask = mask[offset[0] :, offset[1] :]
    h, w = (z.shape[0] // size) * size, (z.shape[1] // size) * size
    zz = np.where(mask, z, np.nan)[:h, :w].reshape(h // size, size, w // size, size)
    valid = mask[:h, :w].reshape(h // size, size, w // size, size).all(axis=(1, 3))
    with np.errstate(invalid="ignore"):
        return np.where(valid, np.nanmean(zz, axis=(1, 3)), np.nan), valid


def high_pass(z: np.ndarray, valid: np.ndarray) -> np.ndarray:
    from scipy.ndimage import uniform_filter

    full = uniform_filter(valid.astype(float), 3, mode="constant") > 1 - 1e-9
    smooth = uniform_filter(np.where(valid, z, 0.0), 3, mode="constant")
    return np.where(full, z - smooth, np.nan)


def m6_correction(v0: Field, v1: Field) -> dict:
    """Spatial scale of the glint correction, offshore."""
    mask = v0.regions["offshore"] & v1.regions["offshore"]
    if v0.sensor == "sentinel2":
        reference = v0.bands["swir"]
        # Align 2 x 2 blocks with SWIR's 20 m pixels: least within-block SWIR variance.
        best = None
        for oy in (0, 1):
            for ox in (0, 1):
                s = reference[oy:, ox:]
                h, w = (s.shape[0] // 2) * 2, (s.shape[1] // 2) * 2
                spread = np.nanmean(np.nanstd(s[:h, :w].reshape(h // 2, 2, w // 2, 2), axis=(1, 3)))
                if best is None or spread < best[0]:
                    best = (spread, (oy, ox))
        offset, fine_size, coarse_size = best[1], 2, 6
    else:
        reference = v0.bands["nir"]
        offset, fine_size, coarse_size = (0, 0), 1, 10
    if fine_size > 1:
        g0, fv = blocks(v0.bands["green"], fine_size, mask, offset)
        g1, _ = blocks(v1.bands["green"], fine_size, mask, offset)
        ref, _ = blocks(reference, fine_size, mask, offset)
    else:
        g0, g1, ref, fv = v0.bands["green"], v1.bands["green"], reference, mask
    hp0, hp1, hpr = high_pass(g0, fv), high_pass(g1, fv), high_pass(ref, fv)
    hpd = hp0 - hp1
    shifts = {}
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            shifted = np.roll(np.roll(hpr, dy, axis=0), dx, axis=1)
            shifts[f"{dy},{dx}"] = correlation(hp0, shifted)
    c0, _ = blocks(v0.bands["green"], coarse_size, mask, offset)
    cr, _ = blocks(reference, coarse_size, mask, offset)
    best_shift = max(shifts, key=lambda k: -1 if shifts[k] != shifts[k] else shifts[k])
    finite = np.isfinite(hp1) & np.isfinite(hpd)
    return {
        "fine_m": fine_size * v0.pixel_m,
        "coarse_m": coarse_size * v0.pixel_m,
        "alignment_offset": list(offset),
        "r_fine": shifts["0,0"],
        "r_coarse": correlation(c0, cr),
        "best_shift": best_shift,
        "r_fine_best_shift": shifts[best_shift],
        "r_fine_by_shift": shifts,
        "correction_fine_variance_ratio": float(np.var(hpd[finite]) / np.var(hp1[finite]))
        if finite.sum() > 20
        else math.nan,
        "r_fine_v1_vs_correction": correlation(hp1, -hpd),
        "fine_sd_v0": robust_sd(hp0),
        "fine_sd_v1": robust_sd(hp1),
    }


# ---------------------------------------------------------------------------------------
# Decision rules
# ---------------------------------------------------------------------------------------


def decide(results: dict, pairs: dict, glinted: set) -> dict:
    out = {
        "D1_noise_floor": {},
        "D2_common_mode": {},
        "D3_seabed_candidate": {},
        "D4_per_pixel_glint": {},
    }
    for name, r in results.items():
        for unit, bins in r.get("m4", {}).items():
            off = bins.get("offshore")
            if not off:
                continue
            for b in BINS:
                if b not in bins:
                    continue
                for role in VISIBLE:
                    ratio = (
                        off[role]["sd"] / bins[b][role]["sd"] if bins[b][role]["sd"] else math.nan
                    )
                    out["D1_noise_floor"].setdefault(f"{name}|{unit}|{b}", {})[role] = {
                        "offshore_over_bin": round(ratio, 3),
                        "met": bool(ratio >= 0.8),
                        "clean_date": name.split("|")[0] in CLEAN,
                    }
                    coherence = r["m2"].get(unit, {}).get(b, {}).get(f"{role}_on_nir", {})
                    shore = r["m3"].get(unit, {}).get("nir", {})
                    met = bool(
                        coherence.get("r2", 0) >= 0.5
                        and shore.get("distance_coef_per_decade", 0) < 0
                        and shore.get("r2_gain", 0) >= 0.05
                    )
                    out["D2_common_mode"].setdefault(f"{name}|{unit}|{b}", {})[role] = {
                        "r2_on_nir": coherence.get("r2"),
                        "nir_distance_coef": shore.get("distance_coef_per_decade"),
                        "nir_r2_gain": shore.get("r2_gain"),
                        "met": met,
                    }
    for key, p in pairs.items():
        if p["days"] > 8:
            continue
        for unit, regions in p["m5"].items():
            off = regions.get("offshore", {})
            for b in BINS:
                if b not in regions:
                    continue
                for role in VISIBLE:
                    r_bin, r_off = regions[b][role]["r"], off.get(role, {}).get("r", math.nan)
                    sizes = [
                        results[s]["m4"][unit]
                        .get(b, {})
                        .get(role, {})
                        .get("sd_after_nir_ratio_to_offshore", math.nan)
                        for s in p["scenes"]
                    ]
                    met = bool(
                        all(x == x and x >= 1.5 for x in sizes)
                        and r_bin == r_bin
                        and r_bin >= 0.3
                        and (r_off != r_off or r_off < 0.2)
                    )
                    out["D3_seabed_candidate"].setdefault(f"{key}|{unit}|{b}", {})[role] = {
                        "r_bin": r_bin,
                        "r_offshore": r_off,
                        "size_ratios": sizes,
                        "met": met,
                    }
    for name, r in results.items():
        m6 = r.get("m6")
        if not m6 or name not in glinted:
            continue
        out["D4_per_pixel_glint"][name] = {
            "r_fine": m6["r_fine"],
            "r_coarse": m6["r_coarse"],
            "smoothing_warranted": bool(m6["r_fine"] < 0.3 and m6["r_coarse"] > 0.7),
            "shift_contributes": bool(m6["r_fine_best_shift"] - m6["r_fine"] >= 0.1),
        }
    return out


# ---------------------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------------------


def run() -> Path:
    from kelpobserve.coastal.provenance import complete_run

    started = provenance()
    spec = rb.load_spec()
    registration = json.loads(REGISTRATION_R5.read_text())
    mask = json.loads((resolve_run(rb.OUTPUT / "r3" / "mask") / "mask.json").read_text())
    bounds = r3.study_bounds(spec)
    glinted = {s for s, row in mask["scenes"].items() if row["glinted"]}
    results, residuals = {}, {}

    def window_for(scene_id, version):
        window = rv.read_window(spec, scene_id, bounds, version)
        if scene_id.startswith("s2"):
            window = window.with_threshold(mask["scenes"][scene_id]["threshold"])
        return window

    v0_also = set(registration["scenes"]["sentinel2"]["V0_also"])
    plan = [
        (s, ("V0", "V1") if s in v0_also else ("V1",))
        for s in registration["scenes"]["sentinel2"]["V1"]
    ]
    plan += [(s, ("V0", "V1")) for s in registration["scenes"]["pleiades_neo"]["dates"]]
    for scene_id, versions in plan:
        pixels = {}
        for version in versions:
            print(f"  {scene_id} {version}", flush=True)
            window = window_for(scene_id, version)
            name = f"{scene_id}|{version}"
            pixel = pixel_field(spec, window, name)
            units = (
                {"s2_10m": pixel}
                if window.sensor == "sentinel2"
                else {
                    "pneo_native": pneo_native_field(spec, window, name),
                    "pneo_10m": grid_field(spec, window, name, bounds),
                }
            )
            entry = {
                "m1": m1_scales(pixel),
                "m2": {"pixel": m2_coherence(pixel)},
                "m3": {"pixel": m3_shore(pixel)},
                "m4": {},
            }
            for unit, f in units.items():
                entry["m2"][unit] = m2_coherence(f)
                entry["m3"][unit] = m3_shore(f)
                entry["m4"][unit], kept = m4_bins(f)
                if unit in ("s2_10m", "pneo_10m"):  # only these enter the cross-date pairs
                    residuals[(name, unit)] = kept
            results[name] = entry
            pixels[version] = pixel
            del window, units
        if len(pixels) == 2:
            results[f"{scene_id}|V1"]["m6"] = m6_correction(pixels["V0"], pixels["V1"])
        del pixels
    pairs = {}
    for sensor, a, b in PAIRS:
        unit = "s2_10m" if sensor == "sentinel2" else "pneo_10m"
        versions = ("V1",) if sensor == "sentinel2" else ("V0", "V1")
        for version in versions:
            key = f"{a}~{b}|{version}"
            days = abs((r3.scene_date(b) - r3.scene_date(a)).days)
            pairs[key] = {
                "days": days,
                "scenes": [f"{a}|{version}", f"{b}|{version}"],
                "m5": {
                    unit: m5_pair(
                        residuals[(f"{a}|{version}", unit)], residuals[(f"{b}|{version}", unit)]
                    )
                },
            }
    decisions = decide(results, pairs, {f"{s}|V1" for s in glinted})
    directory = new_run_directory(OUTPUT / "diagnosis")
    rv.write(
        directory / "diagnosis.json", {"scenes": results, "pairs": pairs, "decisions": decisions}
    )
    complete_run(
        OUTPUT / "diagnosis",
        directory,
        {
            "stage": "r5-diagnosis",
            "created_at": datetime.now(UTC).isoformat(),
            "exploratory": True,
            "outputs": {"diagnosis.json": file_digest(directory / "diagnosis.json")},
            **started,
        },
    )
    return directory


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("run", help="Compute every registered measure and decision.")
    parser.parse_args()
    print(f"r5 written to {run()}")


if __name__ == "__main__":
    main()
