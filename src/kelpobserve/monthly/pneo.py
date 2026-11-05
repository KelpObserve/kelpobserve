"""Extract native six-band PNEO statistics and Lyzenga indices without ML backbones.

The checked-in calibration is frozen from the full, unlabelled scene. Supplying
recalibrate=True recomputes it using the same low-NIR selection and estimator.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform
from rasterio.windows import Window

from kelpobserve.monthly.evaluation import BANDS, PAIRS, digest, write_json


def log_signal(stack, rsw):
    result = np.empty_like(stack, dtype=np.float32)
    for index in range(6):
        result[index] = np.log(np.maximum(stack[index] - rsw[index], 1e-3))
    return result


def calibrate(stack):
    nir = stack[3]
    mask = (nir < np.percentile(nir[nir > 0], 25)) & (nir > 0)
    rsw = np.array([np.percentile(band[mask], 5) for band in stack], dtype=np.float64)
    signal = log_signal(stack, rsw)
    coefficients = {}
    for a, b in PAIRS:
        x, y = signal[BANDS.index(a)][mask], signal[BANDS.index(b)][mask]
        valid = np.isfinite(x) & np.isfinite(y)
        x, y = x[valid], y[valid]
        if len(x) > 200_000:
            indices = np.random.default_rng(0).choice(len(x), size=200_000, replace=False)
            x, y = x[indices], y[indices]
        vx, vy = x.var(ddof=1), y.var(ddof=1)
        covariance = np.cov(x, y, ddof=1)[0, 1]
        denominator = vy - vx + math.sqrt((vy - vx) ** 2 + 4 * covariance**2)
        coefficients[f"DII_{a}_{b}"] = float(2 * covariance / denominator) if denominator else 1.0
    return rsw, coefficients


def patch_features(patch, rsw, coefficients):
    if patch.shape != (6, 5, 5) or not np.isfinite(patch).all() or not (patch > 0).all():
        raise ValueError("Expected a finite, positive six-band 5x5 patch")
    signal = log_signal(patch, rsw)
    records = []
    for window in [1, 3, 5]:
        half = window // 2
        region = slice(2 - half, 3 + half)
        values = patch[:, region, region].astype(np.float64)
        row = {"window_px": window}
        for i, band in enumerate(BANDS):
            row.update(
                {
                    f"{band}_mean": float(values[i].mean()),
                    f"{band}_std": float(values[i].std()),
                    f"{band}_range": float(np.ptp(values[i])),
                }
            )
        for a, b in PAIRS:
            name = f"DII_{a}_{b}"
            indices = (
                signal[BANDS.index(a), region, region] - coefficients[name] * signal[BANDS.index(b), region, region]
            )
            row[name] = float(indices.mean())
        records.append(row)
    return records


def extract(pneo_root: Path, inputs: Path, output: Path, recalibrate: bool = False) -> Path:
    if inputs.resolve() == output.resolve():
        raise ValueError("Extract into a separate folder to preserve the frozen inputs")
    observations = pd.read_csv(inputs / "observations.csv")
    sources = json.loads((inputs / "pneo_sources.json").read_text())
    rows, calibration = [], {}
    for month, campaign in [("7", "july"), ("9", "september")]:
        record = sources[month]
        files = [pneo_root / item["path"] for item in record["files"]]
        for path, item in zip(files, record["files"]):
            if digest(path) != item["sha256"]:
                raise ValueError(f"Source checksum mismatch: {path.name}")
        rgb_path = next(p for p in files if "_RGB_" in p.name)
        ned_path = next(p for p in files if "_NED_" in p.name)
        rsw = np.array(record["rsw"], dtype=np.float64)
        coefficients = record["lyzenga_coefficients"]
        with rasterio.open(rgb_path) as rgb, rasterio.open(ned_path) as ned:
            if (
                rgb.transform != ned.transform
                or rgb.crs != ned.crs
                or rgb.shape != ned.shape
                or rgb.count != 3
                or ned.count != 3
            ):
                raise ValueError("RGB and NED rasters must share a three-band native grid")
            if recalibrate:
                stack = np.concatenate([rgb.read().astype(np.float32), ned.read().astype(np.float32)])
                rsw, coefficients = calibrate(stack)
                del stack
            for obs in observations[observations.campaign == campaign].itertuples():
                xs, ys = transform("EPSG:4326", rgb.crs, [obs.lon], [obs.lat])
                row, col = rgb.index(xs[0], ys[0])
                if not (2 <= row < rgb.height - 2 and 2 <= col < rgb.width - 2):
                    raise ValueError(f"Incomplete native window for {obs.event_id}")
                window = Window(col - 2, row - 2, 5, 5)
                patch = np.concatenate(
                    [rgb.read(window=window).astype(np.float32), ned.read(window=window).astype(np.float32)]
                )
                for values in patch_features(patch, rsw, coefficients):
                    rows.append(dict(event_id=obs.event_id, row=row, column=col, **values))
        calibration[campaign] = {"rsw": rsw.tolist(), "lyzenga_coefficients": coefficients}
    output.mkdir(parents=True, exist_ok=True)
    path = output / "features.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    write_json(
        output / "extraction.json",
        {
            "sources_sha256": digest(inputs / "pneo_sources.json"),
            "observations_sha256": digest(inputs / "observations.csv"),
            "calibration": calibration,
            "recalibrated": recalibrate,
            "windows_px": [1, 3, 5],
            "native_pixel_m": 1.2,
            "labels_used": False,
            "note": "Calibration uses the whole unlabelled scene, including test image locations.",
        },
    )
    manifest = json.loads((inputs / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        if digest(inputs / name) != expected:
            raise ValueError(f"Input checksum mismatch: {name}")
        if name not in ["features.csv", "extraction.json"]:
            shutil.copyfile(inputs / name, output / name)
    manifest["files"] = {name: digest(output / name) for name in manifest["files"]}
    manifest["files"]["extraction.json"] = digest(output / "extraction.json")
    manifest["feature_precision"] = "Native windows re-extracted with this run's recorded calibration."
    write_json(output / "manifest.json", manifest)
    return path
