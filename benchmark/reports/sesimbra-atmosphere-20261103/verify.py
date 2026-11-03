"""Post-run checks for the bounded optical comparison, without cover targets."""

import argparse
import hashlib
import importlib.metadata
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from netCDF4 import Dataset
from rasterio.windows import Window

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))
from kelpobserve import config

PATHS = config.load()
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
parser.add_argument("--comparison-root", type=Path)
parser.add_argument("--cache-root", type=Path, default=PATHS.require("cache") / "optical-atmosphere-20261103")
args = parser.parse_args()
OUT = args.output
CACHE = args.cache_root
AC = PATHS.require("acolite")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


checks = {}
if args.comparison_root:
    count = 0
    for p in args.comparison_root.glob("*.csv*"):
        pd.testing.assert_frame_equal(pd.read_csv(p), pd.read_csv(OUT / p.name), check_exact=True)
        count += 1
    for name in ["calibration.json", "water-states.json", "verification.json"]:
        assert json.loads((args.comparison_root / name).read_text()) == json.loads((OUT / name).read_text())
    checks["style_cleanup_numerical_parity"] = {"csv_files": count, "json_files": 3, "exact": True}
    print("Analysis parity passed", flush=True)

frozen = json.loads((REPO / "benchmark/reports/r7-local-cover-20261102-r3/freeze.json").read_text())
for name, h in frozen["code_sha256"].items():
    assert sha(REPO / name) == h, name
checks["frozen_r7_sources_unchanged"] = len(frozen["code_sha256"])
checks["previous_output_manifests"] = {}
for name in ["r7-local-cover-20261102-r3", "r7-optical-diagnostic-20261102", "pneo-optical-repair-20261103"]:
    p = REPO / "benchmark/reports" / name
    stored = json.loads((p / "outputs.json").read_text())
    manifest = stored.get("sha256", stored)
    for filename, h in manifest.items():
        assert sha(p / filename) == h, (name, filename)
    checks["previous_output_manifests"][name] = len(manifest)
print("Frozen source and previous outputs unchanged", flush=True)

checks["packed_inputs"] = {}
for folder in sorted((CACHE / "inputs").iterdir()):
    provenance = json.loads((folder / "provenance.json").read_text())
    ncf = next(folder.glob("*_L1R.nc"))
    assert sha(ncf) == provenance["l1r_sha256"]
    with Dataset(ncf) as nc:
        for filename, record in provenance["input_files"].items():
            assert sha(filename) == record["sha256"], filename
            with rasterio.open(filename) as src:
                x = src.read(1, window=Window(*provenance["window"])).astype("float32")
                if src.nodata is not None:
                    x[x == src.nodata] = np.nan
            y = np.asarray(nc[record["dataset"]][:].filled(np.nan))
            assert np.array_equal(x, y, equal_nan=True), filename
    checks["packed_inputs"][folder.name] = {
        "source_arrays_exact": len(provenance["input_files"]),
        "hashes_unchanged": True,
    }
print("Packed native inputs passed", flush=True)

checks["exported_rasters"] = {}
checks["glint_applied"] = {}
for branch in ["s2_radcor_fixed", "s2_radcor_auto", "pneo_marine_control", "pneo_continental"]:
    scene = "s2-20260702-S2B-R037" if branch.startswith("s2") else "pneo-20260708"
    pr = json.loads((CACHE / "inputs" / scene / "provenance.json").read_text())
    folders = [CACHE / branch]
    if branch.startswith("s2"):
        folders.append(CACHE / branch / "glint")
    for directory in folders:
        ncf = next(directory.glob("*.nc"))
        rasters = list(directory.glob("*.tif"))
        with Dataset(ncf) as nc:
            for p in rasters:
                name = p.stem[len(ncf.stem) + 1 :]
                x = np.asarray(nc[name][:].filled(np.nan)).squeeze().astype("float32")
                with rasterio.open(p) as src:
                    assert src.crs == rasterio.crs.CRS.from_string(pr["crs"])
                    assert tuple(src.transform) == tuple(pr["transform"])
                    assert np.array_equal(src.read(1), x, equal_nan=True), p
            if directory.name == "glint":
                assert nc.getncattr("glint_correction") == "default"
                assert all(
                    f"{prefix}_{wave}" in nc.variables
                    for prefix in ["rhos", "rhosu"]
                    for wave in [442, 492, 559, 665, 833]
                )
                checks["glint_applied"][branch] = True
        checks["exported_rasters"][str(directory.relative_to(CACHE))] = len(rasters)
print("Exact native export and glint checks passed", flush=True)

checks["pneo_control_native_max_errors"] = {}
for wave in [439, 486, 563, 655, 826]:
    a = next((CACHE / "pneo_marine_control").glob(f"*_rhos_{wave}.tif"))
    b = next((CACHE.parent / "acolite/pneo-20260708").glob(f"*_rhos_{wave}.tif"))
    with rasterio.open(a) as src:
        x = src.read(1)
    with rasterio.open(b) as src:
        y = src.read(1)
    assert np.array_equal(np.isfinite(x), np.isfinite(y))
    delta = float(np.nanmax(np.abs(x - y)))
    assert delta <= 2e-7, (wave, delta)
    checks["pneo_control_native_max_errors"][str(wave)] = delta

frame = pd.read_csv(OUT / "predictions.csv.gz")
bands = ["coastal", "blue", "green", "red"]
assert set(frame.columns).isdisjoint({"cover", "kelp_cover", "observed_cover"})
valid = np.isfinite(frame[["rhos_" + b for b in bands]]).all(axis=1)
assert np.array_equal(valid, frame.all_bands_valid)
common = frame.assign(check_valid=valid).groupby(["scene", "fold", "x", "y"]).check_valid.transform("all")
assert np.array_equal(common, frame.common_support)
old_cols = ["scene", "region", "x", "y"] + [f"rhos_{b}_{v}" for b in bands for v in ["V0", "V1"]]
old = pd.read_csv(REPO / "benchmark/reports/r7-optical-diagnostic-20261102/paired-grid-cells.csv.gz", usecols=old_cols)
keys = ["scene", "region", "x", "y"]
for version in ["V0", "V1"]:
    selected = frame[frame.branch == "archived_" + version]
    paired = selected.merge(old, on=keys, validate="many_to_one")
    assert len(paired) == len(selected)
    for band in bands:
        np.testing.assert_array_equal(paired["rhos_" + band], paired[f"rhos_{band}_{version}"])
calibrations = json.loads((OUT / "calibration.json").read_text())
assert all(x["disjoint"] for k, x in calibrations.items() if "full_reference" not in k)
checks["common_support_independently_recomputed"] = True
checks["archived_controls_exact"] = True
checks["glint_calibration_test_blocks_disjoint"] = True
checks["cover_targets_read"] = False
checks["cover_scores_computed"] = False

head = subprocess.check_output(["git", "-C", str(AC), "rev-parse", "HEAD"], text=True).strip()
assert head == "64a02ff386e2985eef68ae00198b38e04f3c4a1f"
assert not subprocess.check_output(["git", "-C", str(AC), "status", "--porcelain"], text=True).strip()
checks["acolite_source_unchanged"] = head
checks["python"] = sys.version
checks["runtime_packages"] = {
    p: importlib.metadata.version(p) for p in ["numpy", "pandas", "netCDF4", "rasterio", "scipy", "GDAL"]
}
dump(OUT / "post-run-verification.json", checks)
print(json.dumps(checks, indent=2), flush=True)
