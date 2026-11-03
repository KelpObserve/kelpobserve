"""Bounded July atmospheric comparison; cached TOA repacking, no habitat labels."""

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from netCDF4 import Dataset
from pyproj import Transformer
from rasterio.transform import Affine
from rasterio.windows import Window, from_bounds

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from kelpobserve import config

PATHS = config.load()
sys.path.insert(0, str(PATHS.require("acolite")))
import acolite as ac

CACHE = PATHS.require("cache")
ROOT = CACHE / "optical-atmosphere-20261103"
S2 = "s2-20260702-S2B-R037"
PN = "pneo-20260708"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(c)
    return h.hexdigest()


def dump(p, x):
    p.write_text(json.dumps(x, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else str(x)) + "\n")


def export(ncf, scene, datasets):
    """Keep the exact source pixel grid; avoid CF/GDAL half-pixel convention shifts."""
    pr = json.loads((ROOT / "inputs" / scene / "provenance.json").read_text())
    with Dataset(ncf) as nc:
        tags = {"NC_GLOBAL#" + k: str(nc.getncattr(k)) for k in nc.ncattrs()}
        for name in datasets:
            if name not in nc.variables:
                continue
            array = np.asarray(nc[name][:].filled(np.nan)).squeeze()
            assert list(array.shape) == pr["shape"]
            dest = Path(ncf).with_name(Path(ncf).stem + "_" + name + ".tif")
            with rasterio.open(
                dest,
                "w",
                driver="GTiff",
                height=array.shape[0],
                width=array.shape[1],
                count=1,
                dtype="float32",
                nodata=np.nan,
                crs=pr["crs"],
                transform=Affine(*pr["transform"][:6]),
                compress="deflate",
                predictor=3,
                tiled=True,
            ) as dst:
                dst.write(array.astype("float32"), 1)
                dst.update_tags(**tags)


def pack(scene):
    out = ROOT / "inputs" / scene
    out.mkdir(parents=True, exist_ok=False)
    files = sorted((CACHE / "acolite" / scene).glob("*rhot_*.tif"))
    sample = files[0]
    with rasterio.open(sample) as ds:
        w = Window(0, 0, ds.width, ds.height)
        if scene == S2:
            # >5 km of real measured context around every evaluation cell.
            w = from_bounds(475110, 4244730, 491030, 4256690, ds.transform).round_offsets().round_lengths()
        assert (
            w.col_off >= 0 and w.row_off >= 0 and w.col_off + w.width <= ds.width and w.row_off + w.height <= ds.height
        )
        transform = ds.window_transform(w)
        crs = ds.crs
        tags = ds.tags()
        height, width = int(w.height), int(w.width)
    dct = ac.shared.projection_read(str(sample))
    dct.update(
        xrange=(transform.c, transform.c + width * transform.a),
        yrange=(transform.f, transform.f + height * transform.e),
        xdim=width,
        ydim=height,
        dimensions=(width, height),
    )
    ncf = out / (scene + "_L1R.nc")
    gem = ac.gem.gem(str(ncf), new=True)
    meta = {}
    for key in [
        "sensor",
        "isodate",
        "sza",
        "vza",
        "raa",
        "saa",
        "vaa",
        "se_distance",
        "mus",
        "pressure",
        "uoz",
        "uwv",
        "wind",
    ]:
        value = tags.get("NC_GLOBAL#" + key)
        if value is None:
            continue
        if key not in ["sensor", "isodate"]:
            value = float(value)
        meta[key] = value
    meta.update(
        acolite_file_type="L1R",
        scene_pixel_size=(transform.a, transform.e),
        scene_dims=(height, width),
        data_dimensions=(height, width),
        output_name=scene + "_L1R",
        reconstruction="Cached original TOA and geometry arrays, exact subset; no new radiometric conversion",
    )
    gem.gatts = meta
    gem.nc_projection = ac.shared.projection_netcdf(dct)
    records = {}
    for p in files + sorted((CACHE / "acolite" / scene).glob("*L2R_*.tif")):
        name = p.stem.split("_L2R_")[-1]
        if not (name.startswith("rhot_") or name in ["sza", "vza", "raa"]):
            continue
        if str(p) in records:
            continue
        with rasterio.open(p) as ds:
            assert ds.crs == crs and ds.transform.a == transform.a
            data = ds.read(1, window=w).astype("float32")
            if ds.nodata is not None:
                data[data == ds.nodata] = np.nan
            att = {}
            if name.startswith("rhot_"):
                att["wavelength"] = float(name.split("_")[-1])
            gem.write(name, data, ds_att=att)
        records[str(p)] = {"sha256": sha(p), "dataset": name}
    # Geolocation is derived from unchanged pixel centres; explicit projection stored.
    xx, yy = np.meshgrid(
        transform.c + (np.arange(width) + 0.5) * transform.a, transform.f + (np.arange(height) + 0.5) * transform.e
    )
    lon, lat = Transformer.from_crs(crs, 4326, always_xy=True).transform(xx, yy)
    gem.write("lon", lon.astype("float32"))
    gem.write("lat", lat.astype("float32"))
    gem.close()
    with Dataset(ncf) as nc:
        for p, v in records.items():
            with rasterio.open(p) as ds:
                x = ds.read(1, window=w).astype("float32")
                if ds.nodata is not None:
                    x[x == ds.nodata] = np.nan
            y = np.asarray(nc[v["dataset"]][:].filled(np.nan))
            assert np.array_equal(x, y, equal_nan=True), v["dataset"]
    dump(
        out / "provenance.json",
        {
            "scene": scene,
            "input_files": records,
            "window": list(w.flatten()),
            "shape": [height, width],
            "transform": list(transform),
            "crs": str(crs),
            "all_arrays_exact": True,
            "l1r_sha256": sha(ncf),
        },
    )
    print("Packed", ncf, height, width, flush=True)


def process(branch):
    scene = S2 if branch.startswith("s2") else PN
    out = ROOT / branch
    out.mkdir(parents=True, exist_ok=False)
    ncf = ROOT / "inputs" / scene / (scene + "_L1R.nc")
    sensor = "S2B_MSI" if scene == S2 else "PNEO3"
    settings = {
        "output": str(out),
        "ancillary_data": False,
        "dem_pressure": False,
        "dem_pressure_resolved": False,
        "dsf_interface_reflectance": False,
        "dsf_residual_glint_correction": False,
        "l2w_parameters": None,
        "l1r_delete_netcdf": False,
        "l2r_delete_netcdf": False,
        "output_geolocation": True,
        "verbosity": 2,
        "luts": ["ACOLITE-LUT-202110-MOD1", "ACOLITE-LUT-202110-MOD2"],
        "radcor_kernel_radius": 5.0,
        "radcor_edge_extend": False,
        "radcor_crop_centre": False,
        "radcor_write_rhosu": True,
        "radcor_write_rhoe": False,
        "radcor_diagnostic_plots": False,
        "radcor_development": False,
        "radcor_write_rhot": False,
    }
    if branch.startswith("s2"):
        settings["atmospheric_correction_method"] = "radcor"
        if branch == "s2_radcor_fixed":
            settings.update(radcor_force_model="M", radcor_force_aot=0.080870785)
    elif branch == "pneo_marine_control":
        settings.update(
            dsf_aot_estimate="fixed",
            dsf_spectrum_option="percentile",
            dsf_percentile=1.0,
            dsf_fixed_aot=0.078115806,
            dsf_fixed_lut="ACOLITE-LUT-202110-MOD2",
        )
    elif branch == "pneo_continental":
        settings.update(
            dsf_aot_estimate="fixed",
            dsf_spectrum_option="percentile",
            dsf_percentile=1.0,
            luts=["ACOLITE-LUT-202110-MOD1"],
        )
    else:
        raise ValueError(branch)
    ac.settings["user"] = settings.copy()
    ac.settings["run"] = dict(ac.settings["defaults"])
    ac.settings["run"].update(ac.acolite.settings.parse(sensor, settings=settings))
    dump(out / "requested-settings.json", settings)
    t = time.monotonic()
    if branch.startswith("s2"):
        result = ac.adjacency.radcor.radcor(str(ncf))
        files = list(out.glob("*L2R*.nc"))
        assert len(files) == 1, (result, files)
        result = files[0]
    else:
        result, effective = ac.acolite.acolite_l2r(str(ncf))
        result = Path(result)
        dump(out / "effective-settings.json", effective)
    assert result.exists()
    # Export only the bands needed by the predeclared diagnostic.
    wanted = (442, 492, 559, 665, 833, 1610) if scene == S2 else (439, 486, 563, 655, 826)
    ds = [b + "_" + str(w) for b in ("rhos", "rhosu") for w in wanted]
    export(result, scene, ds)
    if scene == S2:
        glint = out / "glint"
        glint.mkdir()
        glint_settings = dict(ac.settings["run"])
        glint_settings.update(output=str(glint), glint_write_rhog_ref=True)
        ac.glint.default(str(result), settings=glint_settings, new_file=True)
        gf = list(glint.glob("*.nc"))
        assert len(gf) == 1
        export(gf[0], scene, ds)
    with Dataset(result) as nc:
        gatts = {
            k: nc.getncattr(k)
            for k in nc.ncattrs()
            if k in ["sensor", "sza", "vza", "raa", "ac_aot_550", "ac_model", "radcor_version", "glint_correction"]
        }
    dump(
        out / "completed.json",
        {
            "branch": branch,
            "seconds": time.monotonic() - t,
            "l2r": str(result),
            "l1r_sha256": sha(ncf),
            "metadata": gatts,
            "processing_source_sha256": sha(Path(__file__)),
        },
    )
    print("Completed", branch, gatts, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["pack", "process"])
    p.add_argument("target")
    p.add_argument("--cache-root", type=Path, default=ROOT)
    a = p.parse_args()
    ROOT = a.cache_root
    (pack if a.action == "pack" else process)(a.target)
