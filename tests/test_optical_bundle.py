import json

import numpy as np
import pytest
import xarray as xr

from kelpobserve.attenuation import optical_bundle as optical


def test_windows_cover_target_days_without_duplicate_downloads():
    assert optical.windows() == [
        ("2025-12-21", "2025-12-27"),
        ("2026-02-19", "2026-02-25"),
        ("2026-06-29", "2026-07-12"),
        ("2026-08-28", "2026-09-06"),
    ]


def ds():
    dims = ("time", "latitude", "longitude")
    return xr.Dataset(
        {
            "BBP": (dims, np.array([[[0.001, 0.002, 0.003]]]), {"units": "m-1"}),
            "BBP_uncertainty": (dims, np.array([[[20.0, 30.0, 40.0]]]), {"units": "%"}),
            "flags": (dims, np.array([[[0.0, 1.0, 0.0]]]), {"flag_masks": 1, "flag_meanings": "LAND"}),
        },
        coords={"time": [np.datetime64("2026-07-08")], "latitude": [38.4], "longitude": [-9.21, -9.20, -9.19]},
    )


def test_product_mask_excludes_land_and_censored_uncertainty():
    d = ds()
    d.BBP_uncertainty.values[0, 0, 2] = 327.67
    assert optical._eligible(d, "BBP").tolist() == [[[True, False, False]]]


def test_product_mask_requires_documented_units_and_flags():
    d = ds()
    d.BBP.attrs["units"] = "mg/m3"
    with pytest.raises(ValueError, match="units"):
        optical._eligible(d, "BBP")
    d = ds()
    d.flags.attrs["flag_masks"] = 2
    with pytest.raises(ValueError, match="schema"):
        optical._eligible(d, "BBP")


def test_bundle_rejects_missing_files_before_scoring(tmp_path):
    p = tmp_path / "manifest.json"
    p.write_text(json.dumps({"files": []}))
    with pytest.raises(ValueError, match="incomplete"):
        optical.audit(p)


def test_fetch_is_pinned_and_does_not_reuse_destination(tmp_path, monkeypatch):
    import copernicusmarine

    requests = []

    def fake_subset(**kwargs):
        requests.append(kwargs)
        (kwargs["output_directory"] / kwargs["output_filename"]).write_bytes(b"test")

    monkeypatch.setattr(copernicusmarine, "subset", fake_subset)
    destination = tmp_path / "bundle"
    manifest = optical.fetch(destination)
    assert len(requests) == 8
    assert all(r["dataset_version"] == "202603" and r["dataset_part"] == "default" for r in requests)
    assert all("flags" in r["variables"] for r in requests)
    assert len(json.loads(manifest.read_text())["files"]) == 8
    with pytest.raises(FileExistsError):
        optical.fetch(destination)


def test_fixed_pixel_audit_does_not_fill_missing_optics_from_another_pixel(tmp_path):
    import pandas as pd

    from kelpobserve.attenuation import cmems

    records = []
    for family, (dataset, variables) in optical.SOURCES.items():
        for start, end in optical.windows():
            times = pd.date_range(start, end)
            shape = (len(times), 2, 2)
            data = {}
            for v in variables:
                data[v] = (
                    ("time", "latitude", "longitude"),
                    np.full(shape, 0.1),
                    {"units": next(iter(optical.UNITS[v]))},
                )
                data[v + "_uncertainty"] = (("time", "latitude", "longitude"), np.full(shape, 25.0), {"units": "%"})
            data["flags"] = (
                ("time", "latitude", "longitude"),
                np.zeros(shape),
                {"flag_masks": 1, "flag_meanings": "LAND"},
            )
            d = xr.Dataset(data, coords={"time": times, "latitude": [38.40, 38.41], "longitude": [-9.21, -9.20]})
            if family == "optics" and "2026-07-08" in times:
                d.BBP.loc[{"time": "2026-07-08", "latitude": 38.40, "longitude": -9.21}] = np.nan
            name = f"{family}-{start}-{end}.nc"
            p = tmp_path / name
            d.to_netcdf(p)
            records.append(
                {
                    "file": name,
                    "sha256": cmems.sha256(p),
                    "request": {
                        "dataset_id": dataset,
                        "dataset_version": "202603",
                        "dataset_part": "default",
                        "start_datetime": start + "T00:00:00",
                        "end_datetime": end + "T23:59:59",
                    },
                }
            )
    p = tmp_path / "manifest.json"
    p.write_text(json.dumps({"files": records}))
    kd = {
        "targets": [
            {
                "target_date": day,
                "candidate": {"date": day, "latitude": 38.40, "longitude": -9.21, "distance_km": 1.0, "kd490": 0.1},
            }
            for day in optical.qualification.TARGETS
        ]
    }
    r = optical.audit(p, kd)
    july = next(t for t in r["targets"] if t["date"] == "2026-07-08")
    assert july["variables"]["BBP"]["value"] is None
    assert july["joint_product_valid_at_fixed_pixel"] is False
    assert july["joint_product_valid_pixels_within_5km"] == 3
    assert july["kd490_relative_change_pct"] == 0
    json.dumps(r, allow_nan=False)
    records[0]["sha256"] = "wrong"
    p.write_text(json.dumps({"files": records}))
    with pytest.raises(ValueError, match="checksum"):
        optical.audit(p, kd)
