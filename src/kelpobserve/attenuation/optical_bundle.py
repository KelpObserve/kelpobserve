"""Version-pinned Atlantic optical inputs for the Sesimbra development audit.

Fixed KD490 candidate coordinates are scored for all optical variables. A missing
value is retained as missing; no variable chooses a different pixel to improve support.
This module does not fit water spectra or score habitat predictions.
"""

from __future__ import annotations

import json
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

import numpy as np

from kelpobserve.attenuation import cmems, qualification

SOURCES = {
    "transp": ("cmems_obs-oc_atl_bgc-transp_my_l3-multi-1km_P1D", ["KD490", "SPM"]),
    "optics": ("cmems_obs-oc_atl_bgc-optics_my_l3-multi-1km_P1D", ["BBP", "CDM"]),
}
UNITS = {"KD490": {"m-1"}, "BBP": {"m-1"}, "CDM": {"m-1"}, "SPM": {"g/m3", "g m-3"}}
WAVELENGTHS = {"KD490": 490, "BBP": 443, "CDM": 443, "SPM": None}


def windows() -> list[tuple[str, str]]:
    ranges = []
    for target in sorted(qualification.TARGETS):
        day = date.fromisoformat(target)
        start, end = day - timedelta(days=3), day + timedelta(days=3)
        if ranges and start <= ranges[-1][1] + timedelta(days=1):
            ranges[-1] = (ranges[-1][0], max(end, ranges[-1][1]))
        else:
            ranges.append((start, end))
    return [(a.isoformat(), b.isoformat()) for a, b in ranges]


def fetch(destination: Path, version: str = "202603", part: str = "default") -> Path:
    """Create a new local bundle; no stream fallback or overwrite of previous inputs."""
    import copernicusmarine

    destination.mkdir(parents=True, exist_ok=False)
    records = []
    for family, (dataset, variables) in SOURCES.items():
        for start, end in windows():
            name = f"{family}-{start}-{end}.nc"
            selected = [item for v in variables for item in (v, v + "_uncertainty")] + ["flags"]
            request = dict(
                dataset_id=dataset,
                dataset_version=version,
                dataset_part=part,
                variables=selected,
                **cmems.REGION,
                start_datetime=start + "T00:00:00",
                end_datetime=end + "T23:59:59",
                coordinates_selection_method="inside",
                disable_progress_bar=True,
            )
            copernicusmarine.subset(**request, output_directory=destination, output_filename=name, overwrite=False)
            path = destination / name
            records.append(
                {
                    "file": name,
                    "sha256": cmems.sha256(path),
                    "bytes": path.stat().st_size,
                    "request": request,
                    "product": cmems.PRODUCT,
                    "toolbox": copernicusmarine.__version__,
                    "fetched_at": datetime.now(UTC).isoformat(timespec="seconds"),
                }
            )
            (destination / "manifest.json").write_text(
                json.dumps({"schema_version": "1.0", "files": records}, indent=2) + "\n"
            )
    return destination / "manifest.json"


def _eligible(ds, variable: str):
    uncertainty = variable + "_uncertainty"
    for name in (variable, uncertainty, "flags"):
        if name not in ds or set(ds[name].dims) != {"time", "latitude", "longitude"}:
            raise ValueError(f"Missing or malformed {name}")
    masks = np.asarray(ds.flags.attrs.get("flag_masks", []), dtype=float).reshape(-1)
    if masks.tolist() != [1.0] or str(ds.flags.attrs.get("flag_meanings", "")).upper().split() != ["LAND"]:
        raise ValueError("Unreviewed flag schema")
    if ds[variable].attrs.get("units") not in UNITS[variable] or ds[uncertainty].attrs.get("units") != "%":
        raise ValueError(f"Unexpected units for {variable}")
    da, ua = ds[variable], ds[uncertainty]
    value, error = da.values, ua.values
    return (
        np.isfinite(value)
        & (value > 0)
        & ~cmems._outside_valid_range(da, value)
        & (ds.flags.values == 0)
        & np.isfinite(error)
        & (error >= 0)
        & (error < 327.665)
        & ~cmems._outside_valid_range(ua, error)
    )


def _number(value):
    return float(value) if np.isfinite(value) else None


def audit(manifest_path: Path, kd_audit: dict | None = None) -> dict:
    import xarray as xr

    manifest = json.loads(manifest_path.read_text())
    expected_requests = {
        (dataset, start + "T00:00:00", end + "T23:59:59") for dataset, _ in SOURCES.values() for start, end in windows()
    }
    observed_requests = []
    grouped = {key: [] for key in SOURCES}
    for record in manifest["files"]:
        request = record["request"]
        observed_requests.append((request["dataset_id"], request["start_datetime"], request["end_datetime"]))
        if not request.get("dataset_version") or not request.get("dataset_part"):
            raise ValueError("Each input requires a pinned dataset version and part")
        path = manifest_path.parent / record["file"]
        if path.resolve().parent != manifest_path.parent.resolve() or cmems.sha256(path) != record["sha256"]:
            raise ValueError(f"Bundle checksum/path mismatch: {record['file']}")
        family = next((key for key, (dataset, _) in SOURCES.items() if request["dataset_id"] == dataset), None)
        if family is None:
            raise ValueError("Unexpected source dataset")
        with xr.open_dataset(path) as ds:
            grouped[family].append(ds.load().transpose("time", "latitude", "longitude"))
    if set(observed_requests) != expected_requests or len(observed_requests) != len(expected_requests):
        raise ValueError("Bundle is incomplete or contains duplicate/unexpected requests")
    datasets = {family: xr.concat(parts, dim="time").sortby("time") for family, parts in grouped.items()}
    for coord in ("time", "latitude", "longitude"):
        if not np.array_equal(datasets["transp"][coord], datasets["optics"][coord]):
            raise ValueError(f"Optics and transparency {coord} grids differ")
    grid = datasets["transp"]
    cmems.validate_subset(grid[["KD490", "KD490_uncertainty", "flags"]])
    valid = {v: _eligible(datasets[family], v) for family, (_, names) in SOURCES.items() for v in names}
    kd_audit = kd_audit or qualification.build()
    distance = cmems.distances_km(grid.latitude.values, grid.longitude.values)
    targets = []
    for target in kd_audit["targets"]:
        candidate = target["candidate"]
        day = target["target_date"]
        if candidate is None or candidate["date"] != day:
            targets.append({"date": day, "status": "no_same_day_fixed_candidate"})
            continue
        ti = int(np.flatnonzero(grid.time.values.astype("datetime64[D]") == np.datetime64(day))[0])
        row = int(np.flatnonzero(grid.latitude.values == candidate["latitude"])[0])
        col = int(np.flatnonzero(grid.longitude.values == candidate["longitude"])[0])
        values = {}
        for family, (_, names) in SOURCES.items():
            ds = datasets[family]
            for variable in names:
                values[variable] = {
                    "value": _number(ds[variable].values[ti, row, col]),
                    "uncertainty_pct": _number(ds[variable + "_uncertainty"].values[ti, row, col]),
                    "flags": _number(ds.flags.values[ti, row, col]),
                    "product_checks_passed": bool(valid[variable][ti, row, col]),
                    "units": ds[variable].attrs["units"],
                    "wavelength_nm": WAVELENGTHS[variable],
                }
        joint = np.logical_and.reduce([mask[ti] for mask in valid.values()])
        original = candidate["kd490"]
        current = values["KD490"]["value"]
        targets.append(
            {
                "date": day,
                "latitude": candidate["latitude"],
                "longitude": candidate["longitude"],
                "distance_km": candidate["distance_km"],
                "variables": values,
                "joint_product_valid_at_fixed_pixel": all(v["product_checks_passed"] for v in values.values()),
                "joint_product_valid_pixels_within_5km": int((joint & (distance <= 5)).sum()),
                "kd490_original_frozen": original,
                "kd490_relative_change_pct": 100 * (current / original - 1) if current is not None else None,
                "status": "candidate_only",
                "reef_representativeness": "unestablished",
            }
        )
    return {
        "created_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "schema_version": "1.0",
        "manifest": manifest,
        "manifest_sha256": cmems.sha256(manifest_path),
        "source_sha256": {
            name: cmems.sha256(path)
            for name, path in {
                "optical_bundle.py": Path(__file__),
                "cmems.py": Path(cmems.__file__),
                "qualification.py": Path(qualification.__file__),
            }.items()
        },
        "selection": "same coordinates and target date as frozen KD490 qualification; no per-variable reselection",
        "uncertainty_interpretation": "provider percentage error bars; no Gaussian or independent-error assumption",
        "targets": targets,
        "scope": "Product checks and source comparison only. No water inversion or habitat scores.",
        "remaining_requirements": [
            "reef representativeness",
            "spectral slopes and phytoplankton absorption",
            "joint product consistency and correlated errors",
            "registered fallback and score rules",
        ],
    }


def markdown(record: dict) -> str:
    lines = [
        "# Sesimbra optical input bundle",
        "",
        f"Generated: {record['created_at']}",
        "",
        (
            "KD490, particulate backscatter (BBP), coloured dissolved and detrital absorption (CDM), and suspended "
            "matter (SPM) were downloaded with flags and provider uncertainties. Dataset versions and parts are "
            "pinned in the manifest; the source files are immutable local inputs. BBP and CDM refer to 443 nm."
        ),
        "",
        (
            "Each row uses the same pixel as the repaired KD490 audit. Values in parentheses are provider "
            "percentage error bars. They have no assigned standard deviation or confidence level here."
        ),
        "",
        "| Date | KD490 (m⁻¹) | BBP443 (m⁻¹) | CDM443 (m⁻¹) | SPM (g/m³) | All four pass product checks |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for target in record["targets"]:
        if "variables" not in target:
            lines.append(f"| {target['date']} | — | — | — | — | {target['status']} |")
            continue
        columns = []
        for variable in ["KD490", "BBP", "CDM", "SPM"]:
            item = target["variables"][variable]
            value, uncertainty = item["value"], item["uncertainty_pct"]
            columns.append(
                f"{value:.6g} ({uncertainty:.2f}%)" if value is not None and uncertainty is not None else "missing"
            )
        lines.append(
            "| " + " | ".join([target["date"], *columns, str(target["joint_product_valid_at_fixed_pixel"])]) + " |"
        )
    lines += [
        "",
        (
            "Passing product checks does not establish that these 1 km pixels describe the reef water. "
            "BBP and CDM constrain parts of the optical model, while spectral slopes and phytoplankton "
            "absorption remain assumptions. Correlated product errors and mutual consistency must be tested "
            "before using these quantities together in inversion. No r7 habitat score was computed."
        ),
        "",
    ]
    return "\n".join(lines)
