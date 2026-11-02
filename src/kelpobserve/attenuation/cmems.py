"""Copernicus Marine KD490 at the Sesimbra reef: freeze a subset, then read the nearest valid pixel.

The Atlantic multi-year ocean-colour product (OCEANCOLOUR_ATL_BGC_L3_MY_009_113) carries daily
KD490 at 1 km from merged sensors, with a per-pixel uncertainty in percent; the Atlantic has no
300 m transparency dataset. Experiments read a frozen local subset and never the live service:
dataset parts are retired (the catalogue gives 24 November 2026 for part "default" of version
202603) and updated in place, so a fetch records the catalogue's dates, the file's SHA-256 and a
digest of its values. A later fetch is checked against that record and never replaces it.

    kelpobserve kd490 fetch     # needs COPERNICUSMARINE_SERVICE_USERNAME and _PASSWORD
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict, dataclass
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

import numpy as np

PRODUCT = "OCEANCOLOUR_ATL_BGC_L3_MY_009_113"
DATASET = "cmems_obs-oc_atl_bgc-transp_my_l3-multi-1km_P1D"
VARIABLES = ("KD490", "KD490_uncertainty", "flags")
SITE = (38.409, -9.202)
REGION = {"minimum_longitude": -9.42, "maximum_longitude": -8.98, "minimum_latitude": 38.24, "maximum_latitude": 38.56}
WINDOW = ("2025-12-01", "2026-09-30")
CREDENTIALS = ("COPERNICUSMARINE_SERVICE_USERNAME", "COPERNICUSMARINE_SERVICE_PASSWORD")


@dataclass(frozen=True)
class Match:
    date: str
    day_offset: int
    distance_km: float
    latitude: float
    longitude: float
    kd490: float
    uncertainty_pct: float
    valid_within_radius: int

    def to_json(self) -> dict:
        return asdict(self)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def catalogue_entry() -> dict:
    """The product's DOI and the dataset's versions and parts, with their dates, as the catalogue lists them."""
    import copernicusmarine

    catalogue = copernicusmarine.describe(dataset_id=DATASET, disable_progress_bar=True).model_dump(mode="json")
    entry: dict = {"doi": None, "versions": []}
    for product in catalogue.get("products", []):
        entry["doi"] = entry["doi"] or product.get("digital_object_identifier")
        for dataset in product.get("datasets", []):
            if dataset.get("dataset_id") != DATASET:
                continue
            for version in dataset.get("versions", []):
                for part in version.get("parts", []):
                    entry["versions"].append({
                        "label": version.get("label"), "part": part.get("name"),
                        **{key: part.get(key) for key in ("released_date", "retired_date", "arco_updated_date")},
                    })
    return entry


def data_digest(path: Path) -> str:
    """SHA-256 of the subset's variables and coordinates, independent of how the file encodes them."""
    import xarray as xr

    digest = hashlib.sha256()
    with xr.open_dataset(path) as dataset:
        for name in sorted([*VARIABLES, *dataset.coords]):
            values = np.asarray(dataset[name].values)
            if values.dtype.kind == "M":
                values = values.astype("datetime64[ns]").astype("int64")
            elif values.dtype.kind == "f":
                values = np.where(np.isnan(values), np.inf, values).astype("float64")
            digest.update(name.encode())
            digest.update(str(values.shape).encode())
            digest.update(np.ascontiguousarray(values).tobytes())
    return digest.hexdigest()


def fetch(destination: Path, start: str = WINDOW[0], end: str = WINDOW[1]) -> dict:
    """Download the subset once and describe it; an existing frozen subset is never replaced."""
    import copernicusmarine

    missing = [name for name in CREDENTIALS if not os.environ.get(name)]
    if missing:
        raise SystemExit(f"Set {' and '.join(missing)} (for example from the .env file the benchmark uses).")
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / f"{DATASET}_{start}_{end}.nc"
    if target.exists():
        raise SystemExit(f"A frozen subset exists already: {target}")
    copernicusmarine.subset(
        dataset_id=DATASET, variables=list(VARIABLES), start_datetime=start, end_datetime=end,
        output_directory=str(destination), output_filename=target.name, disable_progress_bar=True, **REGION,
    )
    catalogue = catalogue_entry()
    return {
        "product": PRODUCT,
        "product_doi": catalogue["doi"],
        "dataset_id": DATASET,
        "dataset_versions_listed": catalogue["versions"],
        "variables": list(VARIABLES),
        "region": REGION,
        "start": start,
        "end": end,
        "file": target.name,
        "sha256": sha256(target),
        "data_sha256": data_digest(target),
        "bytes": target.stat().st_size,
        "fetched_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "toolbox": copernicusmarine.__version__,
        "acknowledgement": "Generated using E.U. Copernicus Marine Service Information",
    }


def matches(path: Path, record: dict) -> bool:
    """True when the file is the recorded subset, byte for byte or value for value."""
    return sha256(path) == record["sha256"] or data_digest(path) == record.get("data_sha256")


def open_subset(path: Path, record: dict | None = None):
    import xarray as xr

    if record and not matches(path, record):
        raise SystemExit(f"The subset at {path} is not the one recorded in the evidence file.")
    return xr.open_dataset(path).load()


def distances_km(latitudes: np.ndarray, longitudes: np.ndarray, site: tuple[float, float] = SITE) -> np.ndarray:
    lat, lon = np.meshgrid(latitudes, longitudes, indexing="ij")
    phi, lam = np.radians(lat), np.radians(lon)
    phi0, lam0 = np.radians(site[0]), np.radians(site[1])
    h = np.sin((phi - phi0) / 2) ** 2 + np.cos(phi0) * np.cos(phi) * np.sin((lam - lam0) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(h))


def nearest_valid(dataset, when: date, max_km: float, max_days: int, site: tuple[float, float] = SITE) -> Match | None:
    """The valid KD490 pixel nearest the site, on the day itself or the closest day within max_days."""
    distance = distances_km(dataset.latitude.values, dataset.longitude.values, site)
    days = dataset.time.values.astype("datetime64[D]")
    for offset in sorted(range(-max_days, max_days + 1), key=lambda d: (abs(d), d)):
        day = np.datetime64(when + timedelta(days=offset), "D")
        index = np.flatnonzero(days == day)
        if index.size == 0:
            continue
        kd = dataset.KD490.isel(time=int(index[0])).values
        usable = np.isfinite(kd) & (distance <= max_km)
        if not usable.any():
            continue
        row, column = np.unravel_index(np.where(usable, distance, np.inf).argmin(), distance.shape)
        uncertainty = dataset.KD490_uncertainty.isel(time=int(index[0])).values[row, column]
        return Match(
            date=str(day), day_offset=offset, distance_km=round(float(distance[row, column]), 3),
            latitude=float(dataset.latitude.values[row]), longitude=float(dataset.longitude.values[column]),
            kd490=float(kd[row, column]), uncertainty_pct=float(uncertainty),
            valid_within_radius=int(usable.sum()),
        )
    return None


def write_provenance(record: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n")
