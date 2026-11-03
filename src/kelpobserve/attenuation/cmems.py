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
    flags: int = 0

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
    with xr.open_dataset(path) as dataset:
        return dataset.load()


def distances_km(latitudes: np.ndarray, longitudes: np.ndarray, site: tuple[float, float] = SITE) -> np.ndarray:
    lat, lon = np.meshgrid(latitudes, longitudes, indexing="ij")
    phi, lam = np.radians(lat), np.radians(lon)
    phi0, lam0 = np.radians(site[0]), np.radians(site[1])
    h = np.sin((phi - phi0) / 2) ** 2 + np.cos(phi0) * np.cos(phi) * np.sin((lam - lam0) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(h))


def validate_subset(dataset) -> None:
    """Require the reviewed L3 LAND flag schema and decoded physical units.

    An unfamiliar flag definition must be reviewed before this selector can use it.
    Unknown flag values are rejected per pixel by quality_masks().
    """
    for name in VARIABLES:
        if name not in dataset or set(dataset[name].dims) != {"time", "latitude", "longitude"}:
            raise ValueError(f"{name} must have time, latitude and longitude dimensions")
    for name, limit in (("latitude", 90), ("longitude", 180)):
        values = np.asarray(dataset[name].values)
        if (values.ndim != 1 or not np.all(np.isfinite(values)) or not values.size
                or len(np.unique(values)) != len(values) or np.any(np.abs(values) > limit)):
            raise ValueError(f"Invalid {name} coordinates")
    days = dataset.time.values.astype("datetime64[D]")
    if days.ndim != 1 or np.isnat(days).any() or len(np.unique(days)) != len(days):
        raise ValueError("Expected at most one daily record per date")
    attrs = dataset.flags.attrs
    try:
        masks = np.asarray(attrs.get("flag_masks", []), dtype=float).reshape(-1)
    except (ValueError, TypeError) as error:
        raise ValueError("Unsupported flag schema; expected LAND mask 1") from error
    if masks.tolist() != [1.0] or str(attrs.get("flag_meanings", "")).upper().split() != ["LAND"]:
        raise ValueError("Unsupported flag schema; expected LAND mask 1")
    if dataset.KD490.attrs.get("units") != "m-1" or dataset.KD490_uncertainty.attrs.get("units") != "%":
        raise ValueError("Expected decoded KD490 in m-1 and uncertainty in percent")


def _outside_valid_range(variable, values: np.ndarray) -> np.ndarray:
    """CF valid_min/max retain packed units after xarray decodes the data."""
    outside = np.zeros(values.shape, dtype=bool)
    scale = variable.encoding.get("scale_factor", 1.0)
    offset = variable.encoding.get("add_offset", 0.0)
    for name, compare in (("valid_min", np.less), ("valid_max", np.greater)):
        if name in variable.attrs:
            outside |= compare(values, float(variable.attrs[name]) * scale + offset)
    return outside


def quality_masks(dataset, index: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Product checks; reason masks may overlap. Caller first validates the schema.

    Uncertainty must be finite, nonnegative and uncensored. The product caps its
    error bars at 327.67%; a capped value cannot provide a quantified upper error.
    No application-specific uncertainty cutoff is imposed here.
    """
    def grid(name):
        return dataset[name].isel(time=index).transpose("latitude", "longitude").values

    kd, uncertainty, flags = (grid(name) for name in VARIABLES)
    reasons = {
        "land": flags == 1,
        "invalid_or_unknown_flags": ~np.isin(flags, [0, 1]),
        "invalid_kd490": ~np.isfinite(kd) | (kd <= 0) | _outside_valid_range(dataset.KD490, kd),
        "invalid_uncertainty": (~np.isfinite(uncertainty) | (uncertainty < 0)
                                | _outside_valid_range(dataset.KD490_uncertainty, uncertainty)),
        "censored_uncertainty": uncertainty >= 327.665,
    }
    usable = ~np.logical_or.reduce(list(reasons.values()))
    return usable, reasons


def search_geometry(dataset, max_km: float, max_days: int, site: tuple[float, float]) -> np.ndarray:
    validate_subset(dataset)
    if not np.isfinite(max_km) or max_km <= 0:
        raise ValueError("max_km must be finite and positive")
    if isinstance(max_days, bool) or not isinstance(max_days, int) or max_days < 0:
        raise ValueError("max_days must be a nonnegative integer")
    if (len(site) != 2 or not np.all(np.isfinite(site)) or abs(site[0]) > 90 or abs(site[1]) > 180):
        raise ValueError("Invalid site latitude/longitude")
    return distances_km(dataset.latitude.values, dataset.longitude.values, site)


def nearest_valid(dataset, when: date, max_km: float, max_days: int, site: tuple[float, float] = SITE) -> Match | None:
    """Nearest product-valid water pixel, preferring the target day then closest day.

    Equidistant dates prefer the earlier day; spatial ties follow latitude/longitude
    array order. This is a candidate for input qualification, not reef validation.
    """
    distance = search_geometry(dataset, max_km, max_days, site)
    days = dataset.time.values.astype("datetime64[D]")
    for offset in sorted(range(-max_days, max_days + 1), key=lambda d: (abs(d), d)):
        day = np.datetime64(when + timedelta(days=offset), "D")
        index = np.flatnonzero(days == day)
        if index.size == 0:
            continue
        usable, _ = quality_masks(dataset, int(index[0]))
        usable &= distance <= max_km
        if not usable.any():
            continue
        row, column = np.unravel_index(np.where(usable, distance, np.inf).argmin(), distance.shape)
        pixel = dataset.isel(time=int(index[0]), latitude=row, longitude=column)
        return Match(
            date=str(day), day_offset=offset, distance_km=round(float(distance[row, column]), 3),
            latitude=float(dataset.latitude.values[row]), longitude=float(dataset.longitude.values[column]),
            kd490=float(pixel.KD490.values), uncertainty_pct=float(pixel.KD490_uncertainty.values),
            valid_within_radius=int(usable.sum()), flags=int(pixel.flags.values),
        )
    return None


def write_provenance(record: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n")
