"""Modelled sea level from Copernicus Marine, on a Lowest Astronomical Tide datum.

The Iberia–Biscay–Ireland (IBI) analysis and forecast gives hourly sea surface
height (``zos``) at about 2.5 km. That is the total level, astronomical tide plus
meteorological surge, on the model's own datum. Bathymetry at the Iberian sites
is referenced to LAT, so the model level has to be re-expressed above LAT before
it can correct a depth.

LAT here is defined as the charts define it: the lowest level of the
*astronomical* tide over a lunar nodal cycle, from a harmonic fit to the model
series (:mod:`kelpobserve.coastal.bathymetry.harmonic`). The lowest value of the
record itself is kept as a labelled alternative. It depends on the record's
length and on whatever storm surge fell inside it, which is why it is only a
proxy.

The nearest model water cell can be kilometres from a site, and a 2.5 km model
does not resolve a lagoon behind barrier islands. For those sites fit
:class:`~kelpobserve.coastal.bathymetry.harmonic.HarmonicConstants` to this
series at the inlet and use a harmonic provider instead.

Time stamps of mean datasets
----------------------------
``copernicusmarine.open_dataset`` serves the ARCO form of a dataset, which aligns
time stamps to the *start* of each averaging interval (Copernicus Marine Toolbox
FAQ). An hourly-mean (``PT1H-m``) value stamped 11:00 is the mean over 11:00-12:00
and belongs at 11:30. :func:`fetch_sea_level` moves the stamps of mean datasets to
the interval centre. Read unshifted, the tide is half an hour early: against the
REFMAR gauge at Le Conquet the fit needed a 30-32 minute lag in every summer
2019-2026, and after the shift it needs none.
"""

from __future__ import annotations

import datetime as dt
import math
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from kelpobserve.coastal.bathymetry.harmonic import (
    DEFAULT_CONSTITUENTS,
    fit_harmonics,
    lowest_astronomical_tide,
)
from kelpobserve.coastal.bathymetry.tide import TideCorrection

IBI_SEA_LEVEL_DATASET = "cmems_mod_ibi_phy_anfc_0.027deg-2D_PT1H-m"
IBI_SEA_LEVEL_VARIABLE = "zos"

#: Uncertainty stated for a modelled level on a harmonic LAT, absent a gauge to
#: measure it: the Sesimbra benchmark registered 0.3 m for its datum proxy.
DEFAULT_UNCERTAINTY_M = 0.3

OpenDataset = Callable[..., Any]

_RESOLUTION = re.compile(r"_P(?P<t>T)?(?P<n>\d+)(?P<unit>[DHM])-(?P<kind>[mi])$")


def interval_centre_shift(dataset_id: str) -> np.timedelta64:
    """How far to move a dataset's time stamps to the centre of its averaging interval.

    Half the interval for a Copernicus Marine mean dataset (``_PT1H-m``, ``_PT15M-m``,
    ``_P1D-m``), whose ARCO stamps mark the start of the interval; zero for an
    instantaneous (``-i``) dataset or an id without a resolution suffix.
    """
    match = _RESOLUTION.search(dataset_id)
    if match is None or match["kind"] == "i":
        return np.timedelta64(0, "s")
    n = int(match["n"])
    if match["t"] and match["unit"] == "H":
        seconds = n * 3600
    elif match["t"] and match["unit"] == "M":
        seconds = n * 60
    elif not match["t"] and match["unit"] == "D":
        seconds = n * 86400
    else:
        raise ValueError(f"No interval-centre rule for the resolution of {dataset_id!r}.")
    return np.timedelta64(seconds // 2, "s")


@dataclass(frozen=True)
class SeaLevelSeries:
    """An hourly model sea-level series at one model cell, on the model datum."""

    times: np.ndarray  # datetime64[s], UTC
    levels_m: np.ndarray
    cell_lon: float
    cell_lat: float
    distance_km: float
    dataset_id: str
    variable: str
    #: Seconds added to the served stamps to put them at the interval centre.
    time_shift_s: float = 0.0

    def __post_init__(self) -> None:
        if self.times.shape != self.levels_m.shape or self.times.size == 0:
            raise ValueError("A sea-level series needs matching, non-empty times and levels.")

    @property
    def station(self) -> str:
        ns = "N" if self.cell_lat >= 0 else "S"
        ew = "E" if self.cell_lon >= 0 else "W"
        return (
            f"{self.dataset_id} cell {abs(self.cell_lat):.4f}{ns} {abs(self.cell_lon):.4f}{ew}, "
            f"{self.distance_km:.1f} km from the site"
        )

    def level_at(self, when: dt.datetime) -> float:
        """Linear interpolation to ``when``; refuses times outside the record."""
        if when.tzinfo is None:
            raise ValueError("Tide times must be timezone-aware; got a naive datetime.")
        stamp = np.datetime64(when.astimezone(dt.UTC).replace(tzinfo=None), "s")
        if stamp < self.times[0] or stamp > self.times[-1]:
            raise ValueError(
                f"{when.isoformat()} is outside the model record "
                f"({self.times[0]} to {self.times[-1]})."
            )
        seconds = (self.times - self.times[0]) / np.timedelta64(1, "s")
        offset = (stamp - self.times[0]) / np.timedelta64(1, "s")
        return float(np.interp(offset, seconds, self.levels_m))


def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 6371.0 * math.asin(math.sqrt(a))


def nearest_water_cell(
    lons: np.ndarray, lats: np.ndarray, field: np.ndarray, lon: float, lat: float
) -> tuple[int, int, float]:
    """Row, column and distance in km of the nearest cell with a finite value."""
    best: tuple[int, int, float] | None = None
    for i, cell_lat in enumerate(lats):
        for j, cell_lon in enumerate(lons):
            if not np.isfinite(field[i, j]):
                continue
            km = haversine_km(lon, lat, float(cell_lon), float(cell_lat))
            if best is None or km < best[2]:
                best = (i, j, km)
    if best is None:
        raise ValueError(f"No model water cell near {lat:.4f}, {lon:.4f}.")
    return best


def fetch_sea_level(
    lon: float,
    lat: float,
    start: dt.datetime,
    end: dt.datetime,
    *,
    dataset_id: str = IBI_SEA_LEVEL_DATASET,
    variable: str = IBI_SEA_LEVEL_VARIABLE,
    search_deg: float = 0.2,
    open_dataset: OpenDataset | None = None,
) -> SeaLevelSeries:
    """Hourly model sea level at the water cell nearest to ``lon``, ``lat``.

    Mean datasets are returned with their stamps at the interval centre
    (:func:`interval_centre_shift`); ``SeaLevelSeries.time_shift_s`` records the move.

    ``open_dataset`` defaults to ``copernicusmarine.open_dataset``, which reads
    credentials from ``COPERNICUSMARINE_SERVICE_USERNAME`` and
    ``COPERNICUSMARINE_SERVICE_PASSWORD``.
    """
    if open_dataset is None:
        import copernicusmarine

        open_dataset = copernicusmarine.open_dataset
    data = open_dataset(
        dataset_id=dataset_id,
        variables=[variable],
        minimum_longitude=lon - search_deg,
        maximum_longitude=lon + search_deg,
        minimum_latitude=lat - search_deg,
        maximum_latitude=lat + search_deg,
        start_datetime=start.astimezone(dt.UTC).strftime("%Y-%m-%dT%H:%M:%S"),
        end_datetime=end.astimezone(dt.UTC).strftime("%Y-%m-%dT%H:%M:%S"),
    )
    grid = data[variable]
    lons, lats = grid.longitude.values, grid.latitude.values
    i, j, km = nearest_water_cell(lons, lats, grid.isel(time=0).values, lon, lat)
    cell = grid.isel(latitude=i, longitude=j)
    levels = np.asarray(cell.values, dtype=np.float64)
    shift = interval_centre_shift(dataset_id)
    times = np.asarray(cell.time.values).astype("datetime64[s]") + shift
    keep = np.isfinite(levels)
    return SeaLevelSeries(
        times=times[keep],
        levels_m=levels[keep],
        cell_lon=float(lons[j]),
        cell_lat=float(lats[i]),
        distance_km=km,
        dataset_id=dataset_id,
        variable=variable,
        time_shift_s=float(shift / np.timedelta64(1, "s")),
    )


@dataclass
class CMEMSTide:
    """Modelled total sea level (tide plus surge) above LAT at acquisition time.

    Build it with :meth:`from_series`, which fixes the datum once for the whole
    series. ``lon`` and ``lat`` in :meth:`water_level` are not used: the series
    is already the one cell nearest the site.
    """

    series: SeaLevelSeries
    lat_on_model_datum_m: float
    datum_method: str
    uncertainty_m: float = DEFAULT_UNCERTAINTY_M
    datum_details: dict[str, Any] = field(default_factory=dict)
    name: str = "cmems"

    @classmethod
    def from_series(
        cls,
        series: SeaLevelSeries,
        *,
        datum_method: str = "harmonic_lat",
        uncertainty_m: float = DEFAULT_UNCERTAINTY_M,
    ) -> CMEMSTide:
        if datum_method == "harmonic_lat":
            constants = fit_harmonics(
                series.times,
                series.levels_m,
                constituents=DEFAULT_CONSTITUENTS,
                datum=f"{series.dataset_id} model datum",
                source=series.station,
            )
            lat = lowest_astronomical_tide(constants)
            details: dict[str, Any] = {
                "definition": "lowest astronomical tide over a nodal cycle, harmonic fit",
                "fit": constants.to_dict(),
            }
        elif datum_method == "record_minimum":
            lat = float(series.levels_m.min())
            details = {
                "definition": "lowest hourly level in the record (proxy; includes surge)",
                "record": [str(series.times[0]), str(series.times[-1])],
            }
        else:
            raise ValueError(
                f"Unknown datum_method {datum_method!r}; use 'harmonic_lat' or 'record_minimum'."
            )
        return cls(series, lat, datum_method, uncertainty_m, details)

    def water_level(self, lon: float, lat: float, when: dt.datetime) -> TideCorrection:
        level = self.series.level_at(when) - self.lat_on_model_datum_m
        return TideCorrection(
            offset_m=level,
            derivation="model",
            source=f"{self.series.dataset_id} {self.series.variable}, {self.series.station}",
            when=when,
            uncertainty_m=self.uncertainty_m,
            metadata={
                "vertical_datum": "LAT",
                "datum_method": self.datum_method,
                "lat_on_model_datum_m": round(self.lat_on_model_datum_m, 4),
                "surge_included": True,
            },
        )
