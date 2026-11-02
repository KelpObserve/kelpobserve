"""Harmonic tide prediction, fitting and Lowest Astronomical Tide.

The astronomical tide at a place is a sum of constituents, each a cosine at a
fixed astronomical speed::

    η(t) = Z0 + Σ f_k(t) · H_k · cos(V_k(t) + u_k(t) − g_k)

``H`` and ``g`` are the site's amplitude and phase lag. ``V`` is the
constituent's astronomical argument; ``f`` and ``u`` are the nodal corrections
for the 18.61-year lunar node cycle, which changes M2 by about ±4 % and K1 and
O1 by 12–19 %.

Why a harmonic predictor rather than the ocean model directly
-------------------------------------------------------------
Regional models such as Copernicus IBI run at about 2.5 km, which cannot resolve
a lagoon like Ria Formosa, and their sea level includes the meteorological surge
of the day. A harmonic set gives the astronomical tide at any time, including
dates the model does not cover, and defines Lowest Astronomical Tide the way the
charts do: the lowest level the predicted tide reaches over a full nodal cycle.
The fitted residual, which includes the surge, becomes the stated uncertainty.

Phase convention
----------------
Phases are Greenwich phase lags against arguments computed from the mean
longitudes of Meeus (1998), with mean lunar time ``τ = 15°·t + h − s`` and ``t``
the hours of UT since midnight. This is the convention of the Instituto
Hidrográfico's published constants: fitted to its 2026 Faro-Olhão predictions,
M2, S2, K1 and O1 come out within 7° of its published phases. Adding 180° to
``τ`` leaves the semidiurnal phases unchanged and turns the diurnal ones by
180°, so only a diurnal constituent can detect that error.

Nodal corrections are the simplified expressions of Pugh (1987), Table 4.3.
"""

from __future__ import annotations

import datetime as dt
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

#: Doodson numbers (τ, s, h, p, N', p1) and the extra phase in degrees.
DOODSON: dict[str, tuple[tuple[int, int, int, int, int, int], float]] = {
    "M2": ((2, 0, 0, 0, 0, 0), 0.0),
    "S2": ((2, 2, -2, 0, 0, 0), 0.0),
    "N2": ((2, -1, 0, 1, 0, 0), 0.0),
    "K2": ((2, 2, 0, 0, 0, 0), 0.0),
    "K1": ((1, 1, 0, 0, 0, 0), 90.0),
    "O1": ((1, -1, 0, 0, 0, 0), -90.0),
    "P1": ((1, 1, -2, 0, 0, 0), -90.0),
    "Q1": ((1, -2, 0, 1, 0, 0), -90.0),
    "M4": ((4, 0, 0, 0, 0, 0), 0.0),
    "MS4": ((4, 2, -2, 0, 0, 0), 0.0),
    "MN4": ((4, -1, 0, 1, 0, 0), 0.0),
}

#: The set fitted by default: the eight major constituents plus the shallow-water
#: overtides that matter in estuaries and lagoons.
DEFAULT_CONSTITUENTS: tuple[str, ...] = (
    "M2",
    "S2",
    "N2",
    "K2",
    "K1",
    "O1",
    "P1",
    "Q1",
    "M4",
    "MS4",
)

#: Rates of the astronomical variables in degrees per hour (τ, s, h, p, N', p1).
_RATES = np.array(
    [
        15.0 + 0.0410686387 - 0.5490165,  # τ = T + h − s
        0.5490165,
        0.0410686387,
        0.0046418,
        0.0022064,
        0.0000020,
    ]
)

#: Hours in the 18.61-year lunar nodal cycle, over which LAT is defined.
NODAL_CYCLE_HOURS = 18.61 * 365.25 * 24


def speed(constituent: str) -> float:
    """Angular speed in degrees per hour."""
    numbers, _ = DOODSON[constituent]
    return float(np.dot(numbers, _RATES))


def _astronomy(times: np.ndarray) -> tuple[np.ndarray, ...]:
    """Mean longitudes τ, s, h, p, N', p1 in degrees at UTC ``times`` (datetime64)."""
    hours = (times - np.datetime64("2000-01-01T12:00:00")) / np.timedelta64(1, "h")
    t_cent = hours / (24 * 36525.0)
    s = 218.3164477 + 481267.88123421 * t_cent
    h = 280.46646 + 36000.76983 * t_cent
    p = 83.3532465 + 4069.0137287 * t_cent
    node = 125.04452 - 1934.136261 * t_cent
    p1 = 282.93735 + 1.71946 * t_cent
    ut_hours = (hours + 12.0) % 24.0
    tau = 15.0 * ut_hours + h - s
    return tau, s, h, p, -node, p1


def _nodal(constituent: str, node_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nodal amplitude factor ``f`` and phase ``u`` (degrees)."""
    n = np.radians(node_deg)
    one = np.ones_like(n)
    m2 = (
        1.0004 - 0.0373 * np.cos(n) + 0.0002 * np.cos(2 * n),
        -2.14 * np.sin(n),
    )
    table = {
        "M2": m2,
        "N2": m2,
        "S2": (one, 0.0 * n),
        "P1": (one, 0.0 * n),
        "K1": (
            1.0060 + 0.1150 * np.cos(n) - 0.0088 * np.cos(2 * n) + 0.0006 * np.cos(3 * n),
            -8.86 * np.sin(n) + 0.68 * np.sin(2 * n) - 0.07 * np.sin(3 * n),
        ),
        "O1": (
            1.0089 + 0.1871 * np.cos(n) - 0.0147 * np.cos(2 * n) + 0.0014 * np.cos(3 * n),
            10.80 * np.sin(n) - 1.34 * np.sin(2 * n) + 0.19 * np.sin(3 * n),
        ),
        "K2": (
            1.0241 + 0.2863 * np.cos(n) + 0.0083 * np.cos(2 * n) - 0.0015 * np.cos(3 * n),
            -17.74 * np.sin(n) + 0.68 * np.sin(2 * n) - 0.04 * np.sin(3 * n),
        ),
    }
    table["Q1"] = table["O1"]
    table["M4"] = (m2[0] ** 2, 2 * m2[1])
    table["MS4"] = m2
    table["MN4"] = (m2[0] ** 2, 2 * m2[1])
    return table[constituent]


def _arguments(constituent: str, times: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``f`` and ``V + u`` in degrees for one constituent at ``times``."""
    numbers, extra = DOODSON[constituent]
    astro = _astronomy(times)
    argument = sum(k * a for k, a in zip(numbers, astro, strict=True)) + extra
    node = -astro[4]
    f, u = _nodal(constituent, node)
    return f, argument + u


def _as_utc64(times: Iterable[dt.datetime] | np.ndarray) -> np.ndarray:
    if isinstance(times, np.ndarray) and np.issubdtype(times.dtype, np.datetime64):
        return times.astype("datetime64[s]")
    values = []
    for when in times:
        if when.tzinfo is None:
            raise ValueError("Tide times must be timezone-aware; got a naive datetime.")
        values.append(when.astimezone(dt.UTC).replace(tzinfo=None))
    return np.array(values, dtype="datetime64[s]")


@dataclass(frozen=True)
class HarmonicConstants:
    """A fitted or published harmonic set for one place.

    ``mean_level_m`` is Z0 on the constants' datum. ``constituents`` maps each
    name to (amplitude in metres, phase lag in degrees).
    """

    constituents: Mapping[str, tuple[float, float]]
    mean_level_m: float
    datum: str
    source: str
    #: Standard deviation of the fit residual, which includes surge and model
    #: error. ``None`` for a published set with no stated error.
    residual_sd_m: float | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        unknown = set(self.constituents) - set(DOODSON)
        if unknown:
            raise ValueError(f"Unknown constituents {sorted(unknown)}; known: {sorted(DOODSON)}.")
        if not self.source.strip():
            raise ValueError("Harmonic constants require their source.")

    def to_dict(self) -> dict[str, object]:
        return {
            "constituents": {
                k: {"amplitude_m": round(a, 5), "phase_deg": round(g % 360.0, 3)}
                for k, (a, g) in self.constituents.items()
            },
            "mean_level_m": round(self.mean_level_m, 4),
            "datum": self.datum,
            "source": self.source,
            "residual_sd_m": None if self.residual_sd_m is None else round(self.residual_sd_m, 4),
            **({"metadata": dict(self.metadata)} if self.metadata else {}),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> HarmonicConstants:
        raw = data["constituents"]
        assert isinstance(raw, Mapping)
        return cls(
            constituents={
                k: (float(v["amplitude_m"]), float(v["phase_deg"])) for k, v in raw.items()
            },
            mean_level_m=float(data["mean_level_m"]),  # type: ignore[arg-type]
            datum=str(data["datum"]),
            source=str(data["source"]),
            residual_sd_m=(
                None if data.get("residual_sd_m") is None else float(data["residual_sd_m"])  # type: ignore[arg-type]
            ),
            metadata=dict(data.get("metadata") or {}),  # type: ignore[call-overload]
        )


def predict(constants: HarmonicConstants, times: Iterable[dt.datetime] | np.ndarray) -> np.ndarray:
    """Astronomical water level on the constants' datum, in metres."""
    stamps = _as_utc64(times)
    level = np.full(stamps.shape, constants.mean_level_m, dtype=np.float64)
    for name, (amplitude, phase) in constants.constituents.items():
        f, argument = _arguments(name, stamps)
        level += f * amplitude * np.cos(np.radians(argument - phase))
    return level


def rayleigh_conflicts(constituents: Sequence[str], duration_hours: float) -> list[str]:
    """Constituent pairs a record of this length cannot separate.

    Two constituents are resolved when the record spans at least one full cycle
    of their speed difference (the Rayleigh criterion).
    """
    conflicts = []
    names = list(constituents)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            gap = abs(speed(a) - speed(b))
            if gap > 0 and duration_hours < 360.0 / gap:
                conflicts.append(f"{a}/{b} needs {360.0 / gap / 24:.0f} days")
    return conflicts


def fit_harmonics(
    times: Iterable[dt.datetime] | np.ndarray,
    levels: np.ndarray,
    *,
    constituents: Sequence[str] = DEFAULT_CONSTITUENTS,
    datum: str,
    source: str,
) -> HarmonicConstants:
    """Least-squares harmonic analysis of a water-level record.

    Refuses a record too short to separate the requested constituents rather
    than returning aliased amplitudes. NaN levels are dropped.
    """
    stamps = _as_utc64(times)
    y = np.asarray(levels, dtype=np.float64)
    if stamps.shape != y.shape:
        raise ValueError(f"times {stamps.shape} and levels {y.shape} differ in shape.")
    keep = np.isfinite(y)
    stamps, y = stamps[keep], y[keep]
    if y.size < 2 * len(constituents) + 1:
        raise ValueError(f"Only {y.size} finite levels for {len(constituents)} constituents.")
    duration = float((stamps.max() - stamps.min()) / np.timedelta64(1, "h"))
    conflicts = rayleigh_conflicts(constituents, duration)
    if conflicts:
        raise ValueError(
            f"A {duration / 24:.0f}-day record cannot separate: {'; '.join(conflicts)}. "
            "Use a longer record or fewer constituents."
        )
    columns = [np.ones_like(y)]
    for name in constituents:
        f, argument = _arguments(name, stamps)
        columns += [f * np.cos(np.radians(argument)), f * np.sin(np.radians(argument))]
    design = np.column_stack(columns)
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    fitted = {}
    for i, name in enumerate(constituents):
        a, b = coefficients[1 + 2 * i], coefficients[2 + 2 * i]
        fitted[name] = (float(np.hypot(a, b)), float(np.degrees(np.arctan2(b, a)) % 360.0))
    return HarmonicConstants(
        constituents=fitted,
        mean_level_m=float(coefficients[0]),
        datum=datum,
        source=source,
        residual_sd_m=float(np.std(residual, ddof=design.shape[1])),
        metadata={
            "n_observations": int(y.size),
            "record_start": str(stamps.min()),
            "record_end": str(stamps.max()),
            "record_days": round(duration / 24, 1),
        },
    )


def lowest_astronomical_tide(
    constants: HarmonicConstants,
    start: dt.datetime = dt.datetime(2020, 1, 1, tzinfo=dt.UTC),
    step_minutes: int = 30,
) -> float:
    """The lowest predicted level over one lunar nodal cycle, on the constants' datum.

    This is the definition of LAT. The half-hour step resolves the minimum of a
    semidiurnal tide to within a centimetre or so.
    """
    n = int(NODAL_CYCLE_HOURS * 60 / step_minutes)
    begin = np.datetime64(start.astimezone(dt.UTC).replace(tzinfo=None), "s")
    stamps = begin + np.arange(n) * np.timedelta64(step_minutes * 60, "s")
    lowest = np.inf
    for chunk in np.array_split(stamps, max(1, n // 50_000)):
        lowest = min(lowest, float(predict(constants, chunk).min()))
    return lowest


def water_level_series(
    constants: HarmonicConstants,
    start: dt.datetime,
    end: dt.datetime,
    step_minutes: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """Predicted levels on a regular grid of times from ``start`` up to ``end``.

    Returns ``(times, levels)``: UTC ``datetime64[s]`` and metres on the constants' datum.
    """
    if end <= start:
        raise ValueError("end must be after start.")
    begin = np.datetime64(start.astimezone(dt.UTC).replace(tzinfo=None), "s")
    stop = np.datetime64(end.astimezone(dt.UTC).replace(tzinfo=None), "s")
    stamps: np.ndarray = np.arange(begin, stop, np.timedelta64(step_minutes * 60, "s"))
    levels = np.concatenate(
        [predict(constants, c) for c in np.array_split(stamps, max(1, stamps.size // 50_000))]
    )
    return stamps, levels


def emersion_fraction(
    depth_datum_m: np.ndarray,
    constants: HarmonicConstants,
    start: dt.datetime,
    end: dt.datetime,
    *,
    step_minutes: int = 10,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Share of time each pixel is out of the water under the predicted tide.

    ``depth_datum_m`` is positive down from the constants' datum, as in
    :func:`kelpobserve.coastal.masks.tidal_exposure`: intertidal ground is negative. A pixel
    is emersed while the level is below its ground, ``h(t) < -depth``. ``weights`` (one per
    step of :func:`water_level_series`) turn the time share into a weighted share, for
    example by daylight. NaN depths give NaN.

    The tide is the astronomical prediction at the constants' station: no surge, and no
    lag or damping across a lagoon.
    """
    _, levels = water_level_series(constants, start, end, step_minutes)
    if weights is None:
        weights = np.ones_like(levels)
    weights = np.asarray(weights, dtype=np.float64)
    if weights.shape != levels.shape:
        raise ValueError(f"weights need one value per time step ({levels.size}).")
    order = np.argsort(levels)
    sorted_levels = levels[order]
    cumulative = np.concatenate([[0.0], np.cumsum(weights[order])])
    total = cumulative[-1]
    if total <= 0:
        raise ValueError("weights sum to zero.")
    depth = np.asarray(depth_datum_m, dtype=np.float64)
    ground = -depth
    out = np.full(depth.shape, np.nan)
    finite = np.isfinite(ground)
    index = np.searchsorted(sorted_levels, ground[finite], side="left")
    out[finite] = cumulative[index] / total
    return out


def referenced_to_lat(constants: HarmonicConstants) -> HarmonicConstants:
    """The same constants with Z0 re-expressed as height above LAT."""
    lat = lowest_astronomical_tide(constants)
    return HarmonicConstants(
        constituents=dict(constants.constituents),
        mean_level_m=constants.mean_level_m - lat,
        datum="LAT",
        source=constants.source,
        residual_sd_m=constants.residual_sd_m,
        metadata={
            **dict(constants.metadata),
            "lat_on_source_datum_m": round(lat, 4),
            "source_datum": constants.datum,
        },
    )
