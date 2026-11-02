"""Measured reflectance spectra: libraries on disk, and their values in sensor bands.

The detectability defaults (:data:`kelpobserve.coastal.detectability.SUBSTRATE_ALBEDO`)
are three band-averaged groups from one Baltic library. Checking or replacing them
needs other libraries, which come in other layouts, measured in other quantities,
and cover other wavelength ranges. This module reads them into one structure and
integrates them over a sensor's band responses. It ships no spectra: every library
is read from a file the caller supplies, so a library whose licence does not allow
redistribution can be used without being copied into this package.

Layouts
-------
* :func:`read_pangaea_wide`: a PANGAEA tab export with one spectrum per row and
  one column per wavelength (``Rrs_350 [1/sr]`` ...), as in PANGAEA 971518.
* :func:`read_pangaea_long`: a PANGAEA tab export with one row per wavelength and
  columns naming the spectrum, as in PANGAEA 824861.
* :func:`read_columns_csv`: one column per spectrum and a wavelength column.
* :func:`read_rows_csv`: one spectrum per row, wavelengths as column names.
* :func:`read_asd`: ASD binary files (version 7 onwards) with their stored white
  reference, as written by the FieldSpec and HandHeld instruments.

Quantities
----------
A library holds either ``"rrs"``, remote-sensing reflectance in sr⁻¹, or
``"reflectance"``, a dimensionless reflectance factor against a white panel. The
only conversion offered is :meth:`SpectralLibrary.as_reflectance`, ``R = pi * Rrs``,
which assumes a Lambertian target; the detectability caveats apply to it.

Band integration
----------------
:func:`resample` weights each spectrum by a band's relative response over the part
of the band the spectrum covers. When the covered part holds less than
``min_coverage`` of the response, the band value is NaN rather than a mean over
the edge that happens to be measured. Coverage is judged on the true response, not
the nominal width: a library starting at 428 nm still covers 99.8 % of Sentinel-2A
band 1 (443 nm), whose response below 428 nm is small.

:func:`read_rsr` reads ACOLITE's response files (ESA's Sentinel-2 functions,
Airbus's Pléiades Neo functions) from an ACOLITE install at run time;
:func:`boxcar_response` is the fallback when no response is available.
"""

from __future__ import annotations

import csv
import re
import struct
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

QUANTITIES = ("rrs", "reflectance")

#: Byte offsets in the ASD file header (ASD File Format, version 8).
_ASD_DATA_TYPE = 186
_ASD_START_WAVELENGTH = 191
_ASD_FORMAT = 199
_ASD_CHANNELS = 204
_ASD_HEADER_BYTES = 484
_ASD_DTYPES = {0: "<f4", 1: "<i4", 2: "<f8"}
#: ASD data types: raw counts need the reference; a stored reflectance does not.
_ASD_RAW, _ASD_REFLECTANCE = 0, 1


def _integrate(y: np.ndarray, x: np.ndarray, axis: int = -1) -> Any:
    """Trapezoid rule; ``np.trapezoid`` from NumPy 2, ``np.trapz`` before it."""
    rule = getattr(np, "trapezoid", None) or np.trapz
    return rule(y, x, axis=axis)


@dataclass(frozen=True)
class SpectralLibrary:
    """Spectra on one wavelength grid, one row each, with per-spectrum metadata."""

    wavelengths_nm: np.ndarray
    values: np.ndarray
    metadata: tuple[dict[str, str], ...]
    quantity: str
    source: str = ""

    def __post_init__(self) -> None:
        wl = np.asarray(self.wavelengths_nm, dtype=np.float64)
        values = np.atleast_2d(np.asarray(self.values, dtype=np.float64))
        if self.quantity not in QUANTITIES:
            raise ValueError(f"quantity must be one of {QUANTITIES}; got {self.quantity!r}.")
        if wl.ndim != 1 or np.any(np.diff(wl) <= 0):
            raise ValueError("wavelengths_nm must be one-dimensional and strictly ascending.")
        if values.shape[1] != wl.size:
            raise ValueError(f"values have {values.shape[1]} columns for {wl.size} wavelengths.")
        if len(self.metadata) != values.shape[0]:
            raise ValueError(
                f"{len(self.metadata)} metadata entries for {values.shape[0]} spectra."
            )
        object.__setattr__(self, "wavelengths_nm", wl)
        object.__setattr__(self, "values", values)

    def __len__(self) -> int:
        return int(self.values.shape[0])

    def select(self, keep: Callable[[dict[str, str]], bool]) -> SpectralLibrary:
        """The spectra whose metadata ``keep`` accepts."""
        rows = [i for i, meta in enumerate(self.metadata) if keep(meta)]
        return SpectralLibrary(
            self.wavelengths_nm,
            self.values[rows] if rows else np.empty((0, self.wavelengths_nm.size)),
            tuple(self.metadata[i] for i in rows),
            self.quantity,
            self.source,
        )

    def where(self, **equals: str) -> SpectralLibrary:
        """The spectra whose metadata match every ``key=value``, ignoring case and spaces."""

        def norm(text: str) -> str:
            return " ".join(str(text).split()).casefold()

        return self.select(
            lambda meta: all(norm(meta.get(k, "")) == norm(v) for k, v in equals.items())
        )

    def as_reflectance(self) -> SpectralLibrary:
        """Reflectance factor; ``pi * Rrs`` for a Lambertian target."""
        if self.quantity == "reflectance":
            return self
        return SpectralLibrary(
            self.wavelengths_nm, np.pi * self.values, self.metadata, "reflectance", self.source
        )


@dataclass(frozen=True)
class BandResponse:
    """Relative spectral response of one sensor band."""

    name: str
    wavelengths_nm: np.ndarray
    response: np.ndarray
    notes: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        wl = np.asarray(self.wavelengths_nm, dtype=np.float64)
        response = np.clip(np.asarray(self.response, dtype=np.float64), 0.0, None)
        if wl.shape != response.shape or wl.size < 2:
            raise ValueError(f"{self.name}: wavelengths and response must match, 2+ samples.")
        order = np.argsort(wl)
        object.__setattr__(self, "wavelengths_nm", wl[order])
        object.__setattr__(self, "response", response[order])

    @property
    def centre_nm(self) -> float:
        """Response-weighted centre wavelength."""
        weight = _integrate(self.response, self.wavelengths_nm)
        return float(_integrate(self.response * self.wavelengths_nm, self.wavelengths_nm) / weight)


def boxcar_response(name: str, centre_nm: float, fwhm_nm: float) -> BandResponse:
    """A flat response over the full width at half maximum, for sensors without one."""
    half = fwhm_nm / 2.0
    wl = np.array(
        [centre_nm - half - 1e-6, centre_nm - half, centre_nm + half, centre_nm + half + 1e-6]
    )
    return BandResponse(name, wl, np.array([0.0, 1.0, 1.0, 0.0]), {"shape": "boxcar"})


def resample(
    library: SpectralLibrary, bands: Sequence[BandResponse], min_coverage: float = 0.95
) -> np.ndarray:
    """Response-weighted band values, ``(n_spectra, n_bands)``.

    Each band is integrated on the library's own wavelength grid. NaN where the
    library covers less than ``min_coverage`` of the band's response, and for a
    spectrum with a gap inside the band.
    """
    out = np.full((len(library), len(bands)), np.nan)
    wl = library.wavelengths_nm
    for j, band in enumerate(bands):
        # Coverage on a fine grid, so neither the band's nor the library's
        # sampling decides where the library's range cuts the response.
        fine = np.linspace(band.wavelengths_nm[0], band.wavelengths_nm[-1], 2001)
        response = np.interp(fine, band.wavelengths_nm, band.response)
        total = _integrate(response, fine)
        covered = _integrate(response * ((fine >= wl[0]) & (fine <= wl[-1])), fine)
        if total <= 0 or covered / total < min_coverage:
            continue
        weight = np.interp(wl, band.wavelengths_nm, band.response, left=0.0, right=0.0)
        inside = weight > 0
        if inside.sum() < 2:
            continue
        rows = library.values[:, inside]
        w = weight[inside]
        integral = _integrate(rows * w, wl[inside], axis=1) / _integrate(w, wl[inside])
        gaps = ~np.all(np.isfinite(rows), axis=1)
        integral[gaps] = np.nan
        out[:, j] = integral
    return out


# ---------------------------------------------------------------------------
# Band responses
# ---------------------------------------------------------------------------

_BAND_HEADER = re.compile(r"band\s+(\w+)", re.IGNORECASE)


def read_rsr(path: str | Path) -> dict[str, BandResponse]:
    """Band responses from an ACOLITE ``data/RSR`` text file, keyed ``B1``, ``B2``...

    A comment line naming ``Band <n>`` starts a band; two-column numeric lines
    follow. Wavelengths are converted from micrometres when the file uses them
    (the Airbus files do; ESA's are in nanometres).
    """
    source = Path(path)
    bands: dict[str, tuple[list[float], list[float]]] = {}
    current: str | None = None
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
        text = line.strip()
        if not text:
            continue
        if text.startswith(("#", ";")):
            match = _BAND_HEADER.search(text)
            if match:
                current = f"B{match.group(1)}"
                bands[current] = ([], [])
            continue
        parts = text.split()
        if current is None or len(parts) < 2:
            continue
        try:
            wavelength, value = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        bands[current][0].append(wavelength)
        bands[current][1].append(value)
    if not bands:
        raise ValueError(f"{source.name}: no 'Band <n>' sections found.")
    out = {}
    for name, (wl, response) in bands.items():
        wl_nm = np.asarray(wl) * (1000.0 if max(wl) < 10.0 else 1.0)
        out[name] = BandResponse(name, wl_nm, np.asarray(response), {"source": source.name})
    return out


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------


def _pangaea_table(path: Path) -> tuple[list[str], list[list[str]]]:
    """Header and rows of a PANGAEA tab export, below its ``*/`` metadata block."""
    text = path.read_text(encoding="utf-8", errors="replace")
    _, marker, body = text.partition("*/")
    lines = (body if marker else text).strip("\n").splitlines()
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    return header, rows


def _number(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        return float("nan")


def read_pangaea_wide(
    path: str | Path, *, prefix: str = "Rrs_", quantity: str = "rrs", source: str = ""
) -> SpectralLibrary:
    """One spectrum per row; wavelength columns named ``<prefix><nm> ...``."""
    header, rows = _pangaea_table(Path(path))
    spectral = [
        (i, float(re.findall(r"[\d.]+", name[len(prefix) :])[0]))
        for i, name in enumerate(header)
        if name.startswith(prefix)
    ]
    if not spectral:
        raise ValueError(f"{Path(path).name}: no columns starting with {prefix!r}.")
    columns = [i for i, _ in spectral]
    meta_columns = [i for i in range(len(header)) if i not in set(columns)]
    values, metadata = [], []
    for raw_row in rows:
        row = raw_row + [""] * (len(header) - len(raw_row))
        values.append([_number(row[i]) for i in columns])
        metadata.append({header[i]: row[i].strip() for i in meta_columns})
    return SpectralLibrary(
        np.array([w for _, w in spectral]), np.array(values), tuple(metadata), quantity, source
    )


def read_pangaea_long(
    path: str | Path,
    *,
    wavelength_column: str,
    value_column: str,
    spectrum_columns: Sequence[str],
    quantity: str = "reflectance",
    source: str = "",
) -> SpectralLibrary:
    """One row per wavelength; ``spectrum_columns`` together name the spectrum.

    Spectra are aligned on the union of their wavelengths, with NaN where a
    spectrum has no value.
    """
    header, rows = _pangaea_table(Path(path))
    index = {name: i for i, name in enumerate(header)}
    for name in (wavelength_column, value_column, *spectrum_columns):
        if name not in index:
            raise KeyError(f"{Path(path).name}: no column {name!r}; columns are {header}.")
    groups: dict[tuple[str, ...], dict[float, float]] = {}
    first_row: dict[tuple[str, ...], list[str]] = {}
    for row in rows:
        key = tuple(row[index[c]].strip() for c in spectrum_columns)
        groups.setdefault(key, {})[_number(row[index[wavelength_column]])] = _number(
            row[index[value_column]]
        )
        first_row.setdefault(key, row)
    grid = np.array(sorted({w for spectrum in groups.values() for w in spectrum if np.isfinite(w)}))
    values = np.array([[groups[key].get(w, np.nan) for w in grid] for key in groups])
    skip = {index[wavelength_column], index[value_column]}
    metadata = tuple(
        {header[i]: first_row[key][i].strip() for i in range(len(header)) if i not in skip}
        for key in groups
    )
    return SpectralLibrary(grid, values, metadata, quantity, source)


def read_columns_csv(
    path: str | Path,
    *,
    wavelength_column: str = "Wavelength",
    quantity: str = "reflectance",
    source: str = "",
    extra_metadata: dict[str, str] | None = None,
    delimiter: str = ",",
    encoding: str = "utf-8-sig",
) -> SpectralLibrary:
    """One column per spectrum; the column name becomes metadata ``name``."""
    with Path(path).open(newline="", encoding=encoding) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        rows = [row for row in reader if row and any(cell.strip() for cell in row)]
    if wavelength_column not in header:
        raise KeyError(f"{Path(path).name}: no column {wavelength_column!r}.")
    w_index = header.index(wavelength_column)
    names = [(i, name.strip()) for i, name in enumerate(header) if i != w_index and name.strip()]
    wl = np.array([_number(row[w_index]) for row in rows])
    values = np.array(
        [[_number(row[i]) if i < len(row) else np.nan for row in rows] for i, _ in names]
    )
    order = np.argsort(wl)
    metadata = tuple(
        {"name": name, "file": Path(path).name, **(extra_metadata or {})} for _, name in names
    )
    return SpectralLibrary(wl[order], values[:, order], metadata, quantity, source)


def read_rows_csv(
    path: str | Path,
    *,
    quantity: str = "reflectance",
    source: str = "",
    delimiter: str = ",",
    encoding: str = "utf-8-sig",
) -> SpectralLibrary:
    """One spectrum per row: numeric column names are wavelengths, the rest metadata."""
    with Path(path).open(newline="", encoding=encoding) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = [name.strip() for name in next(reader)]
        rows = [row for row in reader if row and any(cell.strip() for cell in row)]
    spectral = [(i, _number(name)) for i, name in enumerate(header)]
    spectral = [(i, w) for i, w in spectral if np.isfinite(w)]
    if not spectral:
        raise ValueError(f"{Path(path).name}: no numeric (wavelength) column names.")
    columns = {i for i, _ in spectral}
    order = np.argsort([w for _, w in spectral])
    values = np.array(
        [[_number(row[i]) if i < len(row) else np.nan for i, _ in spectral] for row in rows]
    )
    metadata = tuple(
        {
            header[i]: (row[i].strip() if i < len(row) else "")
            for i in range(len(header))
            if i not in columns
        }
        for row in rows
    )
    wl = np.array([w for _, w in spectral])
    return SpectralLibrary(wl[order], values[:, order], metadata, quantity, source)


def read_asd(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Wavelengths (nm) and reflectance factor from one ASD binary file.

    Raw-count files are divided by the white reference stored after the
    spectrum; files saved as reflectance are returned as stored. Other data
    types (radiance, irradiance) are refused, since they are not reflectance.
    """
    raw = Path(path).read_bytes()
    if raw[:2] != b"as" or raw[2:3] < b"7":
        raise ValueError(f"{Path(path).name}: not an ASD file of version 7 or later.")
    data_type = raw[_ASD_DATA_TYPE]
    start, step = struct.unpack("<ff", raw[_ASD_START_WAVELENGTH : _ASD_START_WAVELENGTH + 8])
    fmt = raw[_ASD_FORMAT]
    channels = struct.unpack("<H", raw[_ASD_CHANNELS : _ASD_CHANNELS + 2])[0]
    if fmt not in _ASD_DTYPES:
        raise ValueError(f"{Path(path).name}: unsupported ASD data format {fmt}.")
    dtype = np.dtype(_ASD_DTYPES[fmt])
    end = _ASD_HEADER_BYTES + channels * dtype.itemsize
    spectrum = np.frombuffer(raw[_ASD_HEADER_BYTES:end], dtype=dtype).astype(np.float64)
    wavelengths = start + step * np.arange(channels)
    # Reference header: flag (2), reference time (8), spectrum time (8), description.
    length = struct.unpack("<H", raw[end + 18 : end + 20])[0]
    ref_start = end + 20 + length
    reference = np.frombuffer(
        raw[ref_start : ref_start + channels * dtype.itemsize], dtype=dtype
    ).astype(np.float64)
    if data_type not in (_ASD_RAW, _ASD_REFLECTANCE):
        raise ValueError(f"{Path(path).name}: ASD data type {data_type} is not reflectance.")
    if reference.size != channels or not np.any(reference > 0):
        raise ValueError(f"{Path(path).name}: no white reference stored after the spectrum.")
    with np.errstate(divide="ignore", invalid="ignore"):
        reflectance = np.where(reference > 0, spectrum / reference, np.nan)
    return wavelengths, reflectance


def read_asd_library(
    paths: Iterable[str | Path],
    metadata: Iterable[dict[str, str]] | None = None,
    *,
    source: str = "",
) -> SpectralLibrary:
    """ASD files on one grid; ``metadata`` in the same order, else the file name only."""
    files = [Path(p) for p in paths]
    meta = list(metadata) if metadata is not None else [{} for _ in files]
    if len(meta) != len(files):
        raise ValueError(f"{len(meta)} metadata entries for {len(files)} files.")
    grids, rows = [], []
    for path in files:
        wl, r = read_asd(path)
        grids.append(wl)
        rows.append(r)
    if any(g.shape != grids[0].shape or not np.allclose(g, grids[0]) for g in grids):
        raise ValueError("ASD files are on different wavelength grids.")
    tagged = tuple({"file": p.name, **m} for p, m in zip(files, meta, strict=True))
    return SpectralLibrary(grids[0], np.vstack(rows), tagged, "reflectance", source)


def describe(library: SpectralLibrary, key: str) -> dict[str, Any]:
    """Count of spectra per value of one metadata field."""
    counts: dict[str, int] = {}
    for meta in library.metadata:
        value = meta.get(key, "")
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))
