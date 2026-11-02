"""Inherent optical properties and Kd at a sensor's bands, from one Kd(490) and Sesimbra's IOP shape.

Only the magnitude comes from Kd(490). The spectral shape is the site's own: the median
non-water absorption and particulate backscatter of E5s, from the forced-560 QAA fits to the
winter Sentinel-2 scenes of 24 December 2025 and 22 February 2026 (the fits that passed E2c's
plausibility check). One factor s scales both, chosen so that the Lee et al. (2005) Kd(490) of
the scaled IOPs equals the given value; a, bb and Kd then follow at every band, with pure water
evaluated exactly at each wavelength. Pure-water optics and the Kd model come from
kelpobserve.coastal.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

E5S_IOP_NM = np.array([442.0, 492.0, 559.0, 665.0, 704.0])


def site_shape(e5s: dict) -> tuple[np.ndarray, np.ndarray]:
    """Median non-water absorption and particulate backscatter at E5S_IOP_NM."""
    from kelpobserve.coastal.optics import water

    median = e5s["iops"]["median"]
    a_nw = np.clip(np.array(median["a"]) - water.a_water(E5S_IOP_NM), 0.0, None)
    bb_p = np.clip(np.array(median["bb"]) - water.bb_water(E5S_IOP_NM), 0.0, None)
    return a_nw, bb_p


def iops_at(wavelengths_nm, scale: float, shape: tuple[np.ndarray, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    from kelpobserve.coastal.optics import water

    wavelengths = np.asarray(wavelengths_nm, dtype=float)
    a_nw = np.interp(wavelengths, E5S_IOP_NM, shape[0])
    bb_p = np.interp(wavelengths, E5S_IOP_NM, shape[1])
    return water.a_water(wavelengths) + scale * a_nw, water.bb_water(wavelengths) + scale * bb_p


def kd(a: np.ndarray, bb: np.ndarray, solar_zenith_deg: float) -> np.ndarray:
    from kelpobserve.coastal.optics.qaa import _kd_lee_2005

    return _kd_lee_2005(np.asarray(a), np.asarray(bb), solar_zenith_deg)


def scale_for(kd490: float, shape, solar_zenith_deg: float) -> float:
    """The factor on the site shape whose Kd(490) equals kd490 (0 when kd490 is at or below pure water)."""

    def gap(scale: float) -> float:
        a, bb = iops_at([490.0], scale, shape)
        return float(kd(a, bb, solar_zenith_deg)[0]) - kd490

    if gap(0.0) >= 0:
        return 0.0
    upper = 1.0
    while gap(upper) < 0:
        upper *= 2.0
        if upper > 1e4:
            raise ValueError(f"Kd(490) {kd490} is beyond any scaling of the site shape.")
    return float(brentq(gap, 0.0, upper, xtol=1e-12))


def at_bands(kd490: float, wavelengths_nm, solar_zenith_deg: float, e5s: dict) -> dict:
    shape = site_shape(e5s)
    scale = scale_for(kd490, shape, solar_zenith_deg)
    a, bb = iops_at(wavelengths_nm, scale, shape)
    return {
        "kd490": kd490,
        "scale_on_site_shape": round(scale, 6),
        "wavelengths_nm": [float(w) for w in wavelengths_nm],
        "a": [round(float(v), 6) for v in a],
        "bb": [round(float(v), 6) for v in bb],
        "kd": [round(float(v), 6) for v in kd(a, bb, solar_zenith_deg)],
    }
