"""Water optics primitives — Rrs conversions, pure-water IOPs, spectral shapes.

Everything in this module operates on numpy arrays. Wavelengths are in
nanometres; reflectances are dimensionless; absorption/backscatter are in m^-1.

Ported verbatim from ``kelp_observe/tools/lee_demo/water.py`` — pure functions
with no I/O, no site coupling.
"""

from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------
# Above-water Rrs ↔ subsurface r_rs (Lee & Carder 1998)
# --------------------------------------------------------------------------
#
# Above-water remote-sensing reflectance:
#     Rrs = pi * rhow / rhow_measurement_geometry
# ACOLITE emits Rrs directly. What Lee et al.'s inversion works with is the
# subsurface value r_rs, related by:
#
#     r_rs = Rrs / (Q * (1 - rho_sky * f_reflectance))
#
# and, empirically (Lee & Carder 1998 eq 3):
#
#     r_rs = Rrs / (0.52 + 1.7 * Rrs)
#
# The inverse:
#
#     Rrs = 0.52 * r_rs / (1 - 1.7 * r_rs)

def rrs_above_to_subsurface(rrs_above: np.ndarray) -> np.ndarray:
    """Lee & Carder 1998, eq 3. Rrs (1/sr) → r_rs (1/sr)."""
    return np.asarray(rrs_above / (0.52 + 1.7 * rrs_above))


def rrs_subsurface_to_above(rrs_sub: np.ndarray) -> np.ndarray:
    """Inverse of `rrs_above_to_subsurface`. r_rs (1/sr) → Rrs (1/sr)."""
    return 0.52 * rrs_sub / np.clip(1.0 - 1.7 * rrs_sub, 1e-6, None)


# --------------------------------------------------------------------------
# Surface reflectance ρ_s → Rrs approximation
# --------------------------------------------------------------------------
#
# ACOLITE's `rhos_*` is the *water-leaving* surface reflectance,
# rhos ≈ pi * Rrs. This is an approximation; the exact conversion needs
# per-geometry BRDF handling. For a demo, pi is close enough — Lee's own
# inversion papers use the same shortcut.

def rhos_to_rrs(rhos: np.ndarray) -> np.ndarray:
    """rhos → Rrs (1/sr). Demo-grade Lambertian approximation."""
    return rhos / np.pi


# --------------------------------------------------------------------------
# Pure-water absorption a_w and backscatter b_bw
# --------------------------------------------------------------------------
#
# Tabulated values from Pope & Fry 1997 (absorption) and Morel 1974
# (backscatter). We interpolate to the sensor centre wavelength once at
# scene setup; the ~1 nm error from linear interpolation is well below the
# noise floor of a 6-band retrieval.

# Wavelength, a_w (m^-1) — Pope & Fry 1997, Table 3 (integrating-cavity
# measurements, 380-700 nm). Values 530-600 nm were previously off by up to
# +20%, which matters because 555/560 nm is the QAA green reference band:
# an inflated a_w there propagates through a(ref) -> b_b -> a(lambda) -> k_b
# and is amplified exponentially by exp(+k_b * H) in the Lee inversion.
# Beyond 700 nm: Röttgers et al., combined pure-water absorption for ESA
# Waterradiance (WOPP v3, R. Röttgers, HZG, November 2016; 20 °C, 0 PSU),
# 10 nm steps to 1000 nm and the Sentinel-2 / Landsat SWIR band centres beyond.
# Pure water absorbs strongly there, so these values serve as floors and QA
# thresholds; the seabed is not seen. Earlier versions extrapolated flat at
# 0.624 m^-1, about 4 times too low at 740 nm and 1,000 times at 1610 nm
# (measurement-definitions audit M9). The NIR temperature coefficient is
# about +0.01 m^-1 per °C at 740-900 nm.
_AW_TABLE = np.array([
    [400, 0.00663],
    [410, 0.00473],
    [420, 0.00454],
    [430, 0.00495],
    [440, 0.00635],
    [450, 0.00922],
    [460, 0.00979],
    [470, 0.01060],
    [480, 0.01270],
    [490, 0.01500],
    [500, 0.02040],
    [510, 0.03250],
    [520, 0.04090],
    [530, 0.04340],
    [540, 0.04740],
    [550, 0.05650],
    [560, 0.06190],
    [570, 0.06950],
    [580, 0.08960],
    [590, 0.13510],
    [600, 0.22240],
    [610, 0.26430],
    [620, 0.27550],
    [630, 0.29160],
    [640, 0.31080],
    [650, 0.34000],
    [660, 0.41000],
    [670, 0.43900],
    [680, 0.46500],
    [690, 0.51600],
    [700, 0.62400],
    [710, 0.85605],
    [720, 1.28344],
    [730, 2.03522],
    [740, 2.4773],
    [750, 2.6125],
    [760, 2.6114],
    [770, 2.47885],
    [780, 2.2706],
    [790, 2.05498],
    [800, 1.964],
    [810, 2.09367],
    [820, 2.4829],
    [830, 3.0905],
    [840, 3.6808],
    [850, 4.38418],
    [860, 4.94008],
    [870, 5.36475],
    [880, 5.83038],
    [890, 6.26645],
    [900, 6.7924],
    [910, 7.8707],
    [920, 11.338],
    [930, 19.38617],
    [940, 29.005],
    [950, 38.08414],
    [960, 44.13729],
    [970, 45.2874],
    [980, 43.851],
    [990, 41.45014],
    [1000, 37.699],
    [1100, 19.8863],
    [1200, 125.4875],
    [1300, 135.80773],
    [1374, 540.93833],
    [1400, 1714.4],
    [1500, 2087.93333],
    [1610, 696.26058],
    [1700, 577.98223],
    [1800, 904.71484],
    [2000, 7444.26399],
    [2190, 1941.4425],
    [2300, 2521.10128],
])
# Pope & Fry 1997 alone, flat at 0.624 m^-1 beyond 700 nm: the tabulation the v1 runs of
# 20 September 2026 used. The v1 replays put it in place of _AW_TABLE (kelpobserve._patch).
_AW_TABLE_TO_700_NM = _AW_TABLE[_AW_TABLE[:, 0] <= 700]

_BBW_550 = 0.00193 / 2.0  # Morel 1974; half of total scattering


def a_water(wavelength_nm: np.ndarray | float) -> np.ndarray | float:
    """Pure-water absorption coefficient a_w(λ), m^-1.

    Linear interpolation of Pope & Fry 1997 (400-700 nm) and Röttgers et al.
    (WOPP v3, 710-2300 nm, 20 °C). Flat outside 400-2300 nm. Beyond about
    700 nm the values are floors and thresholds: pure water hides the seabed.
    """
    return np.interp(wavelength_nm, _AW_TABLE[:, 0], _AW_TABLE[:, 1])


def bb_water(wavelength_nm: np.ndarray | float) -> np.ndarray | float:
    """Pure-water backscatter b_bw(λ), m^-1.

    b_bw(λ) = b_bw(550) × (550/λ)^4.32.
    Morel 1974 pure-water scattering, with the standard b/b_b = 0.5 fraction.
    """
    return _BBW_550 * (550.0 / np.asarray(wavelength_nm)) ** 4.32


# --------------------------------------------------------------------------
# Spectral shapes for particles/CDOM (used by QAA v6 and Lee inversion)
# --------------------------------------------------------------------------

def a_cdm_shape(
    wavelength_nm: np.ndarray | float,
    ref_wavelength_nm: float = 443.0,
    s_cdm: float = 0.015,
) -> np.ndarray | float:
    """Dimensionless spectral shape of CDOM+detritus absorption.

    a_cdm(λ) = a_cdm(λ0) * exp(-S_cdm * (λ - λ0))
    Reference wavelength 443 nm and S_cdm = 0.015 nm^-1 (Bricaud 1981 range).
    Returns the ratio a_cdm(λ) / a_cdm(λ0), so multiply by the retrieved
    a_cdm(λ0) to get absolute values.
    """
    result = np.exp(-s_cdm * (np.asarray(wavelength_nm) - ref_wavelength_nm))
    return np.asarray(result)


def bbp_shape(
    wavelength_nm: np.ndarray | float,
    ref_wavelength_nm: float = 555.0,
    y: float = 1.0,
) -> np.ndarray | float:
    """Dimensionless spectral shape of particle backscatter.

    b_bp(λ) = b_bp(λ0) * (λ0 / λ)^Y
    Y ≈ 1 for typical case-2 waters; 0 for a purely non-selective scatterer;
    2 for pure water. QAA v6 estimates Y from the r_rs ratio itself.
    """
    result = (ref_wavelength_nm / np.asarray(wavelength_nm)) ** y
    return np.asarray(result)
