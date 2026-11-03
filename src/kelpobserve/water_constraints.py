"""Separate absorption/backscatter amplitudes with explicit optical constraints.

This is a conditional water-state construction, without field optical validation.
The winter spectral shapes remain assumptions. No use is made of habitat labels or
image reflectance, and no impossible KD target is silently replaced by pure water.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from kelpobserve.attenuation import bands
from kelpobserve.coastal.optics import water


class WaterConstraintError(ValueError):
    """The inputs do not admit a water state under the declared spectral shape."""


def separate_iops(kd490, bbp443, solar_zenith_deg, wavelengths_nm, e5s, *, cdm443=None, cdm_error_pct=0):
    """Solve non-water absorption after fixing particle backscatter with BBP443.

    The complete interpolation support, including the 443/490 nm constraints, must
    be physical. Unsupported wavelengths and negative reference components raise.
    A CDM conflict is returned explicitly; such a state is not qualified for use.
    Provider percentage error defines a sensitivity bound, not a confidence level.
    """
    values = np.array([kd490, bbp443, solar_zenith_deg], float)
    wave = np.atleast_1d(np.asarray(wavelengths_nm, float))
    if not np.isfinite(values).all() or kd490 <= 0 or bbp443 < 0 or not 0 <= solar_zenith_deg < 90:
        raise WaterConstraintError("Invalid KD490, BBP443 or solar geometry")
    if wave.ndim != 1 or not len(wave) or not np.isfinite(wave).all():
        raise WaterConstraintError("Wavelengths must be a finite nonempty vector")
    ref = bands.E5S_IOP_NM
    if wave.min() < ref.min() or wave.max() > ref.max():
        raise WaterConstraintError("Wavelength extrapolation is not supported")
    an = np.asarray(e5s["iops"]["median"]["a"]) - water.a_water(ref)
    bp = np.asarray(e5s["iops"]["median"]["bb"]) - water.bb_water(ref)
    invalid = (~np.isfinite(an)) | (~np.isfinite(bp)) | (an < -1e-12) | (bp < 0)
    requested = np.r_[wave, 443.0, 490.0]
    required = set()
    for value in requested:
        right = int(np.searchsorted(ref, value))
        if right < len(ref) and ref[right] == value:
            required.add(right)
        else:
            required.update([right - 1, right])
    if any(invalid[i] for i in required):
        raise WaterConstraintError("Requested interpolation uses an unphysical winter reference component")
    b443 = float(np.interp(443, ref, bp))
    if b443 <= 0:
        raise WaterConstraintError("Reference particle backscatter at 443 nm must be positive")
    bscale = bbp443 / b443
    bb490 = water.bb_water(np.array([490.0])) + bscale * np.interp(490, ref, bp)

    def gap(ascale):
        a490 = water.a_water(np.array([490.0])) + ascale * np.interp(490, ref, an)
        return float(bands.kd(a490, bb490, solar_zenith_deg)[0]) - kd490

    if gap(0) > 1e-12 or gap(10000) < 0:
        raise WaterConstraintError("KD490 is infeasible for the specified backscatter")
    ascale = 0.0 if abs(gap(0)) <= 1e-12 else brentq(gap, 0, 10000, xtol=1e-12)
    a = water.a_water(wave) + ascale * np.interp(wave, ref, an)
    bb = water.bb_water(wave) + bscale * np.interp(wave, ref, bp)
    an443 = ascale * float(np.interp(443, ref, an))
    if cdm443 is not None:
        if not np.isfinite(cdm443) or cdm443 < 0 or not np.isfinite(cdm_error_pct) or not 0 <= cdm_error_pct <= 100:
            raise WaterConstraintError("Invalid CDM443 or error percentage")
        cdm_lower = cdm443 * (1 - cdm_error_pct / 100)
        conflict = bool(an443 < cdm_lower)
    else:
        cdm_lower, conflict = None, None
    return {
        "wavelengths_nm": wave.tolist(),
        "a": a.tolist(),
        "bb": bb.tolist(),
        "kd": bands.kd(a, bb, solar_zenith_deg).tolist(),
        "kd490_target": float(kd490),
        "bbp443_target": float(bbp443),
        "a_nonwater443": an443,
        "cdm443_lower": cdm_lower,
        "cdm443_nominal": cdm443,
        "cdm_nominal_conflict": bool(an443 < cdm443) if cdm443 is not None else None,
        "cdm_conflict": conflict,
        "input_constraints_consistent": conflict is False,
        "cdm_check_available": cdm443 is not None,
        "absorption_scale": ascale,
        "backscatter_scale": bscale,
        "invalid_reference_wavelengths_nm": ref[invalid].tolist(),
        "spectral_shape": "winter_E5s_assumed",
        "field_validated": False,
        "status": "input_conflict" if conflict else "conditional_unvalidated",
    }
