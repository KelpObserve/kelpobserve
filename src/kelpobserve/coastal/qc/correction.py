"""Subtract the additive atmospheric-correction residual (experiment E2).

:mod:`.offset` measures the residual; this module decides how much to take
off each band and does it. The hypothesis under test is that the unsubtracted
residual explains three failures at once: QAA switching to its red reference
with an implausible bbp, ``ratio_near_unity``, and red attenuation below the
pure-water floor. Whether any mode fixes them is for the registered
experiment to say, so the default stays ``none``.

Modes, each spectrally explicit so the report shows what was removed:

- ``none``: nothing subtracted.
- ``nir_black``: the NIR black-pixel offset (:func:`.offset.deepwater_additive_offset`),
  flat across bands.
- ``bottom_blind_red``: the median red (~665 nm) reflectance over water deeper
  than ``red_min_depth_m``, flat across bands. At 10 m the red two-way path
  through pure water alone is e^-8, so the seabed is gone; the water's own red
  reflectance is not zero, so this mode overcorrects slightly in turbid water.
- ``spectral``: the NIR offset extrapolated to shorter wavelengths along the
  power law through the NIR and SWIR offsets, as for aerosol path reflectance.
  Falls back to flat, and says so, when SWIR gives no usable slope.

The offset is measured on the deep-water sample chosen before correction and
applied to every pixel; a scalar cannot remove spatially varying glint, which
:mod:`.ac_uncertainty` carries as noise instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

from kelpobserve.coastal.config import OFFSET_MODES, QCConfig
from kelpobserve.coastal.qc.offset import deepwater_additive_offset
from kelpobserve.coastal.scene import Scene

#: Water deeper than this is bottom-blind in the red.
RED_MIN_DEPTH_M = 10.0
RED_NM = 665.0
#: Bands past this carry no water-leaving signal and are not corrected: they
#: are the evidence the correction is judged on.
MAX_CORRECTED_NM = 1000.0
#: Ångström-type exponents outside this range are not aerosol; fall back to flat.
MAX_SPECTRAL_EXPONENT = 4.0


@dataclass(frozen=True)
class OffsetCorrection:
    mode: str
    #: rhos subtracted per band, by wavelength in nm.
    by_band: dict[float, float] = field(default_factory=dict)
    #: Evidence the offset was measured on; empty for ``none``.
    source: dict[str, Any] = field(default_factory=dict)
    flags: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "subtracted_rhos": {f"{wl:g}": round(v, 6) for wl, v in sorted(self.by_band.items())},
            "source": self.source,
            "flags": list(self.flags),
        }


def _median(values: np.ndarray) -> float | None:
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else None


def estimate_offset(
    scene: Scene,
    deep_mask: np.ndarray,
    mode: str,
    *,
    depth_m: np.ndarray | None = None,
    qc_config: QCConfig | None = None,
    red_min_depth_m: float = RED_MIN_DEPTH_M,
) -> OffsetCorrection:
    """How much rhos to subtract from each band under ``mode``."""
    if mode not in OFFSET_MODES:
        raise ValueError(f"Unknown offset mode {mode!r}; expected one of {OFFSET_MODES}.")
    wavelengths = [float(w) for w in scene.wavelengths_nm]
    corrected = [w for w in wavelengths if w < MAX_CORRECTED_NM]
    if mode == "none":
        return OffsetCorrection("none", {w: 0.0 for w in corrected})
    if scene.rhos is None:
        raise ValueError("Offset correction needs surface reflectance (rhos).")
    bands = {w: scene.rhos[i] for i, w in enumerate(wavelengths)}
    nir = deepwater_additive_offset(bands, deep_mask, config=qc_config)
    flags = list(nir.qa_flags)

    if mode == "nir_black":
        return OffsetCorrection(
            mode, {w: nir.value for w in corrected}, {"nir": nir.to_dict()}, tuple(flags)
        )

    if mode == "bottom_blind_red":
        if depth_m is None:
            raise ValueError("bottom_blind_red needs a depth raster.")
        red = min(wavelengths, key=lambda w: abs(w - RED_NM))
        if abs(red - RED_NM) > 20.0:
            raise ValueError(f"No red band near {RED_NM:g} nm.")
        with np.errstate(invalid="ignore"):
            blind = deep_mask & np.isfinite(depth_m) & (depth_m > red_min_depth_m)
        value = _median(bands[red][blind])
        if value is None:
            return OffsetCorrection(
                mode, {w: 0.0 for w in corrected}, {"n_pixels": 0}, ("no_bottom_blind_water",)
            )
        if value < 0:
            flags.append("clamped_negative")
            value = 0.0
        return OffsetCorrection(
            mode,
            {w: value for w in corrected},
            {
                "red_nm": red,
                "red_min_depth_m": red_min_depth_m,
                "n_pixels": int(blind.sum()),
                "median_red_rhos": round(value, 6),
                "nir": nir.to_dict(),
            },
            tuple(flags),
        )

    # spectral
    swir = sorted(w for w in wavelengths if w >= (qc_config or QCConfig()).offset_swir_min_nm)
    swir_value = _median(bands[swir[0]][deep_mask]) if swir else None
    nir_nm = float(np.median(list(nir.by_band))) if nir.by_band else None
    exponent = None
    if nir_nm and nir.value > 0 and swir_value and swir_value > 0:
        exponent = float(np.log(nir.value / swir_value) / np.log(swir[0] / nir_nm))
    if exponent is None or not 0.0 <= exponent <= MAX_SPECTRAL_EXPONENT:
        flags.append("spectral_fallback_flat")
        by_band = {w: nir.value for w in corrected}
    else:
        assert nir_nm is not None
        by_band = {w: nir.value * (w / nir_nm) ** -exponent for w in corrected}
    return OffsetCorrection(
        mode,
        by_band,
        {
            "nir": nir.to_dict(),
            "nir_reference_nm": nir_nm,
            "swir_nm": swir[0] if swir else None,
            "swir_median_rhos": None if swir_value is None else round(swir_value, 6),
            "exponent": None if exponent is None else round(exponent, 4),
        },
        tuple(flags),
    )


def apply_offset(scene: Scene, correction: OffsetCorrection) -> Scene:
    """A copy of ``scene`` with the offset subtracted from rhos and Rrs (offset / π)."""
    if scene.rhos is None:
        raise ValueError("Offset correction needs surface reflectance (rhos).")
    offsets = np.array(
        [correction.by_band.get(float(w), 0.0) for w in scene.wavelengths_nm], dtype=np.float32
    )[:, None, None]
    return replace(
        scene,
        rhos=(scene.rhos - offsets).astype(scene.rhos.dtype),
        rrs_above=(scene.rrs_above - offsets / np.float32(np.pi)).astype(scene.rrs_above.dtype),
        metadata={**scene.metadata, "offset_correction": correction.to_dict()},
    )
