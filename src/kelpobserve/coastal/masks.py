"""Masks — land, cloud, sun-glint, foam, cliff-adjacency buffer, tidal exposure.

Every function returns a boolean array where True means *valid pixel* (keep).
Combine with ``&`` to produce the composite mask that gates every downstream
step. Downstream code must never touch pixels where the composite is False.

Ported from ``kelp_observe/tools/lee_demo/masks.py`` in Phase 1.3. The
prototype's per-function keyword thresholds are now defaults sourced from
:class:`~kelpobserve.coastal.config.MaskConfig`, so one config object
configures the whole chain for a new AOI or sensor.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage

from kelpobserve.coastal.config import MaskConfig

# --------------------------------------------------------------------------
# Individual masks — each returns True where the pixel is valid
# --------------------------------------------------------------------------


def land_mask(nir_reflectance: np.ndarray, threshold: float = 0.10) -> np.ndarray:
    """Water pixels only. NIR is strongly absorbed by clean water.

    The default threshold catches most sun-lit land and optically-shallow
    bright sand; for turbid or whitecap-heavy pixels the glint and foam masks
    pick up the residual. Deliberately generous — tighter thresholds cut into
    optically-shallow water where the seafloor is bright.
    """
    return np.asarray(nir_reflectance < threshold)


def cloud_mask(
    blue_reflectance: np.ndarray,
    nir_reflectance: np.ndarray,
    threshold: float = 0.20,
) -> np.ndarray:
    """Cloudy pixels have high reflectance in *both* blue and NIR.

    Water is only bright in blue; land is only bright in NIR; clouds are
    bright in both. Not a Sen2Cor-grade cloud detector — a fast, opinionated
    threshold for a single-scene retrieval.
    """
    return np.asarray(~((blue_reflectance > threshold) & (nir_reflectance > threshold)))


def glint_mask(
    nir_reflectance: np.ndarray,
    swir_reflectance: np.ndarray | None = None,
    nir_threshold: float = 0.03,
    swir_threshold: float | None = None,
) -> np.ndarray:
    """Reject specular sun-glint.

    Water pixels with elevated NIR are almost always glint-contaminated.
    PNeo has no SWIR — pass ``swir_reflectance=None`` and rely on NIR alone.
    Sentinel-2 has B11/B12; combining them tightens the mask.

    ``swir_threshold`` defaults to ``nir_threshold`` when not supplied.
    """
    valid = nir_reflectance < nir_threshold
    if swir_reflectance is not None:
        thr_swir = nir_threshold if swir_threshold is None else swir_threshold
        valid = valid & (swir_reflectance < thr_swir)
    return np.asarray(valid)


def foam_mask(
    rgb_reflectance: np.ndarray,
    threshold: float = 0.15,
    flatness_max: float = 1.4,
) -> np.ndarray:
    """Whitecaps and surf-line foam — bright, near-flat spectrum.

    ``rgb_reflectance`` has shape (3, H, W) with the visible bands in any
    order. Rejects pixels where all three visible bands exceed ``threshold``
    and their max/min ratio is below ``flatness_max``.
    """
    high = np.all(rgb_reflectance > threshold, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        span = rgb_reflectance.max(axis=0) / np.clip(rgb_reflectance.min(axis=0), 1e-6, None)
    return np.asarray(~(high & (span < flatness_max)))


def adjacency_buffer(land: np.ndarray, buffer_px: int) -> np.ndarray:
    """Drop water pixels within ``buffer_px`` pixels of land.

    ``land`` is the land-*rejection* mask (True over water), matching what
    :func:`land_mask` returns.

    ``buffer_px <= 0`` disables the buffer. The guard is not cosmetic:
    ``scipy.ndimage.binary_dilation`` reads ``iterations=0`` as "dilate until
    nothing changes", so without it a configured zero buffer would flood the
    whole scene with land and silently mask everything.
    """
    if buffer_px <= 0:
        return np.asarray(land)
    dilated_land = ndimage.binary_dilation(~land, iterations=buffer_px, border_value=0)
    return np.asarray(~dilated_land)


# --------------------------------------------------------------------------
# Composite
# --------------------------------------------------------------------------


def composite_mask(
    *,
    blue: np.ndarray,
    red: np.ndarray,
    nir: np.ndarray,
    swir: np.ndarray | None = None,
    config: MaskConfig | None = None,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Combine individual masks. Returns ``(composite, components)``.

    All input arrays share the same shape and are *surface* reflectance
    (ACOLITE ``rhos_*`` or the Hedley-deglinted equivalent).

    The ``components`` dict keeps each individual mask inspectable — when a
    scene fails its rhos-sanity gate, look at the components before touching
    thresholds.
    """
    cfg = config or MaskConfig()
    components: dict[str, np.ndarray] = {
        "land": land_mask(nir, threshold=cfg.land_nir_max),
        "cloud": cloud_mask(blue, nir, threshold=cfg.cloud_reflectance_max),
        "glint": glint_mask(
            nir,
            swir,
            nir_threshold=cfg.glint_nir_max,
            swir_threshold=cfg.glint_swir_max,
        ),
        "foam": foam_mask(
            np.stack([blue, red, nir]),
            threshold=cfg.foam_reflectance_min,
            flatness_max=cfg.foam_flatness_max,
        ),
    }
    components["adjacency"] = adjacency_buffer(
        components["land"], buffer_px=cfg.adjacency_buffer_px
    )

    composite = np.logical_and.reduce(list(components.values()))
    return composite, components


# --------------------------------------------------------------------------
# Deep-water sampling for the QAA fit
# --------------------------------------------------------------------------


def swir_clean_water(
    rhos: np.ndarray,
    wavelengths_nm: np.ndarray,
    domain: np.ndarray,
    config: MaskConfig | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Drop cloud, glint and foam from a deep-water domain using SWIR.

    Pure water absorbs about 670 m⁻¹ at 1610 nm, so water-leaving SWIR
    reflectance is zero and whatever is measured there is contamination. The
    reference domain is chosen on geometry alone, which keeps it independent of
    the result but means it arrives unscreened. Screening on a band that carries
    no benthic signal cannot tune the answer.

    Applies to the domain *before* any capped sampling. Returns the domain
    unchanged, with ``screened: False``, when the sensor has no SWIR band
    within 20 nm of ``config.swir_bands_nm`` (Pléiades Neo).

    Ported from ``kelp_observe/tools/lee_demo/reef_calibration.screen_deepwater``.
    """
    cfg = config or MaskConfig()
    wavelengths = np.asarray(wavelengths_nm, dtype=float)
    bands = [
        int(np.argmin(np.abs(wavelengths - nm)))
        for nm in cfg.swir_bands_nm
        if np.min(np.abs(wavelengths - nm)) <= 20.0
    ]
    if not bands:
        return domain, {"screened": False, "reason": "no SWIR band"}
    clean = domain.copy()
    for index in bands:
        band = rhos[index]
        clean &= np.isfinite(band) & (np.abs(band) <= cfg.swir_max_rhos)
    before = int(domain.sum())
    return clean, {
        "screened": True,
        "swir_bands_nm": [float(wavelengths[i]) for i in bands],
        "swir_max_rhos": cfg.swir_max_rhos,
        "n_before": before,
        "n_after": int(clean.sum()),
        "fraction_kept": round(float(clean.sum() / before), 4) if before else 0.0,
    }


def deep_water_pixels(
    valid: np.ndarray,
    depth: np.ndarray,
    kd_490: np.ndarray | None = None,
    config: MaskConfig | None = None,
) -> np.ndarray:
    """Boolean mask of optically-deep pixels for QAA fitting.

    "Optically deep" means the bottom contribution to ``r_rs`` is below the
    sensor noise floor. Two independent gates, both of which a pixel must
    pass:

    - depth > ``config.deepwater_min_depth_m`` (from EMODnet or Stumpf)
    - z * Kd(490) > ``config.deepwater_min_optical_depth`` (rigorous, but
      requires a first Kd guess — see
      :func:`kelpobserve.coastal.optics.qaa.kd_map`)

    Survivors are subsampled to ``config.deepwater_max_samples`` to keep the
    fit fast; the draw is seeded so a rerun reproduces the same sample.
    """
    cfg = config or MaskConfig()
    deep = valid & (depth > cfg.deepwater_min_depth_m)
    if kd_490 is not None:
        deep = deep & (depth * kd_490 > cfg.deepwater_min_optical_depth)

    if deep.sum() > cfg.deepwater_max_samples:
        rng = np.random.default_rng(cfg.rng_seed)
        idx = np.flatnonzero(deep.ravel())
        keep = rng.choice(idx, size=cfg.deepwater_max_samples, replace=False)
        subsampled = np.zeros_like(deep.ravel())
        subsampled[keep] = True
        deep = subsampled.reshape(deep.shape)
    return np.asarray(deep)


# --------------------------------------------------------------------------
# Tidal exposure: intertidal pixels at acquisition time
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class TidalExposure:
    """Pixels by their water cover at acquisition; the three masks do not overlap.

    ``submerged`` is the valid mask. ``dry`` pixels are above the waterline
    whatever the tide error, and can serve as above-water reference spectra.
    ``awash`` pixels are within the uncertainty of the waterline or under less
    than the minimum water: they are neither, and are excluded from both.
    """

    submerged: np.ndarray
    awash: np.ndarray
    dry: np.ndarray
    stats: dict[str, object] = field(default_factory=dict)


def tidal_exposure(
    depth_datum_m: np.ndarray,
    water_level_m: float,
    uncertainty_m: float = 0.0,
    config: MaskConfig | None = None,
    swir_rhos: np.ndarray | None = None,
) -> TidalExposure:
    """Split pixels into submerged, awash and dry at one water level.

    ``depth_datum_m`` is positive-down depth below the chart datum, negative on
    intertidal ground above it; ``water_level_m`` is the water surface above the
    same datum (a :class:`~kelpobserve.coastal.bathymetry.tide.TideCorrection`
    offset), and ``uncertainty_m`` its standard error, including any datum
    doubt. Water depth at acquisition is their sum.

    With ``config.tidal_dry_swir_min`` set and ``swir_rhos`` given (SWIR surface
    reflectance on the same grid), a tide-dry pixel stays dry only where the
    image agrees; flats that drain late, pools and ground lower than the depth
    grid says go to awash instead.
    """
    cfg = config or MaskConfig()
    if not np.isfinite(water_level_m):
        raise ValueError("The water level must be finite.")
    if not np.isfinite(uncertainty_m) or uncertainty_m < 0:
        raise ValueError("The water-level uncertainty must be finite and nonnegative.")
    margin = cfg.tidal_uncertainty_sigma * uncertainty_m
    depth = np.asarray(depth_datum_m, dtype=float)
    finite = np.isfinite(depth)
    water = np.where(finite, depth + water_level_m, np.nan)
    with np.errstate(invalid="ignore"):
        dry = finite & (water < -margin)
        submerged = finite & (water >= cfg.tidal_min_water_m + margin)
    unconfirmed = np.zeros_like(dry)
    if cfg.tidal_dry_swir_min is not None:
        if swir_rhos is None:
            raise ValueError("tidal_dry_swir_min is set but no SWIR reflectance was given.")
        swir = np.asarray(swir_rhos, dtype=float)
        if swir.shape != depth.shape:
            raise ValueError("swir_rhos must share the depth grid.")
        with np.errstate(invalid="ignore"):
            unconfirmed = dry & ~(swir >= cfg.tidal_dry_swir_min)
        dry = dry & ~unconfirmed
    awash = finite & ~dry & ~submerged
    n = int(finite.sum())
    return TidalExposure(
        submerged=submerged,
        awash=awash,
        dry=dry,
        stats={
            "water_level_m": round(float(water_level_m), 3),
            "uncertainty_m": round(float(uncertainty_m), 3),
            "margin_m": round(float(margin), 3),
            "min_water_m": cfg.tidal_min_water_m,
            "n_pixels": n,
            "n_submerged": int(submerged.sum()),
            "n_awash": int(awash.sum()),
            "n_dry": int(dry.sum()),
            "dry_swir_min": cfg.tidal_dry_swir_min,
            "n_dry_not_confirmed": int(unconfirmed.sum()),
            "fraction_submerged": round(float(submerged.sum() / n), 4) if n else 0.0,
        },
    )
