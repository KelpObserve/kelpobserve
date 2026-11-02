"""Terrain classification from bathymetry geometry alone — slope and local relief.

Promoted from ``sand_control.terrain_classes()`` in the prototype (Phase 1.4),
where it was a private helper inside a one-off diagnostic script. It earns
public API status because it is the only substrate discriminator in the
library that does not read the imagery.

That independence is the whole point. Identifying bare sediment by its
reflectance and then measuring that sediment's reflectance is circular, and
the circularity is invisible in the output. Terrain is orthogonal to optics,
so a terrain-selected class provides a genuine positive control for the Lee
inversion: a substrate whose albedo is predicted in advance, at a location
chosen without looking at the answer.

The classes are ``flat`` and ``rugose``, and they do not partition the scene.
Pixels between the two thresholds belong to neither, deliberately — an
ambiguous pixel contributes nothing to a control and should not be forced
into a class to make the coverage number look better.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import maximum_filter, minimum_filter

from kelpobserve.coastal.config import BathymetryConfig


def terrain_classes(
    depth_m: np.ndarray,
    pixel_size_m: float = 10.0,
    config: BathymetryConfig | None = None,
) -> dict[str, np.ndarray]:
    """Split the seabed into flat and rugose classes using geometry only.

    Parameters
    ----------
    depth_m
        (H, W) bathymetry in metres. Any source — lidar, multibeam, a blended
        Stumpf surface. Non-finite pixels propagate to non-finite slope and
        relief, and are excluded from both classes.
    pixel_size_m
        Ground sample distance, used to scale the gradient into a true slope.
        Passing the wrong value silently rescales every slope in the scene.
    config
        Supplies the four thresholds and the relief window size.

    Returns
    -------
    dict
        ``slope_deg`` and ``relief_m`` continuous rasters, plus the boolean
        ``flat`` and ``rugose`` masks. The continuous rasters are returned
        alongside the masks so a threshold can be re-tuned against the actual
        distribution rather than guessed.
    """
    cfg = config or BathymetryConfig()

    finite = np.where(np.isfinite(depth_m), depth_m, np.nan)
    gy, gx = np.gradient(finite, pixel_size_m, pixel_size_m)
    slope_deg = np.degrees(np.arctan(np.hypot(gx, gy)))
    relief_m = maximum_filter(finite, cfg.relief_window_px) - minimum_filter(
        finite, cfg.relief_window_px
    )

    # The first term is not redundant. np.gradient's central difference never
    # reads the centre pixel, and scipy's min/max filters do not propagate
    # NaN reliably, so a lone data hole in otherwise flat bathymetry comes
    # back with slope 0 and relief 0 — and would be classified as flat. A
    # pixel with no depth measurement must not become a positive control.
    usable = np.isfinite(depth_m) & np.isfinite(slope_deg) & np.isfinite(relief_m)
    return {
        "slope_deg": slope_deg,
        "relief_m": relief_m,
        "flat": usable
        & (slope_deg < cfg.flat_max_slope_deg)
        & (relief_m < cfg.flat_max_relief_m),
        "rugose": usable
        & (slope_deg > cfg.rugose_min_slope_deg)
        & (relief_m > cfg.rugose_min_relief_m),
    }


def topographic_position(surface: np.ndarray, window_px: int) -> np.ndarray:
    """Each pixel minus the NaN-aware mean of a ``window_px`` square around it.

    On a depth raster (positive down) a positive value is deeper than its surroundings;
    on an elevation raster, higher. NaN pixels stay NaN.
    """
    from scipy.ndimage import uniform_filter

    if window_px < 1:
        raise ValueError("window_px must be at least 1.")
    values = np.asarray(surface, dtype=np.float64)
    valid = np.isfinite(values)
    total = uniform_filter(np.where(valid, values, 0.0), size=window_px, mode="constant")
    count = uniform_filter(valid.astype(np.float64), size=window_px, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count
    return np.where(valid, values - mean, np.nan)


def distance_to(mask: np.ndarray, pixel_size_m: float) -> np.ndarray:
    """Straight-line distance in metres from every pixel to the nearest ``True`` pixel."""
    from scipy.ndimage import distance_transform_edt

    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return np.full(mask.shape, np.inf)
    return np.asarray(distance_transform_edt(~mask)) * pixel_size_m


def water_path_distance(
    sources: np.ndarray, passable: np.ndarray, pixel_size_m: float
) -> np.ndarray:
    """Shortest distance in metres from ``sources`` through ``passable`` pixels (8-connected).

    Uses scikit-image's minimum-cost-path search with diagonal steps of √2. Pixels that
    cannot be reached, or are not passable, are ``inf``.
    """
    from skimage.graph import MCP_Geometric

    passable = np.asarray(passable, dtype=bool)
    sources = np.asarray(sources, dtype=bool) & passable
    if not sources.any():
        raise ValueError("No source pixel is passable.")
    costs = np.where(passable, 1.0, np.inf)
    cumulative, _ = MCP_Geometric(costs).find_costs(np.argwhere(sources))
    return np.where(passable, cumulative * pixel_size_m, np.inf)
