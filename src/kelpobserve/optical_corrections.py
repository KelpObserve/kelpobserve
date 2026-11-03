"""Development NIR deglinting with explicit baseline and calibration support.

The output removes an estimated differential surface signal. It does not establish
absolute water-leaving reflectance or correct land adjacency. Negative values are
preserved for quality assessment.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class GlintCalibrationError(ValueError):
    """The supplied reference does not support the requested correction."""


@dataclass(frozen=True)
class GlintCalibration:
    slopes: np.ndarray
    slope_p25: np.ndarray
    slope_p75: np.ndarray
    nir_baseline: float
    nir_p01: float
    nir_p10: float
    nir_p99: float
    block_count: int
    reference_count: int
    baseline_quantile: float = 0.05
    absolute_radiometry_validated: bool = False


def fit_block_glint(nir, visible, blocks, *, min_cells=10, min_span=0.002, min_blocks=6):
    """Fit slopes from deep reference blocks; caller supplies the training split.

    Each qualifying block has equal weight. Block intercepts absorb broad spatial
    background variation. The NIR baseline is fixed at the reference fifth percentile.
    The caller must keep habitat pixels and held-out reference blocks out of this fit.
    """
    nir, visible, blocks = np.asarray(nir, float), np.asarray(visible, float), np.asarray(blocks)
    if nir.ndim != 1 or visible.ndim != 2 or visible.shape[0] != len(nir) or blocks.shape != nir.shape:
        raise GlintCalibrationError("Expected N reference values, N x B spectra and N block ids")
    if min_cells < 3 or min_blocks < 2 or not np.isfinite(min_span) or min_span <= 0:
        raise GlintCalibrationError("Invalid calibration support requirements")
    finite = np.isfinite(nir) & np.isfinite(visible).all(axis=1)
    if not finite.all():
        raise GlintCalibrationError("Reference selection must explicitly exclude nonfinite observations")
    slopes = []
    for block in np.unique(blocks):
        mask = blocks == block
        x, y = nir[mask], visible[mask]
        if len(x) < min_cells or np.quantile(x, 0.95) - np.quantile(x, 0.05) < min_span:
            continue
        dx, dy = x - x.mean(), y - y.mean(axis=0)
        slopes.append((dx[:, None] * dy).sum(axis=0) / (dx * dx).sum())
    if len(slopes) < min_blocks:
        raise GlintCalibrationError(f"Only {len(slopes)} informative blocks; need {min_blocks}")
    p25, median, p75 = np.quantile(slopes, [0.25, 0.5, 0.75], axis=0)
    if not np.isfinite(median).all() or (median <= 0).any() or (median > 3).any():
        raise GlintCalibrationError("Nonpositive or excessive glint slope")
    p01, p05, p10, p99 = np.quantile(nir, [0.01, 0.05, 0.10, 0.99])
    if p05 < 0 or p99 <= p05:
        raise GlintCalibrationError("Reference NIR baseline or dynamic range is invalid")
    return GlintCalibration(median, p25, p75, float(p05), float(p01), float(p10), float(p99), len(slopes), len(nir))


def apply_glint(nir, visible, slopes, *, nir_baseline, support=None):
    """Subtract slope * max(NIR-baseline, 0), retaining diagnostics outside support.

    A baseline must be supplied explicitly. The support flag covers only the NIR
    calibration range; it is not an absolute-radiometry or habitat-validity flag.
    """
    nir, visible, slopes = np.asarray(nir, float), np.asarray(visible, float), np.asarray(slopes, float)
    if slopes.ndim != 1 or not len(slopes) or visible.shape != nir.shape + (len(slopes),):
        raise ValueError("Visible spectra must have shape NIR.shape + (bands,)")
    if not np.isfinite(nir_baseline) or nir_baseline < 0 or not np.isfinite(slopes).all() or (slopes < 0).any():
        raise ValueError("An explicit nonnegative baseline and finite nonnegative slopes are required")
    finite = np.isfinite(nir) & np.isfinite(visible).all(axis=-1)
    amplitude = np.maximum(nir - nir_baseline, 0)
    corrected = visible - amplitude[..., None] * slopes
    corrected = np.where(finite[..., None], corrected, np.nan)
    in_range = finite.copy()
    if support is not None:
        low, high = support
        if not (np.isfinite(low) and np.isfinite(high) and low < high):
            raise ValueError("Invalid NIR calibration range")
        in_range &= (nir >= low) & (nir <= high)
    return {
        "corrected": corrected,
        "removed": amplitude[..., None] * slopes,
        "nir_after": np.where(finite, nir - amplitude, np.nan),
        "within_nir_support": in_range,
        "negative": corrected < 0,
        "absolute_radiometry_validated": False,
    }
