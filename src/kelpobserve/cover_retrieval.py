"""Fixed-endmember Lee cover retrieval, with explicit compatibility diagnostics.

Inputs and tolerances are supplied by the frozen development protocol. Fractions
are spectral mixtures; their relationship to diver cover is assessed separately.
"""

from __future__ import annotations

import numpy as np

from kelpobserve.coastal.inversion import lee
from kelpobserve.coastal.optics import water


def curves(fractions, depth, a, bb, solar_zenith, kelp, background):
    """Above-water dimensionless rhos, shape (fraction, band), under Lee physics."""
    fractions = np.asarray(fractions, dtype=float)
    a, bb, kelp, background = (np.asarray(v, dtype=float) for v in (a, bb, kelp, background))
    if not (np.isfinite(depth) and depth > 0 and np.isfinite(solar_zenith)):
        raise ValueError("Positive finite depth and finite solar zenith are required")
    if a.shape != bb.shape or a.shape != kelp.shape or a.shape != background.shape:
        raise ValueError("Water and spectra must have identical band shapes")
    if not all(np.isfinite(v).all() for v in (a, bb, kelp, background, fractions)) or (a <= 0).any() or (bb < 0).any():
        raise ValueError("Invalid optical inputs")
    u = bb / (a + bb)
    inv_cos = 1 / np.cos(np.deg2rad(lee._subsurface_solar_zenith_deg(solar_zenith)))
    kc, kb = (inv_cos + lee._du_column(u)) * (a + bb), (inv_cos + lee._du_bottom(u)) * (a + bb)
    deep = 0.089 * u + 0.125 * u * u
    rho = background + fractions[..., None] * (kelp - background)
    sub = deep * (1 - np.exp(-kc * depth)) + rho / np.pi * np.exp(-kb * depth)
    return np.pi * water.rrs_subsurface_to_above(sub)


def retrieve(observed, model_curve, fractions, tolerances):
    """Closest fit and feasible cover ranges at fixed rhos tolerances; no probability claim."""
    observed = np.asarray(observed, dtype=float)
    model_curve, fractions = np.asarray(model_curve), np.asarray(fractions)
    if model_curve.shape != (len(fractions), len(observed)):
        raise ValueError("Curve dimensions do not match fractions and observations")
    if not np.isfinite(observed).all() or not np.isfinite(model_curve).all():
        raise ValueError("Nonfinite observation or forward model")
    error = model_curve - observed
    best = int(np.argmin(np.mean(error**2, axis=1)))
    compatible = {}
    for tolerance in tolerances:
        if tolerance <= 0:
            raise ValueError("Tolerance must be positive")
        mask = (np.abs(error) <= tolerance).all(axis=1)
        compatible[str(tolerance)] = {
            "any": bool(mask.any()),
            "best_pct": float(100 * fractions[np.flatnonzero(mask)[np.argmin(np.mean(error[mask] ** 2, axis=1))]])
            if mask.any()
            else None,
            "n_grid_values": int(mask.sum()),
            "min_pct": float(100 * fractions[mask].min()) if mask.any() else None,
            "max_pct": float(100 * fractions[mask].max()) if mask.any() else None,
        }
    return {
        "cover_pct": float(100 * fractions[best]),
        "boundary": bool(best in (0, len(fractions) - 1)),
        "rmse_rhos": float(np.sqrt(np.mean(error[best] ** 2))),
        "max_abs_residual_rhos": float(np.max(np.abs(error[best]))),
        "endpoint_contrast_rhos": (model_curve[-1] - model_curve[0]).tolist(),
        "compatibility": compatible,
    }


def benthic_reflectance(observed, depth, a, bb, solar_zenith):
    """Unclipped diagnostic rho_b; invalid physical values remain visible."""
    a, bb = np.asarray(a), np.asarray(bb)
    u = bb / (a + bb)
    inv_cos = 1 / np.cos(np.deg2rad(lee._subsurface_solar_zenith_deg(solar_zenith)))
    kc, kb = (inv_cos + lee._du_column(u)) * (a + bb), (inv_cos + lee._du_bottom(u)) * (a + bb)
    sub = water.rrs_above_to_subsurface(np.asarray(observed) / np.pi)
    return np.pi * (sub - (0.089 * u + 0.125 * u * u) * (1 - np.exp(-kc * depth))) * np.exp(kb * depth)


def purged_folds(frame, minimum_training=8, minimum_per_stratum=3):
    """Hold a transect number across both strata; purge any shared common-grid cells."""
    folds = []
    for transect in sorted(frame.transect.unique()):
        held = frame[frame.transect == transect]
        rest = frame[frame.transect != transect]
        shared = set().union(*(set(v) for v in held.cell_ids))
        purge = rest.cell_ids.map(lambda values, shared=shared: bool(shared.intersection(values)))
        train = rest[~purge]
        counts = train.groupby("stratum").size()
        sufficient = len(train) >= minimum_training and len(counts) == 2 and bool((counts >= minimum_per_stratum).all())
        folds.append(
            {
                "transect": int(transect),
                "held_indices": held.index.tolist(),
                "train_indices": train.index.tolist(),
                "purged_indices": rest[purge].index.tolist(),
                "sufficient": sufficient,
                "training_by_stratum": {k: int(v) for k, v in counts.items()},
            }
        )
    return folds
