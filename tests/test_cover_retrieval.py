import numpy as np
import pandas as pd
import pytest

from kelpobserve import cover_retrieval as cr
from kelpobserve.coastal.inversion import lee
from kelpobserve.coastal.optics.qaa import SceneIOPs

A = np.array([0.08, 0.10])
BB = np.array([0.01, 0.008])
KELP = np.array([0.04, 0.07])
BASE = np.array([0.02, 0.035])
F = np.linspace(0, 1, 1001)


def test_curve_matches_the_library_forward_model():
    iops = SceneIOPs(np.array([492.0, 560.0]), A, BB, np.zeros(2), 0, 0, 0, 0, 560.0)
    for fraction in [0.0, 0.3, 1.0]:
        rho = BASE + fraction * (KELP - BASE)
        expected = np.pi * lee.forward_model(rho, 7.0, iops, 25.0)
        np.testing.assert_allclose(cr.curves([fraction], 7.0, A, BB, 25.0, KELP, BASE)[0], expected, rtol=1e-12)


@pytest.mark.parametrize("depth", [0.5, 7.0, 17.0])
@pytest.mark.parametrize("fraction", [0.0, 0.17, 0.55, 1.0])
def test_exact_input_cover_recovery_and_unclipped_reflectance(depth, fraction):
    observed = cr.curves([fraction], depth, A, BB, 25.0, KELP, BASE)[0]
    result = cr.retrieve(observed, cr.curves(F, depth, A, BB, 25.0, KELP, BASE), F, [0.001])
    assert result["cover_pct"] == pytest.approx(100 * fraction, abs=0.101)
    assert result["compatibility"]["0.001"]["any"]
    np.testing.assert_allclose(
        cr.benthic_reflectance(observed, depth, A, BB, 25.0), BASE + fraction * (KELP - BASE), atol=1e-11
    )


def test_incompatible_observation_is_not_rescued_by_clipping():
    curve = cr.curves(F, 7.0, A, BB, 25.0, KELP, BASE)
    result = cr.retrieve(curve[-1] + 0.1, curve, F, [0.001, 0.002])
    assert result["cover_pct"] == 100 and result["boundary"]
    assert not result["compatibility"]["0.001"]["any"]
    assert result["compatibility"]["0.001"]["min_pct"] is None


def test_identical_endmembers_retain_full_ambiguity():
    curve = cr.curves(F, 17.0, A, BB, 25.0, BASE, BASE)
    result = cr.retrieve(curve[0], curve, F, [0.0005])
    interval = result["compatibility"]["0.0005"]
    assert interval["min_pct"] == 0 and interval["max_pct"] == 100


def test_missing_input_is_rejected():
    with pytest.raises(ValueError):
        cr.retrieve([np.nan, 0.01], np.ones((1001, 2)), F, [0.001])
    with pytest.raises(ValueError):
        cr.curves(F, -1, A, BB, 25.0, KELP, BASE)


def test_spatial_folds_hold_both_strata_and_purge_shared_cells():
    frame = pd.DataFrame(
        {
            "transect": [1, 1, 2, 2, 3, 3],
            "stratum": ["shallow", "deep"] * 3,
            "cell_ids": [[("a",)], [("b",)], [("a",)], [("c",)], [("d",)], [("e",)]],
        }
    )
    fold = cr.purged_folds(frame, minimum_training=3, minimum_per_stratum=1)[0]
    assert fold["held_indices"] == [0, 1]
    assert fold["purged_indices"] == [2]
    assert fold["train_indices"] == [3, 4, 5] and fold["sufficient"]
    assert not cr.purged_folds(frame, minimum_training=8)[0]["sufficient"]


def test_baselines_do_not_use_held_cover():
    from kelpobserve.local_cover import baselines

    frame = pd.DataFrame(
        [
            {
                "transect": t,
                "stratum": s,
                "cell_ids": [(t, s, j)],
                "eligible": True,
                "cover_pct": 20 * t + j,
                "depth": 5 if s == "shallow" else 15,
                "blue": 0.01 * j,
                "green": 0.02 * j,
            }
            for t in [1, 2, 3]
            for s in ["shallow", "deep"]
            for j in range(4)
        ]
    )
    protocol = {"min_training_units": 8, "min_training_per_stratum": 3, "ridge_alpha": 1}
    original, _ = baselines(frame.copy(), protocol)
    altered = frame.copy()
    altered.loc[altered.transect == 1, "cover_pct"] = 99
    second, _ = baselines(altered, protocol)
    for model in ["mean", "depth", "raw", "raw_depth"]:
        np.testing.assert_equal(
            original.loc[frame.transect == 1, "baseline_" + model].to_numpy(),
            second.loc[frame.transect == 1, "baseline_" + model].to_numpy(),
        )


def test_compatible_estimate_itself_satisfies_the_tolerance():
    fit = cr.retrieve([0.0, 0.0], np.array([[0.0011, 0.0], [0.0009, 0.0009]]), [0.0, 1.0], [0.001])
    assert fit["cover_pct"] == 0
    assert fit["compatibility"]["0.001"]["best_pct"] == 100


def test_score_handles_empty_or_missing_compatible_predictions():
    from kelpobserve.local_cover import paired_summary, score

    assert score([], []) == {"n": 0, "mae_pp": None, "bias_pp": None, "rmse_pp": None}
    assert score(pd.Series([None], dtype=object), [50])["n"] == 0
    frame = pd.DataFrame(
        {"compatible_cover_pct": pd.Series([], dtype=object), "baseline_depth": [], "cover_pct": [], "transect": []}
    )
    assert paired_summary(frame, "compatible_cover_pct", 42, 100) == {"n": 0}


def test_single_transect_does_not_produce_a_bootstrap_interval():
    from kelpobserve.local_cover import paired_summary

    frame = pd.DataFrame(
        {"cover": [0.0, 0.0], "cover_pct": [5.0, 10.0], "baseline_depth": [10.0, 15.0], "transect": [3, 3]}
    )
    result = paired_summary(frame, "cover", 42, 100)
    assert result["descriptive_cluster_interval_95"] is None
    assert result["interval_status"] == "not_estimable_single_transect"
