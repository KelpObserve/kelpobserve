"""Scientific invariants of monthly LOO, not fitted-model performance thresholds."""

import json

import numpy as np
import pandas as pd
import pytest

from kelpobserve.monthly import evaluation as m


def sample():
    return pd.DataFrame(
        {
            "event_id": [f"july-shallow-T1-Q{i}" for i in range(6)],
            "campaign": ["july"] * 6,
            "cover_pct": [0.0, 0.0, 0.0, 0.0, 60.0, 100.0],
            "depth": [1.0, 2.0, 3.0, 4.0, 5.0, 20.0],
            "stratum": ["shallow"] * 6,
            "transect": [1] * 6,
        }
    )


def test_folds_reject_mixed_months_and_exclude_only_one_quadrat():
    q = sample()
    for train, test in m.folds(q):
        assert test not in train and len(train) == len(q) - 1
        assert set(train) | {test} == set(range(len(q)))
    q.loc[1, "campaign"] = "september"
    with pytest.raises(ValueError, match="exactly one month"):
        list(m.folds(q))


def test_held_out_target_cannot_change_its_prediction_or_background():
    q = sample()
    values = np.arange(6 * len(m.STATS), dtype=float).reshape(6, -1) + 10
    features = pd.DataFrame(values, columns=m.STATS).assign(event_id=q.event_id, window_px=5)
    for name in ["mean", "median", "ridge"]:
        for background in ["raw", "global", "local"] if name == "ridge" else ["raw"]:
            experiment = m.Experiment("test", name, "spectral_stats", 5, background)
            before, audits_before = m.predict_experiment(q, features, experiment)
            changed = q.copy()
            changed.loc[0, "cover_pct"] = 95.0  # Change held-out row from empty to occupied.
            after, audits_after = m.predict_experiment(changed, features, experiment)
            assert before[0] == pytest.approx(after[0], abs=1e-10)
            if audits_before:
                assert audits_before[0] == audits_after[0]
                audit = audits_before[0]
                assert q.event_id.iloc[0] not in audit["test_controls"]
                for event, references in audit["training_references"].items():
                    assert event not in references and q.event_id.iloc[0] not in references


def test_depth_standardization_is_training_only():
    q = sample()
    predicted, _ = m.predict_experiment(q, pd.DataFrame(), m.Experiment("depth", "ridge", "depth"))
    # Direct closed-form one-variable ridge solution for the last, outlying depth.
    x = q.depth.to_numpy()[:-1]
    y = q.cover_pct.to_numpy()[:-1]
    z = (x - x.mean()) / x.std()
    coefficient = z @ (y - y.mean()) / (z @ z + 1.0)
    expected = y.mean() + coefficient * ((q.depth.iloc[-1] - x.mean()) / x.std())
    assert predicted[-1] == pytest.approx(expected)


def test_other_month_cannot_change_training_samples_or_predictions():
    q = sample().assign(quality_eligible=True)
    other = q.assign(campaign="september", event_id=lambda x: x.event_id.str.replace("july", "september"))
    combined = pd.concat([q, other])
    before = m.select_samples(combined, "july", "matched_depth")
    combined.loc[combined.campaign == "september", "cover_pct"] = 99
    pd.testing.assert_frame_equal(before, m.select_samples(combined, "july", "matched_depth"))


def test_median_is_recomputed_and_saved_sample_counts_are_explicit():
    observations, features = m.load_inputs()
    expected = {
        ("july", "all_cover"): 30,
        ("september", "all_cover"): 29,
        ("july", "matched_depth"): 30,
        ("september", "matched_depth"): 25,
        ("september", "quality_filtered"): 21,
    }
    for (campaign, subset), n in expected.items():
        q = m.select_samples(observations, campaign, subset)
        assert len(q) == n
        predicted, _ = m.predict_experiment(q, features, m.Experiment("median", "median"))
        manual = [np.median(q.cover_pct.drop(i)) for i in range(len(q))]
        np.testing.assert_array_equal(predicted, manual)
        if campaign == "july":
            assert set(predicted) == {50.0, 60.0}
        if campaign == "september":
            assert set(predicted) == {5.0}


def test_input_tampering_is_rejected(tmp_path):
    (tmp_path / "observations.csv").write_text("changed")
    (tmp_path / "manifest.json").write_text(json.dumps({"files": {"observations.csv": "bad"}}))
    with pytest.raises(ValueError, match="checksum"):
        m.load_inputs(tmp_path)
