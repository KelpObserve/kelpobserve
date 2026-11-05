"""Check that committed metrics, splits and plots describe the same held-out data."""

import json

import numpy as np
import pandas as pd
import pytest

from kelpobserve.monthly import evaluation as m


def test_saved_splits_and_metrics_match_their_inputs():
    observations, _ = m.load_inputs()
    predictions = pd.read_csv(m.RESULTS / "predictions.csv")
    metrics = pd.read_csv(m.RESULTS / "metrics.csv")
    protocol = json.loads((m.RESULTS / "protocol.json").read_text())
    assert protocol["input_manifest_sha256"] == m.digest(m.INPUTS / "manifest.json")
    splits = json.loads((m.RESULTS / "folds.json").read_text())
    for split in splits:
        q = m.select_samples(observations, split["campaign"], split["subset"])
        assert set(split["training"]) == set(q.event_id) - {split["held_out"]}
        assert split["held_out"] in set(q.event_id)
    assert not predictions.duplicated(m.KEYS + ["event_id"]).any()
    for row in metrics.itertuples():
        p = predictions[
            (predictions.campaign == row.campaign)
            & (predictions.subset == row.subset)
            & (predictions.experiment == row.experiment)
        ]
        q = m.select_samples(observations, row.campaign, row.subset)
        assert len(p) == len(q) == row.n
        np.testing.assert_array_equal(p.set_index("event_id").loc[q.event_id].observed, q.cover_pct)
        assert p.train_n.eq(len(q) - 1).all()
        expected = m.score(p.observed.to_numpy(), p.predicted.to_numpy())
        assert row.mae_pp == pytest.approx(expected["mae_pp"])
        assert row.rmse_pp == pytest.approx(expected["rmse_pp"])
        assert row.r2 == pytest.approx(expected["r2"])
