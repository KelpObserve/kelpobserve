"""Revision r6: date combination, matched windows and the registered condition sets."""

import numpy as np
import pytest
from rasterio.transform import from_origin

from benchmark import revision, revision_r5, revision_r6


def window_with(values, valid):
    shape = values.shape
    return revision.Window(
        scene_id="s2",
        sensor="sentinel2",
        transform=from_origin(0, 10 * shape[0], 10, 10),
        crs="EPSG:32629",
        bands={r: values.copy() for r in revision.ROLES},
        mask_ok=valid,
        water=np.ones(shape, bool),
        lidar=np.full(shape, 10.0),
        solar_zenith_deg=21.0,
        mask_rhot=np.zeros(shape),
    )


def test_combination_averages_valid_dates_and_respects_the_cell_set():
    a, b = np.full((2, 2), 0.02), np.full((2, 2), 0.04)
    va = np.array([[True, True], [False, True]])
    vb = np.array([[True, False], [True, True]])
    wa, wb = window_with(a, va), window_with(b, vb)
    common, mask_common = revision_r6.combine([wa, wb], [[a] * 3, [b] * 3], require_all=True)
    union, mask_union = revision_r6.combine([wa, wb], [[a] * 3, [b] * 3], require_all=False)
    assert mask_common.tolist() == [[True, False], [False, True]]
    assert common[0][0, 0] == pytest.approx(0.03) and np.isnan(common[0][0, 1])
    assert union[0][0, 1] == pytest.approx(0.02) and union[0][1, 0] == pytest.approx(0.04)
    assert mask_union.all()


def test_restrict_only_narrows_validity():
    w = window_with(np.full((2, 2), 0.02), np.array([[True, False], [True, True]]))
    narrowed = revision_r6.restrict(w, np.array([[True, True], [False, True]]))
    assert narrowed.valid().tolist() == [[True, False], [False, True]]


def test_window_stats_separate_a_plane_from_local_noise():
    rng = np.random.default_rng(0)
    rows, cols = np.indices((300, 300))
    plane = 0.00002 * cols  # 0.002 per 100 m along x
    noise = rng.normal(0, 0.001, (300, 300))
    green = 0.02 + plane + noise
    f = revision_r5.Field(
        "x",
        "sentinel2",
        {"blue": green, "green": green, "nir": noise * 0},
        np.ones((300, 300), bool),
        np.full((300, 300), 60.0),
        1.0,
    )
    f.regions = {"offshore": np.ones((300, 300), bool)} | {
        b: np.zeros((300, 300), bool) for b in revision_r5.BINS
    }
    stats = revision_r6.window_stats(f, 100.0)["offshore"]
    assert stats["windows"] == 9
    assert stats["local"]["median"] == pytest.approx(0.001, rel=0.05)
    assert stats["broad"]["median"] == pytest.approx(0.00002 * np.std(np.arange(100)), rel=0.05)
    assert stats["local_after_nir"] is None  # NIR without variance is not used


def test_condition_sets_match_the_registration():
    pneo = revision_r6.pneo_conditions()
    assert len(pneo) == 2 * 2 * (2 * 25 + 1) + 2 * 2 * 5
    s2 = revision_r6.s2_conditions()
    assert len(s2) == 4 * (2 * (25 + 3) + 1)
    assert len({c["id"] for c in pneo}) == len(pneo) and len({c["id"] for c in s2}) == len(s2)


def test_valid_area_aggregation_ignores_undefined_cells():
    from types import SimpleNamespace

    fine_t = from_origin(0.0, 30.0, 1.0, 1.0)
    rows, cols = np.indices((30, 30))
    xs, ys = fine_t * (cols + 0.5, rows + 0.5)
    fine = revision_r6.r4.Fine(fine_t, (30, 30), np.column_stack([xs.ravel(), ys.ravel()]))
    values = np.full((30, 30), 0.02)
    values[:, :5] = np.nan  # no lidar along the western edge
    window = SimpleNamespace(
        transform=from_origin(0.0, 30.0, 10.0, 10.0), shape=(3, 3), crs="EPSG:32629"
    )
    out = revision_r6.to_sensor_valid(values, fine, window)
    assert np.allclose(out, 0.02)  # a constant stays constant, edges included


def test_offset_is_the_stratum_median_and_keeps_between_unit_differences():
    import pandas as pd

    real = pd.DataFrame(
        {
            "unit_id": ["a", "b", "c", "d"],
            "stratum": ["shallow", "shallow", "shallow", "deep"],
            "blue": [0.010, 0.012, 0.020, 0.030],
            "green": [0.010, 0.012, 0.020, 0.030],
            "red": [0.010, 0.012, 0.020, 0.030],
        }
    )
    model = real.assign(blue=0.015, green=0.015, red=0.015)
    offsets = revision_r6.offsets_from(real, model)
    assert offsets["shallow"]["green"] == pytest.approx(-0.003)  # median of -0.005, -0.003, +0.005
    assert offsets["deep"]["green"] == pytest.approx(0.015)
    adjusted = revision_r6.apply_offsets(real, offsets)
    assert adjusted.green.tolist() == pytest.approx([0.013, 0.015, 0.023, 0.015])
    shallow = adjusted[adjusted.stratum == "shallow"].green
    assert np.diff(shallow).tolist() == pytest.approx(np.diff(real.green[:3]).tolist())


def test_outcome_rules_follow_the_amendment():
    setting, contrast, model = "pneo|V0|native", "bare", "spectral_depth_ridge"
    keys = [f"{setting}|{contrast}|{contrast}|{r}|0.0" for r in (5.0, 10.0, 25.0)]

    def entries(shares, recovered=False):
        return {
            k: {model: {"share_ge_2pp": s, "recovered": recovered or s >= 0.8}}
            for k, s in zip(keys, shares)
        }

    null = {f"{setting}|None|{contrast}|None|0.0": {model: {"share_ge_2pp": 0.05}}}
    before = entries([0.40, 0.45, 0.45])
    out = revision_r6.classify(entries([0.85, 0.5, 0.5]), before, null)[
        f"{setting}|{contrast}|{model}"
    ]
    assert out["outcome"] == "recovery_appears"
    out = revision_r6.classify(entries([0.55, 0.6, 0.6]), before, null)[
        f"{setting}|{contrast}|{model}"
    ]
    assert out["outcome"] == "improves_but_fails"
    out = revision_r6.classify(entries([0.42, 0.47, 0.44]), before, null)[
        f"{setting}|{contrast}|{model}"
    ]
    assert out["outcome"] == "barely_changes"
    noisy = {k: {model: {"share_ge_2pp": 0.3}} for k in null}
    out = revision_r6.classify(entries([0.85, 0.5, 0.5]), before, noisy)[
        f"{setting}|{contrast}|{model}"
    ]
    assert out["outcome"] != "recovery_appears"  # the null must stay below 0.2
