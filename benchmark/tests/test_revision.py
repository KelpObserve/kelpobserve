"""Revision r1: QC propagation, flag classes, scenarios and the injected signal."""

import json
import math

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from rasterio.transform import from_origin

from benchmark import revision


def make_window(shape=(12, 12), pixel=1.2):
    bands = {role: np.full(shape, 0.02) for role in revision.ROLES}
    return revision.Window(
        scene_id="pneo-test",
        sensor="pleiades_neo",
        transform=from_origin(0.0, shape[0] * pixel, pixel, pixel),
        crs="EPSG:32629",
        bands=bands,
        mask_ok=np.ones(shape, bool),
        water=np.ones(shape, bool),
        lidar=np.full(shape, 10.0),
        solar_zenith_deg=21.0,
    )


def make_events(points):
    return gpd.GeoDataFrame(
        {
            "event_id": [f"july-shallow-T1-Q{i + 1}" for i in range(len(points))],
            "transect": [1] * len(points),
            "stratum": ["shallow"] * len(points),
        },
        geometry=gpd.points_from_xy(*zip(*points)),
        crs=32629,
    )


@pytest.fixture
def flat_lidar(monkeypatch):
    monkeypatch.setattr(revision, "lidar_at", lambda spec, xy: np.full(len(xy), 8.0))


def test_flag_classes_follow_the_registration():
    classes = revision.classify_flags(
        [
            "a_cdm_443_not_computable(no_411nm_band)",
            "a_below_pure_water_floor(826nm: 0.41 < 2.85)",
            "a_below_pure_water_floor(655nm: 0.27 < 0.38)",
            "bbp_555_implausible(0.06217)",
            "red_reference_selected(655nm)",
        ]
    )
    assert classes["informational"] == ["a_cdm_443_not_computable(no_411nm_band)"]
    assert classes["band_limited"] == ["a_below_pure_water_floor(826nm: 0.41 < 2.85)"]
    assert len(classes["fatal"]) == 3


def test_informational_flag_alone_does_not_block():
    classes = revision.classify_flags(["a_cdm_443_not_computable(no_411nm_band)"])
    assert classes["fatal"] == []


def test_native_unit_with_a_masked_pixel_is_excluded(flat_lidar):
    window = make_window()
    window.mask_ok[5, 6] = False  # one pixel of the 3 x 3 around (6, 6)
    events = make_events([(6 * 1.2 + 0.6, (12 - 6) * 1.2 - 0.6), (2 * 1.2 + 0.6, 10 * 1.2 - 0.6)])
    units = revision.native_units({}, events, window, 1.0, "july")
    assert list(units.reason) == ["masked", "ok"]
    assert math.isnan(units.green[0]) and units.green[1] == pytest.approx(0.02)
    assert units.depth[1] == pytest.approx(9.0)
    # v1's finite-only rule keeps it: the asymmetry the revision removes.
    as_v1 = revision.native_units({}, events, window, 1.0, "july", rule="v1")
    assert list(as_v1.reason) == ["ok", "ok"]


def test_native_unit_touching_land_is_excluded(flat_lidar):
    window = make_window()
    window.water[:, :3] = False
    events = make_events([(2 * 1.2 + 0.6, 6 * 1.2 - 0.6)])
    assert revision.native_units({}, events, window, 1.0, "july").reason[0] == "nonwater"


def test_window_size_widens_the_footprint(flat_lidar):
    window = make_window()
    window.mask_ok[4, 3] = False  # inside the 5 x 5 around (6, 5), outside the 3 x 3
    events = make_events([(5 * 1.2 + 0.6, (12 - 6) * 1.2 - 0.6)])
    assert revision.native_units({}, events, window, 0.0, "july", size=3).reason[0] == "ok"
    assert revision.native_units({}, events, window, 0.0, "july", size=5).reason[0] == "masked"


def test_shift_moves_the_extraction(flat_lidar):
    window = make_window()
    window.mask_ok[:, 8:] = False
    events = make_events([(4 * 1.2 + 0.6, 6 * 1.2 - 0.6)])
    assert revision.native_units({}, events, window, 0.0, "july").reason[0] == "ok"
    shift = np.array([[6.0, 0.0]])  # 5 pixels east
    assert revision.native_units({}, events, window, 0.0, "july", shift=shift).reason[0] == "masked"


def test_grid_unit_needs_valid_not_merely_finite(monkeypatch):
    window = make_window(shape=(25, 25), pixel=1.2)
    window.water[:, :12] = False  # the western half of the window is land
    grid = from_origin(0.0, 30.0, 10.0, 10.0)
    monkeypatch.setattr(
        revision.rb, "common_grid", lambda spec, bounds, ref: (grid, (3, 3), "EPSG:32629")
    )
    monkeypatch.setattr(revision, "raster_on", lambda *a, **k: np.full((3, 3), 9.0, "float32"))
    spec = {
        "common_grid": {"min_valid_fraction": 0.8},
        "selected_grid_scene": "x",
        "sources": {"lidar": "unused"},
    }
    events = make_events([(5.0, 15.0), (25.0, 15.0)])
    units = revision.grid_units(spec, events, window, 1.0, "july")
    by_column = units.set_index("column")
    assert by_column.loc[0, "reason"] == "nonwater" and math.isnan(by_column.loc[0, "green"])
    assert by_column.loc[2, "reason"] == "ok" and by_column.loc[2, "depth"] == pytest.approx(10.0)
    as_v1 = revision.grid_units(spec, events, window, 1.0, "july", rule="v1")
    assert set(as_v1.reason) == {"ok"}


def test_reversed_order_swaps_within_transect(monkeypatch):
    events = pd.DataFrame(
        {
            "event_id": [f"july-deep-T1-Q{q}" for q in range(1, 6)] + ["july-shallow-T1-Q1"],
            "campaign": ["july"] * 6,
            "stratum": ["deep"] * 5 + ["shallow"],
            "transect": [1] * 6,
            "quadrat": [1, 2, 3, 4, 5, 1],
        }
    )
    monkeypatch.setattr(gpd, "read_file", lambda *a, **k: events)
    targets = pd.DataFrame(
        {
            "event_id": events.event_id,
            "cover_pct": [10.0, 20.0, 30.0, 40.0, 50.0, 99.0],
            "cover_status": ["measured"] * 6,
        }
    )
    out = revision.reversed_targets(targets).set_index("event_id").cover_pct
    assert out["july-deep-T1-Q1"] == 50.0 and out["july-deep-T1-Q2"] == 40.0
    assert out["july-deep-T1-Q3"] == 30.0
    assert math.isnan(out["july-shallow-T1-Q1"])  # its partner Q5 is not in the data


def test_signal_is_zero_without_kelp_and_fades_with_depth():
    kelp, understory = np.array([0.042, 0.073, 0.066]), np.array([0.020, 0.034, 0.041])
    iops = {"a": np.array([0.09, 0.095, 0.42]), "bb": np.array([0.015, 0.013, 0.0098])}
    depth = np.array([[5.0, 15.0]])
    none = revision.signal_rhos(np.zeros_like(depth), depth, kelp, understory, iops, 21.0)
    full = revision.signal_rhos(np.ones_like(depth), depth, kelp, understory, iops, 21.0)
    assert all(np.allclose(b, 0.0) for b in none)
    green = full[1][0]
    assert green[0] > green[1] > 0  # kelp brighter than understory, fading with depth
    half = revision.signal_rhos(np.full_like(depth, 0.5), depth, kelp, understory, iops, 21.0)
    assert half[1][0, 0] == pytest.approx(green[0] / 2, rel=0.02)


def test_signal_skips_pixels_without_depth():
    depth = np.array([[np.nan, 6.0]])
    out = revision.signal_rhos(
        np.ones_like(depth),
        depth,
        np.array([0.04, 0.07, 0.06]),
        np.array([0.3, 0.3, 0.3]),
        {"a": np.array([0.09, 0.095, 0.42]), "bb": np.array([0.015, 0.013, 0.0098])},
        21.0,
    )
    assert out[0][0, 0] == 0.0 and out[0][0, 1] < 0  # kelp darker than bare substrate


def test_loto_scores_every_family_on_the_same_units():
    spec = {
        "evaluation": {"min_training_units": 4, "min_training_units_per_stratum": 2},
        "models": {"ridge": {"alpha": 1.0}},
    }
    rng = np.random.default_rng(0)
    rows = []
    for transect in (1, 2, 3):
        for stratum in ("shallow", "deep"):
            for k in range(2):
                rows.append(
                    {
                        "unit_id": f"u{transect}{stratum}{k}",
                        "transect": transect,
                        "stratum": stratum,
                        "depth": 6.0 if stratum == "shallow" else 16.0,
                        **{r: rng.uniform(0.01, 0.03) for r in revision.ROLES},
                        "cover_pct": rng.uniform(0, 100),
                    }
                )
    units = pd.DataFrame(rows)
    units.loc[0, "nir"] = np.nan  # this unit must leave every family, not only the NIR ones
    family_list = [
        {"id": "training_mean", "kind": "mean", "predictors": []},
        {"id": "depth_ridge", "kind": "ridge", "predictors": ["depth"]},
        {"id": "nir_ridge", "kind": "ridge", "predictors": ["nir"]},
    ]
    frame, folds = revision.loto(units, spec, family_list)
    counts = frame.groupby("family").unit_id.nunique()
    assert set(counts) == {11}
    assert "u1shallow0" not in set(frame.unit_id)
    assert all(f["sufficient"] for f in folds)


def test_clean_turns_nan_into_null():
    value = {"a": float("nan"), "b": [np.float32(1.5), np.inf], "c": np.int64(3)}
    assert json.dumps(revision.clean(value)) == '{"a": null, "b": [1.5, null], "c": 3}'


def test_clipped_read_fills_the_overhang_with_nan(tmp_path):
    import rasterio
    from rasterio.windows import Window as RasterWindow

    path = tmp_path / "band.tif"
    data = np.arange(16, dtype="float32").reshape(4, 4)
    profile = {
        "driver": "GTiff",
        "height": 4,
        "width": 4,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:32629",
        "transform": from_origin(0, 4, 1, 1),
        "nodata": -1.0,
    }
    data[0, 0] = -1.0
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data, 1)
    out = revision.read_clipped(path, RasterWindow(-1, -1, 3, 3))
    assert np.isnan(out[0]).all() and np.isnan(out[:, 0]).all()
    assert np.isnan(out[1, 1])  # nodata becomes NaN
    assert out[2, 2] == 5.0


def test_improvement_rule_and_missing_measures():
    v0 = {
        "M1_offshore": {"nir_median": 0.010},
        "M3_common_mode": {"ratio": 0.5},
        "M4_contrast": {"ratio": 3.0},
    }
    better = {
        "M1_offshore": {"nir_median": 0.003},
        "M3_common_mode": {"ratio": 0.1},
        "M4_contrast": {"ratio": 2.8},
    }
    assert revision.improvement(v0, better, 0.001)["improves"] is True
    assert revision.improvement(v0, better, 0.004)["improves"] is False
    blind = better | {"M1_offshore": {"nir_median": None}}
    assert revision.improvement(v0, blind, 0.001)["improves"] is None
