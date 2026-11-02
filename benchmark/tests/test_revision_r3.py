"""Revision r3: the glint-tolerant mask, acceptance rule, field aggregation and ranking."""

from datetime import UTC, datetime

import numpy as np
import pytest
from rasterio.transform import from_origin

from benchmark import revision, revision_r3

CLEAN = {
    "A1_nir_median": 0.008,
    "A3_ratio": 0.2,
    "A4_green_mad": 0.001,
    "A5_green_slope": 0.00015,
}


def make_window(shape=(30, 30), pixel=10.0, rhot=0.004):
    return revision.Window(
        scene_id="s2-test",
        sensor="sentinel2",
        transform=from_origin(0.0, shape[0] * pixel, pixel, pixel),
        crs="EPSG:32629",
        bands={role: np.full(shape, 0.02) for role in revision.ROLES},
        mask_ok=np.ones(shape, bool),
        water=np.ones(shape, bool),
        lidar=np.full(shape, 15.0),
        solar_zenith_deg=21.0,
        mask_rhot=np.full(shape, rhot),
    )


def test_threshold_is_v1_on_clean_scenes_and_rises_with_glint():
    assert revision_r3.tolerant_threshold(0.0215, 0.004, 0.005) == 0.0215
    assert revision_r3.tolerant_threshold(0.0215, 0.005, 0.005) == 0.0215
    assert revision_r3.tolerant_threshold(0.0215, 0.035, 0.005) == pytest.approx(0.0515)


def test_glinted_water_survives_the_tolerant_mask_but_a_bright_target_does_not():
    window = make_window(rhot=0.035)
    window.mask_rhot[5, 5] = 0.09  # cloud, boat or shoreline mixing above the glint field
    v1 = window.with_threshold(0.0215)
    assert not v1.valid().any()
    tolerant = window.with_threshold(revision_r3.tolerant_threshold(0.0215, 0.035, 0.005))
    assert tolerant.valid().sum() == 30 * 30 - 1 and not tolerant.valid()[5, 5]


def test_background_water_keeps_off_land_and_out_of_the_shallows(monkeypatch):
    window = make_window()
    window.water[:, :5] = False  # land along the western edge
    window.lidar[:, 20:] = 6.0  # shallow band in the east
    monkeypatch.setattr(revision, "raster_on", lambda *a, **k: np.full((30, 30), np.nan))
    background = revision_r3.background_water(
        {"sources": {"dtm": "x"}}, window, np.ones((30, 30), bool)
    )
    assert not background[:, :14].any()  # closer than 100 m to land (column 4)
    assert background[:, 14:20].all()  # 100 m and more counts
    assert not background[:, 20:].any()  # shallower than 10 m


def test_acceptance_needs_every_check_and_handles_numpy_booleans():
    clean = [CLEAN, CLEAN | {"A1_nir_median": 0.004}]
    good = {
        "A1_nir_median": np.float64(0.009),
        "A2_glint_slope": np.float64(0.05),
        "A3_ratio": np.float64(0.15),
        "A4_green_mad": np.float64(0.0015),
        "A5_green_slope": np.float64(0.0002),
    }
    assert revision_r3.acceptability(good, clean, glinted=True)["accepted"] is True
    leaky = good | {"A2_glint_slope": 0.6}
    assert revision_r3.acceptability(leaky, clean, glinted=True)["accepted"] is False
    assert (
        revision_r3.acceptability(leaky, clean, glinted=False)["checks"]["A2"] == "not applicable"
    )
    missing = good | {"A1_nir_median": float("nan")}
    assert revision_r3.acceptability(missing, clean, glinted=True)["accepted"] is False
    flipped = good | {"A5_green_slope": -0.0002}
    assert revision_r3.acceptability(flipped, clean, glinted=True)["accepted"] is False


def test_fine_field_weights_a_disc_by_its_area_in_the_pixel():
    window = make_window(shape=(3, 3))
    xy = np.array([[15.0, 15.0]])  # centre of the middle 10 m pixel
    field = revision_r3.fine_field(window, xy, np.array([1.0]), np.array(["shallow"]))
    # The background is the stratum mean, 1.0 here, so the pixel is 1.0 whatever the disc.
    assert field[1, 1] == pytest.approx(1.0)
    field = revision_r3.fine_field(
        window, np.vstack([xy, [[25.0, 25.0]]]), np.array([1.0, 0.0]), np.array(["s", "s"])
    )
    # On the 1 m grid a 2.5 m disc covers 16 of 100 m2 (19.6 exactly): 0.16 + 0.84 * 0.5.
    assert field[1, 1] == pytest.approx(0.58)


def test_ranking_prefers_the_date_nearest_the_survey():
    spec = {
        "campaigns": {
            "july": {"plausible_survey_dates": ["2026-07-07", "2026-07-08"], "pneo": None}
        }
    }
    rows = [
        {"scene": "s2-20260702-S2B-R037", "valid_water_fraction": 0.99},
        {"scene": "s2-20260709-S2A-R037", "valid_water_fraction": 0.85},
    ]
    for row in rows:
        row["sensing_utc"] = datetime(2026, 7, 1, tzinfo=UTC)
    ranked = sorted(rows, key=lambda r: revision_r3.rank_key(r, "july", spec))
    assert ranked[0]["scene"] == "s2-20260709-S2A-R037"


def test_c1_matches_e2c():
    assert revision_r3.c1({"reference_nm": 559.0, "bbp_555": 0.009, "kd_490": 0.15}) is True
    assert revision_r3.c1({"reference_nm": 665.0, "bbp_555": 0.009, "kd_490": 0.15}) is False
    assert revision_r3.c1({"pixels": 0, "available": False}) is None
