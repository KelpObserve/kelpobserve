"""Revision r5: scale, coherence and correction diagnostics on synthetic fields."""

import math

import numpy as np
import pytest

from benchmark import revision_r5 as r5


def test_semivariance_recovers_white_noise_variance():
    rng = np.random.default_rng(0)
    z = rng.normal(0, 0.002, (300, 300))
    mask = np.ones_like(z, bool)
    assert math.sqrt(r5.semivariance(z, mask, 1)) == pytest.approx(0.002, rel=0.05)
    assert math.sqrt(r5.semivariance(z, mask, 10)) == pytest.approx(0.002, rel=0.05)


def test_semivariance_grows_with_lag_for_a_gradient():
    z = np.tile(np.arange(200) * 1e-4, (200, 1))
    mask = np.ones_like(z, bool)
    assert r5.semivariance(z, mask, 10) > r5.semivariance(z, mask, 1)


def test_robust_sd_ignores_outliers():
    rng = np.random.default_rng(1)
    v = rng.normal(0, 0.001, 5000)
    v[:50] = 1.0
    assert r5.robust_sd(v) == pytest.approx(0.001, rel=0.1)


def test_blocks_require_every_pixel_valid():
    z = np.arange(16, dtype=float).reshape(4, 4)
    mask = np.ones((4, 4), bool)
    mask[0, 0] = False
    means, valid = r5.blocks(z, 2, mask)
    assert not valid[0, 0] and np.isnan(means[0, 0])
    assert means[1, 1] == pytest.approx(np.mean([10, 11, 14, 15]))


def test_high_pass_removes_a_constant_and_keeps_a_spike():
    z = np.full((9, 9), 0.02)
    z[4, 4] += 0.009
    hp = r5.high_pass(z, np.ones_like(z, bool))
    assert hp[1, 1] == pytest.approx(0.0, abs=1e-12)
    assert hp[4, 4] == pytest.approx(0.008)


def make_field(green, nir, swir=None, sensor="pleiades_neo"):
    shape = green.shape
    bands = {"blue": green.copy(), "green": green, "nir": nir}
    if swir is not None:
        bands["swir"] = swir
    f = r5.Field(
        "x",
        sensor,
        bands,
        np.ones(shape, bool),
        np.full(shape, 60.0),
        1.2 if sensor == "pleiades_neo" else 10.0,
    )
    f.regions = {"offshore": np.ones(shape, bool)}
    return f


def test_correction_diagnostic_flags_fine_scale_noise_from_the_reference():
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(2)
    shape = (240, 240)
    smooth = gaussian_filter(rng.normal(0, 1, shape), 8)
    glint = 0.004 * smooth / smooth.std()  # smooth glint, shared by both bands
    reference = glint + rng.normal(0, 0.002, shape)  # plus the reference band's own noise
    green = 0.02 + glint + rng.normal(0, 0.0005, shape)  # plus green's own, smaller noise
    v0 = make_field(green, reference)
    v1 = make_field(green - reference, reference)
    m6 = r5.m6_correction(v0, v1)
    assert m6["r_coarse"] > 0.7 and m6["r_fine"] < 0.3
    assert m6["r_fine_v1_vs_correction"] > 0.9  # V1's fine variation is the subtracted noise


def test_correction_diagnostic_detects_a_one_pixel_shift():
    rng = np.random.default_rng(3)
    base = rng.normal(0, 0.003, (200, 200))
    reference = np.roll(base, 1, axis=1)
    m6 = r5.m6_correction(
        make_field(0.02 + base, reference), make_field(0.02 + base - reference, reference)
    )
    assert m6["best_shift"] == "0,-1" and m6["r_fine_best_shift"] - m6["r_fine"] >= 0.1


def test_decision_rules_follow_the_registration():
    results = {
        "s2-20260709-S2A-R037|V1": {
            "m4": {
                "s2_10m": {
                    "offshore": {"blue": {"sd": 0.001}, "green": {"sd": 0.001}},
                    "5-10": {
                        "blue": {"sd": 0.002, "sd_after_nir_ratio_to_offshore": 2.0},
                        "green": {"sd": 0.0011, "sd_after_nir_ratio_to_offshore": 2.0},
                    },
                }
            },
            "m2": {"s2_10m": {"5-10": {"blue_on_nir": {"r2": 0.6}, "green_on_nir": {"r2": 0.1}}}},
            "m3": {"s2_10m": {"nir": {"distance_coef_per_decade": -0.001, "r2_gain": 0.1}}},
            "m6": {"r_fine": 0.2, "r_coarse": 0.8, "r_fine_best_shift": 0.25},
        },
        "s2-20260702-S2B-R037|V1": {
            "m4": {
                "s2_10m": {
                    "5-10": {
                        "blue": {"sd_after_nir_ratio_to_offshore": 1.6},
                        "green": {"sd_after_nir_ratio_to_offshore": 1.6},
                    }
                }
            },
            "m2": {},
            "m3": {},
        },
    }
    pairs = {
        "a~b|V1": {
            "days": 7,
            "scenes": list(results),
            "m5": {
                "s2_10m": {
                    "5-10": {"blue": {"r": 0.4}, "green": {"r": 0.1}},
                    "offshore": {"blue": {"r": 0.05}, "green": {"r": 0.05}},
                }
            },
        }
    }
    d = r5.decide(results, pairs, {"s2-20260709-S2A-R037|V1"})
    key = "s2-20260709-S2A-R037|V1|s2_10m|5-10"
    assert (
        d["D1_noise_floor"][key]["blue"]["met"] is False
        and d["D1_noise_floor"][key]["green"]["met"] is True
    )
    assert (
        d["D2_common_mode"][key]["blue"]["met"] is True
        and d["D2_common_mode"][key]["green"]["met"] is False
    )
    seabed = d["D3_seabed_candidate"]["a~b|V1|s2_10m|5-10"]
    assert seabed["blue"]["met"] is True and seabed["green"]["met"] is False
    glint = d["D4_per_pixel_glint"]["s2-20260709-S2A-R037|V1"]
    assert glint["smoothing_warranted"] is True and glint["shift_contributes"] is False
