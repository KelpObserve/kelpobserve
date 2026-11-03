import json

import numpy as np
from kelpobserve import config, cover_retrieval
from kelpobserve.attenuation import bands
from kelpobserve.coastal.optics import water
from rasterio.transform import from_origin

import diagnose as d


def test_vector_signal_matches_frozen_forward_model():
    a, bb = np.array([0.06, 0.08]), np.array([0.004, 0.003])
    for depth in (0.5, 10, 60):
        for rho in (0, 0.6):
            expected = cover_retrieval.curves([0], depth, a, bb, 23, [rho, rho], [rho, rho])[0]
            assert np.allclose(d.optical_signal(a, bb, 23, depth, rho), expected, atol=1e-12)


def test_horizon_sees_wall_and_clear_sea():
    terrain = np.zeros((800, 800))
    transform = from_origin(480000, 4255000, 10, 10)
    point = np.array([[484000, 4251000]])
    # A 100 m high wall starting about 100 m north of the water target.
    terrain[350:391, 300:500] = 100
    angles, coverage = d.ray_horizons(terrain, transform, point, [0, 180], 0)
    assert 43 < angles[0, 0] < 47
    assert angles[0, 1] == 0
    assert coverage.min() == 1


def test_independent_water_constraints_are_reproduced():
    spectra = json.loads((config.OPTICS / "reports/e5s-kelp-spectra/results.json").read_text())
    audit = json.loads((config.BENCHMARK / "registration/evidence/cmems-optical-inputs-20261102.json").read_text())
    match = next(t for t in audit["targets"] if t["date"] == "2026-07-08")
    rows = d.candidate_states(spectra, match, [443, 490], 21)
    for row in rows:
        if row["candidate"] and row["valid_solution"]:
            assert np.isclose(row["bb"][0] - water.bb_water(np.array([443]))[0], row["bbp443_target"])
            assert np.isclose(bands.kd(np.array(row["a"]), np.array(row["bb"]), 21)[1], row["kd490_target"])
            assert row["cdm_conflict"] == (row["a_nonwater443"] < row["cdm443_lower"])


def test_impossible_kd_bbp_is_reported():
    spectra = json.loads((config.OPTICS / "reports/e5s-kelp-spectra/results.json").read_text())
    match = {
        "variables": {k: {"value": v, "uncertainty_pct": 0} for k, v in [("KD490", 0.02), ("BBP", 1.0), ("CDM", 0.02)]}
    }
    rows = d.candidate_states(spectra, match, [492, 559], 21)
    assert all(not r["valid_solution"] for r in rows if r["candidate"])
