"""Revision r4: forward and inverse physics, oracle, cover field, aggregation, conditions."""

from types import SimpleNamespace

import numpy as np
import pytest
from rasterio.transform import from_origin

from benchmark import revision, revision_r4

A, BB, SZA = 0.09, 0.0154, 21.0


def test_surface_rhos_is_pi_times_the_library_forward_model():
    from kelpobserve.coastal.inversion.lee import forward_model

    rho_b, depth = np.array([[0.05, 0.3]]), np.array([[6.0, 17.0]])
    iops = SimpleNamespace(a=np.array([A]), bb=np.array([BB]))
    expected = np.pi * forward_model(rho_b, depth, iops, SZA)
    assert np.allclose(revision_r4.surface_rhos(rho_b, depth, A, BB, SZA), expected)


def test_bottom_reflectance_inverts_the_forward_model_at_both_depths():
    rho_b, depth = np.array([0.02, 0.07, 0.3]), np.array([6.0, 17.0, 17.0])
    rhos = revision_r4.surface_rhos(rho_b, depth, A, BB, SZA)
    assert np.allclose(revision_r4.bottom_reflectance(rhos, depth, A, BB, SZA), rho_b)


def test_bottom_reflectance_is_not_clipped_where_the_library_would_drop_it():
    from kelpobserve.coastal.inversion.lee import invert_scene

    depth = np.array([[17.0]])
    rhos = revision_r4.surface_rhos(np.array([[0.05]]), depth, A, BB, SZA) + 0.02
    mine = revision_r4.bottom_reflectance(rhos, depth, A, BB, SZA)
    iops = SimpleNamespace(wavelengths_nm=np.array([490.0]), a=np.array([A]), bb=np.array([BB]))
    library = invert_scene(
        (rhos / np.pi)[None], np.array([490.0]), depth, iops, np.ones((1, 1), bool), SZA
    )
    assert np.isfinite(mine[0, 0]) and mine[0, 0] > 1.5
    assert np.isnan(library[0, 0, 0])


def test_oracle_recovers_a_known_mixture_and_clips():
    kelp, bg = np.array([0.04, 0.07, 0.06]), np.array([0.02, 0.034, 0.04])
    mix = [0.3 * kelp[b] + 0.7 * bg[b] for b in range(2)]
    assert revision_r4.unmix([np.array([m]) for m in mix], kelp, bg)[0] == pytest.approx(30.0)
    assert revision_r4.unmix([np.array([1.0]), np.array([1.0])], kelp, bg)[0] == 100.0


def test_cover_field_respects_the_radius_and_the_stratum_background():
    transform = from_origin(0.0, 40.0, 1.0, 1.0)
    rows, cols = np.indices((40, 40))
    xs, ys = transform * (cols + 0.5, rows + 0.5)
    fine = revision_r4.Fine(transform, (40, 40), np.column_stack([xs.ravel(), ys.ravel()]))
    xy = np.array([[10.0, 20.0], [30.0, 20.0]])
    field = revision_r4.cover_field(fine, xy, np.array([1.0, 0.0]), np.array(["s", "s"]), 5.0)
    assert field[20, 10] == 1.0 and field[20, 30] == 0.0
    assert field[5, 5] == 0.5  # outside both discs: the stratum mean
    assert (field == 1.0).sum() == pytest.approx(np.pi * 25, rel=0.1)


def test_aggregation_conserves_the_signal():
    fine_t = from_origin(0.0, 100.0, 1.0, 1.0)
    rows, cols = np.indices((100, 100))
    xs, ys = fine_t * (cols + 0.5, rows + 0.5)
    fine = revision_r4.Fine(fine_t, (100, 100), np.column_stack([xs.ravel(), ys.ravel()]))
    values = np.zeros((100, 100))
    values[40:50, 40:50] = 1.0
    window = SimpleNamespace(
        transform=from_origin(0.0, 100.0, 10.0, 10.0), shape=(10, 10), crs="EPSG:32629"
    )
    box = revision_r4.to_sensor(values, fine, window, psf=False)
    # Fine rows 40-49 (y 50-60) and columns 40-49 fall in sensor row 4, column 4.
    assert box[4, 4] == pytest.approx(1.0) and box.sum() == pytest.approx(1.0)
    blurred = revision_r4.to_sensor(values, fine, window, psf=True)
    assert blurred.sum() == pytest.approx(1.0, rel=1e-3) and blurred[4, 4] < 1.0


def test_conditions_match_the_registration():
    conditions = revision_r4.conditions()
    blocks = {}
    for c in conditions:
        blocks[c["block"]] = blocks.get(c["block"], 0) + 1
    assert blocks == {"A": 40, "B": 200, "C_iop": 8, "B_null": 4, "C_psf": 8}
    assert len({c["id"] for c in conditions}) == len(conditions)


def test_evaluate_scores_the_oracle_on_ideal_units():
    import pandas as pd

    kelp, bg = np.array([0.04, 0.07, 0.06]), np.array([0.02, 0.034, 0.04])
    iops = {"median": {"a": np.array([A, A, 0.42]), "bb": np.array([BB, BB, 0.01])}}
    rng = np.random.default_rng(0)
    rows = []
    for t in (1, 2, 3):
        for stratum, depth in (("shallow", 7.0), ("deep", 17.0)):
            for k in range(3):
                cover = rng.uniform(0, 1)
                rhos = [
                    float(
                        revision_r4.surface_rhos(
                            cover * kelp[b] + (1 - cover) * bg[b],
                            depth,
                            iops["median"]["a"][b],
                            iops["median"]["bb"][b],
                            SZA,
                        )
                    )
                    for b in range(3)
                ]
                rows.append(
                    {
                        "unit_id": f"{t}{stratum}{k}",
                        "transect": t,
                        "stratum": stratum,
                        "depth": depth,
                        "blue": rhos[0],
                        "green": rhos[1],
                        "red": rhos[2],
                        "nir": 0.0,
                        "cover_pct": 100 * cover,
                    }
                )
    spec = {
        "evaluation": {"min_training_units": 8, "min_training_units_per_stratum": 3},
        "models": {"ridge": {"alpha": 1.0}},
    }
    result = revision_r4.evaluate(pd.DataFrame(rows), spec, kelp, bg, iops, SZA)
    assert result["mae_oracle_unmixing"] == pytest.approx(0.0, abs=1e-6)
    assert result["mae_corrected_depth_ridge"] < result["mae_depth_ridge"]
    assert result["units"] == 18


def test_geometries_cover_both_sentinel2_supports():
    assert set(revision_r4.GEOMETRIES) == {
        "G1_pneo_native",
        "G2_pneo_10m",
        "G3_s2_2jul",
        "G4_s2_9jul",
    }
    assert revision.ROLES == ("blue", "green", "red", "nir")
