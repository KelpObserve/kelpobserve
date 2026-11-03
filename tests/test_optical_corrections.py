import numpy as np
import pytest

from kelpobserve.optical_corrections import GlintCalibrationError, apply_glint, fit_block_glint


def reference():
    blocks = np.repeat(np.arange(8), 40)
    nir = np.tile(np.linspace(0.005, 0.035, 40), 8)
    slopes = np.array([0.75, 0.9])
    # Block-dependent water backgrounds must not bias the glint slopes.
    background = np.column_stack([0.01 + blocks * 0.002, 0.008 + blocks * 0.003])
    visible = background + (nir - 0.005)[:, None] * slopes
    return nir, visible, blocks, slopes


def test_block_slopes_ignore_between_block_water_variation():
    x, y, blocks, truth = reference()
    result = fit_block_glint(x, y, blocks)
    np.testing.assert_allclose(result.slopes, truth, atol=1e-12)
    assert result.block_count == 8
    assert not result.absolute_radiometry_validated


def test_explicit_baseline_preserves_background_and_negative_diagnostics():
    result = apply_glint(
        [0.005, 0.025, 0.001], [[0.01, 0.02], [0.025, 0.038], [-0.001, 0.004]], [0.75, 0.9], nir_baseline=0.005
    )
    np.testing.assert_allclose(result["corrected"], [[0.01, 0.02], [0.01, 0.02], [-0.001, 0.004]])
    assert result["negative"][2, 0]
    np.testing.assert_allclose(result["nir_after"], [0.005, 0.005, 0.001])


def test_support_flags_extrapolation_without_hiding_diagnostics():
    result = apply_glint([0.001, 0.02, 0.2], np.ones((3, 2)), [0.8, 0.9], nir_baseline=0.004, support=(0.003, 0.05))
    assert result["within_nir_support"].tolist() == [False, True, False]
    assert np.isfinite(result["corrected"]).all()


def test_low_variation_reference_is_rejected():
    with pytest.raises(GlintCalibrationError, match="informative blocks"):
        fit_block_glint(np.ones(100) * 0.01, np.ones((100, 2)) * 0.02, np.repeat(np.arange(10), 10))


def test_nonfinite_reference_requires_explicit_selection():
    x, y, b, _ = reference()
    y[0, 0] = np.nan
    with pytest.raises(GlintCalibrationError, match="nonfinite"):
        fit_block_glint(x, y, b)


def test_invalid_or_implicit_baseline_is_not_accepted():
    with pytest.raises(TypeError):
        apply_glint([0.01], [[0.02]], [0.9])
    with pytest.raises(ValueError):
        apply_glint([0.01], [[0.02]], [0.9], nir_baseline=-0.001)
