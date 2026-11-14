import numpy as np
import pytest

from kelpobserve.monthly.evaluation import PAIRS
from kelpobserve.monthly.pneo import patch_features


def test_nested_native_windows_use_containing_pixel_without_resampling():
    patch = np.arange(150, dtype=np.float32).reshape(6, 5, 5) + 10
    coefficients = {f"DII_{a}_{b}": 1.0 for a, b in PAIRS}
    rows = patch_features(patch, np.zeros(6), coefficients)
    assert [r["window_px"] for r in rows] == [1, 3, 5]
    assert rows[0]["R_mean"] == patch[0, 2, 2]
    assert rows[0]["R_std"] == rows[0]["R_range"] == 0
    assert rows[1]["R_mean"] == pytest.approx(patch[0, 1:4, 1:4].mean())
    assert rows[2]["R_range"] == 24
    expected = np.log(patch[5, 2, 2]) - np.log(patch[2, 2, 2])
    assert rows[0]["DII_DB_B"] == pytest.approx(expected)
    patch[0, 0, 0] = 0
    with pytest.raises(ValueError, match="positive"):
        patch_features(patch, np.zeros(6), coefficients)
