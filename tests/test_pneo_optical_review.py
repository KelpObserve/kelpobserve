from types import SimpleNamespace

import numpy as np

from kelpobserve.pneo_optical_review import ROLES, NativeCorrection


def test_negative_native_nir_is_handled_before_area_averaging():
    correction = object.__new__(NativeCorrection)
    correction.nir = np.array([[-0.002, 0.004]])
    correction.valid = np.ones((1, 2), dtype=bool)
    correction.data = {band: np.full((1, 2), 0.01) for band in ROLES}
    correction.row = correction.col = np.array([0])
    correction.transform = None
    correction.shape = (1, 1)
    correction.win = SimpleNamespace(crs="EPSG:32629")
    correction.rv = SimpleNamespace(_average=lambda values, *_: np.array([[np.nanmean(values)]]))
    # Native corrected values are .010 and .006, whose mean is .008.
    # Applying the correction to mean NIR (.001) would incorrectly return .009.
    np.testing.assert_allclose(correction.apply(np.ones(4), 0), [[0.008] * 4])
    # An explicit .003 baseline changes only the positive excess on the second pixel.
    np.testing.assert_allclose(correction.apply(np.ones(4), 0.003), [[0.0095] * 4])
