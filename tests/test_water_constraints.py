import json

import numpy as np
import pytest

from kelpobserve import config
from kelpobserve.attenuation import bands
from kelpobserve.coastal.optics import water
from kelpobserve.water_constraints import WaterConstraintError, separate_iops


@pytest.fixture
def reference():
    return json.loads((config.OPTICS / "reports/e5s-kelp-spectra/results.json").read_text())


def test_kd_and_backscatter_are_separately_reproduced(reference):
    s = separate_iops(0.0784245, 0.000601314, 21, [443, 490, 563], reference, cdm443=0.0795, cdm_error_pct=57)
    np.testing.assert_allclose(s["bb"][0] - water.bb_water(np.array([443]))[0], 0.000601314, atol=1e-12)
    np.testing.assert_allclose(bands.kd(np.array(s["a"]), np.array(s["bb"]), 21)[1], 0.0784245, atol=1e-12)
    assert s["absorption_scale"] != s["backscatter_scale"]
    assert s["cdm_nominal_conflict"]
    assert not s["cdm_conflict"]
    assert s["input_constraints_consistent"]
    assert not s["field_validated"]


def test_kd_below_feasible_water_is_rejected(reference):
    with pytest.raises(WaterConstraintError, match="infeasible"):
        separate_iops(0.01, 0.1, 21, [486, 563], reference)


def test_cdm_conflict_is_not_qualified(reference):
    s = separate_iops(0.08, 0.001, 21, [486, 563], reference, cdm443=1)
    assert s["cdm_conflict"]
    assert not s["input_constraints_consistent"]
    assert s["status"] == "input_conflict"


def test_missing_cdm_check_is_reported(reference):
    s = separate_iops(0.08, 0.001, 21, [486, 563], reference)
    assert not s["cdm_check_available"]
    assert not s["input_constraints_consistent"]


@pytest.mark.parametrize("wave", [[439], [700], [704], [800]])
def test_unphysical_reference_or_extrapolation_is_rejected(reference, wave):
    with pytest.raises(WaterConstraintError):
        separate_iops(0.08, 0.001, 21, wave, reference)


def test_outside_requested_support_is_disclosed(reference):
    s = separate_iops(0.08, 0.001, 21, [486, 563], reference)
    assert s["invalid_reference_wavelengths_nm"] == [704.0]
