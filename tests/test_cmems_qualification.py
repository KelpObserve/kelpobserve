import json
from datetime import date

import numpy as np
import pytest
import xarray as xr

from kelpobserve.attenuation import cmems, qualification

DAY = date(2026, 7, 8)
SITE = (0.0, 0.0)


def data():
    shape = (3, 3, 3)
    dims = ("time", "latitude", "longitude")
    return xr.Dataset({
        "KD490": (dims, np.full(shape, .1), {"units": "m-1", "valid_min": 0., "valid_max": 1000.}),
        "KD490_uncertainty": (dims, np.full(shape, 30.), {"units": "%"}),
        "flags": (dims, np.zeros(shape), {"flag_masks": 1, "flag_meanings": "LAND"}),
    }, coords={"time": np.array(["2026-07-07", "2026-07-08", "2026-07-09"], dtype="datetime64[ns]"),
               "latitude": [-.01, 0., .01], "longitude": [-.01, 0., .01]})


def select(ds, **kwargs):
    return cmems.nearest_valid(ds, DAY, site=SITE, **{"max_km": 5, "max_days": 1, **kwargs})


def test_finite_land_is_rejected_and_rejection_is_auditable():
    ds = data()
    ds.flags.values[:, 1, 1] = 1
    result = qualification.audit_match(ds, DAY, site=SITE)
    assert result["previous_finite_only_candidate"]["flags"] == 1
    assert result["previous_finite_only_candidate"]["rejection_reasons"] == ["land"]
    assert result["candidate"]["flags"] == 0
    assert result["candidate"]["distance_km"] == pytest.approx(1.112, abs=.001)
    assert result["candidate"]["valid_within_radius"] == 8
    day = result["search_days"][0]
    assert day["rejected_unique_within_radius"] == 1
    assert day["eligible_within_radius"] == 8
    assert day["rejection_reasons_overlapping"]["land"] == 1
    assert result["qualification"]["r7_fallback_selected"] is False


@pytest.mark.parametrize("variable,value", [
    ("KD490", 0), ("KD490", -1), ("KD490", np.nan), ("KD490", np.inf), ("KD490", 1001),
    ("KD490_uncertainty", np.nan), ("KD490_uncertainty", -1), ("KD490_uncertainty", np.inf),
    ("KD490_uncertainty", 327.67), ("flags", np.nan), ("flags", 2), ("flags", .5), ("flags", -1),
])
def test_invalid_values_cannot_be_selected(variable, value):
    ds = data()
    ds[variable].values[:, 1, 1] = value
    assert select(ds).distance_km > 1
    assert select(ds, max_km=.5) is None


def test_packed_valid_range_is_applied_after_decoding():
    ds = data()
    ds.KD490_uncertainty.attrs["valid_max"] = 5000
    ds.KD490_uncertainty.encoding["scale_factor"] = .01
    ds.KD490_uncertainty.values[:, 1, 1] = 51
    assert select(ds).distance_km > 1
    ds.KD490_uncertainty.values[:, 1, 1] = 50
    assert select(ds).distance_km == 0


@pytest.mark.parametrize("attrs", [{}, {"flag_masks": 2, "flag_meanings": "LAND"},
                                   {"flag_masks": [1, 2], "flag_meanings": "LAND INTERPOLATED"}])
def test_unreviewed_flag_schema_fails_closed(attrs):
    ds = data()
    ds.flags.attrs = attrs
    with pytest.raises(ValueError, match="flag schema"):
        select(ds)


def test_string_flag_mask_and_transposed_variables_are_supported():
    ds = data().transpose("longitude", "time", "latitude")
    ds.flags.attrs["flag_masks"] = "1"
    assert select(ds).kd490 == .1


@pytest.mark.parametrize("variable", ["flags", "KD490_uncertainty"])
def test_missing_quality_variable_fails_closed(variable):
    with pytest.raises(ValueError, match="dimensions"):
        select(data().drop_vars(variable))


def test_unknown_units_and_duplicate_days_fail_closed():
    ds = data()
    ds.KD490_uncertainty.attrs["units"] = "fraction"
    with pytest.raises(ValueError, match="percent"):
        select(ds)
    ds = data().isel(time=[0, 1, 1])
    with pytest.raises(ValueError, match="daily record"):
        select(ds)


@pytest.mark.parametrize("kwargs", [{"max_km": 0}, {"max_km": np.nan},
                                    {"max_days": -1}, {"max_days": 1.5}, {"max_days": True}])
def test_invalid_search_limits_raise(kwargs):
    with pytest.raises(ValueError):
        select(data(), **kwargs)


def test_quality_rejection_triggers_date_fallback_within_limits():
    ds = data()
    ds.flags.values[1] = 1
    assert select(ds).day_offset == -1
    assert select(ds, max_days=0) is None
    ds.flags.values[0] = 1
    assert select(ds).day_offset == 1
    ds.flags.values[2] = 1
    result = qualification.audit_match(ds, DAY, max_days=1, site=SITE)
    assert result["candidate"] is None
    assert result["qualification"]["status"] == "unavailable"
    json.dumps(result, allow_nan=False)


def test_temporal_check_keeps_coordinate_and_does_not_substitute_neighbours():
    ds = data()
    ds.KD490.values[0, 1, 1] = .05
    ds.KD490.values[2, 1, 1] = .2
    ds.flags.values[2, 1, 1] = 1
    result = qualification.audit_match(ds, DAY, site=SITE)
    temporal = result["temporal_support"]
    assert temporal["kd490"]["count"] == 2
    assert temporal["kd490"]["min"] == .05
    assert temporal["kd490"]["max"] == .1
    assert sum(d["status"] == "missing_day" for d in temporal["days"]) == 4
    assert sum(d["status"] == "rejected" for d in temporal["days"]) == 1


def test_neighbours_exclude_selected_and_land_pixels():
    ds = data()
    ds.KD490.values[1, 1, 1] = .5
    ds.KD490.values[1, 0, 0] = 10
    ds.flags.values[1, 0, 0] = 1
    result = qualification.audit_match(ds, DAY, site=SITE)
    neighbours = result["spatial_support"]["adjacent_3x3_excluding_selected"]
    assert neighbours["eligible_count"] == 7
    assert neighbours["kd490"]["max"] == .1
    assert result["spatial_support"]["within_1_km"]["kd490"]["count"] == 1
    assert result["spatial_support"]["ring_1_3_km"]["kd490"]["count"] == 7


def test_uncertainty_sensitivity_never_reselects_a_lower_uncertainty_pixel():
    ds = data()
    ds.KD490_uncertainty.values[1, 1, 1] = 60
    result = qualification.audit_match(ds, DAY, site=SITE)
    assert result["candidate"]["distance_km"] == 0
    assert result["candidate"]["uncertainty_pct"] == 60
    assert not any(item["fixed_candidate_meets_cap"] for item in result["uncertainty_sensitivity"])
    assert result["uncertainty_sensitivity"][-1]["eligible_pixels_within_radius_meeting_cap"] == 8
    assert result["qualification"]["status"] == "candidate_only"


def test_hash_mismatch_is_rejected_before_dataset_is_opened(tmp_path):
    path = tmp_path / "subset.nc"
    data().to_netcdf(path)
    with pytest.raises(SystemExit, match="not the one recorded"):
        cmems.open_subset(path, {"sha256": "incorrect", "data_sha256": "also-incorrect"})


def test_audit_writes_strict_json_and_report_without_overwriting(tmp_path):
    target = {"role": "synthetic test", **qualification.audit_match(data(), DAY, site=SITE)}
    record = {"created_at": "test", "subset": {"sha256": "test"}, "targets": [target]}
    folder = qualification.write(record, tmp_path / "audit")
    saved = json.loads((folder / "qualification.json").read_text())
    assert saved["targets"][0]["candidate"]["kd490"] == .1
    assert "candidate_only" in (folder / "qualification.md").read_text()
    with pytest.raises(FileExistsError):
        qualification.write(record, folder)
