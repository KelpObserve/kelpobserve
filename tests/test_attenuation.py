from datetime import date

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from kelpobserve.attenuation import bands, cmems

E5S = {"iops": {"median": {"a": [0.10204, 0.08979, 0.09484, 0.4245, 0.624],
                           "bb": [0.0184, 0.01542, 0.01265, 0.00982, 0.00906]}}}


def subset(values: dict[str, np.ndarray]) -> xr.Dataset:
    latitudes = np.array([38.40, 38.41, 38.42])
    longitudes = np.array([-9.21, -9.20, -9.19])
    times = pd.to_datetime(list(values))
    kd = np.stack(list(values.values()))
    return xr.Dataset(
        {"KD490": (("time", "latitude", "longitude"), kd, {"units": "m-1"}),
         "KD490_uncertainty": (("time", "latitude", "longitude"), np.full_like(kd, 30.0), {"units": "%"}),
         "flags": (("time", "latitude", "longitude"), np.zeros_like(kd),
                   {"flag_masks": 1, "flag_meanings": "LAND"})},
        coords={"time": times, "latitude": latitudes, "longitude": longitudes},
    )


def test_nearest_valid_prefers_the_day_itself_then_the_closest_day():
    nan = np.full((3, 3), np.nan)
    near = nan.copy()
    near[1, 1] = 0.09
    far = nan.copy()
    far[0, 0] = 0.20
    data = subset({"2026-07-07": near, "2026-07-08": far, "2026-07-09": near})
    same_day = cmems.nearest_valid(data, date(2026, 7, 8), max_km=5.0, max_days=1)
    assert (same_day.day_offset, same_day.kd490) == (0, 0.20)
    within_metres = cmems.nearest_valid(data, date(2026, 7, 8), max_km=0.5, max_days=1)
    assert (within_metres.day_offset, within_metres.kd490) == (-1, 0.09)
    assert cmems.nearest_valid(data, date(2026, 7, 8), max_km=0.5, max_days=0) is None


def test_the_site_shape_is_scaled_to_the_given_kd490():
    shape = bands.site_shape(E5S)
    a, bb = bands.iops_at([490.0], bands.scale_for(0.0861, shape, 20.83), shape)
    assert bands.kd(a, bb, 20.83)[0] == pytest.approx(0.0861, abs=1e-9)
    record = bands.at_bands(0.0861, [486.0, 563.0, 655.0], 20.83, E5S)
    assert record["kd"][0] == pytest.approx(0.0861, abs=0.002)
    assert record["kd"][2] > record["kd"][1] > record["kd"][0]


def test_kd490_at_or_below_pure_water_leaves_pure_water():
    assert bands.scale_for(0.01, bands.site_shape(E5S), 30.0) == 0.0
