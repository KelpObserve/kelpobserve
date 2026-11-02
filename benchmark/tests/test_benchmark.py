"""Sesimbra habitat benchmark: reference-dataset, leakage and selection invariants."""

import importlib.util
import json
import math
import sys
from datetime import UTC, date, datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "benchmark"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("sesimbra_benchmark", SCRIPTS / "run_benchmark.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)
sys.path.remove(str(SCRIPTS))

COLUMNS = [
    "Date",
    "Depth",
    "Transect",
    "Quadrat",
    "Species",
    module.JUVENILE,
    module.ADULT,
    "Size",
    "Erosion",
    "Epibiont",
    module.COVER,
    "Notes",
]


def make_spec(**overrides):
    spec = {
        "experiment_id": "sesimbra-habitat-v1",
        "year": 2026,
        "taxon_mapping": {"target": "saccorhiza polyschides", "absence": "EMPTY"},
        "campaigns": {
            "july": {
                "month": 7,
                "role": "development",
                "plausible_survey_dates": ["2026-07-06", "2026-07-07"],
                "strata": {
                    "workbook": {"7": "shallow", "16": "deep"},
                    "transects": {"7m": "shallow", "17m": "deep"},
                },
                "pneo": None,
            },
            "september": {
                "month": 9,
                "role": "evaluation",
                "plausible_survey_dates": ["2026-09-01", "2026-09-02"],
                "strata": {
                    "workbook": {"6": "shallow", "15": "deep"},
                    "transects": {"7m": "shallow", "17m": "deep"},
                },
                "pneo": None,
            },
        },
        "ccmar_confirmation": {"confirmed_campaigns": []},
    }
    spec.update(overrides)
    return spec


def workbook(*records):
    rows = pd.DataFrame(
        [
            {
                "Date": pd.Timestamp(r[0]),
                "Depth": r[1],
                "Transect": r[2],
                "Quadrat": r[3],
                "Species": r[4],
                module.JUVENILE: r[5] if len(r) > 5 else None,
                module.ADULT: r[6] if len(r) > 6 else None,
                "Size": 10.0,
                "Erosion": None,
                "Epibiont": None,
                module.COVER: r[7] if len(r) > 7 else None,
                "Notes": None,
            }
            for r in records
        ],
        columns=COLUMNS,
    )
    rows.insert(0, "source_row", rows.index + 2)
    return rows


SAC = "saccorhiza polyschides"


def test_plant_rows_collapse_to_one_event_without_multiplying_cover():
    rows = workbook(
        ("2026-07-06", 7, 1, 1, SAC, 1, None, 60.0),
        ("2026-07-06", 7, 1, 1, SAC, 1),
        ("2026-07-06", 7, 1, 1, SAC, None, 1),
    )
    events = module.aggregate_events(rows, make_spec())
    assert len(events) == 1
    event = events.iloc[0]
    assert event.cover_pct == 60.0 and event.plant_rows == 3
    assert event.cover_status == "measured" and json.loads(event.source_rows) == [2, 3, 4]


def test_explicit_absence_survives_blank_cover_and_positive_missing_stays_missing():
    rows = workbook(
        ("2026-09-02", 6, 1, 1, "EMPTY"),
        ("2026-09-02", 15, 2, 3, SAC),
        ("2026-07-07", 16, 1, 1, "EMPTY", None, None, 0.0),
    )
    events = module.aggregate_events(rows, make_spec()).set_index(["campaign", "stratum"])
    derived = events.loc[("september", "shallow")]
    assert derived.cover_pct == 0.0 and derived.cover_status == "absence_derived_zero"
    assert not derived.occurrence
    positive = events.loc[("september", "deep")]
    assert math.isnan(positive.cover_pct) and positive.cover_status == "missing_positive"
    assert positive.occurrence
    assert events.loc[("july", "deep")].cover_status == "absence_recorded_zero"


@pytest.mark.parametrize(
    "records, message",
    [
        (
            [
                ("2026-07-06", 7, 1, 1, SAC, 1, None, 60.0),
                ("2026-07-06", 7, 1, 1, SAC, 1, None, 50.0),
            ],
            "Repeated cover",
        ),
        ([("2026-07-06", 7, 1, 1, SAC), ("2026-07-06", 7, 1, 1, "EMPTY")], "both EMPTY"),
        ([("2026-07-06", 7, 1, 1, "EMPTY", None, None, 5.0)], "nonzero cover"),
        ([("2026-07-06", 7, 1, 1, "laminaria ochroleuca")], "Unmapped taxa"),
        ([("2026-07-06", 9, 1, 1, SAC)], "Unmapped workbook depth"),
    ],
)
def test_ambiguous_quadrat_records_raise(records, message):
    with pytest.raises(ValueError, match=message):
        module.aggregate_events(workbook(*records), make_spec())


def test_transect_dates_use_portuguese_months():
    assert module.parse_transect_date("01/set", 2026) == date(2026, 9, 1)
    assert module.parse_transect_date("07/jul", 2026) == date(2026, 7, 7)
    with pytest.raises(ValueError, match="Unrecognised"):
        module.parse_transect_date("01/sep", 2026)


def points(*keys):
    import geopandas as gpd
    from shapely.geometry import Point

    return gpd.GeoDataFrame(
        [
            {
                "campaign": c,
                "stratum": s,
                "transect": t,
                "quadrat": q,
                "transects_date": "2026-07-07",
                "transects_depth_label": "7m",
                "transects_fid": i,
                "geometry": Point(-9.2 + i * 1e-4, 38.4),
            }
            for i, (c, s, t, q) in enumerate(keys)
        ],
        geometry="geometry",
        crs="EPSG:4326",
    )


def test_duplicate_and_unmatched_keys_are_exposed():
    events = module.aggregate_events(
        workbook(
            ("2026-07-06", 7, 1, 1, SAC, 1, None, 60.0), ("2026-07-06", 7, 1, 2, SAC, 1, None, 5.0)
        ),
        make_spec(),
    )
    joined = module.join_events(
        events,
        points(("july", "shallow", 1, 1), ("july", "shallow", 1, 1), ("july", "shallow", 1, 3)),
    ).set_index("event_id")
    assert set(joined.loc["july-shallow-T1-Q1", "join_issues"].iloc[0]) == {
        "duplicate_transect_key"
    }
    assert joined.loc["july-shallow-T1-Q2", "join_issues"] == ["no_transect_point"]
    assert joined.loc["july-shallow-T1-Q3", "join_issues"] == ["no_workbook_event"]


def test_unregistered_survey_dates_raise_and_discrepant_registered_dates_survive():
    spec = make_spec()
    joined = pd.DataFrame(
        [
            {
                "event_id": "e",
                "campaign": "september",
                "workbook_date": "2026-09-02",
                "transects_date": "2026-09-01",
            }
        ]
    )
    module.check_dates(joined, spec)
    joined.loc[0, "transects_date"] = "2026-08-30"
    with pytest.raises(ValueError, match="unregistered survey date"):
        module.check_dates(joined, spec)


def test_depth_qc_flags_displaced_event_and_passes_consistent_one():
    levels = (0.5, 3.5)
    status, low, high, neighbour = module.depth_consistency(15, 13.1, 11.0, 15.0, levels, 3.0)
    assert status == "pass" and low <= 3.0 and high >= -3.0
    status, low, _, neighbour = module.depth_consistency(15, 5.2, 5.2, 12.6, levels, 3.0)
    assert status == "flag" and low > 3.0 and neighbour is True
    assert module.depth_consistency(15, 5.2, 5.2, 6.0, levels, 3.0)[3] is False
    assert module.depth_consistency(15, 13.1, 11.0, 15.0, None, 3.0)[0] == "pending_water_level"
    assert module.depth_consistency(6, math.nan, math.nan, math.nan, levels, 3.0)[0] == "no_lidar"
    assert module.depth_consistency(math.nan, 5.0, 5.0, 5.0, levels, 3.0)[0] == "no_label"


def test_join_classes_follow_join_issues_and_depth_qc():
    spec = make_spec()
    assert module.classify_join("july", [], "pass", spec) == ("reconciled_by_rule", [])
    assert module.classify_join("july", [], "flag", spec)[0] == "provisional"
    assert module.classify_join("july", [], "no_lidar", spec)[0] == "provisional"
    assert module.classify_join("july", [], "pending_water_level", spec)[0] == "pending"
    assert module.classify_join("july", ["no_transect_point"], "pass", spec)[0] == "provisional"
    spec["ccmar_confirmation"]["confirmed_campaigns"] = ["july"]
    assert module.classify_join("july", [], "pass", spec) == ("confirmed", [])
    assert module.classify_join("july", [], "flag", spec)[0] == "provisional"


def test_water_levels_need_lat_and_every_plausible_date():
    spec = make_spec()
    levels = pd.DataFrame(
        {
            "campaign": ["september"] * 2,
            "date": ["2026-09-01", "2026-09-02"],
            "min_m": [0.6, 0.4],
            "max_m": [3.2, 3.4],
            "datum": ["LAT", "LAT"],
        }
    )
    assert module.water_level_interval("september", spec, levels) == (0.4, 3.4)
    assert module.water_level_interval("september", spec, levels.iloc[:1]) is None
    with pytest.raises(ValueError, match="LAT"):
        module.water_level_interval("september", spec, levels.assign(datum="MSL"))


def test_survey_window_must_hold_for_every_plausible_date():
    assert module.survey_window(["2026-07-06", "2026-07-07", "2026-07-08"], 4) == (
        date(2026, 7, 4),
        date(2026, 7, 10),
    )
    with pytest.raises(ValueError, match="No acquisition date"):
        module.survey_window(["2026-07-01", "2026-07-20"], 4)


def scene(name, when):
    return SimpleNamespace(name=name, sensing_datetime=when, tile="T29SMC")


def test_scene_ranking_is_deterministic_under_date_uncertainty():
    pneo = datetime(2026, 9, 3, 11, 34, tzinfo=UTC)
    candidates = [
        scene("S2A_0831", datetime(2026, 8, 31, 11, 33)),
        scene("S2C_0905", datetime(2026, 9, 5, 11, 21)),
        scene("S2B_0903", datetime(2026, 9, 3, 11, 33)),
        scene("S2B_0831", datetime(2026, 8, 31, 11, 21)),
    ]
    ranked = module.rank_candidates(candidates, ["2026-09-01", "2026-09-02"], pneo)
    assert [r["candidate"].name for r in ranked] == ["S2B_0903", "S2A_0831", "S2B_0831", "S2C_0905"]
    no_pneo = module.rank_candidates(candidates[::-1], ["2026-09-01", "2026-09-02"], None)
    assert [r["candidate"].name for r in no_pneo][:3] == ["S2A_0831", "S2B_0831", "S2B_0903"]


def test_newest_processing_baseline_wins_per_datatake():
    when = datetime(2026, 6, 2, 11, 21)
    kept = module.latest_baseline(
        [
            scene("S2B_MSIL1C_20260602T112109_N0511_R037_T29SMC_X", when),
            scene("S2B_MSIL1C_20260602T112109_N0512_R037_T29SMC_Y", when),
        ]
    )
    assert [c.name for c in kept] == ["S2B_MSIL1C_20260602T112109_N0512_R037_T29SMC_Y"]


def test_reference_mask_applies_depth_distance_and_footprint_rules():
    shape = (60, 60)
    depth = np.full(shape, 80.0)
    depth[:, :10] = np.nan  # land on the west
    depth[:, 50:] = 30.0  # too shallow in the east
    valid = np.zeros(shape, bool)
    valid[5:55, 5:58] = True
    habitat = np.zeros(shape, bool)
    habitat[0:3, 20:23] = True
    rules = {
        "min_dtm_depth_m": 50.0,
        "min_distance_from_land_m": 300.0,
        "min_distance_from_habitat_m": 300.0,
        "min_inside_footprint_m": 90.0,
    }
    mask = module.reference_mask(valid, depth, habitat, 30.0, rules)
    rows, cols = np.nonzero(mask)
    # Distances are centre to centre: column 19 is exactly 10 pixels (300 m) from land.
    assert cols.min() == 19 and cols.max() == 49  # land rule west, depth rule east
    assert rows.min() == 7 and rows.max() == 52  # three pixels (90 m) inside the footprint
    assert not mask[:12, 20:23].any() and mask[12, 21]  # 300 m from habitat


def test_workbook_revisions_expose_recoded_counts_without_cover_changes():
    spec = make_spec()
    previous = workbook(("2026-07-06", 7, 1, 1, SAC, 2, None, 60.0), ("2026-07-06", 7, 1, 1, SAC))
    current = workbook(
        ("2026-07-06", 7, 1, 1, SAC, 1, None, 60.0),
        ("2026-07-06", 7, 1, 1, SAC, 1),
        ("2026-09-02", 6, 1, 1, "EMPTY"),
    )
    report = module.workbook_revisions(previous, current, spec)
    assert report["cells_differing_by_column"] == {module.JUVENILE: 2}
    assert report["events_with_changed_cover_pct"] == []
    assert report["events_with_changed_juvenile_total"] == []  # 2 == 1 + 1


@pytest.fixture
def synthetic_sources(tmp_path):
    import geopandas as gpd
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import Point

    rows = workbook(
        ("2026-07-06", 7, 1, 1, SAC, 1, None, 60.0),
        ("2026-07-06", 7, 1, 1, SAC, 1),
        ("2026-07-07", 16, 1, 1, "EMPTY", None, None, 0.0),
        ("2026-09-02", 6, 1, 1, "EMPTY"),
        ("2026-09-02", 15, 1, 1, SAC),
    ).drop(columns="source_row")
    kelp = tmp_path / "kelp.xlsx"
    rows.to_excel(kelp, sheet_name="Folha1", index=False)
    rows[rows.Date.dt.month == 7].to_excel(
        tmp_path / "previous.xlsx", sheet_name="Folha1", index=False
    )
    coords = [(-9.2020, 38.4096), (-9.2020, 38.4086), (-9.2010, 38.4096), (-9.2010, 38.4086)]
    labels = [("07/jul", "7m"), ("07/jul", "17m"), ("01/set", "7m"), ("01/set", "17m")]
    gpd.GeoDataFrame(
        [
            {
                "Date": d,
                "Depth": z,
                "Transect number": "T1",
                "Quadrat": "Q1",
                "GPS Lat": y,
                "GPS Long": x,
                "geometry": Point(x, y),
            }
            for (d, z), (x, y) in zip(labels, coords)
        ],
        geometry="geometry",
        crs="EPSG:4326",
    ).to_file(tmp_path / "transects.gpkg", layer="transects", driver="GPKG")
    previous = tmp_path / "previous.gpkg"
    common = {
        "Transect number": "T1",
        "GPS Lat": coords[0][1],
        "GPS Long": coords[0][0],
        "geometry": Point(coords[0]),
    }
    gpd.GeoDataFrame(
        [
            {
                "Date": "07/jul",
                "Depth": "7m",
                "Quadrat": "#20",
                "Label": "7m#20",
                "Time": "13:32:00",
                **common,
            }
        ],
        crs="EPSG:4326",
    ).to_file(previous, layer="timed", driver="GPKG")
    gpd.GeoDataFrame(
        [{"Date": "07/jul", "Depth": "7m", "Quadrat": "Q1", "Label": "JulT1", **common}],
        crs="EPSG:4326",
    ).to_file(previous, layer="numbered", driver="GPKG")
    lidar = tmp_path / "lidar.tif"
    with rasterio.open(
        lidar,
        "w",
        driver="GTiff",
        height=20,
        width=20,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=from_origin(-9.205, 38.412, 0.0005, 0.0005),
        nodata=float("nan"),
    ) as dst:
        dst.write(np.full((1, 20, 20), 5.0, dtype="float32"))
    spec = make_spec(
        sources={
            "kelp_workbook": str(kelp),
            "kelp_sheet": "Folha1",
            "transects": {"path": str(tmp_path / "transects.gpkg"), "layer": "transects"},
            "previous_workbook": str(tmp_path / "previous.xlsx"),
            "previous_transects": {
                "path": str(previous),
                "timed_layer": "timed",
                "numbered_layer": "numbered",
            },
            "lidar": str(lidar),
        },
        expected_events=4,
        quadrat_order={"assumption": "survey order", "status": "assumed_pending_ccmar"},
        water_levels={"path": None, "status": "missing"},
        depth_consistency={"tolerance_m": 3.0, "label_source": "workbook"},
    )
    path = tmp_path / "experiment.json"
    path.write_text(json.dumps(spec))
    return path


def test_inventory_keeps_targets_out_of_the_event_table(synthetic_sources, tmp_path):
    directory = module.inventory(synthetic_sources, tmp_path / "out")
    events = pd.read_csv(directory / "events.csv")
    assert not set(module.TARGET_FIELDS) & set(events.columns)
    development = pd.read_csv(directory / "targets_development.csv")
    september = pd.read_csv(directory / "targets_september.csv")
    assert set(development.campaign) == {"july"} and set(september.campaign) == {"september"}
    assert dict(zip(september.event_id, september.cover_status)) == {
        "september-shallow-T1-Q1": "absence_derived_zero",
        "september-deep-T1-Q1": "missing_positive",
    }
    assert set(events.depth_qc_status) == {"pending_water_level"}
    assert set(events.join_class) == {"pending"} and not events.eligible_registered.any()
    manifest = json.loads((directory / "manifest.json").read_text())
    assert manifest["complete"] and set(manifest["outputs"]) >= {
        "targets_september.csv",
        "events.gpkg",
    }
    pointer = json.loads((tmp_path / "out/inventory/latest.json").read_text())
    assert (tmp_path / "out/inventory" / pointer["run_dir"]).resolve() == directory


def test_inventory_withholds_completion_when_events_are_missing(synthetic_sources, tmp_path):
    spec = json.loads(synthetic_sources.read_text())
    spec["expected_events"] = 5
    synthetic_sources.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="Expected 5 events"):
        module.inventory(synthetic_sources, tmp_path / "out")
    assert not (tmp_path / "out/inventory/latest.json").exists()


REGISTRATION = Path(__file__).resolve().parents[1] / "registration"


def test_registered_acolite_settings_exclude_pansharpening_and_inert_keys():
    common = module.read_settings_file(REGISTRATION / "acolite/common.txt")
    pneo = module.read_settings_file(REGISTRATION / "acolite/pleiades_neo.txt")
    s2 = module.read_settings_file(REGISTRATION / "acolite/sentinel2.txt")
    assert common["pans"] == "False" and pneo["pleiades_skip_pan"] == "True"
    assert "glint_correction" not in common | pneo | s2
    assert common["dsf_residual_glint_correction"] == "False"
    assert s2["dsf_aot_estimate"] == "tiled" and pneo["dsf_aot_estimate"] == "fixed"
    assert not {"inputfile", "output", "limit"} & set(common | pneo | s2)


def test_settings_parser_rejects_malformed_lines(tmp_path):
    path = tmp_path / "settings.txt"
    path.write_text("# comment\nlimit=1,2,3,4  # trailing\n\n")
    assert module.read_settings_file(path) == {"limit": "1,2,3,4"}
    path.write_text("dsf_aot_estimate fixed\n")
    with pytest.raises(ValueError, match="Malformed"):
        module.read_settings_file(path)


def test_sentinel2_scene_ids_carry_date_platform_and_orbit():
    name = "S2B_MSIL1C_20260903T113319_N0512_R080_T29SMC_20260903T153653.SAFE.zip"
    assert module.s2_scene_id(name) == "s2-20260903-S2B-R080"
    with pytest.raises(ValueError, match="Unrecognised"):
        module.s2_scene_id("LC09_L1TP_204033_20260903.tar")


def test_acolite_limit_is_south_west_north_east():
    south, west, north, east = module.geographic_limit(
        (480000, 4249000, 486000, 4253000), "EPSG:32629"
    )
    assert south < north and west < east
    assert 38.3 < south < 38.5 and -9.3 < west < -9.1


def test_other_acolite_revisions_are_refused(monkeypatch):
    config = {"path": "/nonexistent", "commit": "abc"}
    monkeypatch.setattr(
        module,
        "acolite_revision",
        lambda path: {"commit": "abc", "describe": "", "uncommitted_changes": True},
    )
    with pytest.raises(ValueError, match="differs from the registered revision"):
        module.verify_acolite(config)
    monkeypatch.setattr(
        module,
        "acolite_revision",
        lambda path: {"commit": "abc", "describe": "", "uncommitted_changes": False},
    )
    assert module.verify_acolite(config)["commit"] == "abc"


CHECKS = {
    "valid_water_min_fraction": 0.8,
    "quadrat_cell_min_fraction": 0.8,
    "sentinel2_max_median_rhos_1614": 0.01,
}


def test_qc_metrics_count_valid_water_and_distinct_quadrat_cells():
    shape = (10, 10)
    blue, green, red = (np.full(shape, 0.02) for _ in range(3))
    red[0, :] = np.nan  # incomplete band on one row
    mask_rhot = np.full(shape, 0.01)
    mask_rhot[1, :] = 0.05  # cloud or glint above the ACOLITE threshold
    water = np.ones(shape, bool)
    water[:, 9] = False  # land column
    glint = np.full(shape, 0.004)
    glint[1, :] = 0.2
    metrics = module.qc_metrics(
        [blue, green, red],
        mask_rhot,
        0.0215,
        water,
        [(0, 0), (5, 5), (5, 5), (6, 6), (1, 3)],
        glint=glint,
        nir=np.full(shape, 0.003),
    )
    assert metrics["water_pixels"] == 90 and metrics["valid_water_pixels"] == 72
    assert metrics["valid_water_fraction"] == pytest.approx(0.8)
    assert metrics["quadrat_cells"] == 4 and metrics["valid_quadrat_cells"] == 2
    assert metrics["median_rhos_glint_band"] == pytest.approx(0.004)  # before the cloud mask
    assert metrics["median_rhos_nir_quadrat_cells"] == pytest.approx(0.003)


def test_qc_verdict_applies_thresholds_and_glint_only_to_sentinel2():
    good = {
        "valid_water_fraction": 0.9,
        "quadrat_cell_fraction": 1.0,
        "median_rhos_glint_band": 0.004,
    }
    assert module.qc_verdict(good, "sentinel2", CHECKS) == (True, [])
    glinted = good | {"median_rhos_glint_band": 0.02}
    assert module.qc_verdict(glinted, "sentinel2", CHECKS) == (False, ["glint_above_maximum"])
    assert module.qc_verdict(glinted, "pleiades_neo", CHECKS) == (True, [])
    sparse = good | {"quadrat_cell_fraction": 0.5, "valid_water_fraction": 0.7}
    assert module.qc_verdict(sparse, "sentinel2", CHECKS)[1] == [
        "valid_water_fraction_below_minimum",
        "quadrat_cell_fraction_below_minimum",
    ]
    assert module.qc_verdict({**good, "median_rhos_glint_band": None}, "sentinel2", CHECKS)[1] == [
        "glint_band_missing"
    ]


def candidate(name, offset, pneo, processed=True, passed=True, water=0.9):
    return {
        "name": name,
        "max_offset_days": offset,
        "pneo_offset_s": pneo,
        "processed": processed,
        "passed": passed,
        "valid_water_fraction": water,
    }


def test_scene_selection_breaks_ties_by_valid_water_and_waits_for_unprocessed_scenes():
    tie = [candidate("S2A_a", 0, None, water=0.85), candidate("S2B_b", 0, None, water=0.95)]
    assert module.select_scene(tie) == {"status": "selected", "name": "S2B_b"}
    july = [
        candidate("jul09", 3, 85000, passed=False),
        candidate("jul05", 3, 259000, passed=False),
        candidate("jul10", 4, 172000, processed=False),
    ]
    assert module.select_scene(july) == {
        "status": "undetermined",
        "reason": "not processed",
        "waiting_for": ["jul10"],
    }
    ahead_unprocessed = [candidate("sep03", 2, 60), candidate("aug31", 2, 259000, processed=False)]
    assert module.select_scene(ahead_unprocessed)["name"] == "sep03"
    assert module.select_scene([candidate("x", 1, None, passed=False)])["status"] == "unavailable"


def test_nearest_band_respects_tolerance():
    bands = {490.0: Path("b490"), 1614.0: Path("b1614")}
    assert module.nearest_band(bands, 489)[0] == 490.0
    assert module.nearest_band(bands, 560) == (None, None)
    assert module.nearest_band(bands, 1600, tolerance_nm=None)[0] == 1614.0


def test_amended_window_overrides_the_registered_one():
    spec = {"sentinel2": {"window_days": 4}}
    assert module.window_days(spec) == (4, False)
    spec["sentinel2"]["window_days_amended"] = 7
    assert module.window_days(spec) == (7, True)


def test_screen_matches_candidates_by_datatake_minute_and_platform():
    when = datetime(2026, 7, 2, 11, 21, 19, tzinfo=UTC)
    name = "S2B_MSIL1C_20260702T112119_N0512_R037_T29SMC_20260702T140000.SAFE"
    assert module.platform_of(name) == "Sentinel-2B"
    screen = {module.screen_key(when, "Sentinel-2B"): {"b11_median": 0.0122, "screened_out": False}}
    candidate = SimpleNamespace(name=name, sensing_datetime=when, tile="T29SMC")
    assert module.screen_fields(candidate, screen, 0.02) == {
        "screen_b11_median": 0.0122,
        "screened_out": False,
    }
    glinted = {module.screen_key(when, "Sentinel-2B"): {"b11_median": 0.21, "screened_out": True}}
    assert module.screen_fields(candidate, glinted, 0.02)["screened_out"] is True
    # An unscreened candidate is never silently treated as screened out.
    assert module.screen_fields(candidate, {}, 0.02) == {
        "screen_b11_median": None,
        "screened_out": False,
    }


def test_screened_out_candidates_are_determined_failures_not_pending():
    ranked = [
        candidate("glinted", 3, 1000, processed=True, passed=False),
        candidate("clean", 6, 5000, processed=True, passed=True, water=0.9),
    ]
    assert module.select_scene(ranked) == {"status": "selected", "name": "clean"}


def test_sensitivity_reports_alternatives_without_changing_the_verdict():
    valid_swir = np.ones((4, 4), bool)
    valid_swir[0, 0] = False  # one masked cell, which happens to hold three quadrats
    valid_nir = np.ones((4, 4), bool)
    cells = {(0, 0): 3, (1, 1): 1, (2, 2): 1, (3, 3): 1}
    out = module.qc_sensitivity(valid_swir, valid_nir, cells)
    assert out["registered_cell_fraction"] == 0.75
    assert out["would_pass_at_relaxed_cell_fraction"] is True  # at 0.75
    assert out["per_quadrat_fraction"] == 0.5 and out["would_pass_per_quadrat"] is False
    assert out["nir_mask_cell_fraction"] == 1.0 and out["would_pass_with_nir_mask"] is True
    assert "Labelled sensitivity" in out["note"]


def test_sentinel2_sensing_time_comes_from_the_product_name():
    when = module.s2_sensing_time(
        "S2B_MSIL1C_20260831T112109_N0512_R037_T29SMC_20260831T134650.SAFE.zip"
    )
    assert when == datetime(2026, 8, 31, 11, 21, 9, tzinfo=UTC)
    with pytest.raises(ValueError, match="No sensing time"):
        module.s2_sensing_time("PNEO3_STD_ORT.zip")


def test_fit_stage_cannot_reach_the_holdout_targets(tmp_path):
    assert (
        module.targets_path(tmp_path, "development", stage="fit").name == "targets_development.csv"
    )
    assert (
        module.targets_path(tmp_path, "holdout", stage="evaluate").name == "targets_september.csv"
    )
    with pytest.raises(ValueError, match="cannot open the held-out targets"):
        module.targets_path(tmp_path, "holdout", stage="fit")


def test_unit_targets_average_members_and_never_impute_missing_cover():
    units = pd.DataFrame(
        {
            "unit_id": ["a", "b", "c"],
            "event_ids": [json.dumps(["e1", "e2"]), json.dumps(["e3"]), json.dumps(["e4"])],
        }
    )
    targets = pd.DataFrame(
        {
            "event_id": ["e1", "e2", "e3", "e4"],
            "cover_pct": [60.0, 40.0, 0.0, float("nan")],
            "cover_status": ["measured", "measured", "absence_derived_zero", "missing_positive"],
        }
    )
    out = module.attach_targets(units, targets)
    assert out.cover_pct.tolist()[:2] == [50.0, 0.0]
    assert math.isnan(out.cover_pct.iloc[2]) and out.n_scored_quadrats.tolist() == [2, 1, 0]


def test_ridge_standardises_on_training_data_only_and_scoring_reports_clipping():
    train = pd.DataFrame({"depth": [5.0, 10.0, 15.0, 20.0], "cover_pct": [80.0, 60.0, 40.0, 20.0]})
    spec = {"models": {"ridge": {"alpha": 1.0}}}
    fitted = module.fit_family({"kind": "ridge", "predictors": ["depth"]}, train, spec)
    assert fitted["centre"] == [12.5] and fitted["n_train"] == 4
    predicted = module.predict_family(fitted, pd.DataFrame({"depth": [5.0, 20.0]}))
    assert predicted[0] > predicted[1]  # shallower predicts more cover
    scored = module.score(np.array([-10.0, 50.0]), np.array([0.0, 50.0]))
    assert scored["bias_pp"] == -5.0 and scored["clipping_frequency"] == 0.5
    assert scored["mae_clipped_pp"] == 0.0 and scored["mae_pp"] == 5.0
