"""Habitat context: classification of the park map and the photoquadrat index."""

import json

import pandas as pd
import pytest

from benchmark import habitat_map as hm
from benchmark import report_maps as maps


def test_every_park_code_prefix_has_a_group():
    assert hm.group_of("A3.31_PT16") == "saccorhiza"
    assert hm.group_of("A3.15_PT12") == "infralittoral_rock"
    assert hm.group_of("A5.23_PT18") == "sand"
    assert hm.group_of("A5.37_PT24*") == "mud"
    assert hm.group_of("A1.2_PT9_PT10") == "littoral_rock"
    with pytest.raises(ValueError):
        hm.group_of("B1.1")


def test_community_names_come_from_the_proposed_unit_or_eunis():
    note = (
        "A3.31_PT16 - Proposed new EUNIS unit : [Saccorhiza polyschides] on infralittoral rock."
        " PMLS = Parque Marinho Luiz Saldanha"
    )
    assert hm.community_name("A3.31_PT16", note) == "Saccorhiza polyschides on infralittoral rock"
    replacement = (
        "A5.13_PT17 - Proposed replacement of EUNIS unit A5.246 : [Ervilia castanea] beds in"
        " infralittoral coarse sand. PMLS = Parque"
    )
    assert hm.community_name("A5.13_PT17", replacement).startswith("Ervilia castanea beds")
    assert hm.community_name("A5.37", None) == "Deep circalittoral mud"
    assert hm.community_name("A5.37*", "") == "Deep circalittoral mud"
    assert hm.community_name("A5.37", float("nan")) == "Deep circalittoral mud"


def test_tracks_are_named_by_month_or_by_recording_date():
    assert hm.track_campaign("June 2026") == "june"
    assert hm.track_campaign("July 2026") == "july"
    assert hm.track_campaign("01/09/26, 15:33:37") == "september"


def test_photoquadrats_are_indexed_by_survey_event(tmp_path):
    source = tmp_path / "published"
    for size in ("full", "thumb"):
        (source / size).mkdir(parents=True)
        (source / size / "16m_T1_Q1.jpg").write_bytes(b"jpeg")
    (source / "photoquadrats.json").write_text(
        json.dumps(
            {
                "campaign": "july_2026",
                "items": [
                    {
                        "depth_bucket_m": 16,
                        "transect_int": 1,
                        "quadrat_int": 1,
                        "source_filename": "Q1_DSC05354.JPG",
                        "full_relpath": "full/16m_T1_Q1.jpg",
                        "thumb_relpath": "thumb/16m_T1_Q1.jpg",
                        "exif_datetime": "2026:07:08 13:27:35",
                    }
                ],
            }
        )
    )
    spec = {"campaigns": {"july": {"strata": {"workbook": {"7": "shallow", "16": "deep"}}}}}
    frame = pd.DataFrame({"event_id": ["july-deep-T1-Q1", "july-shallow-T1-Q1"]})
    site = tmp_path / "site"
    index = maps.photoquadrats(spec, frame, [source], site)
    assert index["june"]["status"] == "pending" and index["september"]["items"] == {}
    record = index["july"]["items"]["july-deep-T1-Q1"]
    assert record["taken"] == "2026:07:08 13:27:35"
    assert (site / record["thumb"]).read_bytes() == b"jpeg"
    assert record["full"] == "photoquadrats/july/full/july-deep-T1-Q1.jpg"


def test_unknown_photoquadrat_is_refused(tmp_path):
    source = tmp_path / "published"
    source.mkdir()
    (source / "photoquadrats.json").write_text(
        json.dumps(
            {
                "campaign": "july_2026",
                "items": [
                    {
                        "depth_bucket_m": 16,
                        "transect_int": 9,
                        "quadrat_int": 1,
                        "source_filename": "Q1_X.JPG",
                        "full_relpath": "full/x.jpg",
                        "thumb_relpath": "thumb/x.jpg",
                        "exif_datetime": None,
                    }
                ],
            }
        )
    )
    spec = {"campaigns": {"july": {"strata": {"workbook": {"16": "deep"}}}}}
    with pytest.raises(ValueError, match="matches no survey event"):
        maps.photoquadrats(
            spec, pd.DataFrame({"event_id": ["july-deep-T1-Q1"]}), [source], tmp_path
        )
