"""Attenuation table for the scenes whose in-scene QAA reference failed (technical report, section 4.2).

For each Pleiades Neo scene of the v1 optics-reference run it lists the in-scene QAA result, the
valid CMEMS KD490 pixel nearest the reef, the fallback rule's choice, and the IOPs and Kd at the
scene's bands that the CMEMS value implies with Sesimbra's IOP shape.

    kelpobserve kd490 table
"""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import UTC, date, datetime

from kelpobserve.attenuation import bands, cmems, fallback, qualification
from kelpobserve.config import ARCHIVE, BENCHMARK, OPTICS, REPO, Paths

SEARCH_KM = 5.0
SEARCH_DAYS = 3
SUBSETS = REPO / ".cache" / "cmems-kd490"
EVIDENCE = BENCHMARK / "registration" / "evidence" / "cmems-kd490-subset.json"
REFERENCE = ARCHIVE / "optics-reference" / "runs" / "20260920T115500-a5823dc0e7cc" / "reference_iops.json"
E5S = OPTICS / "reports" / "e5s-kelp-spectra" / "results.json"


def scene_date(scene_id: str) -> date:
    stamp = scene_id.split("-")[1]
    return date(int(stamp[:4]), int(stamp[4:6]), int(stamp[6:8]))


def build(paths: Paths) -> dict:
    evidence = json.loads(EVIDENCE.read_text())
    subset = cmems.open_subset(SUBSETS / evidence["file"], evidence)
    reference = json.loads(REFERENCE.read_text())
    e5s = json.loads(E5S.read_text())

    scenes = []
    for scene_id, qaa in reference["scenes"].items():
        match = cmems.nearest_valid(subset, scene_date(scene_id), SEARCH_KM, SEARCH_DAYS)
        audit = qualification.audit_match(subset, scene_date(scene_id), SEARCH_KM, SEARCH_DAYS)
        try:
            choice = asdict(fallback.choose(qaa, match))
        except NotImplementedError as missing:
            choice = {"source": None, "reason": str(missing)}
        scenes.append({
            "scene": scene_id,
            "qaa": {"wavelengths_nm": qaa["wavelengths_nm"], "kd": qaa["kd"], "qa_flags": qaa["qa_flags"],
                    "corrected_branch_available": qaa["corrected_branch_available"]},
            "cmems": match.to_json() if match else None,
            "cmems_qualification": audit,
            "choice": choice,
            "cmems_at_bands_status": "conditional diagnostic; candidate has not been accepted as reef water",
            "cmems_at_bands": bands.at_bands(match.kd490, qaa["wavelengths_nm"], qaa["solar_zenith_deg"], e5s)
            if match else None,
        })
    return {
        "created_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "subset": evidence,
        "matchup_policy": qualification.POLICY,
        "search": {"max_km": SEARCH_KM, "max_days": SEARCH_DAYS},
        "spectral_shape": "E5s median IOPs (forced-560 QAA, Sentinel-2, 24 December 2025 and 22 February 2026)",
        "scenes": scenes,
    }


def write(record: dict) -> str:
    folder = REPO / "out" / "kd490" / datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    folder.mkdir(parents=True)
    (folder / "attenuation.json").write_text(json.dumps(record, indent=2) + "\n")
    return str(folder / "attenuation.json")
