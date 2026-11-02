"""Which attenuation the water-column correction uses for one scene: in-scene QAA, CMEMS KD490, or none.

The registered v1 rule (reference_qc, approved 20 September 2026) leaves the corrected branch
unavailable when the QAA fit over the deep-water reference raises any plausibility flag; both
Pleiades Neo scenes failed it. This function decides what stands in for them. Fix it, and record
it in benchmark/registration/revision-r7.json, before any corrected-branch
result with CMEMS attenuation is computed.
"""

from __future__ import annotations

from dataclasses import dataclass

from kelpobserve.attenuation.cmems import Match


@dataclass(frozen=True)
class Choice:
    source: str
    reason: str


def choose(qaa: dict, match: Match | None) -> Choice:
    """qaa is one scene of the v1 reference_iops.json (keys kd, wavelengths_nm, qa_flags,
    corrected_branch_available, reference_pixels, solar_zenith_deg). match is the valid CMEMS
    KD490 pixel nearest the reef within 5 km and 3 days, or None. Return Choice("qaa" | "cmems" |
    "none", reason)."""
    raise NotImplementedError("Fix the fallback rule before using CMEMS attenuation.")
