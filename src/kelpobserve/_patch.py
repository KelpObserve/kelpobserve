"""Rebind an experiment's anchors for a replay, in every process of the replay.

The experiment's modules are used unchanged; only their anchors are rebound, before any stage
module is imported, so that module constants computed from them point into the workspace.

Benchmark experiments: the runner's ROOT, REGISTRATION and OUTPUT; load_spec relocates the
recorded paths of the registration; scene_directories relocates the recorded ACOLITE output
folders; verify_lock, when the experiment has a lock, returns the check already made against
this checkout. For the v1 experiments, the pure-water absorption is the tabulation the v1 runs
used, flat beyond 700 nm.

Optics experiments (the coastal_validation modules): common.ROOT, from which every module
takes its folders; the stage module's verify_lock returns the check already made; and
e5s.ACOLITE_RSR follows the configured ACOLITE checkout.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path

from kelpobserve.locks import sha256
from kelpobserve.relocation import Relocation

_INSTALLED: set[str] = set()


def installed(job_path: str) -> bool:
    return job_path in _INSTALLED


def lock_result(job: dict) -> dict:
    lock = Path(job["lock"])
    frozen = json.loads(lock.read_text())
    return {"path": str(lock), "sha256": sha256(lock), "frozen_at": frozen.get("frozen_at"),
            "n_files": len(frozen["files"]), "replay_check": job["checks"].get("lock")}


def benchmark(job: dict, relocate: Relocation, workspace: Path) -> None:
    from benchmark import run_benchmark as rb

    recorded_scenes = rb.scene_directories

    def load_spec(path: Path | None = None) -> dict:
        return relocate.json(json.loads(Path(path or rb.REGISTRATION / "experiment.json").read_text()))

    def scene_directories(spec: dict) -> dict[str, dict]:
        return {key: {**record, "output_dir": str(relocate(record["output_dir"]))}
                for key, record in recorded_scenes(spec).items()}

    rb.ROOT, rb.REGISTRATION, rb.OUTPUT = workspace, workspace / "registration", workspace / "out"
    rb.load_spec = load_spec
    rb.scene_directories = scene_directories
    if job.get("lock"):
        rb.verify_lock = lambda path: lock_result(job)
    if job.get("pure_water_to_700_nm"):
        from kelpobserve.coastal.optics import water

        water._AW_TABLE = water._AW_TABLE_TO_700_NM


def optics(job: dict, workspace: Path) -> None:
    from coastal_validation import common

    common.ROOT = workspace
    for name in job["modules"]:
        module = importlib.import_module(name)
        module.verify_lock = lambda: lock_result(job)
        if hasattr(module, "ACOLITE_RSR") and job.get("acolite"):
            module.ACOLITE_RSR = Path(job["acolite"]) / "data" / "RSR"


def install(job_path: str) -> None:
    if job_path in _INSTALLED:
        return
    job = json.loads(Path(job_path).read_text())
    workspace = Path(job["workspace"])
    if job.get("kind") == "optics":
        optics(job, workspace)
    else:
        benchmark(job, Relocation.from_json(job["relocation"]), workspace)
    _INSTALLED.add(job_path)
