"""Replay registered experiments in an isolated workspace and compare them with the archived runs.

A replay builds a workspace from the experiment's registration and the archived runs it reads,
checks the inputs it can check against their lock, runs each stage's own command line in a
separate interpreter, and compares every replayed output with the archived one. Archived runs
are never modified.

    kelpobserve replay v1
    kelpobserve replay r4
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from kelpobserve import compare, locks
from kelpobserve.config import ARCHIVE, BENCHMARK, LIBRARY, OPTICS, REPO, Paths
from kelpobserve.relocation import ACOLITE, CACHE, DATA, E5S, LIB, OUT, REG, ROOT, ROOTS, Relocation


@dataclass(frozen=True)
class Stage:
    module: str
    argv: tuple[str, ...]
    output: str
    call: str = "main"
    archived: str | None = None

    @property
    def label(self) -> str:
        name = self.module.removeprefix("benchmark.").removeprefix("coastal_validation.")
        return " ".join((name, *self.argv)) if self.call == "main" else f"{name} {self.call}"


def v1(*argv: str) -> Stage:
    return Stage("benchmark.run_benchmark", argv, argv[0])


def revision(module: str, folder: str, *argv: str) -> Stage:
    return Stage(f"benchmark.{module}", argv, folder)


@dataclass(frozen=True)
class Experiment:
    name: str
    title: str
    registration: Path
    stages: tuple[Stage, ...]
    seed: tuple[str, ...] | None = None
    pointers: dict[str, str] = field(default_factory=dict)
    lock: Path | None = None
    verify_prepared: bool = False
    not_replayed: dict[str, str] = field(default_factory=dict)
    kind: str = "benchmark"
    renames: tuple[tuple[str, str], ...] = ()
    optional: tuple[Stage, ...] = ()
    # The v1 runs used pure-water absorption flat beyond 700 nm (kelpobserve.coastal.optics.water).
    pure_water_to_700_nm: bool = False
    # Configured paths besides the cache that the stages or the lock check read (kelpobserve.config.KEYS).
    needs: tuple[str, ...] = ("acolite",)


CURRENT = BENCHMARK / "registration"
# The EMODnet rasters the registration records under $DATA/external, as this repository carries them.
EMODNET = BENCHMARK / "inputs" / "emodnet"
WRITES_CACHE = "writes into the shared imagery cache; replay it only against a separate copy of the cache"
NETWORK = "reads Sentinel-2 assets from the Planetary Computer, whose holdings can change; name it in --stages to replay it"

EXPERIMENTS = {
    "v1": Experiment(
        name="v1",
        title="Registered benchmark: July leave-one-transect-out fit and the sealed September evaluation (report section 5)",
        registration=CURRENT,
        stages=(v1("fit"), v1("evaluate")),
        seed=("inventory", "prepare"),
        pointers={"inventory": "runs/20260920T115122-1f59c78deecc", "prepare": "runs/20260919T223226-fb8361d31a85"},
        lock=CURRENT / "input-lock.json",
        pure_water_to_700_nm=True,
        needs=("data",),
    ),
    "v1-readiness": Experiment(
        name="v1-readiness",
        title="Registered scene checks, water-column readiness and the cross-sensor comparison (report section 4)",
        registration=CURRENT,
        stages=(v1("scene-qc"), v1("optics-reference"), v1("cross-sensor")),
        seed=("inventory", "prepare", "scenes"),
        pointers={"inventory": "runs/20260919T190733-31c3ba3e0e11", "prepare": "runs/20260919T223226-fb8361d31a85"},
        verify_prepared=True,
        pure_water_to_700_nm=True,
    ),
    "r1": Experiment(
        name="r1",
        title="r1: QC-consistent support, July models with controls and sensitivities, radiometry, injection (report 6.1)",
        registration=CURRENT,
        stages=(
            revision("revision", "r1/qc", "qc"),
            revision("revision", "r1/models", "models"),
            revision("revision", "r1/models-V1", "models", "--variant", "V1"),
            revision("revision", "r1/radiometry", "radiometry"),
            revision("revision", "r1/injection", "injection"),
        ),
        optional=(revision("revision", "r1/screen", "screen"), revision("revision", "r1/windows", "windows")),
        not_replayed={"revision screen": NETWORK, "revision windows": NETWORK, "r1 acolite": WRITES_CACHE},
    ),
    "r3": Experiment(
        name="r3",
        title="r3: glint-tolerant Sentinel-2 mask, models on the selected scenes, injection, windows (report 6.2)",
        registration=CURRENT,
        stages=(
            revision("revision_r3", "r3/mask", "mask"),
            revision("revision_r3", "r3/models", "models"),
            revision("revision_r3", "r3/injection", "injection"),
            revision("revision_r3", "r3/windows", "windows"),
        ),
        not_replayed={"r3 fetch": "downloads from CDSE", "r3 acolite": WRITES_CACHE},
        renames=(("r2-", "r3-"), ("/r2/", "/r3/"), ("registration_r2", "registration_r3")),
    ),
    "r4": Experiment(
        name="r4",
        title="r4: detectability envelope, 100 draws of synthetic cover through the real sampling (report 6.3)",
        registration=CURRENT,
        stages=(revision("revision_r4", "r4/envelope", "run", "--draws", "100"),),
    ),
    "r5": Experiment(
        name="r5",
        title="r5: label-free diagnosis of the unexplained variation (report 6.4)",
        registration=CURRENT,
        stages=(revision("revision_r5", "r5/diagnosis", "run"),),
    ),
    "r6": Experiment(
        name="r6",
        title="r6: background bias, windows, paired Pleiades Neo and Sentinel-2 recovery, offset control (report 6.5)",
        registration=CURRENT,
        stages=(
            revision("revision_r6", "r6/bias", "bias"),
            revision("revision_r6", "r6/windows", "windows"),
            revision("revision_r6", "r6/pneo", "pneo"),
            revision("revision_r6", "r6/s2", "s2"),
            revision("revision_r6", "r6/offset", "offset"),
        ),
        not_replayed={"r6 check": WRITES_CACHE},
    ),
    "e5s": Experiment(
        name="e5s",
        title="E5s: kelp spectral libraries through Sesimbra's water, separability and detection depth (supporting report)",
        registration=OPTICS / "registration",
        stages=(Stage("coastal_validation.e5s", (), "reports/e5s-kelp-spectra", call="run"),),
        lock=OPTICS / "registration" / "e5s-kelp-spectra-lock.json",
        kind="optics",
    ),
}


def optics_relocation(experiment: Experiment, paths: Paths, workspace: Path, archive: bool = False) -> Relocation:
    project = OPTICS if archive else workspace
    rules: dict[str, Path | str] = {
        CACHE: paths.require("optics_cache"),
        REG: project / "registration",
        OUT: project / "out",
        f"{ROOT}/coastal_validation": OPTICS / "coastal_validation",
        ROOT: project,
        LIB: LIBRARY,
        ACOLITE: paths.require("acolite"),
    }
    return Relocation.of(rules)


def relocation(experiment: Experiment, paths: Paths, workspace: Path, archive: bool = False) -> Relocation:
    """Recorded paths mapped into the workspace, or into the archive when checking the lock."""
    out = ARCHIVE if archive else workspace / "out"
    registration = experiment.registration if archive else workspace / "registration"
    rules: dict[str, Path | str] = {
        f"{DATA}/external": EMODNET,
        CACHE: paths.require("cache"),
        OUT: out,
        REG: registration,
        LIB: LIBRARY,
        f"{ROOT}/inputs": BENCHMARK / "inputs",
        f"{ROOT}/benchmark": BENCHMARK / "benchmark",
        f"{ROOT}/docs": BENCHMARK / "docs",
        ROOT: BENCHMARK if archive else workspace,
        E5S: OPTICS / "reports" / "e5s-kelp-spectra" / "results.json",
    }
    for key, token in (("data", DATA), ("acolite", ACOLITE)):
        if key in experiment.needs:
            rules[token] = paths.require(key)
        elif getattr(paths, key):
            rules[token] = getattr(paths, key)
    return Relocation.of(rules)


def seed_runs(experiment: Experiment, workspace: Path, relocate: Relocation) -> None:
    folders = experiment.seed or tuple(sorted(p.name for p in ARCHIVE.iterdir() if p.is_dir()))
    for folder in folders:
        shutil.copytree(ARCHIVE / folder, workspace / "out" / folder)
    for path in (workspace / "out").rglob("*.json"):
        text = path.read_text()
        if any(token in text for token in ROOTS):
            path.write_text(json.dumps(relocate.json(json.loads(text)), indent=2) + "\n")
    for folder, run in experiment.pointers.items():
        pointer = {"schema_version": "2.0", "run_dir": run}
        (workspace / "out" / folder / "latest.json").write_text(json.dumps(pointer, indent=2) + "\n")


def new_optics_workspace(experiment: Experiment, paths: Paths, workspace: Path) -> tuple[Path, Relocation]:
    relocate = optics_relocation(experiment, paths, workspace)
    shutil.copytree(OPTICS / "registration", workspace / "registration")
    shutil.copytree(OPTICS / "reports", workspace / "reports")
    (workspace / ".cache").symlink_to(paths.require("optics_cache"), target_is_directory=True)
    return workspace, relocate


def new_workspace(experiment: Experiment, paths: Paths) -> tuple[Path, Relocation]:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    workspace = paths.workspace / f"{experiment.name}-{stamp}"
    if experiment.kind == "optics":
        return new_optics_workspace(experiment, paths, workspace)
    relocate = relocation(experiment, paths, workspace)
    shutil.copytree(experiment.registration, workspace / "registration")
    seed_runs(experiment, workspace, relocate)
    (workspace / ".cache").symlink_to(paths.require("cache"), target_is_directory=True)
    (workspace / "inputs").symlink_to(BENCHMARK / "inputs", target_is_directory=True)
    (workspace / "benchmark").symlink_to(BENCHMARK / "benchmark", target_is_directory=True)
    return workspace, relocate


def check_prepared(workspace: Path, relocate: Relocation, cache: Path, digests: Path) -> dict:
    """Every Sentinel-2 and Pleiades Neo product the prepare run recorded, hashed where it is now."""
    pointer = json.loads((workspace / "out" / "prepare" / "latest.json").read_text())["run_dir"]
    manifest = json.loads((workspace / "out" / "prepare" / pointer / "manifest.json").read_text())
    table = json.loads(digests.read_text()) if digests.is_file() else {}
    verified, changed, missing = 0, [], []
    for scene in manifest["scenes"]:
        folder = relocate(scene["output_dir"])
        for name, expected in scene.get("outputs", {}).items():
            path = folder / name
            if not path.is_file():
                missing.append(f"{scene['id']}/{name}")
            elif locks.sha256(path, table) == expected:
                verified += 1
            else:
                changed.append(f"{scene['id']}/{name}")
    digests.write_text(json.dumps(table))
    return {"prepared_products": {"verified": verified, "changed": changed, "missing": missing}}


def archived_run(folder: str) -> Path:
    pointer = json.loads((ARCHIVE / folder / "latest.json").read_text())["run_dir"]
    return ARCHIVE / folder / pointer


def uncopied(run: Path) -> dict[str, str]:
    index = json.loads((ARCHIVE.parent / "index.json").read_text())["files"]
    prefix = run.relative_to(ARCHIVE).as_posix() + "/"
    return {f["path"][len(prefix):]: f["sha256"] for f in index if f["path"].startswith(prefix) and not f["copied"]}


def run(name: str, paths: Paths, only: tuple[str, ...] | None = None) -> dict:
    experiment = EXPERIMENTS[name]
    workspace, relocate = new_workspace(experiment, paths)
    digests = paths.workspace / "digests.json"
    report = {"experiment": name, "title": experiment.title, "workspace": str(workspace),
              "python": sys.version.split()[0],
              "not_replayed": experiment.not_replayed, "checks": {}, "stages": {}}
    optics = experiment.kind == "optics"
    if experiment.lock:
        recorded = (optics_relocation if optics else relocation)(experiment, paths, workspace, archive=True)
        checked = locks.check(experiment.lock, recorded, digests)
        report["checks"]["lock"] = checked.summary()
        if not checked.ok:
            report["stopped"] = "the input lock does not verify; nothing was run"
            return finish(workspace, report)
    if experiment.verify_prepared:
        report["checks"].update(check_prepared(workspace, relocate, paths.require("cache"), digests))
        if report["checks"]["prepared_products"]["changed"] or report["checks"]["prepared_products"]["missing"]:
            report["stopped"] = "a prepared image product no longer matches its recorded hash; nothing was run"
            return finish(workspace, report)
    named = lambda s: bool(only) and (s.label in only or s.output in only)
    stages = [s for s in (*experiment.stages, *experiment.optional) if named(s) or (not only and s in experiment.stages)]
    job = {"workspace": str(workspace), "kind": experiment.kind, "relocation": relocate.to_json(),
           "modules": sorted({s.module for s in stages}),
           "lock": str(experiment.lock) if experiment.lock else None,
           "acolite": str(paths.acolite) if paths.acolite else None, "checks": report["checks"],
           "pure_water_to_700_nm": experiment.pure_water_to_700_nm}
    (workspace / "job.json").write_text(json.dumps(job, indent=2) + "\n")
    environment = {
        **os.environ,
        "KELPOBSERVE_REPLAY_JOB": str(workspace / "job.json"),
        "MPLCONFIGDIR": str(workspace / ".matplotlib"),
        "PYTHONPATH": os.pathsep.join([
            str(REPO / "src" / "kelpobserve" / "replay_site"),
            str(OPTICS if optics else BENCHMARK),
            str(REPO / "src"),
        ]),
    }
    # r3 ran under the name r2; registration/revision-r3.json records the rename.
    table = compare.tokens(relocate, workspace, paths) + [list(pair) for pair in experiment.renames]
    for stage in stages:
        log = workspace / "logs" / f"{stage.output.replace('/', '-')}.log"
        log.parent.mkdir(exist_ok=True)
        with log.open("w") as handle:
            done = subprocess.run([sys.executable, "-m", "kelpobserve._driver", stage.module, stage.call, *stage.argv],
                                  cwd=workspace, env=environment, stdout=handle, stderr=subprocess.STDOUT, check=False)
        if done.returncode:
            report["stages"][stage.label] = {"reproduced": False, "failed": str(log)}
            report["stopped"] = f"{stage.label} failed; see {log}"
            return finish(workspace, report)
        if optics:
            archived, replayed, by_hash = OPTICS / (stage.archived or stage.output), workspace / stage.output, {}
        else:
            archived = archived_run(stage.output)
            pointer = json.loads((workspace / "out" / stage.output / "latest.json").read_text())["run_dir"]
            replayed, by_hash = workspace / "out" / stage.output / pointer, uncopied(archived)
        report["stages"][stage.label] = compare.runs(archived, replayed, by_hash, table)
    return finish(workspace, report)


def recompare(workspace: Path, paths: Paths) -> dict:
    """Compare a finished replay again, without rerunning it (after a change to the comparison)."""
    report = json.loads((workspace / "replay.json").read_text())
    experiment = EXPERIMENTS[report["experiment"]]
    relocate = Relocation.from_json(json.loads((workspace / "job.json").read_text())["relocation"])
    table = compare.tokens(relocate, workspace, paths) + [list(pair) for pair in experiment.renames]
    for stage in (*experiment.stages, *experiment.optional):
        if stage.label not in report["stages"] or "failed" in report["stages"][stage.label]:
            continue
        if experiment.kind == "optics":
            archived, replayed, by_hash = OPTICS / (stage.archived or stage.output), workspace / stage.output, {}
        else:
            archived = Path(report["stages"][stage.label]["archived"])
            replayed, by_hash = Path(report["stages"][stage.label]["replayed"]), uncopied(archived)
        report["stages"][stage.label] = compare.runs(archived, replayed, by_hash, table)
    return finish(workspace, report)


def finish(workspace: Path, report: dict) -> dict:
    (workspace / "replay.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
