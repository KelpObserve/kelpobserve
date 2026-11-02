"""Compare a replayed run directory with its archived counterpart, file by file.

A file is identical when its bytes are, and equal when its values are once two kinds of field
are set aside: paths, which are compared after mapping the replay's local prefixes to the root
tokens the archived runs record ($OUT, $CACHE, ...), and provenance fields that describe the run
rather than its result (creation times, git state, runner and registration checksums, and the
digests of inputs this repository stores with root tokens). Provenance differences are listed
separately. Archived files too large to copy are compared by the SHA-256 in the archive index.
"""

from __future__ import annotations

import functools
import hashlib
import io
import json
import math
from pathlib import Path

import pandas as pd

from kelpobserve.config import ARCHIVE, BENCHMARK, LIBRARY, OPTICS, Paths
from kelpobserve.relocation import ACOLITE, CACHE, DATA, E5S, LIB, OUT, REG, ROOT, ROOTS, Relocation

PROVENANCE = {
    "created_at", "started_at", "finished_at", "elapsed_s", "duration_s", "lock", "runner_sha256", "git",
    "fit_run", "targets_used", "inventory_run", "run_dir", "registration", "revision_runner_sha256",
    "v1_runner_sha256", "benchmark_git", "library_git", "code", "ran_at", "input_lock",
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def is_provenance(key: str) -> bool:
    return key in PROVENANCE or key.startswith("registration") or key.endswith("_runner_sha256")


@functools.cache
def tokenised_inputs() -> frozenset[str]:
    """SHA-256 of the files of this repository that a stage reads and that hold recorded paths as
    root tokens. A replay records the digest of such a file as it is stored here, which differs
    from the digest the archived run recorded only because of those paths."""
    folders = (BENCHMARK / "registration", BENCHMARK / "inputs", OPTICS / "registration", OPTICS / "reports")
    found = set()
    for path in (p for folder in folders for p in folder.rglob("*") if p.is_file()):
        data = path.read_bytes()
        if any(token.encode() in data for token in ROOTS):
            found.add(digest(data))
    return frozenset(found)


def settle_manifest(files: dict) -> None:
    """A manifest's recorded hash of an output differs whenever that output differs only in its
    recorded paths; such an output compares as equal, and so does the manifest entry naming it."""
    entry = files.get("manifest.json")
    if not entry or entry["status"] != "different":
        return
    kept, moved = [], []
    for difference in entry["differences"]:
        name = difference.split(":")[0].removeprefix("/outputs/")
        if difference.startswith("/outputs/") and files.get(name, {}).get("status") == "equal":
            moved.append(f"{difference.split(':')[0]}: hash of an output that differs only in recorded paths")
        else:
            kept.append(difference)
    if moved:
        entry["provenance_differences"] = entry.get("provenance_differences", []) + moved
        entry["differences"] = kept
        if not kept:
            entry["status"] = "equal"
            del entry["differences"]


def tokens(relocate: Relocation, workspace: Path, paths: Paths) -> list[list[str]]:
    pairs = {
        str(workspace / "out"): OUT, str(ARCHIVE): OUT,
        str(workspace / ".cache"): CACHE,
        str(workspace / "registration"): REG, str(BENCHMARK / "registration"): REG, str(OPTICS / "registration"): REG,
        str(LIBRARY): LIB,
        str(OPTICS / "reports" / "e5s-kelp-spectra" / "results.json"): E5S,
        str(workspace / "inputs" / "emodnet"): f"{DATA}/external",
        str(BENCHMARK / "inputs" / "emodnet"): f"{DATA}/external",
        str(workspace / "inputs"): f"{ROOT}/inputs", str(BENCHMARK / "inputs"): f"{ROOT}/inputs",
        str(workspace): ROOT, str(BENCHMARK): ROOT, str(OPTICS): ROOT,
    }
    for key, token in (("cache", CACHE), ("optics_cache", CACHE), ("data", DATA), ("acolite", ACOLITE)):
        value = getattr(paths, key)
        if value:
            pairs[str(value)] = token
    return sorted(([k, v] for k, v in pairs.items()), key=lambda kv: -len(kv[0]))


def normalise(value: str, table: list[list[str]]) -> str:
    for prefix, token in table:
        if prefix in value:
            value = value.replace(prefix, token)
    return value


def json_differences(left, right, table, path: str = "", out: list | None = None, provenance: list | None = None):
    out = [] if out is None else out
    provenance = [] if provenance is None else provenance
    if isinstance(left, dict) and isinstance(right, dict):
        for key in sorted(set(left) | set(right)):
            target = provenance if is_provenance(key) else out
            if key not in left or key not in right:
                target.append(f"{path}/{key}: only in {'replay' if key in right else 'archive'}")
            elif is_provenance(key):
                inner: list[str] = []
                json_differences(left[key], right[key], table, f"{path}/{key}", inner, inner)
                provenance.extend(inner)
            else:
                json_differences(left[key], right[key], table, f"{path}/{key}", out, provenance)
    elif isinstance(left, list) and isinstance(right, list) and len(left) == len(right):
        for index, (a, b) in enumerate(zip(left, right)):
            json_differences(a, b, table, f"{path}[{index}]", out, provenance)
    elif isinstance(left, str) and isinstance(right, str):
        if normalise(left, table) == normalise(right, table):
            pass
        elif right in tokenised_inputs():
            provenance.append(f"{path}: digest of an input stored with root tokens")
        else:
            out.append(f"{path}: archive {left[:90]} | replay {right[:90]}")
    elif isinstance(left, float) and isinstance(right, float) and math.isnan(left) and math.isnan(right):
        pass
    elif left != right:
        out.append(f"{path}: archive {str(left)[:90]} | replay {str(right)[:90]}")
    return out, provenance


def csv_differences(left: bytes, right: bytes) -> list[str]:
    a, b = pd.read_csv(io.BytesIO(left)), pd.read_csv(io.BytesIO(right))
    if list(a.columns) != list(b.columns) or a.shape != b.shape:
        return [f"shape {a.shape} | {b.shape}; columns only on one side: {sorted(set(a.columns) ^ set(b.columns))}"]
    out = []
    for column in a.columns:
        x, y = a[column], b[column]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            gap = (x - y).abs().max(skipna=True)
            if not ((x.isna() == y.isna()).all() and (pd.isna(gap) or gap == 0)):
                out.append(f"{column}: largest absolute difference {gap}")
        elif not x.fillna("<NA>").astype(str).equals(y.fillna("<NA>").astype(str)):
            out.append(f"{column}: values differ")
    return out


def runs(archived: Path, replayed: Path, by_hash: dict[str, str] | None = None,
         table: list[list[str]] | None = None) -> dict:
    files, table = {}, table or []
    for source in sorted(p for p in archived.rglob("*") if p.is_file()):
        rel = source.relative_to(archived).as_posix()
        target = replayed / rel
        if not target.is_file():
            files[rel] = {"status": "missing from the replay"}
            continue
        left, right = source.read_bytes(), target.read_bytes()
        if digest(left) == digest(right):
            files[rel] = {"status": "identical"}
            continue
        provenance: list[str] = []
        if source.suffix == ".json":
            differences, provenance = json_differences(json.loads(left), json.loads(right), table)
        elif source.suffix == ".csv":
            differences = csv_differences(left, right)
        else:
            differences = ["bytes differ"]
        entry = {"status": "different", "differences": differences[:12]} if differences else {"status": "equal"}
        if provenance:
            entry["provenance_differences"] = provenance[:12]
        files[rel] = entry
    for rel, expected in sorted((by_hash or {}).items()):
        target = replayed / rel
        if not target.is_file():
            files[rel] = {"status": "missing from the replay"}
        else:
            files[rel] = {"status": "identical" if digest(target.read_bytes()) == expected else "different",
                          "compared": "SHA-256 in the archive index"}
    settle_manifest(files)
    unmatched = [rel for rel, entry in files.items() if entry["status"] not in ("identical", "equal")]
    return {"archived": str(archived), "replayed": str(replayed), "reproduced": not unmatched, "files": files}
