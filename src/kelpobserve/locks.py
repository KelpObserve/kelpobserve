"""Check an input lock against the files of this checkout.

Each recorded path is relocated and hashed. A file counts as verified when its SHA-256 equals
the recorded one; any other file is reported as changed or missing.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

from kelpobserve.relocation import Relocation


def sha256(path: Path, cache: dict | None = None) -> str:
    stat = path.stat()
    key = str(path.resolve())
    if cache is not None and cache.get(key, [None, None])[:2] == [stat.st_size, stat.st_mtime_ns]:
        return cache[key][2]
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    if cache is not None:
        cache[key] = [stat.st_size, stat.st_mtime_ns, digest.hexdigest()]
    return digest.hexdigest()


@dataclass
class LockCheck:
    lock: Path
    verified: list[str] = field(default_factory=list)
    changed: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.changed and not self.missing

    def summary(self) -> dict:
        return {
            "lock": str(self.lock),
            "verified": len(self.verified),
            "changed": self.changed,
            "missing": self.missing,
        }


def check(lock_path: Path, relocation: Relocation, cache_path: Path | None = None) -> LockCheck:
    lock = json.loads(lock_path.read_text())
    cache = json.loads(cache_path.read_text()) if cache_path and cache_path.is_file() else {}
    result = LockCheck(lock_path)
    for recorded, expected in lock["files"].items():
        local = relocation(recorded)
        if not local.is_file():
            result.missing.append(recorded)
        elif sha256(local, cache) == expected:
            result.verified.append(recorded)
        else:
            result.changed.append(recorded)
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache))
    return result
