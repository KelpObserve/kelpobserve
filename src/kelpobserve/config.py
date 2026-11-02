"""Local paths of the inputs this repository does not carry.

The paths are read from kelpobserve.toml at the repository root (start from
kelpobserve.example.toml), or from the file named by KELPOBSERVE_CONFIG. Relative paths
resolve against the repository root. Every key is optional: a replay that needs a missing
path stops and names it.
"""

from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BENCHMARK = REPO / "benchmark"
OPTICS = REPO / "optics"
LIBRARY = REPO / "src" / "kelpobserve"
ARCHIVE = BENCHMARK / "archive" / "out"

KEYS = {
    "data": "the CCMAR deliveries and Pleiades Neo archives that v1's input lock pins, as ccmar-diving/ and pneo/",
    "cache": "the benchmark imagery cache: acolite/, acolite-r1/, l1c/, pneo/ and the other folders its stages write",
    "optics_cache": "the optics cache holding the spectral libraries that E5s reads",
    "acolite": "an ACOLITE checkout at revision 64a02ff, whose sensor defaults and response functions the stages read",
    "workspace": "where replays write their runs (default: out/replay in the repository)",
}


@dataclass(frozen=True)
class Paths:
    data: Path | None = None
    cache: Path | None = None
    optics_cache: Path | None = None
    acolite: Path | None = None
    workspace: Path = REPO / "out" / "replay"

    def require(self, key: str) -> Path:
        value = getattr(self, key)
        if value is None:
            raise SystemExit(f"Set paths.{key} in kelpobserve.toml: {KEYS[key]}.")
        if not value.exists():
            raise SystemExit(f"paths.{key} does not exist: {value}")
        return value


def config_file() -> Path:
    return Path(os.environ.get("KELPOBSERVE_CONFIG", REPO / "kelpobserve.toml"))


def load(path: Path | None = None) -> Paths:
    path = path or config_file()
    table = tomllib.loads(path.read_text()).get("paths", {}) if path.is_file() else {}
    unknown = set(table) - set(KEYS)
    if unknown:
        raise SystemExit(f"Unknown keys in {path}: {sorted(unknown)}")
    resolved = {key: (REPO / Path(value).expanduser()).resolve() for key, value in table.items()}
    return Paths(**resolved)
