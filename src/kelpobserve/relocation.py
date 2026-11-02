"""Recorded paths of the experiments, mapped to this checkout.

The registrations, locks and archived runs write each path they record relative to a named
root: $DATA (the field and reference data), $CACHE (the imagery cache, or for the optics
analyses the optics cache), $OUT (the run outputs), $REG (the registration), $LIB (the
package holding kelpobserve.coastal), $ROOT (the experiment's folder), $E5S (the E5s results)
or $ACOLITE (the ACOLITE checkout). A Relocation rewrites the longest matching prefix of such a
path to its local equivalent and leaves other paths unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

DATA, CACHE, OUT, REG, LIB, ROOT, E5S, ACOLITE = "$DATA", "$CACHE", "$OUT", "$REG", "$LIB", "$ROOT", "$E5S", "$ACOLITE"
ROOTS = (DATA, CACHE, OUT, REG, LIB, ROOT, E5S, ACOLITE)


@dataclass(frozen=True)
class Relocation:
    rules: tuple[tuple[str, str], ...]

    @classmethod
    def of(cls, rules: dict[str, Path | str]) -> Relocation:
        ordered = sorted(((old.rstrip("/"), str(new)) for old, new in rules.items()), key=lambda r: -len(r[0]))
        return cls(tuple(ordered))

    def text(self, value: str) -> str:
        for old, new in self.rules:
            if value == old or value.startswith(old + "/"):
                return new + value[len(old):]
        return value

    def __call__(self, path: str | Path) -> Path:
        return Path(self.text(str(path)))

    def json(self, value: Any) -> Any:
        if isinstance(value, str):
            return self.text(value)
        if isinstance(value, list):
            return [self.json(item) for item in value]
        if isinstance(value, dict):
            return {key: self.json(item) for key, item in value.items()}
        return value

    def to_json(self) -> list[list[str]]:
        return [list(rule) for rule in self.rules]

    @classmethod
    def from_json(cls, rules: list[list[str]]) -> Relocation:
        return cls(tuple((old, new) for old, new in rules))
