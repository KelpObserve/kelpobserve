import hashlib
import json

import pytest

from kelpobserve import compare, locks, replay
from kelpobserve.config import Paths
from kelpobserve.relocation import OUT, ROOT, Relocation


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_relocation_takes_the_longest_prefix_and_leaves_other_paths():
    relocate = Relocation.of({"/a": "/x", "/a/b": "/y"})
    assert relocate.text("/a/b/c.json") == "/y/c.json"
    assert relocate.text("/a/c.json") == "/x/c.json"
    assert relocate.text("/ab/c.json") == "/ab/c.json"
    assert relocate.json({"k": ["/a/b", 3, {"p": "/a"}]}) == {"k": ["/y", 3, {"p": "/x"}]}


def test_root_tokens_resolve_only_at_a_path_boundary():
    relocate = Relocation.of({OUT: "/w/out", ROOT: "/w", f"{ROOT}/inputs": "/repo/inputs"})
    assert relocate.text(f"{OUT}/fit/runs/x") == "/w/out/fit/runs/x"
    assert relocate.text(f"{ROOT}/inputs/habitat_aoi.geojson") == "/repo/inputs/habitat_aoi.geojson"
    assert relocate.text(f"{ROOT}/benchmark/run_benchmark.py") == "/w/benchmark/run_benchmark.py"
    assert relocate.text(f"{OUT}PUT/x") == f"{OUT}PUT/x"


def test_lock_check_classifies_every_entry(tmp_path):
    same, edited = tmp_path / "same.csv", tmp_path / "edited.csv"
    same.write_text("a,b\n1,2\n")
    edited.write_text("changed\n")
    lock = tmp_path / "lock.json"
    lock.write_text(json.dumps({"frozen_at": "2026-09-20", "files": {
        "$DATA/same.csv": digest(same.read_bytes()),
        "$DATA/edited.csv": digest(b"original\n"),
        "$DATA/gone.csv": digest(b""),
    }}))
    result = locks.check(lock, Relocation.of({"$DATA": tmp_path}))
    assert result.verified == ["$DATA/same.csv"]
    assert result.changed == ["$DATA/edited.csv"]
    assert result.missing == ["$DATA/gone.csv"]
    assert not result.ok


def test_runs_compare_values_paths_and_provenance_separately(tmp_path):
    archived, replayed = tmp_path / "archived", tmp_path / "replayed"
    for folder in (archived, replayed):
        folder.mkdir()
    (archived / "same.csv").write_text("x\n1\n")
    (replayed / "same.csv").write_text("x\n1\n")
    (archived / "result.json").write_text(json.dumps({"mae": 26.3, "path": "$OUT/r", "created_at": "t0"}))
    (replayed / "result.json").write_text(json.dumps({"path": "/machine-b/out/r", "mae": 26.3, "created_at": "t1"}))
    (archived / "worse.csv").write_text("x\n1.0\n")
    (replayed / "worse.csv").write_text("x\n1.5\n")
    (replayed / "big.csv").write_text("draws\n")
    table = [["/machine-b/out", "$OUT"]]
    result = compare.runs(archived, replayed, {"big.csv": digest(b"draws\n")}, table)
    files = result["files"]
    assert files["same.csv"]["status"] == "identical"
    assert files["result.json"]["status"] == "equal"
    assert files["result.json"]["provenance_differences"] == ["/created_at: archive t0 | replay t1"]
    assert files["worse.csv"]["status"] == "different"
    assert files["big.csv"] == {"status": "identical", "compared": "SHA-256 in the archive index"}
    assert result["reproduced"] is False


def test_the_emodnet_rasters_resolve_into_the_repository_and_match_the_v1_lock(tmp_path):
    relocate = replay.relocation(replay.EXPERIMENTS["r3"], Paths(cache=tmp_path, acolite=tmp_path), tmp_path)
    pinned = json.loads(replay.EXPERIMENTS["v1"].lock.read_text())["files"]
    rasters = {path: expected for path, expected in pinned.items() if path.startswith("$DATA/external/")}
    assert len(rasters) == 2
    for path, expected in rasters.items():
        assert relocate(path).parent == replay.EMODNET
        assert locks.sha256(relocate(path)) == expected


def test_a_replay_stops_on_an_unset_path_it_needs_and_only_then(tmp_path):
    with pytest.raises(SystemExit, match="paths.acolite"):
        replay.relocation(replay.EXPERIMENTS["r3"], Paths(cache=tmp_path), tmp_path)
    with pytest.raises(SystemExit, match="paths.data"):
        replay.relocation(replay.EXPERIMENTS["v1"], Paths(cache=tmp_path), tmp_path)
    replay.relocation(replay.EXPERIMENTS["v1"], Paths(cache=tmp_path, data=tmp_path), tmp_path)
