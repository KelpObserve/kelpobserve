from pathlib import Path

import pandas as pd
import pytest

from kelpobserve.cli import main
from kelpobserve.monthly import evaluation


def test_monthly_namespace_preserves_existing_commands(capsys):
    with pytest.raises(SystemExit) as stopped:
        main(["--help"])
    assert stopped.value.code == 0
    help_text = capsys.readouterr().out
    for name in [
        "experiments",
        "replay",
        "recompare",
        "kd490",
        "optics-inputs",
        "local-cover",
        "pneo-optical-review",
        "monthly",
    ]:
        assert name in help_text
    assert main(["experiments"]) == 0
    listing = capsys.readouterr().out
    for name in ["v1:", "r1:", "r3:", "r4:", "r5:", "r6:", "e5s:"]:
        assert name in listing


def test_monthly_evaluation_dispatches_explicit_paths(monkeypatch, tmp_path):
    seen = []

    def evaluate(inputs, output):
        seen.append((inputs, output))
        return pd.DataFrame(columns=["subset", "campaign", "experiment", "n", "mae_pp", "rmse_pp"])

    monkeypatch.setattr(evaluation, "evaluate", evaluate)
    destination = tmp_path / "result"
    assert main(["monthly", "evaluate", "--inputs", "local-inputs", "--output", str(destination)]) == 0
    assert seen == [(Path("local-inputs"), destination)]
