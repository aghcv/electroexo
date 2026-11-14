from __future__ import annotations

import csv
import json
import runpy
from pathlib import Path

import pytest
import yaml

CHECKER = Path(__file__).resolve().parents[1] / "examples" / "slurm" / "check_demo.py"
CHECK_DEMO = runpy.run_path(str(CHECKER))["check_demo"]


def _demo(path: Path, *, value: float = 1.0, generated_at: str = "first") -> Path:
    path.mkdir()
    metadata = {"scenario_name": "release_demo", "n_timepoints": 7,
                "warnings": ["Exploratory defaults"], "generated_at": generated_at}
    (path / "run_metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    (path / "summary.json").write_text(
        json.dumps({"scenario_name": "release_demo", "yield": value, "warnings": metadata["warnings"]}),
        encoding="utf-8",
    )
    (path / "parameters_used.yaml").write_text(
        yaml.safe_dump({"scenario": {"name": "release_demo"}, "parameters_nested": {"value": value}}),
        encoding="utf-8",
    )
    for name in ("state_timeseries.csv", "ev_outputs.csv"):
        with (path / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time_s", "value", "quality_pass"])
            writer.writerows((time, value, False) for time in range(0, 31, 5))
    for suffix in ("json", "md", "tex"):
        (path / f"abbreviations.{suffix}").write_text("{}", encoding="utf-8")
    return path


def test_demo_checker_accepts_warnings_and_tolerated_reference(tmp_path: Path) -> None:
    reference = _demo(tmp_path / "reference")
    actual = _demo(tmp_path / "actual", value=1.0 + 1e-8, generated_at="second")
    assert CHECK_DEMO(actual, reference) == ["Exploratory defaults"]


@pytest.mark.parametrize("corruption", ["missing", "time", "nonfinite", "summary", "parameters", "failed"])
def test_demo_checker_rejects_invalid_artifacts(tmp_path: Path, corruption: str) -> None:
    run = _demo(tmp_path / "run")
    if corruption == "missing":
        (run / "abbreviations.tex").unlink()
    elif corruption == "time":
        path = run / "ev_outputs.csv"
        path.write_text(path.read_text(encoding="utf-8").replace("30,", "29,"), encoding="utf-8")
    elif corruption == "nonfinite":
        path = run / "state_timeseries.csv"
        path.write_text(path.read_text(encoding="utf-8").replace("1.0", "nan"), encoding="utf-8")
    elif corruption == "summary":
        (run / "summary.json").write_text('{"value": Infinity}', encoding="utf-8")
    elif corruption == "parameters":
        (run / "parameters_used.yaml").write_text("[]", encoding="utf-8")
    else:
        (run / "run_provenance.json").write_text('{"exit_code": 2}', encoding="utf-8")
    with pytest.raises(ValueError):
        CHECK_DEMO(run)


def test_demo_checker_rejects_numerical_reference_difference(tmp_path: Path) -> None:
    reference = _demo(tmp_path / "reference")
    actual = _demo(tmp_path / "actual", value=1.01)
    with pytest.raises(ValueError, match="differs from reference"):
        CHECK_DEMO(actual, reference)


@pytest.mark.parametrize("tolerance", [float("nan"), float("inf"), -1.0])
def test_demo_checker_rejects_invalid_tolerance(tmp_path: Path, tolerance: float) -> None:
    with pytest.raises(ValueError, match="rtol must"):
        CHECK_DEMO(tmp_path, rtol=tolerance)
