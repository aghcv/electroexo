"""Validate a completed release demo and optionally compare a reference run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import yaml

ARTIFACTS = (
    "summary.json", "run_metadata.json", "state_timeseries.csv", "ev_outputs.csv",
    "parameters_used.yaml", "abbreviations.json", "abbreviations.md", "abbreviations.tex",
)
TIME_GRID = [float(value) for value in range(0, 31, 5)]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _finite(value, location: str) -> None:
    if isinstance(value, float):
        _require(math.isfinite(value), f"{location}: nonfinite number")
    elif isinstance(value, dict):
        for key, item in value.items():
            _finite(item, f"{location}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _finite(item, f"{location}[{index}]")


def _table(path: Path) -> tuple[list[str], list[list[float | bool]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        records = list(csv.reader(handle))
    _require(bool(records) and "time_s" in records[0], f"{path}: missing time_s column")
    headers, rows = records[0], records[1:]
    _require(len(rows) == len(TIME_GRID), f"{path}: expected seven rows")
    parsed = []
    for index, row in enumerate(rows):
        _require(len(row) == len(headers), f"{path}: row {index} has wrong column count")
        values = []
        for column, cell in zip(headers, row):
            if cell in ("True", "False"):
                values.append(cell == "True")
            else:
                try:
                    value = float(cell)
                except ValueError as exc:
                    raise ValueError(f"{path}: row {index}, {column}: invalid number") from exc
                _finite(value, f"{path}: row {index}, {column}")
                values.append(value)
        parsed.append(values)
    times = [row[headers.index("time_s")] for row in parsed]
    _require(not any(isinstance(time, bool) for time in times) and times == TIME_GRID,
             f"{path}: expected time_s grid 0..30 in steps of 5")
    return headers, parsed


def _load(run_dir: Path) -> dict:
    for name in ARTIFACTS:
        path = run_dir / name
        _require(path.is_file() and path.stat().st_size > 0, f"Missing or empty artifact: {path}")
    metadata = json.loads((run_dir / "run_metadata.json").read_text(encoding="utf-8"))
    _require(isinstance(metadata, dict), f"{run_dir}: metadata must be an object")
    _require(metadata.get("scenario_name") == "release_demo", f"{run_dir}: wrong scenario name")
    _require(metadata.get("n_timepoints") == 7, f"{run_dir}: expected seven time points")
    _require(isinstance(metadata.get("warnings"), list), f"{run_dir}: warnings must be a list")
    _require(all(isinstance(item, str) for item in metadata["warnings"]),
             f"{run_dir}: warnings must contain strings")
    metadata.pop("generated_at", None)
    summary = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    _require(isinstance(summary, dict) and bool(summary), f"{run_dir}: summary must be an object")
    _finite(summary, "summary")
    parameters = yaml.safe_load((run_dir / "parameters_used.yaml").read_text(encoding="utf-8"))
    _require(isinstance(parameters, dict), f"{run_dir}: parameters must be a mapping")
    for key in ("scenario", "parameters_nested"):
        _require(isinstance(parameters.get(key), dict), f"{run_dir}: missing parameter section {key}")
    _finite(parameters, "parameters")
    provenance_path = run_dir / "run_provenance.json"
    if provenance_path.exists():
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        _require(isinstance(provenance, dict), f"{run_dir}: provenance must be an object")
        _require(provenance.get("exit_code") == 0, f"{run_dir}: provenance reports an unsuccessful run")
    return {
        "metadata": metadata, "summary": summary, "parameters": parameters,
        **{name: _table(run_dir / name) for name in ("state_timeseries.csv", "ev_outputs.csv")},
    }


def _compare(actual, expected, location: str, rtol: float, atol: float) -> None:
    if (
        isinstance(actual, (int, float)) and not isinstance(actual, bool)
        and isinstance(expected, (int, float)) and not isinstance(expected, bool)
    ):
        _require(math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol),
                 f"{location}: {actual!r} differs from reference {expected!r}")
    elif isinstance(actual, dict) and isinstance(expected, dict):
        _require(actual.keys() == expected.keys(), f"{location}: keys differ from reference")
        for key in actual:
            _compare(actual[key], expected[key], f"{location}.{key}", rtol, atol)
    elif isinstance(actual, (list, tuple)) and isinstance(expected, (list, tuple)):
        _require(len(actual) == len(expected), f"{location}: length differs from reference")
        for index, (item, reference) in enumerate(zip(actual, expected)):
            _compare(item, reference, f"{location}[{index}]", rtol, atol)
    else:
        _require(type(actual) is type(expected) and actual == expected,
                 f"{location}: {actual!r} differs from reference {expected!r}")


def check_demo(
    run_dir: Path, reference: Path | None = None, rtol: float = 1e-7, atol: float = 1e-10,
) -> list[str]:
    """Raise ValueError on incomplete, nonfinite, or numerically different results."""
    _require(math.isfinite(rtol) and rtol >= 0, "rtol must be finite and nonnegative")
    _require(math.isfinite(atol) and atol >= 0, "atol must be finite and nonnegative")
    actual = _load(run_dir)
    if reference is not None:
        _compare(actual, _load(reference), "run", rtol, atol)
    return actual["metadata"]["warnings"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-10)
    arguments = parser.parse_args()
    try:
        warnings = check_demo(**vars(arguments))
    except (OSError, ValueError, yaml.YAMLError) as exc:
        parser.exit(1, f"Demo verification failed: {exc}\n")
    for warning in warnings:
        print(f"Model warning: {warning}")
    print(f"Demo verified: {arguments.run_dir}")


if __name__ == "__main__":
    main()
