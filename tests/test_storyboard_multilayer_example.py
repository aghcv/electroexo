from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

from examples.storyboard_multilayer_overview import build_storyboard_dataset, write_outputs


def test_storyboard_multilayer_example_builds_summary(tmp_path: Path) -> None:
    summary, outputs = build_storyboard_dataset()

    assert len(outputs) == 3
    assert not summary.empty
    assert {
        "scenario",
        "dose_index",
        "peak_ca_i",
        "peak_repair_state",
        "cumulative_small_EV",
        "cumulative_apoptotic_body",
    }.issubset(summary.columns)

    mild = summary[summary["scenario"] == "mild_reversible_window"].iloc[0]
    productive = summary[summary["scenario"] == "productive_secretory_window"].iloc[0]
    injury = summary[summary["scenario"] == "injury_apoptotic_window"].iloc[0]

    assert productive["cumulative_small_EV"] > mild["cumulative_small_EV"]
    assert injury["cumulative_apoptotic_body"] > productive["cumulative_apoptotic_body"]
    assert injury["viability_fraction"] < mild["viability_fraction"]

    write_outputs(summary, outputs, tmp_path, make_plots=False)
    assert (tmp_path / "multilayer_storyboard_summary.csv").exists()


def test_storyboard_multilayer_example_cli(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "examples" / "storyboard_multilayer_overview.py"

    subprocess.run(
        [
            sys.executable,
            str(script),
            "--out",
            str(tmp_path),
            "--no-plots",
        ],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )

    summary = pd.read_csv(tmp_path / "multilayer_storyboard_summary.csv")
    assert set(summary["scenario"]) == {
        "mild_reversible_window",
        "productive_secretory_window",
        "injury_apoptotic_window",
    }
