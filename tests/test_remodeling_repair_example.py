from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

from examples.compare_remodeling_repair import build_response_table, write_outputs


def test_remodeling_repair_example_builds_summary(tmp_path: Path) -> None:
    summary, timeseries = build_response_table()

    assert not summary.empty
    assert not timeseries.empty
    assert {
        "scenario",
        "peak_ca_submembrane_uM",
        "peak_ps_exposure",
        "peak_calpain_activity",
        "peak_annexin_activity",
        "peak_repair_state",
        "cumulative_repair_shedding",
    }.issubset(summary.columns)
    buffered = summary[summary["scenario"] == "buffered_microdomain"].iloc[0]
    baseline = summary[summary["scenario"] == "baseline_repair"].iloc[0]
    calpain_inhibited = summary[summary["scenario"] == "calpain_inhibited"].iloc[0]

    assert buffered["peak_ca_submembrane_uM"] < baseline["peak_ca_submembrane_uM"]
    assert calpain_inhibited["peak_calpain_activity"] < baseline["peak_calpain_activity"]

    write_outputs(summary, timeseries, tmp_path, make_plots=False)
    assert (tmp_path / "remodeling_repair_summary.csv").exists()
    assert (tmp_path / "remodeling_repair_timeseries.csv").exists()


def test_remodeling_repair_example_cli(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "examples" / "compare_remodeling_repair.py"

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

    summary = pd.read_csv(tmp_path / "remodeling_repair_summary.csv")
    assert set(summary["scenario"]) == {
        "baseline_repair",
        "buffered_microdomain",
        "calpain_inhibited",
        "annexin_lysosome_strong",
    }
