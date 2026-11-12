from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

from examples.compare_ev_biogenesis_release import build_response_table, write_outputs


def test_ev_biogenesis_example_builds_summary(tmp_path: Path) -> None:
    summary, timeseries = build_response_table()

    assert not summary.empty
    assert not timeseries.empty
    assert {
        "scenario",
        "peak_secretory_bias",
        "peak_lysosomal_routing",
        "peak_docked_MVB_pool",
        "cumulative_small_EV",
        "cumulative_apoptotic_body",
    }.issubset(summary.columns)

    baseline = summary[summary["scenario"] == "baseline_ev_biogenesis"].iloc[0]
    limited = summary[summary["scenario"] == "rab_escrt_limited"].iloc[0]
    ceramide = summary[summary["scenario"] == "ceramide_secretory_bias"].iloc[0]
    apoptotic = summary[summary["scenario"] == "injury_apoptotic_shift"].iloc[0]

    assert limited["cumulative_small_EV"] < baseline["cumulative_small_EV"]
    assert ceramide["peak_secretory_bias"] > baseline["peak_secretory_bias"]
    assert ceramide["cumulative_small_EV"] > baseline["cumulative_small_EV"]
    assert apoptotic["cumulative_apoptotic_body"] > baseline["cumulative_apoptotic_body"]

    write_outputs(summary, timeseries, tmp_path, make_plots=False)
    assert (tmp_path / "ev_biogenesis_summary.csv").exists()
    assert (tmp_path / "ev_biogenesis_timeseries.csv").exists()


def test_ev_biogenesis_example_cli(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "examples" / "compare_ev_biogenesis_release.py"

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

    summary = pd.read_csv(tmp_path / "ev_biogenesis_summary.csv")
    assert set(summary["scenario"]) == {
        "baseline_ev_biogenesis",
        "rab_escrt_limited",
        "ceramide_secretory_bias",
        "injury_apoptotic_shift",
    }
