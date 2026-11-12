from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

from examples.compare_ion_transport_bioenergetics import build_response_table, write_outputs


def test_ion_transport_bioenergetics_example_builds_summary(tmp_path: Path) -> None:
    summary, timeseries = build_response_table()

    assert not summary.empty
    assert not timeseries.empty
    assert {
        "scenario",
        "peak_ca_i_uM",
        "peak_ca_mito_uM",
        "peak_ros",
        "min_atp",
        "min_mitochondrial_potential",
        "peak_osmotic_stress",
        "peak_J_Ca_pore_uM_s",
    }.issubset(summary.columns)
    calcium_limited = summary[summary["scenario"] == "calcium_limited_medium"].iloc[0]
    baseline = summary[summary["scenario"] == "baseline_100ns"].iloc[0]
    assert calcium_limited["peak_J_Ca_pore_uM_s"] < baseline["peak_J_Ca_pore_uM_s"]

    write_outputs(summary, timeseries, tmp_path, make_plots=False)
    assert (tmp_path / "ion_transport_bioenergetics_summary.csv").exists()
    assert (tmp_path / "ion_transport_bioenergetics_timeseries.csv").exists()


def test_ion_transport_bioenergetics_example_cli(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "examples" / "compare_ion_transport_bioenergetics.py"

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

    summary = pd.read_csv(tmp_path / "ion_transport_bioenergetics_summary.csv")
    assert set(summary["scenario"]) == {
        "baseline_100ns",
        "wide_pulse_high_dose",
        "calcium_limited_medium",
        "osmotic_recovery_slow",
    }
