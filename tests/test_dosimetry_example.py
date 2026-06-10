from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd


def test_compare_dosimetry_models_example(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "examples" / "compare_dosimetry_models.py"

    subprocess.run(
        [
            sys.executable,
            str(script),
            "--out",
            str(tmp_path),
            "--no-plots",
            "--profile-points",
            "12",
        ],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )

    summary = pd.read_csv(tmp_path / "dosimetry_model_summary.csv")
    profiles = pd.read_csv(tmp_path / "dosimetry_temperature_profiles.csv")

    assert set(summary["dosimetry_model"]) == {"legacy", "joule_adiabatic", "joule_lumped_thermal"}
    assert len(profiles) == 4 * 3 * 12

    low_prr = _end_temp_rise(summary, "hbss_like_high_conductivity_20Hz", "joule_lumped_thermal")
    high_prr = _end_temp_rise(summary, "hbss_like_high_conductivity_200Hz", "joule_lumped_thermal")
    assert high_prr > low_prr

    legacy_dish_energy = _absorbed_energy(summary, "dish_exponential_high_prr", "legacy")
    joule_dish_energy = _absorbed_energy(summary, "dish_exponential_high_prr", "joule_adiabatic")
    assert joule_dish_energy < legacy_dish_energy


def _end_temp_rise(summary: pd.DataFrame, scenario: str, model: str) -> float:
    row = summary[(summary["scenario"] == scenario) & (summary["dosimetry_model"] == model)]
    return float(row["end_temp_rise_C"].iloc[0])


def _absorbed_energy(summary: pd.DataFrame, scenario: str, model: str) -> float:
    row = summary[(summary["scenario"] == scenario) & (summary["dosimetry_model"] == model)]
    return float(row["absorbed_energy_density_mJ_mm3"].iloc[0])
