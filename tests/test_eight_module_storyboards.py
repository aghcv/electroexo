from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from examples.storyboard_eight_modules import (
    OUTPUT_DT_S,
    build_storyboards,
    plot_study,
    write_outputs,
)


@pytest.fixture(scope="module")
def storyboards():
    return build_storyboards()


def test_storyboard_outputs_are_reproducible_and_complete(storyboards, tmp_path: Path) -> None:
    write_outputs(storyboards, tmp_path, make_plots=False)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert len(manifest["studies"]) == 5
    assert len(manifest["modules"]) == 8
    assert len(list((tmp_path / "timeseries").glob("*.csv"))) == 15
    assert len(list((tmp_path / "scenarios").glob("*.json"))) == 15
    assert len(list((tmp_path / "snapshots").glob("*.json"))) == 15
    assert len(pd.read_csv(tmp_path / "storyboard_summary.csv")) == 15

    for study in manifest["studies"]:
        assert len(study["rows"]) == 3
        for row in study["rows"]:
            assert len(row["panels"]) == 8
            assert row["parameters"]["manufacturing_qc"]["harvest_time_h"] == 0.5
            for module, panel in zip(manifest["modules"], row["panels"], strict=True):
                assert len(panel["x"]) == len(panel["y"])
                assert np.isfinite(panel["y"]).all()
                assert min(panel["y"]) >= module["y_min"] - 1e-12
                assert max(panel["y"]) <= module["y_max"] + 1e-12
                if module["kind"] == "bar":
                    assert len(panel["y"]) == 1
                    assert panel["value"] == panel["y"][0]

    first_name = next(iter(storyboards.outputs))
    raw = pd.read_csv(tmp_path / "timeseries" / f"{first_name}.csv")
    assert len(raw) == 7201
    assert np.allclose(np.diff(raw["time_s"]), OUTPUT_DT_S)
    assert raw["time_s"].iloc[-1] == 1800.0
    calcium_curve = manifest["studies"][0]["rows"][0]["panels"][2]
    acute_x = np.array(calcium_curve["x"])
    assert np.allclose(np.diff(acute_x[acute_x <= 10.0]), 0.25)

    # Exercise real figure export once at low raster resolution; production uses 600 dpi.
    figure_stem = tmp_path / "figure_smoke"
    plot_study(manifest["studies"][0], manifest["modules"], figure_stem, dpi=72)
    for suffix in ("png", "svg"):
        assert figure_stem.with_suffix(f".{suffix}").stat().st_size > 1000


def test_isolation_changes_only_recovered_yield(storyboards) -> None:
    study = next(
        study for study in storyboards.manifest["studies"]
        if study["id"] == "05_isolation_recovery"
    )
    baseline = study["rows"][0]
    for row in study["rows"][1:]:
        for column in range(7):
            np.testing.assert_array_equal(
                row["panels"][column]["y"], baseline["panels"][column]["y"]
            )
    recovered = np.array([row["panels"][7]["value"] for row in study["rows"]])
    np.testing.assert_allclose(recovered / recovered[0], [1.0, 3.0, 6.0])


def test_higher_pump_capacity_shortens_calcium_response(storyboards) -> None:
    rows = storyboards.summary[storyboards.summary["study"] == "03_calcium_clearance"]
    assert len(rows) == 3
    duration = rows["calcium_above_1uM_duration_s"].to_numpy()
    assert duration[0] > duration[1] > duration[2]
    viability = rows["viability_score_percent"].to_numpy()
    assert viability[0] < viability[1] < viability[2]
    # This is a clearance intervention, so electrical exposure must be unchanged.
    assert rows["deposited_energy_J_ml"].nunique() == 1
    assert rows["delta_Vm_V"].nunique() == 1
