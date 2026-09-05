from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from electro_exocytosis.config import SimulationScenario
from electro_exocytosis.io.readers import load_parameter_overrides, load_scenario



def test_load_baseline_scenario() -> None:
    scenario = load_scenario(Path("examples/scenario_baseline.yaml"))
    assert scenario.scenario.name == "baseline_nsPEF_EV_release"
    assert scenario.pulse.amplitude_kV_cm == 10
    assert scenario.simulation.output_dt_s == 60



def test_invalid_negative_amplitude() -> None:
    with pytest.raises(ValidationError):
        SimulationScenario.model_validate(
            {
                "scenario": {"name": "bad"},
                "pulse": {
                    "amplitude_kV_cm": -1,
                    "pulse_width_ns": 100,
                    "pulse_number": 1,
                    "repetition_rate_Hz": 1,
                },
                "exposure": {},
            }
        )



def test_defaults_are_applied() -> None:
    scenario = SimulationScenario.model_validate(
        {
            "scenario": {"name": "defaults"},
            "pulse": {
                "amplitude_kV_cm": 5,
                "pulse_width_ns": 50,
                "pulse_number": 2,
                "repetition_rate_Hz": 10,
            },
            "exposure": {},
        }
    )
    assert scenario.exposure.geometry == "cuvette"
    assert scenario.exposure.dosimetry_model == "legacy"
    assert scenario.cell_state.cell_type == "generic"
    assert scenario.simulation.t_end_s == 7200.0
    assert scenario.extracellular_medium.initial_volume_ml == 1.0
    assert scenario.extracellular_medium.sampling_events == []


def test_extracellular_medium_and_sampling_event_are_validated() -> None:
    scenario = SimulationScenario.model_validate(
        {
            "scenario": {"name": "sampling"},
            "pulse": {
                "amplitude_kV_cm": 5,
                "pulse_width_ns": 50,
                "pulse_number": 2,
                "repetition_rate_Hz": 10,
            },
            "exposure": {"cell_density_per_ml": 2.0e6},
            "extracellular_medium": {
                "initial_volume_ml": 5.0,
                "sampling_events": [
                    {
                        "time_s": 1800.0,
                        "sampled_volume_ml": 0.5,
                        "replacement_volume_ml": 0.5,
                    }
                ],
            },
        }
    )
    assert scenario.extracellular_medium.initial_volume_ml == pytest.approx(5.0)
    assert scenario.extracellular_medium.sampling_events[0].time_s == pytest.approx(
        1800.0
    )


def test_parameter_override_file_must_be_a_mapping(tmp_path: Path) -> None:
    valid_file = tmp_path / "parameters.yml"
    valid_file.write_text(
        "extracellular_kinetics:\n  effective_loss_rate_s: 0.001\n",
        encoding="utf-8",
    )
    assert load_parameter_overrides(valid_file)["extracellular_kinetics"][
        "effective_loss_rate_s"
    ] == pytest.approx(0.001)

    invalid_file = tmp_path / "invalid.yml"
    invalid_file.write_text("- not\n- a\n- mapping\n", encoding="utf-8")
    with pytest.raises(ValueError, match="YAML mapping"):
        load_parameter_overrides(invalid_file)
