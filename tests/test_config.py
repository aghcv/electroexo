from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from electro_exocytosis.config import SimulationScenario
from electro_exocytosis.io.readers import load_scenario



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
    assert scenario.cell_state.cell_type == "generic"
    assert scenario.simulation.t_end_s == 7200.0
