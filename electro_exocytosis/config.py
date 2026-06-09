from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


class PulseConfig(BaseModel):
    amplitude_kV_cm: float = Field(gt=0, description="Peak electric field, kV/cm")
    pulse_width_ns: float = Field(gt=0)
    pulse_number: int = Field(ge=1)
    repetition_rate_Hz: float = Field(gt=0)
    waveform: Literal["square", "bipolar", "exponential"] = "square"


class ExposureConfig(BaseModel):
    geometry: Literal["cuvette", "dish", "flow"] = "cuvette"
    medium_conductivity_S_m: float = Field(default=1.5, gt=0)
    temperature_C: float = Field(default=37.0)
    cell_density_per_ml: float = Field(default=1e6, gt=0)


class CellStateConfig(BaseModel):
    cell_type: str = "generic"
    membrane_modifier: float = 1.0
    calcium_handling_modifier: float = 1.0
    baseline_EV_release_modifier: float = 1.0
    stress_sensitivity_modifier: float = 1.0


class SimulationConfig(BaseModel):
    t_start_s: float = 0.0
    t_end_s: float = 7200.0
    output_dt_s: float = 1.0
    numerical_method: Literal["solve_ivp", "euler"] = "solve_ivp"


class ScenarioConfig(BaseModel):
    name: str
    mode: Literal["cell_based_electro_exocytosis", "direct_EV_engineering"] = "cell_based_electro_exocytosis"


class SimulationScenario(BaseModel):
    scenario: ScenarioConfig
    pulse: PulseConfig
    exposure: ExposureConfig
    cell_state: CellStateConfig = Field(default_factory=CellStateConfig)
    simulation: SimulationConfig = Field(default_factory=SimulationConfig)
