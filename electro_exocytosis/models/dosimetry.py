from __future__ import annotations

from dataclasses import dataclass
from math import exp

from electro_exocytosis.config import ExposureConfig
from electro_exocytosis.models.pulse import PulseDescriptors

GEOMETRY_FIELD_FACTORS = {"cuvette": 1.0, "dish": 0.85, "flow": 0.90}


@dataclass(slots=True)
class DosimetryResult:
    field_uniformity_factor: float
    effective_E_V_m: float
    joule_heat_density_J_m3: float
    temperature_rise_K: float
    mean_E2_factor: float = 1.0
    adiabatic_temperature_rise_K: float = 0.0
    thermal_retention_factor: float = 1.0
    dosimetry_model: str = "legacy"


# TODO-literature-review: update geometry factors for cuvette, dish, and flow exposures with experimental/EM references.
def compute_dosimetry(descriptors: PulseDescriptors, exposure: ExposureConfig) -> DosimetryResult:
    """Compute geometry-corrected Layer 1 dosimetry outputs."""
    geometry_factor = GEOMETRY_FIELD_FACTORS[exposure.geometry]
    effective_e = descriptors.E_peak_V_m * geometry_factor
    if exposure.dosimetry_model == "legacy":
        mean_e2_factor = geometry_factor
        joule_heat_density = descriptors.energy_density_J_m3 * geometry_factor
        adiabatic_temperature_rise = descriptors.temperature_rise_K * geometry_factor
        thermal_retention_factor = 1.0
    else:
        mean_e2_factor = geometry_factor ** 2
        joule_heat_density = descriptors.energy_density_J_m3 * mean_e2_factor
        adiabatic_temperature_rise = joule_heat_density / exposure.volumetric_heat_capacity_J_m3_K
        thermal_retention_factor = _thermal_retention_factor(descriptors, exposure)

    temperature_rise = adiabatic_temperature_rise * thermal_retention_factor
    return DosimetryResult(
        field_uniformity_factor=geometry_factor,
        mean_E2_factor=mean_e2_factor,
        effective_E_V_m=effective_e,
        joule_heat_density_J_m3=joule_heat_density,
        temperature_rise_K=temperature_rise,
        adiabatic_temperature_rise_K=adiabatic_temperature_rise,
        thermal_retention_factor=thermal_retention_factor,
        dosimetry_model=exposure.dosimetry_model,
    )


def _thermal_retention_factor(descriptors: PulseDescriptors, exposure: ExposureConfig) -> float:
    """Fraction of adiabatic heating retained at the end of a finite pulse train."""
    if exposure.dosimetry_model != "joule_lumped_thermal":
        return 1.0

    train_duration_s = max(descriptors.train_duration_s, descriptors.pulse_width_s)
    tau_loss_s = exposure.thermal_relaxation_time_s
    retention = (tau_loss_s / train_duration_s) * (1.0 - exp(-train_duration_s / tau_loss_s))
    return exposure.thermal_efficiency * retention
