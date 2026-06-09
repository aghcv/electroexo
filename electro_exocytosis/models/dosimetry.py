from __future__ import annotations

from dataclasses import dataclass

from electro_exocytosis.config import ExposureConfig
from electro_exocytosis.models.pulse import PulseDescriptors


@dataclass(slots=True)
class DosimetryResult:
    field_uniformity_factor: float
    effective_E_V_m: float
    joule_heat_density_J_m3: float
    temperature_rise_K: float


# TODO-literature-review: update geometry factors for cuvette, dish, and flow exposures with experimental/EM references.
def compute_dosimetry(descriptors: PulseDescriptors, exposure: ExposureConfig) -> DosimetryResult:
    """Compute geometry-corrected Layer 1 dosimetry outputs."""
    geometry_factor = {"cuvette": 1.0, "dish": 0.85, "flow": 0.90}[exposure.geometry]
    effective_e = descriptors.E_peak_V_m * geometry_factor
    joule_heat_density = descriptors.energy_density_J_m3 * geometry_factor
    temperature_rise = descriptors.temperature_rise_K * geometry_factor
    return DosimetryResult(
        field_uniformity_factor=geometry_factor,
        effective_E_V_m=effective_e,
        joule_heat_density_J_m3=joule_heat_density,
        temperature_rise_K=temperature_rise,
    )
