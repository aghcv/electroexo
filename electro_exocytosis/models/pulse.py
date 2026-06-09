from __future__ import annotations

from dataclasses import dataclass

from electro_exocytosis.config import ExposureConfig, PulseConfig
from electro_exocytosis.units import kV_cm_to_V_m, ns_to_s


@dataclass(slots=True)
class PulseDescriptors:
    E_peak_V_m: float
    pulse_width_s: float
    pulse_number: int
    repetition_rate_Hz: float
    total_duration_s: float
    energy_density_J_m3: float
    dose_index: float
    temperature_rise_K: float
    waveform: str


# TODO-literature-review: replace placeholder energy and heating formulas with calibrated nsPEF dosimetry relations.
def compute_pulse_descriptors(pulse: PulseConfig, exposure: ExposureConfig) -> PulseDescriptors:
    """Compute Layer 1 / Tables A1-A2 pulse descriptors from scenario inputs."""
    e_peak = kV_cm_to_V_m(pulse.amplitude_kV_cm)
    pulse_width_s = ns_to_s(pulse.pulse_width_ns)
    total_duration_s = pulse_width_s * pulse.pulse_number
    energy_density = (e_peak ** 2) * exposure.medium_conductivity_S_m * pulse_width_s * pulse.pulse_number
    dose_index = e_peak * pulse_width_s * pulse.pulse_number
    temperature_rise = energy_density / 4.0e6
    return PulseDescriptors(
        E_peak_V_m=e_peak,
        pulse_width_s=pulse_width_s,
        pulse_number=pulse.pulse_number,
        repetition_rate_Hz=pulse.repetition_rate_Hz,
        total_duration_s=total_duration_s,
        energy_density_J_m3=energy_density,
        dose_index=dose_index,
        temperature_rise_K=temperature_rise,
        waveform=pulse.waveform,
    )
