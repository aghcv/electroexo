from __future__ import annotations

from dataclasses import dataclass

from electro_exocytosis.config import ExposureConfig, PulseConfig
from electro_exocytosis.units import kV_cm_to_V_m, ns_to_s

DEFAULT_WAVEFORM_ENERGY_FACTORS = {
    "square": 1.0,
    "bipolar": 1.0,
    "exponential": 0.5,
}


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
    train_duration_s: float = 0.0
    adiabatic_temperature_rise_K: float = 0.0
    waveform_energy_factor: float = 1.0
    dosimetry_model: str = "legacy"


def compute_pulse_descriptors(pulse: PulseConfig, exposure: ExposureConfig) -> PulseDescriptors:
    """Compute Layer 1 / Tables A1-A2 pulse descriptors from scenario inputs.

    ``legacy`` preserves the original v0.1 energy and heating estimates.
    The Joule models use the same absorbed-energy skeleton with an optional
    waveform factor: u = sigma * E_peak^2 * pulse_width * pulse_number * chi.
    """
    e_peak = kV_cm_to_V_m(pulse.amplitude_kV_cm)
    pulse_width_s = ns_to_s(pulse.pulse_width_ns)
    total_duration_s = pulse_width_s * pulse.pulse_number
    train_duration_s = pulse.pulse_number / pulse.repetition_rate_Hz
    waveform_energy_factor = _waveform_energy_factor(pulse, exposure)
    if exposure.dosimetry_model == "legacy":
        waveform_energy_factor = 1.0

    energy_density = (
        (e_peak ** 2)
        * exposure.medium_conductivity_S_m
        * pulse_width_s
        * pulse.pulse_number
        * waveform_energy_factor
    )
    dose_index = e_peak * pulse_width_s * pulse.pulse_number
    temperature_capacity = (
        4.0e6 if exposure.dosimetry_model == "legacy" else exposure.volumetric_heat_capacity_J_m3_K
    )
    temperature_rise = energy_density / temperature_capacity
    return PulseDescriptors(
        E_peak_V_m=e_peak,
        pulse_width_s=pulse_width_s,
        pulse_number=pulse.pulse_number,
        repetition_rate_Hz=pulse.repetition_rate_Hz,
        total_duration_s=total_duration_s,
        train_duration_s=train_duration_s,
        energy_density_J_m3=energy_density,
        dose_index=dose_index,
        temperature_rise_K=temperature_rise,
        adiabatic_temperature_rise_K=temperature_rise,
        waveform_energy_factor=waveform_energy_factor,
        dosimetry_model=exposure.dosimetry_model,
        waveform=pulse.waveform,
    )


def _waveform_energy_factor(pulse: PulseConfig, exposure: ExposureConfig) -> float:
    """Return normalized integral of squared waveform over pulse width."""
    if exposure.waveform_energy_factor is not None:
        return exposure.waveform_energy_factor
    return DEFAULT_WAVEFORM_ENERGY_FACTORS[pulse.waveform]
