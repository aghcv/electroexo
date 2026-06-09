from __future__ import annotations

from electro_exocytosis.config import ExposureConfig, PulseConfig
from electro_exocytosis.models.pulse import compute_pulse_descriptors



def test_compute_pulse_descriptors_units() -> None:
    pulse = PulseConfig(amplitude_kV_cm=10, pulse_width_ns=100, pulse_number=10, repetition_rate_Hz=1)
    exposure = ExposureConfig()
    desc = compute_pulse_descriptors(pulse, exposure)
    assert desc.E_peak_V_m == 1e6
    assert desc.pulse_width_s == 100e-9
    assert desc.dose_index > 0
    assert desc.energy_density_J_m3 > 0
