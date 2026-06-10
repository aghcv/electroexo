from __future__ import annotations

from math import exp

import pytest

from electro_exocytosis.config import ExposureConfig, PulseConfig
from electro_exocytosis.models.dosimetry import compute_dosimetry
from electro_exocytosis.models.pulse import compute_pulse_descriptors


def test_compute_pulse_descriptors_units() -> None:
    pulse = PulseConfig(amplitude_kV_cm=10, pulse_width_ns=100, pulse_number=10, repetition_rate_Hz=1)
    exposure = ExposureConfig()
    desc = compute_pulse_descriptors(pulse, exposure)
    assert desc.E_peak_V_m == 1e6
    assert desc.pulse_width_s == 100e-9
    assert desc.dose_index > 0
    assert desc.energy_density_J_m3 > 0


def test_legacy_dosimetry_preserves_original_energy_and_heating() -> None:
    pulse = PulseConfig(amplitude_kV_cm=10, pulse_width_ns=100, pulse_number=10, repetition_rate_Hz=1)
    exposure = ExposureConfig()
    desc = compute_pulse_descriptors(pulse, exposure)
    dosimetry = compute_dosimetry(desc, exposure)

    assert desc.dosimetry_model == "legacy"
    assert desc.waveform_energy_factor == 1.0
    assert desc.train_duration_s == 10.0
    assert desc.energy_density_J_m3 == pytest.approx(1.5e6)
    assert desc.temperature_rise_K == pytest.approx(1.5e6 / 4.0e6)
    assert dosimetry.mean_E2_factor == 1.0
    assert dosimetry.joule_heat_density_J_m3 == pytest.approx(desc.energy_density_J_m3)
    assert dosimetry.temperature_rise_K == pytest.approx(desc.temperature_rise_K)


def test_joule_adiabatic_uses_waveform_energy_factor() -> None:
    pulse = PulseConfig(
        amplitude_kV_cm=10,
        pulse_width_ns=100,
        pulse_number=10,
        repetition_rate_Hz=1,
        waveform="exponential",
    )
    exposure = ExposureConfig(dosimetry_model="joule_adiabatic")
    desc = compute_pulse_descriptors(pulse, exposure)
    dosimetry = compute_dosimetry(desc, exposure)

    assert desc.waveform_energy_factor == 0.5
    assert desc.energy_density_J_m3 == pytest.approx(0.5 * 1.5e6)
    assert dosimetry.joule_heat_density_J_m3 == pytest.approx(desc.energy_density_J_m3)
    assert dosimetry.temperature_rise_K == pytest.approx(
        desc.energy_density_J_m3 / exposure.volumetric_heat_capacity_J_m3_K
    )


def test_joule_dosimetry_scales_energy_with_mean_e_squared_geometry() -> None:
    pulse = PulseConfig(amplitude_kV_cm=10, pulse_width_ns=100, pulse_number=10, repetition_rate_Hz=1)
    exposure = ExposureConfig(geometry="dish", dosimetry_model="joule_adiabatic")
    desc = compute_pulse_descriptors(pulse, exposure)
    dosimetry = compute_dosimetry(desc, exposure)

    assert dosimetry.field_uniformity_factor == pytest.approx(0.85)
    assert dosimetry.mean_E2_factor == pytest.approx(0.85 ** 2)
    assert dosimetry.effective_E_V_m == pytest.approx(0.85 * desc.E_peak_V_m)
    assert dosimetry.joule_heat_density_J_m3 == pytest.approx(desc.energy_density_J_m3 * 0.85 ** 2)


def test_lumped_thermal_model_reduces_end_of_train_temperature() -> None:
    pulse = PulseConfig(amplitude_kV_cm=100, pulse_width_ns=10, pulse_number=500, repetition_rate_Hz=20)
    exposure = ExposureConfig(
        dosimetry_model="joule_lumped_thermal",
        medium_conductivity_S_m=1.4,
        thermal_relaxation_time_s=5.0,
    )
    desc = compute_pulse_descriptors(pulse, exposure)
    dosimetry = compute_dosimetry(desc, exposure)

    expected_retention = (5.0 / 25.0) * (1.0 - exp(-25.0 / 5.0))
    assert desc.train_duration_s == pytest.approx(25.0)
    assert dosimetry.thermal_retention_factor == pytest.approx(expected_retention)
    assert dosimetry.temperature_rise_K < dosimetry.adiabatic_temperature_rise_K
    assert dosimetry.temperature_rise_K == pytest.approx(
        dosimetry.adiabatic_temperature_rise_K * expected_retention
    )
