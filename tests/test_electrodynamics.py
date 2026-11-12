from __future__ import annotations

from math import exp

import pytest

from electro_exocytosis.config import CellStateConfig, ExposureConfig, PulseConfig
from electro_exocytosis.models.dosimetry import compute_dosimetry
from electro_exocytosis.models.electrodynamics import (
    ElectrodynamicsParams,
    compute_electrodynamics_state,
)
from electro_exocytosis.models.pulse import compute_pulse_descriptors


def _electro_state(
    *,
    amplitude_kV_cm: float = 1.0,
    pulse_width_ns: float = 100.0,
    params: ElectrodynamicsParams | dict[str, object] | None = None,
):
    pulse = PulseConfig(
        amplitude_kV_cm=amplitude_kV_cm,
        pulse_width_ns=pulse_width_ns,
        pulse_number=1,
        repetition_rate_Hz=1,
    )
    exposure = ExposureConfig()
    descriptors = compute_pulse_descriptors(pulse, exposure)
    dosimetry = compute_dosimetry(descriptors, exposure)
    return compute_electrodynamics_state(descriptors, dosimetry, CellStateConfig(), params)


def test_cell_radius_parameter_scales_membrane_voltage() -> None:
    small_cell = _electro_state(params={"cell_radius_m": 5e-6})
    large_cell = _electro_state(params={"cell_radius_m": 10e-6})
    nested_cell = _electro_state(params={"electrodynamics": {"cell_radius_m": 20e-6}})

    assert small_cell.cell_radius_m == pytest.approx(5e-6)
    assert large_cell.cell_radius_m == pytest.approx(10e-6)
    assert large_cell.delta_Vm == pytest.approx(2.0 * small_cell.delta_Vm)
    assert nested_cell.delta_Vm == pytest.approx(2.0 * large_cell.delta_Vm)


def test_membrane_charging_factor_depends_on_pulse_width() -> None:
    params = ElectrodynamicsParams(membrane_charging_tau_s=100e-9)
    short_pulse = _electro_state(pulse_width_ns=10.0, params=params)
    long_pulse = _electro_state(pulse_width_ns=300.0, params=params)

    assert short_pulse.membrane_charging_factor == pytest.approx(1.0 - exp(-0.1))
    assert long_pulse.membrane_charging_factor > 0.9
    assert short_pulse.delta_Vm < long_pulse.delta_Vm


def test_custom_permeability_and_pore_parameters_are_used() -> None:
    baseline = _electro_state(amplitude_kV_cm=0.1)
    permissive = _electro_state(
        amplitude_kV_cm=0.1,
        params={
            "permeability_threshold_V": 0.05,
            "pore_density_max_m2": 2.0e12,
            "pore_density_half_voltage_V": 0.1,
        },
    )

    assert permissive.membrane_permeability > baseline.membrane_permeability
    assert permissive.pore_density > baseline.pore_density
