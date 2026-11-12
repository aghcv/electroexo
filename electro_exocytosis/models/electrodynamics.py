from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
from typing import Any

import numpy as np

from electro_exocytosis.config import CellStateConfig
from electro_exocytosis.models.dosimetry import DosimetryResult
from electro_exocytosis.models.pulse import PulseDescriptors


@dataclass(slots=True)
class ElectrodynamicsState:
    delta_Vm: float
    delta_V_ER: float
    delta_V_mito: float
    delta_V_MVB: float
    pore_density: float
    membrane_permeability: float
    cell_radius_m: float = 10e-6
    schwan_limit_V: float = 0.0
    membrane_charging_factor: float = 1.0
    membrane_charging_tau_s: float = 1e-7
    electrodynamics_model: str = "reduced_schwan_permeability"


@dataclass(slots=True)
class ElectrodynamicsParams:
    cell_radius_m: float = 10e-6
    schwan_factor: float = 1.5
    membrane_charging_tau_s: float = 1e-7
    delta_V_ER_fraction: float = 0.3
    delta_V_mito_fraction: float = 0.2
    delta_V_MVB_fraction: float = 0.15
    permeability_threshold_V: float = 0.25
    permeability_slope: float = 10.0
    pore_density_max_m2: float = 1e12
    pore_density_half_voltage_V: float = 0.5
    pore_density_hill_coefficient: float = 3.0


def compute_electrodynamics_state(
    descriptors: PulseDescriptors,
    dosimetry: DosimetryResult,
    cell_state: CellStateConfig,
    params: Mapping[str, Any] | ElectrodynamicsParams | None = None,
) -> ElectrodynamicsState:
    """Compute Layer 2, Tables A3-A4 electrodynamic perturbation state.

    The reduced model keeps the Schwan-style steady-state voltage amplitude but
    adds a first-order charging factor so nanosecond pulse width affects the
    membrane voltage reached during each pulse. Pore density and permeability
    remain phenomenological proxies pending calibrated pore-kinetic models.
    """
    p = _coerce_electrodynamics_params(params)
    schwan_limit = (
        dosimetry.effective_E_V_m
        * p.cell_radius_m
        * p.schwan_factor
        * cell_state.membrane_modifier
    )
    charging_factor = _charging_factor(descriptors.pulse_width_s, p.membrane_charging_tau_s)
    delta_vm = schwan_limit * charging_factor
    delta_v_er = p.delta_V_ER_fraction * delta_vm
    delta_v_mito = p.delta_V_mito_fraction * delta_vm
    delta_v_mvb = p.delta_V_MVB_fraction * delta_vm
    membrane_permeability = _logistic_permeability(
        delta_vm,
        threshold_v=p.permeability_threshold_V,
        slope=p.permeability_slope,
    )
    pore_density = _hill_pore_density(
        delta_vm,
        max_density_m2=p.pore_density_max_m2,
        half_voltage_v=p.pore_density_half_voltage_V,
        hill_coefficient=p.pore_density_hill_coefficient,
    )
    return ElectrodynamicsState(
        delta_Vm=delta_vm,
        delta_V_ER=delta_v_er,
        delta_V_mito=delta_v_mito,
        delta_V_MVB=delta_v_mvb,
        pore_density=pore_density,
        membrane_permeability=membrane_permeability,
        cell_radius_m=p.cell_radius_m,
        schwan_limit_V=schwan_limit,
        membrane_charging_factor=charging_factor,
        membrane_charging_tau_s=p.membrane_charging_tau_s,
    )


def _coerce_electrodynamics_params(
    params: Mapping[str, Any] | ElectrodynamicsParams | None,
) -> ElectrodynamicsParams:
    if params is None:
        return ElectrodynamicsParams()
    if isinstance(params, ElectrodynamicsParams):
        return params
    nested_params = params.get("electrodynamics")
    if isinstance(nested_params, Mapping):
        params = nested_params
    defaults = ElectrodynamicsParams()
    values = {
        field.name: float(params.get(field.name, getattr(defaults, field.name)))
        for field in fields(ElectrodynamicsParams)
    }
    return ElectrodynamicsParams(**values)


def _charging_factor(pulse_width_s: float, tau_s: float) -> float:
    tau_s = max(tau_s, 1e-15)
    return float(np.clip(1.0 - np.exp(-pulse_width_s / tau_s), 0.0, 1.0))


def _logistic_permeability(delta_vm: float, threshold_v: float, slope: float) -> float:
    exponent = -slope * (delta_vm - threshold_v)
    permeability = 1.0 / (1.0 + np.exp(np.clip(exponent, -700.0, 700.0)))
    return float(np.clip(permeability, 0.0, 1.0))


def _hill_pore_density(
    delta_vm: float,
    max_density_m2: float,
    half_voltage_v: float,
    hill_coefficient: float,
) -> float:
    delta_vm = max(delta_vm, 0.0)
    half_voltage_v = max(half_voltage_v, 1e-15)
    hill_coefficient = max(hill_coefficient, 1e-12)
    numerator = delta_vm ** hill_coefficient
    denominator = half_voltage_v ** hill_coefficient + numerator
    return float(max_density_m2 * numerator / denominator)
