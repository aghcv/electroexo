from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
from typing import Any, Callable

import numpy as np

from electro_exocytosis.models.electrodynamics import ElectrodynamicsState


ION_STATE_NAMES = [
    "Ca_i",
    "Ca_ER",
    "Ca_mito",
    "mitochondrial_potential",
    "ROS",
    "ATP",
    "Na_i",
    "K_i",
    "Cl_i",
    "osmotic_stress",
]


@dataclass(slots=True)
class IonTransportParams:
    Ca_baseline_uM: float = 0.1
    Ca_ext_uM: float = 1200.0
    Ca_ER_uM: float = 500.0
    Ca_mito_baseline_uM: float = 0.2
    Ca_max_uM: float = 10.0
    tau_pore_reseal_s: float = 0.5
    tau_Ca_homeostasis_s: float = 600.0
    J_Ca_pore_factor: float = 1.5
    J_ER_release_factor: float = 1.5
    ER_activation_threshold_V: float = 0.05
    ER_activation_slope: float = 30.0
    SERCA_Vmax_uM_s: float = 0.15
    SERCA_K_uM: float = 0.4
    SERCA_ATP_K: float = 0.3
    PMCA_Vmax_uM_s: float = 0.20
    PMCA_K_uM: float = 0.5
    NCX_Vmax_uM_s: float = 0.06
    NCX_K_uM: float = 0.8
    calcium_hill_n: float = 2.0
    mitochondrial_Ca_K_uM: float = 0.8
    mitochondrial_uptake_Vmax_uM_s: float = 0.08
    mitochondrial_release_rate_s: float = 0.01
    mPTP_release_rate_s: float = 0.08
    tau_mito_recovery_s: float = 600.0
    mitochondrial_Ca_depolarization_factor_s: float = 0.001
    mitochondrial_ROS_depolarization_factor_s: float = 0.001
    mitochondrial_electro_depolarization_factor_s: float = 0.001
    Na_baseline_mM: float = 12.0
    Na_ext_mM: float = 145.0
    K_baseline_mM: float = 140.0
    K_ext_mM: float = 4.0
    Cl_baseline_mM: float = 10.0
    Cl_ext_mM: float = 110.0
    J_Na_pore_factor_mM_s: float = 1.5
    J_K_pore_factor_mM_s: float = 1.5
    J_Cl_pore_factor_mM_s: float = 1.0
    tau_ion_recovery_s: float = 300.0
    tau_osmotic_s: float = 30.0
    ROS_baseline: float = 0.1
    ROS_production_factor: float = 0.0001
    ROS_mito_factor_s: float = 0.0002
    ROS_depolarization_factor_s: float = 0.0005
    ROS_osmotic_factor_s: float = 0.0002
    ROS_pulse_factor_s: float = 0.0001
    tau_ROS_s: float = 120.0
    ATP_baseline: float = 1.0
    ATP_depletion_factor: float = 0.05
    ATP_pump_cost_factor: float = 0.0002
    ATP_osmotic_cost_factor: float = 0.0002
    ATP_ROS_cost_factor: float = 0.0001
    tau_ATP_s: float = 600.0
    ATP_damage_threshold: float = 0.5
    osmotic_damage_scale: float = 2.0
    mitochondrial_damage_scale: float = 0.5
    ATP_damage_scale: float = 0.5


def build_ion_transport_rhs(
    params: IonTransportParams,
    electro_state: ElectrodynamicsState,
    t_pulse_end: float,
) -> Callable[[float, list[float] | np.ndarray], list[float]]:
    """Build the Layer 3 RHS for ion, organelle, ROS, and ATP dynamics."""

    def rhs(t: float, y: list[float] | np.ndarray) -> list[float]:
        fluxes = compute_ion_transport_fluxes(params, electro_state, t, y, t_pulse_end)
        state = ion_state_to_dict(y, params)

        d_ca_i = (
            fluxes["J_Ca_pore"]
            + fluxes["J_ER_release"]
            + fluxes["J_Ca_homeostasis"]
            - fluxes["J_SERCA"]
            - fluxes["J_PMCA"]
            - fluxes["J_NCX"]
            - fluxes["J_mito_uptake"]
            + fluxes["J_mito_release"]
        )
        d_ca_er = fluxes["J_SERCA"] - fluxes["J_ER_release"]
        d_ca_mito = fluxes["J_mito_uptake"] - fluxes["J_mito_release"]
        d_mito_potential = (
            (1.0 - state["mitochondrial_potential"]) / max(params.tau_mito_recovery_s, 1e-9)
            - fluxes["mitochondrial_depolarization"]
        )
        d_ros = (
            fluxes["ROS_source"]
            - (state["ROS"] - params.ROS_baseline) / max(params.tau_ROS_s, 1e-9)
        )
        d_atp = fluxes["ATP_production"] - fluxes["ATP_consumption"]
        d_na_i = fluxes["J_Na_pore"] + (params.Na_baseline_mM - state["Na_i"]) / max(params.tau_ion_recovery_s, 1e-9)
        d_k_i = -fluxes["J_K_pore"] + (params.K_baseline_mM - state["K_i"]) / max(params.tau_ion_recovery_s, 1e-9)
        d_cl_i = fluxes["J_Cl_pore"] + (params.Cl_baseline_mM - state["Cl_i"]) / max(params.tau_ion_recovery_s, 1e-9)
        d_osmotic = (fluxes["osmotic_target"] - state["osmotic_stress"]) / max(params.tau_osmotic_s, 1e-9)

        return [
            float(d_ca_i),
            float(d_ca_er),
            float(d_ca_mito),
            float(d_mito_potential),
            float(d_ros),
            float(d_atp),
            float(d_na_i),
            float(d_k_i),
            float(d_cl_i),
            float(d_osmotic),
        ]

    return rhs


def compute_ion_transport_fluxes(
    params: IonTransportParams,
    electro_state: ElectrodynamicsState,
    t: float,
    y: list[float] | np.ndarray,
    t_pulse_end: float,
) -> dict[str, float]:
    """Compute named Layer 3 fluxes for diagnostics and coupling."""
    state = ion_state_to_dict(y, params)
    local_t = max(float(t) - float(t_pulse_end), 0.0)
    pore_decay = float(np.exp(-local_t / max(params.tau_pore_reseal_s, 1e-9)))
    pore_activation = float(np.clip(electro_state.membrane_permeability * pore_decay, 0.0, 1.0))
    er_activation = float(
        _logistic_activation(
            electro_state.delta_V_ER,
            threshold=params.ER_activation_threshold_V,
            slope=params.ER_activation_slope,
        )
        * pore_decay
    )

    ca_i = state["Ca_i"]
    ca_er = state["Ca_ER"]
    ca_mito = state["Ca_mito"]
    mito_potential = state["mitochondrial_potential"]
    ros = state["ROS"]
    atp = state["ATP"]
    osmotic_stress = state["osmotic_stress"]

    ca_drive = max(params.Ca_ext_uM - ca_i, 0.0) / max(params.Ca_ext_uM, 1e-9)
    j_ca_pore = params.J_Ca_pore_factor * pore_activation * ca_drive
    er_drive = max(ca_er - ca_i, 0.0) / max(params.Ca_ER_uM, 1e-9)
    j_er_release = params.J_ER_release_factor * er_activation * er_drive

    ca_excess_uM = max(ca_i - params.Ca_baseline_uM, 0.0)
    ca_hill = _hill(ca_excess_uM, params.SERCA_K_uM, params.calcium_hill_n)
    atp_gate = atp / (params.SERCA_ATP_K + atp) if atp > 0 else 0.0
    er_capacity_gate = max(params.Ca_ER_uM - ca_er, 0.0) / max(params.Ca_ER_uM, 1e-9)
    j_serca = params.SERCA_Vmax_uM_s * ca_hill * atp_gate * er_capacity_gate
    j_pmca = params.PMCA_Vmax_uM_s * _hill(ca_excess_uM, params.PMCA_K_uM, params.calcium_hill_n) * atp_gate
    j_ncx = params.NCX_Vmax_uM_s * _hill(ca_excess_uM, params.NCX_K_uM, params.calcium_hill_n)
    j_ca_homeostasis = max(params.Ca_baseline_uM - ca_i, 0.0) / max(
        params.tau_Ca_homeostasis_s, 1e-9
    )

    mito_ca_hill = _hill(ca_excess_uM, params.mitochondrial_Ca_K_uM, params.calcium_hill_n)
    j_mito_uptake = params.mitochondrial_uptake_Vmax_uM_s * mito_ca_hill * mito_potential
    depolarization_gate = max(1.0 - mito_potential, 0.0)
    j_mito_release = (
        params.mitochondrial_release_rate_s * max(ca_mito - params.Ca_mito_baseline_uM, 0.0)
        + params.mPTP_release_rate_s * depolarization_gate * ca_mito
    )

    j_na_pore = params.J_Na_pore_factor_mM_s * pore_activation * max(params.Na_ext_mM - state["Na_i"], 0.0) / max(params.Na_ext_mM, 1e-9)
    j_k_pore = params.J_K_pore_factor_mM_s * pore_activation * max(state["K_i"] - params.K_ext_mM, 0.0) / max(params.K_baseline_mM, 1e-9)
    j_cl_pore = params.J_Cl_pore_factor_mM_s * pore_activation * max(params.Cl_ext_mM - state["Cl_i"], 0.0) / max(params.Cl_ext_mM, 1e-9)

    baseline_osmoles = params.Na_baseline_mM + params.K_baseline_mM + params.Cl_baseline_mM
    current_osmoles = state["Na_i"] + state["K_i"] + state["Cl_i"]
    osmotic_target = max((current_osmoles - baseline_osmoles) / max(baseline_osmoles, 1e-9), 0.0)

    ca_excess = max(ca_i / max(params.Ca_baseline_uM, 1e-9) - 1.0, 0.0)
    mito_ca_excess = max(ca_mito / max(params.Ca_mito_baseline_uM, 1e-9) - 1.0, 0.0)
    mitochondrial_depolarization = (
        params.mitochondrial_Ca_depolarization_factor_s * mito_ca_excess
        + params.mitochondrial_ROS_depolarization_factor_s * max(ros - params.ROS_baseline, 0.0)
        + params.mitochondrial_electro_depolarization_factor_s * pore_decay * max(electro_state.delta_V_mito, 0.0)
    ) * mito_potential
    ros_source = (
        params.ROS_production_factor * ca_excess
        + params.ROS_mito_factor_s * mito_ca_excess
        + params.ROS_depolarization_factor_s * depolarization_gate
        + params.ROS_osmotic_factor_s * osmotic_stress
        + params.ROS_pulse_factor_s * pore_activation
    )
    pump_load = j_serca + j_pmca + j_ncx
    atp_production = max(params.ATP_baseline - atp, 0.0) * mito_potential / max(params.tau_ATP_s, 1e-9)
    atp_consumption = (
        params.ATP_pump_cost_factor * pump_load
        + params.ATP_osmotic_cost_factor * osmotic_stress
        + params.ATP_ROS_cost_factor * max(ros - params.ROS_baseline, 0.0)
        + params.ATP_depletion_factor * max(ros - params.ROS_baseline, 0.0) / max(params.tau_ATP_s, 1e-9)
    )

    return {
        "pore_activation": float(pore_activation),
        "ER_activation": float(er_activation),
        "J_Ca_pore": float(j_ca_pore),
        "J_ER_release": float(j_er_release),
        "J_Ca_homeostasis": float(j_ca_homeostasis),
        "J_SERCA": float(j_serca),
        "J_PMCA": float(j_pmca),
        "J_NCX": float(j_ncx),
        "J_mito_uptake": float(j_mito_uptake),
        "J_mito_release": float(j_mito_release),
        "J_Na_pore": float(j_na_pore),
        "J_K_pore": float(j_k_pore),
        "J_Cl_pore": float(j_cl_pore),
        "osmotic_target": float(osmotic_target),
        "mitochondrial_depolarization": float(mitochondrial_depolarization),
        "ROS_source": float(ros_source),
        "ATP_production": float(atp_production),
        "ATP_consumption": float(atp_consumption),
    }


def ion_state_to_dict(y: list[float] | np.ndarray, params: IonTransportParams) -> dict[str, float]:
    """Return a clipped state mapping for Layer 3 variables."""
    values = [float(value) for value in y[: len(ION_STATE_NAMES)]]
    state = dict(zip(ION_STATE_NAMES, values, strict=True))
    state["Ca_i"] = max(state["Ca_i"], 0.0)
    state["Ca_ER"] = float(np.clip(state["Ca_ER"], 0.0, params.Ca_ER_uM))
    state["Ca_mito"] = max(state["Ca_mito"], 0.0)
    state["mitochondrial_potential"] = float(np.clip(state["mitochondrial_potential"], 0.0, 1.0))
    state["ROS"] = max(state["ROS"], 0.0)
    state["ATP"] = float(np.clip(state["ATP"], 0.0, params.ATP_baseline))
    state["Na_i"] = max(state["Na_i"], 0.0)
    state["K_i"] = max(state["K_i"], 0.0)
    state["Cl_i"] = max(state["Cl_i"], 0.0)
    state["osmotic_stress"] = max(state["osmotic_stress"], 0.0)
    return state


def get_ion_state_names() -> list[str]:
    """Return the state order used by Layer 3 ion transport."""
    return list(ION_STATE_NAMES)


def get_ion_initial_conditions(params: IonTransportParams) -> list[float]:
    """Return initial conditions in ``get_ion_state_names()`` order."""
    return [
        params.Ca_baseline_uM,
        params.Ca_ER_uM,
        params.Ca_mito_baseline_uM,
        1.0,
        params.ROS_baseline,
        params.ATP_baseline,
        params.Na_baseline_mM,
        params.K_baseline_mM,
        params.Cl_baseline_mM,
        0.0,
    ]


def coerce_ion_transport_params(
    params: Mapping[str, Any] | IonTransportParams | None,
) -> IonTransportParams:
    """Build ion-transport parameters from flat or nested configuration."""
    if params is None:
        return IonTransportParams()
    if isinstance(params, IonTransportParams):
        return params
    nested_params = params.get("ion_transport")
    if isinstance(nested_params, Mapping):
        params = nested_params
    defaults = IonTransportParams()
    values = {
        field.name: float(params.get(field.name, getattr(defaults, field.name)))
        for field in fields(IonTransportParams)
    }
    return IonTransportParams(**values)


def ion_transport_defaults() -> dict[str, float]:
    """Return ion transport defaults as a plain dict."""
    return asdict(IonTransportParams())


def _hill(value: float, half_value: float, hill_coefficient: float) -> float:
    value = max(value, 0.0)
    half_value = max(half_value, 1e-12)
    hill_coefficient = max(hill_coefficient, 1e-12)
    numerator = value**hill_coefficient
    return float(numerator / (half_value**hill_coefficient + numerator))


def _logistic_activation(value: float, threshold: float, slope: float) -> float:
    exponent = -slope * (value - threshold)
    activation = 1.0 / (1.0 + np.exp(np.clip(exponent, -700.0, 700.0)))
    return float(np.clip(activation, 0.0, 1.0))
