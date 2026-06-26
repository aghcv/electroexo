from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
from typing import Any, Callable


EV_STATE_NAMES = [
    "MVB_pool",
    "ILV_load",
    "docked_MVB_pool",
    "budding_pool",
    "apoptotic_commitment",
]


@dataclass(slots=True)
class EVReleaseParams:
    MVB_pool_baseline: float = 1.0
    ILV_load_baseline: float = 0.7
    docked_MVB_baseline: float = 0.2
    budding_pool_baseline: float = 0.15
    apoptotic_commitment_baseline: float = 0.02

    baseline_sEV_rate: float = 1.0
    baseline_mlEV_rate: float = 0.12
    baseline_AB_rate: float = 0.01

    k_MVB_maturation_s: float = 0.0020
    k_MVB_docking_s: float = 0.0010
    k_MVB_lysosomal_s: float = 0.0008
    k_ILV_ESCRT_dependent_s: float = 0.0024
    k_ILV_ESCRT_independent_s: float = 0.0016
    k_ILV_release_s: float = 0.0013
    k_ILV_lysosomal_s: float = 0.0010
    k_docked_MVB_turnover_s: float = 0.0011
    k_budding_s: float = 0.0012
    k_budding_turnover_s: float = 0.0010
    k_apoptotic_commitment_s: float = 0.0008
    k_apoptotic_resolution_s: float = 0.0004

    K_rab_Ca_uM: float = 0.8
    n_rab_Ca: float = 2.0
    K_fusion_Ca_uM: float = 1.0
    n_fusion_Ca: float = 2.5
    K_synaptotagmin7_Ca_uM: float = 0.7
    n_synaptotagmin7_Ca: float = 2.0
    K_PS_budding: float = 0.35
    n_PS_budding: float = 2.0
    K_damage_apoptotic: float = 0.7
    n_damage_apoptotic: float = 3.0
    K_lysosomal_routing: float = 0.45
    n_lysosomal_routing: float = 2.0
    K_ILV_release: float = 0.5
    MVB_voltage_threshold_V: float = 0.04
    MVB_voltage_slope: float = 30.0
    atp_half_max: float = 0.4
    ros_reference: float = 0.1
    ros_stress_scale: float = 0.4

    rab_conversion_baseline: float = 0.35
    rab_conversion_Ca_weight: float = 0.30
    rab_conversion_voltage_weight: float = 0.25

    escrt_baseline: float = 0.40
    escrt_rab_weight: float = 0.35
    escrt_repair_weight: float = 0.15
    escrt_damage_penalty: float = 0.25

    ceramide_baseline: float = 0.25
    ceramide_PS_weight: float = 0.35
    ceramide_repair_weight: float = 0.25
    ceramide_ROS_weight: float = 0.20

    acidification_baseline: float = 0.40
    acidification_ATP_weight: float = 0.35
    acidification_ceramide_relief: float = 0.45
    lysosomal_damage_weight: float = 0.20

    rab_docking_baseline: float = 0.30
    rab27_weight: float = 0.35
    rab11_weight: float = 0.25
    rab35_weight: float = 0.20
    actin_barrier_penalty: float = 0.20

    fusion_baseline: float = 0.25
    fusion_SNARE_weight: float = 0.35
    fusion_Munc13_weight: float = 0.25
    fusion_synaptotagmin_weight: float = 0.25

    budding_baseline: float = 0.20
    budding_arf6_weight: float = 0.30
    budding_ps_weight: float = 0.30
    budding_tension_weight: float = 0.20
    budding_calpain_weight: float = 0.20

    scission_baseline: float = 0.25
    scission_rhoa_weight: float = 0.35
    scission_calpain_weight: float = 0.25
    scission_ps_weight: float = 0.20
    repair_shedding_mlEV_weight: float = 2.0

    apoptosis_baseline: float = 0.05
    apoptosis_damage_weight: float = 0.45
    apoptosis_ATP_loss_weight: float = 0.20
    apoptosis_ROS_weight: float = 0.15
    apoptosis_tension_weight: float = 0.15


def build_ev_release_rhs(
    params: EVReleaseParams,
) -> Callable[[float, list[float]], list[float]]:
    """Build a reduced Layer 5 RHS callable.

    The returned function expects keyword-only upstream drivers via ``**kwargs``:
    ``Ca_i``, ``Ca_submembrane``, ``ROS``, ``ATP``, ``damage_state``,
    ``delta_V_MVB``, ``pore_activation``, ``PS_exposure``, ``calpain_activity``,
    ``annexin_activity``, ``actomyosin_tension``, ``actin_disruption``,
    ``repair_state``, and ``repair_shedding_rate``.
    """

    def rhs(t: float, y: list[float], **kwargs: float) -> list[float]:
        _ = t
        fluxes = compute_ev_release_fluxes(params, y, **kwargs)
        return compute_ev_release_derivatives(params, y, fluxes)

    return rhs


def compute_ev_release_fluxes(
    params: EVReleaseParams,
    y: list[float],
    *,
    Ca_i: float,
    Ca_submembrane: float,
    ROS: float,
    ATP: float,
    damage_state: float,
    delta_V_MVB: float,
    pore_activation: float,
    PS_exposure: float,
    calpain_activity: float,
    annexin_activity: float,
    actomyosin_tension: float,
    actin_disruption: float,
    repair_state: float,
    repair_shedding_rate: float,
) -> dict[str, float]:
    """Compute Layer 5 mechanistic signals and subtype release rates."""
    state = ev_state_to_dict(y, params)
    ca_i = max(float(Ca_i), 0.0)
    ca_submembrane = max(float(Ca_submembrane), 0.0)
    ros = max(float(ROS), 0.0)
    atp = max(float(ATP), 0.0)
    damage_state = max(float(damage_state), 0.0)
    pore_activation = _clip01(float(pore_activation))
    ps_exposure = _clip01(float(PS_exposure))
    calpain_activity = _clip01(float(calpain_activity))
    annexin_activity = _clip01(float(annexin_activity))
    actomyosin_tension = _clip01(float(actomyosin_tension))
    actin_disruption = _clip01(float(actin_disruption))
    repair_state = _clip01(float(repair_state))
    repair_shedding_rate = max(float(repair_shedding_rate), 0.0)

    atp_gate = atp / (params.atp_half_max + atp) if atp > 0 else 0.0
    ros_gate = _clip01(max(ros - params.ros_reference, 0.0) / max(params.ros_stress_scale, 1e-12))
    voltage_gate = _logistic_activation(
        float(delta_V_MVB),
        threshold=params.MVB_voltage_threshold_V,
        slope=params.MVB_voltage_slope,
    )
    ca_rab_gate = _hill(ca_i, params.K_rab_Ca_uM, params.n_rab_Ca)
    ca_fusion_gate = _hill(ca_submembrane, params.K_fusion_Ca_uM, params.n_fusion_Ca)
    syt7_gate = _hill(ca_submembrane, params.K_synaptotagmin7_Ca_uM, params.n_synaptotagmin7_Ca)
    ps_gate = _hill(ps_exposure, params.K_PS_budding, params.n_PS_budding)
    damage_gate = _hill(damage_state, params.K_damage_apoptotic, params.n_damage_apoptotic)
    munc13_gate = _clip01(0.65 * ca_fusion_gate + 0.20 * atp_gate + 0.15 * pore_activation)

    rab_conversion_signal = _clip01(
        params.rab_conversion_baseline
        + params.rab_conversion_Ca_weight * ca_rab_gate
        + params.rab_conversion_voltage_weight * voltage_gate
    )
    escrt_dependent_signal = _clip01(
        params.escrt_baseline
        + params.escrt_rab_weight * rab_conversion_signal
        + params.escrt_repair_weight * repair_state
        - params.escrt_damage_penalty * damage_gate
    )
    ceramide_signal = _clip01(
        params.ceramide_baseline
        + params.ceramide_PS_weight * ps_gate
        + params.ceramide_repair_weight * repair_state
        + params.ceramide_ROS_weight * ros_gate
    )
    acidification_signal = _clip01(
        params.acidification_baseline
        + params.acidification_ATP_weight * atp_gate
        - params.acidification_ceramide_relief * ceramide_signal
    )
    lysosomal_routing = _hill(
        acidification_signal + params.lysosomal_damage_weight * (damage_gate + ros_gate),
        params.K_lysosomal_routing,
        params.n_lysosomal_routing,
    )
    secretory_bias = _clip01(1.0 - lysosomal_routing)
    rab_docking_signal = _clip01(
        params.rab_docking_baseline
        + params.rab27_weight * rab_conversion_signal
        + params.rab11_weight * munc13_gate
        + params.rab35_weight * actin_disruption
        - params.actin_barrier_penalty * actomyosin_tension
    )
    fusion_signal = (
        _clip01(
            params.fusion_baseline
            + params.fusion_SNARE_weight * ca_fusion_gate
            + params.fusion_Munc13_weight * munc13_gate
            + params.fusion_synaptotagmin_weight * syt7_gate
            + 0.05 * annexin_activity
        )
        * atp_gate
        * (0.30 + 0.70 * secretory_bias)
    )
    budding_signal = _clip01(
        params.budding_baseline
        + params.budding_arf6_weight * repair_state
        + params.budding_ps_weight * ps_gate
        + params.budding_tension_weight * actomyosin_tension
        + params.budding_calpain_weight * calpain_activity
    )
    scission_signal = _clip01(
        params.scission_baseline
        + params.scission_rhoa_weight * actomyosin_tension
        + params.scission_calpain_weight * calpain_activity
        + params.scission_ps_weight * ps_gate
    )
    apoptotic_blebbing_signal = _clip01(
        params.apoptosis_baseline
        + params.apoptosis_damage_weight * damage_gate
        + params.apoptosis_ATP_loss_weight * max(1.0 - atp_gate, 0.0)
        + params.apoptosis_ROS_weight * ros_gate
        + params.apoptosis_tension_weight * actomyosin_tension
    )

    ilv_availability = state["ILV_load"] / (params.K_ILV_release + state["ILV_load"])
    release_propensity = (
        params.k_ILV_release_s
        * fusion_signal
        * state["docked_MVB_pool"]
        * ilv_availability
        * (0.50 + 0.25 * escrt_dependent_signal + 0.25 * ceramide_signal)
        * (0.50 + 0.50 * secretory_bias)
    )
    mlEV_propensity = (
        scission_signal
        * state["budding_pool"]
        * (1.0 + params.repair_shedding_mlEV_weight * repair_shedding_rate)
    )
    ab_propensity = apoptotic_blebbing_signal * state["apoptotic_commitment"]

    return {
        "rab_conversion_signal": float(rab_conversion_signal),
        "escrt_dependent_signal": float(escrt_dependent_signal),
        "ceramide_signal": float(ceramide_signal),
        "acidification_signal": float(acidification_signal),
        "lysosomal_routing": float(lysosomal_routing),
        "secretory_bias": float(secretory_bias),
        "rab_docking_signal": float(rab_docking_signal),
        "fusion_signal": float(fusion_signal),
        "budding_signal": float(budding_signal),
        "scission_signal": float(scission_signal),
        "apoptotic_blebbing_signal": float(apoptotic_blebbing_signal),
        "release_propensity": float(release_propensity),
        "sEV_rate": float(params.baseline_sEV_rate * release_propensity),
        "mlEV_rate": float(params.baseline_mlEV_rate * mlEV_propensity),
        "AB_rate": float(params.baseline_AB_rate * ab_propensity),
    }


def compute_ev_release_derivatives(
    params: EVReleaseParams,
    y: list[float],
    fluxes: Mapping[str, float],
) -> list[float]:
    """Return Layer 5 state derivatives from the current state and fluxes."""
    state = ev_state_to_dict(y, params)
    d_mvb = (
        params.k_MVB_maturation_s * fluxes["rab_conversion_signal"] * fluxes["secretory_bias"]
        - params.k_MVB_docking_s * fluxes["rab_docking_signal"] * state["MVB_pool"]
        - params.k_MVB_lysosomal_s * fluxes["lysosomal_routing"] * state["MVB_pool"]
    )
    d_ilv = (
        (
            params.k_ILV_ESCRT_dependent_s * fluxes["escrt_dependent_signal"]
            + params.k_ILV_ESCRT_independent_s * fluxes["ceramide_signal"]
        )
        * state["MVB_pool"]
        - fluxes["release_propensity"]
        - params.k_ILV_lysosomal_s * fluxes["lysosomal_routing"] * state["ILV_load"]
    )
    d_docked = (
        params.k_MVB_docking_s * fluxes["rab_docking_signal"] * state["MVB_pool"]
        - params.k_docked_MVB_turnover_s * fluxes["fusion_signal"] * state["docked_MVB_pool"]
    )
    d_budding = (
        params.k_budding_s * fluxes["budding_signal"]
        - params.k_budding_turnover_s * fluxes["scission_signal"] * state["budding_pool"]
    )
    d_apoptotic = (
        params.k_apoptotic_commitment_s * fluxes["apoptotic_blebbing_signal"]
        - params.k_apoptotic_resolution_s * fluxes["secretory_bias"] * state["apoptotic_commitment"]
    )
    return [
        float(d_mvb),
        float(d_ilv),
        float(d_docked),
        float(d_budding),
        float(d_apoptotic),
    ]


def coerce_ev_release_params(
    params: Mapping[str, Any] | EVReleaseParams | None,
) -> EVReleaseParams:
    """Build EV release parameters from flat or nested configuration."""
    if params is None:
        return EVReleaseParams()
    if isinstance(params, EVReleaseParams):
        return params
    nested_params = params.get("ev_release")
    if isinstance(nested_params, Mapping):
        params = nested_params
    defaults = EVReleaseParams()
    values = {
        field.name: float(params.get(field.name, getattr(defaults, field.name)))
        for field in fields(EVReleaseParams)
    }
    return EVReleaseParams(**values)


def get_ev_state_names() -> list[str]:
    """Return the state order used by Layer 5 EV release."""
    return list(EV_STATE_NAMES)


def get_ev_initial_conditions(params: EVReleaseParams) -> list[float]:
    """Return initial conditions in ``get_ev_state_names()`` order."""
    return [
        float(params.MVB_pool_baseline),
        float(params.ILV_load_baseline),
        float(params.docked_MVB_baseline),
        float(params.budding_pool_baseline),
        float(params.apoptotic_commitment_baseline),
    ]


def ev_state_to_dict(
    y: list[float],
    params: EVReleaseParams,
) -> dict[str, float]:
    """Return a clipped state mapping for Layer 5 variables."""
    values = [max(float(value), 0.0) for value in y[: len(EV_STATE_NAMES)]]
    state = dict(zip(EV_STATE_NAMES, values, strict=True))
    state["MVB_pool"] = max(state["MVB_pool"], 0.0)
    state["ILV_load"] = max(state["ILV_load"], 0.0)
    state["docked_MVB_pool"] = max(state["docked_MVB_pool"], 0.0)
    state["budding_pool"] = max(state["budding_pool"], 0.0)
    state["apoptotic_commitment"] = max(state["apoptotic_commitment"], 0.0)
    return state


def ev_release_defaults() -> dict[str, float]:
    """Return EV release defaults as a dict."""
    return asdict(EVReleaseParams())


def _hill(value: float, half_value: float, hill_coefficient: float) -> float:
    value = max(float(value), 0.0)
    half_value = max(float(half_value), 1e-12)
    hill_coefficient = max(float(hill_coefficient), 1e-12)
    numerator = value**hill_coefficient
    return float(numerator / (half_value**hill_coefficient + numerator))


def _logistic_activation(value: float, threshold: float, slope: float) -> float:
    exponent = -float(slope) * (float(value) - float(threshold))
    return float(1.0 / (1.0 + pow(2.718281828459045, exponent)))


def _clip01(value: float) -> float:
    return float(max(0.0, min(1.0, value)))
