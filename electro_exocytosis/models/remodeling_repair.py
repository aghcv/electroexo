from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
from typing import Any


@dataclass(slots=True)
class RemodelingParams:
    K_PS_uM: float = 1.0
    n_PS: float = 2.0
    PS_max: float = 1.0
    K_calpain_uM: float = 2.0
    n_calpain: float = 2.0
    K_annex_uM: float = 0.5
    n_annex: float = 1.5
    tau_repair_s: float = 300.0
    microdomain_gain: float = 2.0
    microdomain_pore_gain: float = 4.0
    microdomain_osmotic_gain: float = 1.0
    K_scramblase_uM: float = 0.8
    n_scramblase: float = 2.0
    flippase_Ca_K_uM: float = 0.8
    flippase_min_activity: float = 0.2
    K_lysosomal_repair_uM: float = 0.6
    n_lysosomal_repair: float = 2.0
    K_actomyosin_uM: float = 1.5
    n_actomyosin: float = 2.0
    actin_calpain_weight: float = 0.55
    actin_ps_weight: float = 0.25
    actin_osmotic_weight: float = 0.20
    resealing_annexin_weight: float = 0.50
    resealing_lysosome_weight: float = 0.35
    resealing_ps_weight: float = 0.15
    resealing_calpain_penalty: float = 0.25
    shedding_rate_scale: float = 0.05


def compute_remodeling_state(
    Ca_i: float,
    params: RemodelingParams,
    *,
    osmotic_stress: float = 0.0,
    mitochondrial_potential: float = 1.0,
    pore_activation: float = 0.0,
) -> dict[str, float]:
    """Compute reduced Layer 4 remodeling and repair observables."""
    ca_i = max(float(Ca_i), 0.0)
    osmotic_stress = max(float(osmotic_stress), 0.0)
    mitochondrial_stress = max(1.0 - float(mitochondrial_potential), 0.0)
    pore_activation = _clip01(float(pore_activation))

    local_ca = ca_i * (
        1.0
        + params.microdomain_gain * pore_activation
        + params.microdomain_pore_gain * pore_activation**2
        + params.microdomain_osmotic_gain * osmotic_stress
    )
    scramblase = _hill(local_ca, params.K_scramblase_uM, params.n_scramblase)
    flippase = params.flippase_min_activity + (1.0 - params.flippase_min_activity) * (
        1.0 - _hill(local_ca, params.flippase_Ca_K_uM, params.n_PS)
    )
    ps_exposure = _clip01(params.PS_max * scramblase * (1.0 - 0.5 * flippase))

    calpain = _hill(local_ca, params.K_calpain_uM, params.n_calpain)
    annexin = _hill(local_ca, params.K_annex_uM, params.n_annex)
    lysosomal_repair = _hill(local_ca, params.K_lysosomal_repair_uM, params.n_lysosomal_repair)
    actomyosin = _clip01(
        _hill(local_ca, params.K_actomyosin_uM, params.n_actomyosin)
        + 0.5 * osmotic_stress
        + 0.25 * mitochondrial_stress
    )
    actin_disruption = _clip01(
        params.actin_calpain_weight * calpain
        + params.actin_ps_weight * ps_exposure
        + params.actin_osmotic_weight * osmotic_stress
    )
    repair_state = _clip01(
        params.resealing_annexin_weight * annexin
        + params.resealing_lysosome_weight * lysosomal_repair
        + params.resealing_ps_weight * ps_exposure
        - params.resealing_calpain_penalty * calpain
    )
    repair_shedding_rate = (
        params.shedding_rate_scale
        * repair_state
        * ps_exposure
        * (0.5 + 0.5 * actomyosin)
        * (1.0 + osmotic_stress)
    )
    return {
        "Ca_submembrane": float(local_ca),
        "PS_exposure": float(ps_exposure),
        "scramblase_activity": float(scramblase),
        "flippase_activity": float(flippase),
        "calpain_activity": float(calpain),
        "annexin_activity": float(annexin),
        "lysosomal_repair_activity": float(lysosomal_repair),
        "actomyosin_tension": float(actomyosin),
        "actin_disruption": float(actin_disruption),
        "repair_state": float(repair_state),
        "repair_shedding_rate": float(repair_shedding_rate),
    }


def coerce_remodeling_params(
    params: Mapping[str, Any] | RemodelingParams | None,
) -> RemodelingParams:
    """Build remodeling parameters from flat or nested configuration."""
    if params is None:
        return RemodelingParams()
    if isinstance(params, RemodelingParams):
        return params
    nested_params = params.get("remodeling_repair")
    if isinstance(nested_params, Mapping):
        params = nested_params
    defaults = RemodelingParams()
    values = {
        field.name: float(params.get(field.name, getattr(defaults, field.name)))
        for field in fields(RemodelingParams)
    }
    return RemodelingParams(**values)


def remodeling_defaults() -> dict[str, float]:
    """Return remodeling defaults as a dict."""
    return asdict(RemodelingParams())


def _hill(value: float, half_value: float, hill_coefficient: float) -> float:
    value = max(float(value), 0.0)
    half_value = max(float(half_value), 1e-12)
    hill_coefficient = max(float(hill_coefficient), 1e-12)
    numerator = value**hill_coefficient
    return float(numerator / (half_value**hill_coefficient + numerator))


def _clip01(value: float) -> float:
    return float(max(0.0, min(1.0, value)))
