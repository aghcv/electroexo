from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable


@dataclass(slots=True)
class InjuryParams:
    K_apoptosis_damage: float = 0.5
    K_necrosis_damage: float = 0.9
    K_stress_damage: float = 0.25
    n_apoptosis: float = 3.0
    n_necrosis: float = 5.0
    debris_fraction_scale: float = 0.2
    aggregate_fraction_scale: float = 0.08
    apoptotic_body_contamination_weight: float = 0.75
    necrotic_debris_weight: float = 0.85
    marker_panel_weight: float = 0.70
    contamination_threshold: float = 0.3
    viability_threshold: float = 0.5
    damage_rate: float = 0.01
    repair_rate: float = 0.005


# TODO-literature-review: replace phenomenological damage accumulation with calibrated injury pathways.
def build_injury_rhs(params: InjuryParams) -> Callable[[float, list[float], float], list[float]]:
    """Build placeholder Layer 7 damage accumulation dynamics."""

    def rhs(t: float, y: list[float], damage_input: float) -> list[float]:
        (damage,) = y
        d_damage = params.damage_rate * damage_input - params.repair_rate * damage / (1.0 + damage)
        return [float(d_damage)]

    return rhs



def compute_quality_gate(
    damage: float,
    debris: float,
    params: InjuryParams,
    *,
    cumulative_sEV: float = 0.0,
    cumulative_mlEV: float = 0.0,
    cumulative_AB: float = 0.0,
) -> dict[str, float | bool]:
    """Compute Layer 7 fate fractions, measured-particle mixture, and quality gate."""
    damage = max(float(damage), 0.0)
    stress = _hill(damage, params.K_stress_damage, 2.0)
    apoptosis = _hill(damage, params.K_apoptosis_damage, params.n_apoptosis)
    necrosis = _hill(damage, params.K_necrosis_damage, params.n_necrosis)
    apoptosis = min(apoptosis, 0.8)
    necrosis = min(necrosis, 0.9)
    stressed_viable = max(0.0, min(1.0, stress * (1.0 - apoptosis) * (1.0 - necrosis)))
    viability = max(0.0, 1.0 - 0.25 * stressed_viable - 0.6 * apoptosis - 0.9 * necrosis)
    ev_particles = max(float(cumulative_sEV) + float(cumulative_mlEV), 0.0)
    apoptotic_particles = max(float(cumulative_AB), 0.0)
    debris_particles = max(float(debris), 0.0) + params.necrotic_debris_weight * necrosis
    aggregate_particles = params.aggregate_fraction_scale * (ev_particles + apoptotic_particles + debris_particles)
    measured_particles = ev_particles + apoptotic_particles + debris_particles + aggregate_particles
    bona_fide_fraction = ev_particles / measured_particles if measured_particles > 0 else max(0.0, 1.0 - debris)
    marker_score = _clip01(
        params.marker_panel_weight * bona_fide_fraction
        + (1.0 - params.marker_panel_weight) * (1.0 - apoptosis)
    )
    purity = _clip01(
        bona_fide_fraction
        * marker_score
        * (1.0 - params.apoptotic_body_contamination_weight * _fraction(apoptotic_particles, measured_particles))
    )
    quality_pass = bool(
        purity >= (1.0 - params.contamination_threshold)
        and viability >= params.viability_threshold
    )
    return {
        "viability_fraction": float(viability),
        "stressed_viable_fraction": float(stressed_viable),
        "apoptosis_fraction": float(apoptosis),
        "necrosis_fraction": float(necrosis),
        "debris_particles": float(debris_particles),
        "aggregate_particles": float(aggregate_particles),
        "measured_particles": float(measured_particles),
        "bona_fide_EV_fraction": float(bona_fide_fraction),
        "marker_panel_score": float(marker_score),
        "purity_score": float(purity),
        "quality_pass": quality_pass,
    }


def _hill(value: float, half_max: float, coefficient: float) -> float:
    value = max(float(value), 0.0)
    half_max = max(float(half_max), 1e-12)
    coefficient = max(float(coefficient), 1e-12)
    numerator = value**coefficient
    return numerator / (half_max**coefficient + numerator) if value > 0 else 0.0


def _fraction(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if denominator > 0 else 0.0


def _clip01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))



def injury_defaults() -> dict:
    """Return injury defaults as a dict."""
    return asdict(InjuryParams())
