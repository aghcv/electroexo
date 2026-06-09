from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable


@dataclass(slots=True)
class InjuryParams:
    K_apoptosis_damage: float = 0.5
    K_necrosis_damage: float = 0.9
    n_apoptosis: float = 3.0
    n_necrosis: float = 5.0
    debris_fraction_scale: float = 0.2
    contamination_threshold: float = 0.3
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



def compute_quality_gate(damage: float, debris: float, params: InjuryParams) -> dict[str, float | bool]:
    """Compute placeholder viability, injury fractions, and pass/fail quality metrics."""
    apoptosis = (damage ** params.n_apoptosis) / (params.K_apoptosis_damage ** params.n_apoptosis + damage ** params.n_apoptosis) if damage > 0 else 0.0
    necrosis = (damage ** params.n_necrosis) / (params.K_necrosis_damage ** params.n_necrosis + damage ** params.n_necrosis) if damage > 0 else 0.0
    apoptosis = min(apoptosis, 0.8)
    necrosis = min(necrosis, 0.9)
    viability = max(0.0, 1.0 - 0.6 * apoptosis - 0.9 * necrosis)
    purity = max(0.0, 1.0 - debris)
    quality_pass = bool(purity >= (1.0 - params.contamination_threshold) and viability >= 0.5)
    return {
        "viability_fraction": float(viability),
        "apoptosis_fraction": float(apoptosis),
        "necrosis_fraction": float(necrosis),
        "purity_score": float(purity),
        "quality_pass": quality_pass,
    }



def injury_defaults() -> dict:
    """Return injury defaults as a dict."""
    return asdict(InjuryParams())
