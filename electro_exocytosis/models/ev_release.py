from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable

from electro_exocytosis.config import CellStateConfig


@dataclass(slots=True)
class EVReleaseParams:
    baseline_sEV_rate: float = 1.0
    baseline_mlEV_rate: float = 0.1
    baseline_AB_rate: float = 0.01
    K_sEV_Ca_uM: float = 0.5
    n_sEV_Ca: float = 1.5
    K_mlEV_PS: float = 0.3
    n_mlEV_PS: float = 2.0
    K_AB_damage: float = 0.7
    n_AB_damage: float = 3.0
    sEV_Ca_scale: float = 5.0
    mlEV_PS_scale: float = 3.0
    AB_damage_scale: float = 10.0


# TODO-literature-review: replace placeholder EV subtype release kinetics with curated literature-supported rate laws.
def build_ev_release_rhs(params: EVReleaseParams, cell_state: CellStateConfig) -> Callable[[float, list[float], float, float, float], list[float]]:
    """Build a placeholder Layer 5 EV release RHS using algebraic Ca/PS/damage couplings."""

    def rhs(t: float, y: list[float], Ca_i: float, PS_exposure: float, damage_state: float) -> list[float]:
        _, _, _ = y
        s_hill = (Ca_i ** params.n_sEV_Ca) / (params.K_sEV_Ca_uM ** params.n_sEV_Ca + Ca_i ** params.n_sEV_Ca)
        m_hill = (PS_exposure ** params.n_mlEV_PS) / (params.K_mlEV_PS ** params.n_mlEV_PS + PS_exposure ** params.n_mlEV_PS)
        d_hill = (damage_state ** params.n_AB_damage) / (params.K_AB_damage ** params.n_AB_damage + damage_state ** params.n_AB_damage)
        sev = params.baseline_sEV_rate * cell_state.baseline_EV_release_modifier * (1.0 + params.sEV_Ca_scale * s_hill)
        mlev = params.baseline_mlEV_rate * cell_state.baseline_EV_release_modifier * (1.0 + params.mlEV_PS_scale * m_hill)
        ab = params.baseline_AB_rate * (1.0 + params.AB_damage_scale * cell_state.stress_sensitivity_modifier * d_hill)
        return [float(sev), float(mlev), float(ab)]

    return rhs



def ev_release_defaults() -> dict:
    """Return EV release defaults as a dict."""
    return asdict(EVReleaseParams())
