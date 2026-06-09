from __future__ import annotations

from copy import deepcopy

from electro_exocytosis.config import CellStateConfig


# TODO-literature-review: replace broad modifier scaling with cell-type and disease-state specific evidence.
def apply_cell_state_modifiers(base_params: dict, cell_state: CellStateConfig) -> dict:
    """Apply cross-cutting cell-state modifiers to a flat parameter dictionary."""
    params = deepcopy(base_params)
    for key in ["J_Ca_pore_factor", "tau_Ca_release_s", "tau_Ca_uptake_s", "Ca_max_uM"]:
        if key in params:
            if key == "tau_Ca_uptake_s":
                params[key] = params[key] / max(cell_state.calcium_handling_modifier, 1e-6)
            else:
                params[key] = params[key] * cell_state.calcium_handling_modifier
    for key in ["baseline_sEV_rate", "baseline_mlEV_rate", "baseline_AB_rate"]:
        if key in params:
            params[key] = params[key] * cell_state.baseline_EV_release_modifier
    for key in ["damage_rate", "AB_damage_scale", "ROS_production_factor"]:
        if key in params:
            params[key] = params[key] * cell_state.stress_sensitivity_modifier
    return params
