from __future__ import annotations

from copy import deepcopy

from electro_exocytosis.config import CellStateConfig


# TODO-literature-review: replace broad modifier scaling with cell-type and disease-state specific evidence.
def apply_cell_state_modifiers(base_params: dict, cell_state: CellStateConfig) -> dict:
    """Apply cross-cutting cell-state modifiers to a flat parameter dictionary."""
    params = deepcopy(base_params)
    for key in [
        "J_Ca_pore_factor",
        "J_ER_release_factor",
        "SERCA_Vmax_uM_s",
        "PMCA_Vmax_uM_s",
        "NCX_Vmax_uM_s",
        "mitochondrial_uptake_Vmax_uM_s",
        "Ca_max_uM",
    ]:
        if key in params:
            params[key] = float(params[key]) * cell_state.calcium_handling_modifier
    for key in ["baseline_sEV_rate", "baseline_mlEV_rate", "baseline_AB_rate"]:
        if key in params:
            params[key] = float(params[key]) * cell_state.baseline_EV_release_modifier
    for key in [
        "damage_rate",
        "AB_damage_scale",
        "ROS_production_factor",
        "ROS_mito_factor_s",
        "ROS_depolarization_factor_s",
        "ROS_osmotic_factor_s",
        "ATP_depletion_factor",
    ]:
        if key in params:
            params[key] = float(params[key]) * cell_state.stress_sensitivity_modifier
    return params
