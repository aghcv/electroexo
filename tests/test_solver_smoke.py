from __future__ import annotations

from pathlib import Path

from electro_exocytosis.io.readers import load_scenario
from electro_exocytosis.simulation import Simulation



def test_solver_smoke() -> None:
    scenario = load_scenario(Path("examples/scenario_baseline.yaml"))
    result = Simulation(scenario).run()
    assert len(result.t_array) > 0
    assert not result.state_timeseries.empty
    assert not result.ev_timeseries.empty
    for key in [
        "scenario_name",
        "mode",
        "dose_index",
        "peak_ca_i",
        "peak_ca_mito",
        "peak_ros",
        "min_atp",
        "min_mitochondrial_potential",
        "peak_osmotic_stress",
        "peak_ps_exposure",
        "peak_calpain_activity",
        "peak_annexin_activity",
        "peak_repair_state",
        "peak_secretory_bias",
        "peak_docked_MVB_pool",
        "cumulative_repair_shedding",
        "cumulative_small_EV",
        "purity_score",
        "protein_enrichment",
        "RNA_enrichment",
        "lipid_enrichment",
        "antigen_enrichment",
        "bona_fide_EV_fraction",
        "cell_normalized_yield",
        "optimization_objective",
        "warnings",
    ]:
        assert key in result.summary
    for column in [
        "Ca_mito",
        "mitochondrial_potential",
        "Na_i",
        "K_i",
        "Cl_i",
        "osmotic_stress",
        "J_Ca_pore",
        "J_ER_release",
        "J_Na_pore",
        "J_K_pore",
        "J_Cl_pore",
        "Ca_submembrane",
        "PS_exposure",
        "scramblase_activity",
        "calpain_activity",
        "annexin_activity",
        "lysosomal_repair_activity",
        "repair_state",
        "repair_shedding_rate",
        "MVB_pool",
        "ILV_load",
        "docked_MVB_pool",
        "budding_pool",
        "apoptotic_commitment",
        "secretory_bias",
        "lysosomal_routing",
        "fusion_signal",
    ]:
        assert column in result.state_timeseries
