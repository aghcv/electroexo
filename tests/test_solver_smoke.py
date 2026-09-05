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
        "terminal_extracellular_total_concentration_particles_per_ml",
        "peak_extracellular_total_concentration_particles_per_ml",
        "extracellular_decline_from_peak_fraction",
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

    for column in [
        "sEV_extracellular_concentration_particles_per_ml",
        "mlEV_extracellular_concentration_particles_per_ml",
        "AB_extracellular_concentration_particles_per_ml",
        "total_extracellular_concentration_particles_per_ml",
        "measured_particle_concentration_particles_per_ml",
        "extracellular_medium_volume_ml",
        "viable_producer_fraction",
    ]:
        assert column in result.ev_timeseries


def test_simulation_extracellular_stock_can_decline() -> None:
    scenario = load_scenario(Path("examples/scenario_baseline.yaml"))
    scenario.extracellular_medium.use_time_varying_viability = False
    result = Simulation(
        scenario,
        params_override={
            "extracellular_kinetics": {
                "source_scale_particles_per_model_unit": 0.0,
                "initial_sEV_concentration_particles_per_ml": 1.0e9,
                "degradation_rate_s": 0.0002,
            }
        },
    ).run()

    concentration = result.ev_timeseries[
        "total_extracellular_concentration_particles_per_ml"
    ].to_numpy()
    assert concentration[-1] < concentration[0]
    assert result.summary["extracellular_decline_from_peak_fraction"] > 0.5
