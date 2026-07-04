from __future__ import annotations

import pytest

from electro_exocytosis.models.cargo_potency import CargoPotencyParams, compute_cargo_state
from electro_exocytosis.models.injury_quality import InjuryParams, compute_quality_gate
from electro_exocytosis.models.manufacturing_qc import ManufacturingParams, compute_manufacturing_outputs


def test_cargo_state_tracks_stress_sorting_and_direct_loading() -> None:
    params = CargoPotencyParams()
    low = compute_cargo_state(
        0.1,
        0.1,
        1.0,
        params,
        cumulative_sEV=1.0,
        cumulative_mlEV=0.1,
        cumulative_AB=0.0,
        escrt_signal=0.1,
        ceramide_signal=0.1,
        secretory_bias=0.2,
    )
    high = compute_cargo_state(
        1.2,
        0.8,
        0.4,
        params,
        cumulative_sEV=2.0,
        cumulative_mlEV=0.3,
        cumulative_AB=0.05,
        escrt_signal=0.8,
        ceramide_signal=0.8,
        secretory_bias=0.9,
        direct_loading_drive=1.0,
    )

    assert high["protein_enrichment"] > low["protein_enrichment"]
    assert high["RNA_enrichment"] > low["RNA_enrichment"]
    assert high["lipid_enrichment"] > low["lipid_enrichment"]
    assert high["antigen_enrichment"] > low["antigen_enrichment"]
    assert high["direct_loaded_cargo"] > low["direct_loaded_cargo"]
    assert high["potency_score"] > low["potency_score"]


def test_quality_gate_penalizes_apoptotic_and_necrotic_mixtures() -> None:
    params = InjuryParams()
    clean = compute_quality_gate(
        0.05,
        0.01,
        params,
        cumulative_sEV=2.0,
        cumulative_mlEV=0.2,
        cumulative_AB=0.01,
    )
    damaged = compute_quality_gate(
        1.2,
        0.3,
        params,
        cumulative_sEV=2.0,
        cumulative_mlEV=0.2,
        cumulative_AB=0.8,
    )

    assert damaged["apoptosis_fraction"] > clean["apoptosis_fraction"]
    assert damaged["necrosis_fraction"] > clean["necrosis_fraction"]
    assert damaged["measured_particles"] > clean["measured_particles"]
    assert damaged["purity_score"] < clean["purity_score"]
    assert damaged["viability_fraction"] < clean["viability_fraction"]


def test_manufacturing_outputs_reflect_recovery_purity_and_objective() -> None:
    params = ManufacturingParams(isolation_efficiency=0.5, isolation_method_factor=0.8)
    baseline = compute_manufacturing_outputs(
        1.0,
        0.2,
        0.05,
        0.9,
        params,
        potency_score=0.5,
        purity_score=0.8,
    )
    stronger = compute_manufacturing_outputs(
        2.0,
        0.4,
        0.05,
        0.95,
        params,
        potency_score=0.9,
        purity_score=0.9,
    )

    assert baseline["process_recovery"] == pytest.approx(0.4)
    assert stronger["total_measured_particles"] > baseline["total_measured_particles"]
    assert stronger["cell_normalized_yield"] > baseline["cell_normalized_yield"]
    assert stronger["optimization_objective"] > baseline["optimization_objective"]
