from pathlib import Path

import numpy as np

from electro_exocytosis.experimental_bridge import ExperimentalObservationBridge
from tools.fit_ffrci_normalized_ev_kinetics import (
    CONTROL_ORDER,
    load_normalized_observations,
)


DATA_DIR = (
    Path(__file__).resolve().parents[1] / "data" / "experimental" / "ffrci_data_sharing"
)


def test_normalized_loader_keeps_controls_separate() -> None:
    observations, controls = load_normalized_observations(
        DATA_DIR, max_particle_diameter_nm=200.0
    )

    assert set(observations["control"]) == set(CONTROL_ORDER)
    assert len(observations) == 18
    assert len(controls) == 2
    assert observations["observed_fold_of_sham"].gt(0.0).all()
    assert observations.groupby("control")["sham_concentration"].nunique().eq(1).all()
    assert not np.isclose(
        observations.loc[
            observations["control"] == "sham2", "sham_concentration"
        ].iloc[0],
        observations.loc[
            observations["control"] == "sham_media", "sham_concentration"
        ].iloc[0],
    )


def test_single_cell_conversion_preserves_fold_change_with_shared_bridge() -> None:
    bridge = ExperimentalObservationBridge(
        initial_cell_count=5.0e6, medium_volume_ml=5.0
    )
    raw, _ = load_normalized_observations(
        DATA_DIR, max_particle_diameter_nm=200.0
    )
    per_cell, _ = load_normalized_observations(
        DATA_DIR,
        max_particle_diameter_nm=200.0,
        observation_bridge=bridge,
    )

    assert np.allclose(
        raw["observed_fold_of_sham"], per_cell["observed_fold_of_sham"]
    )
    assert np.allclose(
        per_cell["stimulated_particles_per_cell"],
        per_cell["stimulated_concentration"] / 1.0e6,
    )
