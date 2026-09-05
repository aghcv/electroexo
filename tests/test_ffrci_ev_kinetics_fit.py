from pathlib import Path

import numpy as np

from electro_exocytosis.experimental_bridge import ExperimentalObservationBridge
from tools.fit_ffrci_ev_kinetics import (
    CurrentModelPredictor,
    load_longitudinal_observations,
    load_replicate_observations,
)


DATA_DIR = (
    Path(__file__).resolve().parents[1] / "data" / "experimental" / "ffrci_data_sharing"
)


def test_longitudinal_loader_retains_nominal_sample_sets() -> None:
    replicates = load_replicate_observations(DATA_DIR)
    joint = load_longitudinal_observations(DATA_DIR)

    assert len(replicates) == 26
    assert len(joint) == 9
    group_sizes = replicates.groupby(["condition", "replicate"]).size()
    assert len(group_sizes) == 9
    assert (group_sizes == 3).sum() == 8
    assert group_sizes.loc[("3p40kV", 3)] == 2


def test_longitudinal_loader_can_filter_to_particles_below_200_nm() -> None:
    all_sizes = load_replicate_observations(DATA_DIR)
    below_200 = load_replicate_observations(
        DATA_DIR, max_particle_diameter_nm=200.0
    )

    assert len(below_200) == len(all_sizes) == 26
    assert (below_200["observed_concentration"] < all_sizes["observed_concentration"]).all()
    assert (below_200["particle_diameter_filter"] == "diameter_lt_200_nm").all()

    joint = load_longitudinal_observations(
        DATA_DIR, max_particle_diameter_nm=200.0
    )
    assert len(joint) == 9
    assert (joint["particle_diameter_filter"] == "diameter_lt_200_nm").all()


def test_loader_bridges_particles_per_ml_to_per_cell() -> None:
    bridge = ExperimentalObservationBridge(
        initial_cell_count=5.0e6, medium_volume_ml=5.0
    )
    observations = load_replicate_observations(
        DATA_DIR,
        max_particle_diameter_nm=200.0,
        observation_bridge=bridge,
    )

    expected = observations["measured_concentration_particles_per_ml"] / 1.0e6
    assert np.allclose(observations["observed_concentration"], expected)
    assert (observations["observation_unit"] == "particle_equivalents_per_cell").all()


def test_baseline_release_rates_scale_cached_current_model_prediction() -> None:
    observations = load_longitudinal_observations(DATA_DIR)
    observations = observations[observations["condition"] == "1p20kV"].reset_index(
        drop=True
    )
    predictor = CurrentModelPredictor(
        observations,
        pulse_width_ns=60.0,
        repetition_rate_hz=1.0,
    )

    default_prediction = predictor.predict(predictor.defaults)
    runs_after_default = predictor.simulator_runs
    doubled_rates = dict(predictor.defaults)
    for name in ("baseline_sEV_rate", "baseline_mlEV_rate", "baseline_AB_rate"):
        doubled_rates[name] *= 2.0
    doubled_prediction = predictor.predict(doubled_rates)

    assert np.all(default_prediction > 0.0)
    assert np.allclose(doubled_prediction, 2.0 * default_prediction)
    assert predictor.simulator_runs == runs_after_default
