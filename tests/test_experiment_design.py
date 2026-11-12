from __future__ import annotations

import numpy as np
import pandas as pd

from electro_exocytosis.experiment_design import (
    build_parameter_inventory,
    build_state_inventory,
    local_fisher_diagnostics,
    morris_screen,
    priority_scores,
    sobol_saltelli,
)


def test_parameter_and_state_inventories_remain_separate() -> None:
    parameters = build_parameter_inventory()
    states = build_state_inventory()

    assert len(parameters) == 241
    assert parameters["parameter_id"].is_unique
    assert states["state_id"].is_unique
    assert "state.ion_transport.Ca_i" in set(states["state_id"])
    assert "parameter.ion_transport.Ca_baseline_uM" in set(parameters["parameter_id"])
    assert set(parameters["parameter_id"]).isdisjoint(states["state_id"])
    inactive = parameters.set_index("qualified_name").loc[
        "remodeling_repair.tau_repair_s"
    ]
    assert inactive["evidence_confidence"] == "inactive in current equations"

    fitted = build_parameter_inventory(
        fit_parameters=pd.DataFrame(
            [
                {
                    "variant": "selected",
                    "module": "Source conversion",
                    "parameter": "scale",
                    "fitted": 2.0,
                }
            ]
        ),
        selected_fit_variant="selected",
    )
    assert "parameter.ev_size_observation.scale" in set(fitted["parameter_id"])


def test_morris_identifies_linear_effects_reproducibly() -> None:
    def model(x: np.ndarray) -> np.ndarray:
        return np.asarray([2.0 * x[0] + 0.5 * x[1]])

    first, raw_first = morris_screen(
        model,
        [0.0, 0.0],
        [1.0, 1.0],
        trajectories=6,
        seed=12,
        parameter_ids=["a", "b"],
        target_ids=["y"],
    )
    second, raw_second = morris_screen(
        model,
        [0.0, 0.0],
        [1.0, 1.0],
        trajectories=6,
        seed=12,
        parameter_ids=["a", "b"],
        target_ids=["y"],
    )

    assert raw_first.equals(raw_second)
    effects = first.set_index("parameter_id")
    assert effects.loc["a", "morris_mu_star"] == 2.0
    assert effects.loc["b", "morris_mu_star"] == 0.5
    assert np.allclose(effects["morris_sigma"], 0.0)


def test_sobol_distinguishes_dominant_additive_parameter() -> None:
    result = sobol_saltelli(
        lambda x: np.asarray([x[0] + 0.1 * x[1]]),
        [0.0, 0.0],
        [1.0, 1.0],
        base_samples=1024,
        seed=9,
        parameter_ids=["a", "b"],
        target_ids=["y"],
    ).set_index("parameter_id")

    assert result.loc["a", "sobol_ST"] > 0.9
    assert result.loc["b", "sobol_ST"] < 0.03
    assert "Monte Carlo" in result.loc["a", "estimator"]


def test_local_fisher_and_priority_components_are_explicit() -> None:
    diagnostics, correlation, summary = local_fisher_diagnostics(
        lambda x: np.asarray([x[0] + x[1], x[0] - x[1]]),
        [1.0, 1.0],
        [0.0, 0.0],
        [2.0, 2.0],
        parameter_ids=["a", "b"],
    )
    assert summary["effective_rank"] == 2
    assert np.allclose(np.diag(correlation), 1.0)
    assert set(diagnostics["identifiability_status"]) == {"weak"}

    ranked = priority_scores(
        diagnostics.assign(
            sensitivity=[1.0, 0.5],
            uncertainty=1.0,
            identifiability_gap=[0.8, 0.2],
            measurability=0.8,
            feasibility=0.9,
        )
    )
    assert ranked.iloc[0]["parameter_id"] == "a"
    assert "geometric" in ranked.iloc[0]["decision_metric"]
