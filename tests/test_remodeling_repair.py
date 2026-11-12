from __future__ import annotations

import pytest

from electro_exocytosis.models.remodeling_repair import (
    RemodelingParams,
    coerce_remodeling_params,
    compute_remodeling_state,
)


def test_remodeling_state_tracks_microdomain_and_repair_submodules() -> None:
    low = compute_remodeling_state(0.2, RemodelingParams(), pore_activation=0.0)
    high = compute_remodeling_state(2.0, RemodelingParams(), pore_activation=1.0, osmotic_stress=0.05)

    assert high["Ca_submembrane"] > high["Ca_i"] if "Ca_i" in high else high["Ca_submembrane"] > 2.0
    assert high["PS_exposure"] > low["PS_exposure"]
    assert high["calpain_activity"] > low["calpain_activity"]
    assert high["annexin_activity"] > low["annexin_activity"]
    assert high["lysosomal_repair_activity"] > low["lysosomal_repair_activity"]
    assert high["repair_shedding_rate"] > low["repair_shedding_rate"]


def test_calpain_inhibition_reduces_actin_disruption() -> None:
    baseline = compute_remodeling_state(2.0, RemodelingParams(), pore_activation=1.0)
    inhibited = compute_remodeling_state(
        2.0,
        RemodelingParams(K_calpain_uM=8.0, actin_calpain_weight=0.2),
        pore_activation=1.0,
    )

    assert inhibited["calpain_activity"] < baseline["calpain_activity"]
    assert inhibited["actin_disruption"] < baseline["actin_disruption"]


def test_nested_remodeling_params_are_coerced_to_float() -> None:
    params = coerce_remodeling_params({"remodeling_repair": {"microdomain_gain": "0.5", "K_calpain_uM": "8.0"}})

    assert params.microdomain_gain == pytest.approx(0.5)
    assert params.K_calpain_uM == pytest.approx(8.0)
