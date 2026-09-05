from __future__ import annotations

import pytest

from electro_exocytosis.models.ev_release import (
    EVReleaseParams,
    coerce_ev_release_params,
    compute_ev_release_fluxes,
    compute_ev_release_derivatives,
    get_ev_initial_conditions,
)


def test_calcium_and_rab_drivers_raise_small_ev_release() -> None:
    params = EVReleaseParams()
    state = get_ev_initial_conditions(params)

    low = compute_ev_release_fluxes(
        params,
        state,
        Ca_i=0.1,
        Ca_submembrane=0.2,
        ROS=0.1,
        ATP=0.8,
        damage_state=0.05,
        delta_V_MVB=0.01,
        pore_activation=0.1,
        PS_exposure=0.1,
        calpain_activity=0.1,
        annexin_activity=0.1,
        actomyosin_tension=0.2,
        actin_disruption=0.1,
        repair_state=0.1,
        repair_shedding_rate=0.0,
    )
    high = compute_ev_release_fluxes(
        params,
        state,
        Ca_i=1.5,
        Ca_submembrane=2.0,
        ROS=0.15,
        ATP=1.0,
        damage_state=0.05,
        delta_V_MVB=0.08,
        pore_activation=0.8,
        PS_exposure=0.2,
        calpain_activity=0.2,
        annexin_activity=0.3,
        actomyosin_tension=0.15,
        actin_disruption=0.5,
        repair_state=0.4,
        repair_shedding_rate=0.0,
    )

    assert high["rab_conversion_signal"] > low["rab_conversion_signal"]
    assert high["rab_docking_signal"] > low["rab_docking_signal"]
    assert high["fusion_signal"] > low["fusion_signal"]
    assert high["sEV_rate"] > low["sEV_rate"]


def test_ps_tension_and_calpain_raise_medium_large_ev_release() -> None:
    params = EVReleaseParams()
    state = get_ev_initial_conditions(params)

    low = compute_ev_release_fluxes(
        params,
        state,
        Ca_i=0.3,
        Ca_submembrane=0.5,
        ROS=0.1,
        ATP=0.9,
        damage_state=0.05,
        delta_V_MVB=0.02,
        pore_activation=0.2,
        PS_exposure=0.05,
        calpain_activity=0.05,
        annexin_activity=0.1,
        actomyosin_tension=0.1,
        actin_disruption=0.1,
        repair_state=0.05,
        repair_shedding_rate=0.0,
    )
    high = compute_ev_release_fluxes(
        params,
        state,
        Ca_i=0.8,
        Ca_submembrane=1.0,
        ROS=0.2,
        ATP=0.9,
        damage_state=0.2,
        delta_V_MVB=0.03,
        pore_activation=0.4,
        PS_exposure=0.8,
        calpain_activity=0.9,
        annexin_activity=0.3,
        actomyosin_tension=0.9,
        actin_disruption=0.6,
        repair_state=0.6,
        repair_shedding_rate=0.05,
    )

    assert high["budding_signal"] > low["budding_signal"]
    assert high["scission_signal"] > low["scission_signal"]
    assert high["mlEV_rate"] > low["mlEV_rate"]


def test_damage_ros_and_atp_loss_raise_apoptotic_and_lysosomal_signals() -> None:
    params = EVReleaseParams()
    state = get_ev_initial_conditions(params)

    low = compute_ev_release_fluxes(
        params,
        state,
        Ca_i=0.4,
        Ca_submembrane=0.6,
        ROS=0.1,
        ATP=1.0,
        damage_state=0.05,
        delta_V_MVB=0.02,
        pore_activation=0.2,
        PS_exposure=0.1,
        calpain_activity=0.1,
        annexin_activity=0.1,
        actomyosin_tension=0.1,
        actin_disruption=0.1,
        repair_state=0.2,
        repair_shedding_rate=0.0,
    )
    high = compute_ev_release_fluxes(
        params,
        state,
        Ca_i=1.0,
        Ca_submembrane=1.2,
        ROS=0.8,
        ATP=0.15,
        damage_state=1.5,
        delta_V_MVB=0.05,
        pore_activation=0.7,
        PS_exposure=0.6,
        calpain_activity=0.5,
        annexin_activity=0.2,
        actomyosin_tension=0.7,
        actin_disruption=0.3,
        repair_state=0.1,
        repair_shedding_rate=0.02,
    )

    assert high["lysosomal_routing"] > low["lysosomal_routing"]
    assert high["apoptotic_blebbing_signal"] > low["apoptotic_blebbing_signal"]
    assert high["AB_rate"] > low["AB_rate"]


def test_nested_ev_release_params_are_coerced_to_float() -> None:
    params = coerce_ev_release_params(
        {
            "ev_release": {
                "k_MVB_maturation_s": "0.004",
                "ceramide_baseline": "0.55",
            }
        }
    )

    assert params.k_MVB_maturation_s == pytest.approx(0.004)
    assert params.ceramide_baseline == pytest.approx(0.55)


def test_voltage_stimulus_is_not_sustained_after_pore_resealing() -> None:
    params = EVReleaseParams()
    state = get_ev_initial_conditions(params)
    common = {
        "Ca_i": 0.1,
        "Ca_submembrane": 0.1,
        "ROS": 0.1,
        "ATP": 1.0,
        "damage_state": 0.0,
        "PS_exposure": 0.0,
        "calpain_activity": 0.0,
        "annexin_activity": 0.0,
        "actomyosin_tension": 0.0,
        "actin_disruption": 0.0,
        "repair_state": 0.0,
        "repair_shedding_rate": 0.0,
    }
    resealed = compute_ev_release_fluxes(
        params,
        state,
        delta_V_MVB=0.2,
        pore_activation=0.0,
        **common,
    )
    no_voltage = compute_ev_release_fluxes(
        params,
        state,
        delta_V_MVB=0.0,
        pore_activation=0.0,
        **common,
    )

    assert resealed["rab_conversion_signal"] == pytest.approx(
        no_voltage["rab_conversion_signal"]
    )


def test_homeostatic_reference_balances_baseline_ev_pools() -> None:
    params = EVReleaseParams()
    state = get_ev_initial_conditions(params)
    reference = compute_ev_release_fluxes(
        params,
        state,
        Ca_i=0.1,
        Ca_submembrane=0.1,
        ROS=0.1,
        ATP=1.0,
        damage_state=0.0,
        delta_V_MVB=0.0,
        pore_activation=0.0,
        PS_exposure=0.0,
        calpain_activity=0.0,
        annexin_activity=0.0,
        actomyosin_tension=0.0,
        actin_disruption=0.0,
        repair_state=0.0,
        repair_shedding_rate=0.0,
    )
    derivatives = compute_ev_release_derivatives(
        params,
        state,
        reference,
        homeostatic_reference_fluxes=reference,
    )

    assert derivatives == pytest.approx([0.0] * len(state), abs=1.0e-15)
