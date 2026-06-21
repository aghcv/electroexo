from __future__ import annotations

import pytest

from electro_exocytosis.models.electrodynamics import ElectrodynamicsState
from electro_exocytosis.models.ion_transport import (
    IonTransportParams,
    build_ion_transport_rhs,
    coerce_ion_transport_params,
    compute_ion_transport_fluxes,
    get_ion_initial_conditions,
    get_ion_state_names,
)


def _electro_state(permeability: float = 1.0) -> ElectrodynamicsState:
    return ElectrodynamicsState(
        delta_Vm=1.0,
        delta_V_ER=0.3,
        delta_V_mito=0.2,
        delta_V_MVB=0.1,
        pore_density=1.0e12,
        membrane_permeability=permeability,
    )


def test_ion_transport_state_names_cover_table_a5_submodules() -> None:
    names = get_ion_state_names()

    assert {
        "Ca_i",
        "Ca_ER",
        "Ca_mito",
        "mitochondrial_potential",
        "ROS",
        "ATP",
        "Na_i",
        "K_i",
        "Cl_i",
        "osmotic_stress",
    }.issubset(names)


def test_pore_fluxes_follow_permeability_and_extracellular_calcium() -> None:
    params = IonTransportParams()
    low_ca_params = IonTransportParams(Ca_ext_uM=0.1)
    y0 = get_ion_initial_conditions(params)

    open_fluxes = compute_ion_transport_fluxes(params, _electro_state(1.0), 0.0, y0, t_pulse_end=1.0)
    closed_fluxes = compute_ion_transport_fluxes(params, _electro_state(0.0), 0.0, y0, t_pulse_end=1.0)
    low_ca_fluxes = compute_ion_transport_fluxes(low_ca_params, _electro_state(1.0), 0.0, y0, t_pulse_end=1.0)

    assert open_fluxes["J_Ca_pore"] > 0.0
    assert open_fluxes["J_Na_pore"] > 0.0
    assert open_fluxes["J_K_pore"] > 0.0
    assert closed_fluxes["J_Ca_pore"] == pytest.approx(0.0)
    assert low_ca_fluxes["J_Ca_pore"] < open_fluxes["J_Ca_pore"]


def test_rhs_advances_ca_ions_and_mitochondrial_stress() -> None:
    params = IonTransportParams()
    rhs = build_ion_transport_rhs(params, _electro_state(1.0), t_pulse_end=1.0)
    dydt = rhs(0.0, get_ion_initial_conditions(params))
    derivatives = dict(zip(get_ion_state_names(), dydt, strict=True))

    assert derivatives["Ca_i"] > 0.0
    assert derivatives["Ca_ER"] < 0.0
    assert derivatives["Na_i"] > 0.0
    assert derivatives["K_i"] < 0.0
    assert derivatives["mitochondrial_potential"] < 0.0


def test_nested_ion_transport_params_are_coerced_to_float() -> None:
    params = coerce_ion_transport_params({"ion_transport": {"Ca_ext_uM": "5.0", "tau_ion_recovery_s": "1200.0"}})

    assert params.Ca_ext_uM == pytest.approx(5.0)
    assert params.tau_ion_recovery_s == pytest.approx(1200.0)
