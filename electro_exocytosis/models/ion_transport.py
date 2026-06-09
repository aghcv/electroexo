from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable

import numpy as np

from electro_exocytosis.models.electrodynamics import ElectrodynamicsState


@dataclass(slots=True)
class IonTransportParams:
    tau_Ca_release_s: float = 0.5
    tau_Ca_uptake_s: float = 30.0
    Ca_baseline_uM: float = 0.1
    Ca_ER_uM: float = 500.0
    Ca_max_uM: float = 10.0
    J_Ca_pore_factor: float = 0.01
    ROS_baseline: float = 0.1
    ROS_production_factor: float = 0.001  # TODO-literature-review: set to keep ROS in ~1–10x baseline range during nsPEF stress
    tau_ROS_s: float = 120.0
    ATP_baseline: float = 1.0
    ATP_depletion_factor: float = 0.3
    tau_ATP_s: float = 600.0


# TODO-literature-review: replace release, uptake, ROS, and ATP couplings with literature-supported flux laws.
def build_ion_transport_rhs(
    params: IonTransportParams,
    electro_state: ElectrodynamicsState,
    t_pulse_end: float,
) -> Callable[[float, list[float]], list[float]]:
    """Build a placeholder Layer 3 RHS for Ca2+, ROS, and ATP dynamics."""

    def rhs(t: float, y: list[float]) -> list[float]:
        ca_i, ca_er, ros, atp = y
        local_t = max(t - t_pulse_end, 0.0)
        j_ca_in = (
            params.J_Ca_pore_factor
            * electro_state.membrane_permeability
            * max(ca_er - ca_i, 0.0)
            * np.exp(-local_t / max(params.tau_Ca_release_s, 1e-9))
        )
        j_ca_out = max(ca_i - params.Ca_baseline_uM, 0.0) / params.tau_Ca_uptake_s
        d_ca_i = j_ca_in - j_ca_out
        d_ca_er = -j_ca_in
        d_ros = params.ROS_production_factor * max(ca_i / params.Ca_baseline_uM - 1.0, 0.0) - (ros - params.ROS_baseline) / params.tau_ROS_s
        d_atp = -params.ATP_depletion_factor * max(ros - params.ROS_baseline, 0.0) - (atp - params.ATP_baseline) / params.tau_ATP_s
        return [float(d_ca_i), float(d_ca_er), float(d_ros), float(d_atp)]

    return rhs



def get_ion_initial_conditions(params: IonTransportParams) -> list[float]:
    """Return initial conditions [Ca_i, Ca_ER, ROS, ATP]."""
    return [params.Ca_baseline_uM, params.Ca_ER_uM, params.ROS_baseline, params.ATP_baseline]



def ion_transport_defaults() -> dict:
    """Return ion transport defaults as a plain dict."""
    return asdict(IonTransportParams())
