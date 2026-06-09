from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electro_exocytosis.config import CellStateConfig
from electro_exocytosis.models.dosimetry import DosimetryResult
from electro_exocytosis.models.pulse import PulseDescriptors


@dataclass(slots=True)
class ElectrodynamicsState:
    delta_Vm: float
    delta_V_ER: float
    delta_V_mito: float
    delta_V_MVB: float
    pore_density: float
    membrane_permeability: float


# TODO-literature-review: replace simplified Schwan scaling and permeability equations with validated membrane/organelle models.
def compute_electrodynamics_state(
    descriptors: PulseDescriptors,
    dosimetry: DosimetryResult,
    cell_state: CellStateConfig,
) -> ElectrodynamicsState:
    """Compute Layer 2, Tables A3-A4 electrodynamic perturbation state."""
    # TODO-literature-review: cell_radius_m is currently a fixed placeholder (10 µm, typical mammalian cell).
    # Future versions should accept cell_radius_m from the parameter set (see default_parameters.yaml key
    # 'cell_radius_m') so cell-type-specific sizes can be configured without code changes.
    cell_radius_m = 10e-6
    delta_vm = dosimetry.effective_E_V_m * cell_radius_m * 1.5 * cell_state.membrane_modifier
    delta_v_er = 0.3 * delta_vm
    delta_v_mito = 0.2 * delta_vm
    delta_v_mvb = 0.15 * delta_vm
    threshold_v = 0.25
    slope = 10.0
    membrane_permeability = float(1.0 / (1.0 + np.exp(-slope * (delta_vm - threshold_v))))
    pore_density = 1e12 * (delta_vm ** 3) / (0.5 ** 3 + delta_vm ** 3)
    return ElectrodynamicsState(
        delta_Vm=delta_vm,
        delta_V_ER=delta_v_er,
        delta_V_mito=delta_v_mito,
        delta_V_MVB=delta_v_mvb,
        pore_density=float(pore_density),
        membrane_permeability=float(np.clip(membrane_permeability, 0.0, 1.0)),
    )
