"""
electroexo.results
===================
Container for simulation output.

:class:`SimulationResults` stores the full time-course trajectories of all
state variables together with derived quantities and convenience properties.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

__all__ = ["SimulationResults"]


@dataclass
class SimulationResults:
    """Output from a completed :class:`~electroexo.simulation.SimulationEngine` run.

    Attributes
    ----------
    t : ndarray, shape (n,)
        Simulation time points (ms).
    v : ndarray, shape (n,)
        Membrane potential trajectory (mV).
    ca : ndarray, shape (n,)
        Cytosolic Ca²⁺ concentration trajectory (µM).
    docked : ndarray, shape (n,)
        Docked vesicle pool trajectory (vesicle count).
    reserve : ndarray, shape (n,)
        Reserve pool trajectory (vesicle count).
    released : ndarray, shape (n,)
        Cumulative exocytosis events trajectory (vesicle count).
    m : ndarray, shape (n,)
        Na⁺ channel activation gate.
    h : ndarray, shape (n,)
        Na⁺ channel inactivation gate.
    n : ndarray, shape (n,)
        K⁺ channel activation gate.
    d : ndarray, shape (n,)
        Ca²⁺ channel activation gate.
    f : ndarray, shape (n,)
        Ca²⁺ channel inactivation gate.
    i_stim : ndarray, shape (n,)
        Applied stimulus current trajectory (µA/cm²).
    exo_rate : ndarray, shape (n,)
        Exocytosis rate trajectory (ms⁻¹).
    dt : float
        Integration time step (ms).
    """

    t: np.ndarray
    v: np.ndarray
    ca: np.ndarray
    docked: np.ndarray
    reserve: np.ndarray
    released: np.ndarray
    m: np.ndarray
    h: np.ndarray
    n: np.ndarray
    d: np.ndarray
    f: np.ndarray
    i_stim: np.ndarray
    exo_rate: np.ndarray
    dt: float

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def total_release_events(self) -> float:
        """Total number of exocytosis events over the entire simulation."""
        return float(self.released[-1])

    @property
    def peak_voltage(self) -> float:
        """Maximum membrane potential reached (mV)."""
        return float(np.max(self.v))

    @property
    def peak_calcium(self) -> float:
        """Maximum cytosolic Ca²⁺ concentration reached (µM)."""
        return float(np.max(self.ca))

    @property
    def n_spikes(self) -> int:
        """Approximate number of action potentials detected (threshold = 0 mV)."""
        threshold = 0.0
        above = self.v > threshold
        crossings = np.diff(above.astype(int))
        return int(np.sum(crossings == 1))

    def release_rate(self) -> np.ndarray:
        """Instantaneous release rate estimated from cumulative releases (ms⁻¹).

        Returns
        -------
        ndarray, shape (n-1,)
            Finite-difference approximation of ``d(released)/dt``.
        """
        return np.diff(self.released) / self.dt
