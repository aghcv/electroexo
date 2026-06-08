"""
electroexo.simulation
======================
Main simulation engine for electro-exocytosis modelling.

:class:`SimulationEngine` wires together the biophysical sub-models and a
numerical solver to produce a full :class:`~electroexo.results.SimulationResults`
trajectory for a given electrical stimulus.

Quick start
-----------
>>> from electroexo import SimulationEngine
>>> from electroexo.stimulus import StepStimulus
>>> engine = SimulationEngine()
>>> results = engine.run(t_end=200.0, stimulus=StepStimulus(amplitude=10.0,
...                                                         t_start=50.0,
...                                                         t_end=150.0))
>>> print(f"Spikes: {results.n_spikes}, Release events: {results.total_release_events:.2f}")
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from electroexo.models.calcium import CalciumDynamics
from electroexo.models.exocytosis import CooperativeExocytosis
from electroexo.models.membrane import Membrane
from electroexo.models.vesicle import VesiclePool
from electroexo.results import SimulationResults
from electroexo.solver import RK4Solver
from electroexo.stimulus import BaseStimulus, ConstantStimulus

__all__ = ["SimulationEngine"]

# State vector indices
_IDX_V = 0          # membrane potential
_IDX_M = 1          # Na+ activation
_IDX_H = 2          # Na+ inactivation
_IDX_N = 3          # K+ activation
_IDX_D = 4          # Ca2+ activation
_IDX_F = 5          # Ca2+ inactivation
_IDX_CA = 6         # cytosolic Ca2+
_IDX_RESERVE = 7    # reserve pool
_IDX_DOCKED = 8     # docked (RRP) pool
_IDX_RELEASED = 9   # cumulative release events


class SimulationEngine:
    """Full electro-exocytosis simulation engine.

    Combines:
    * Hodgkin-Huxley + L-type Ca²⁺ membrane model
    * Single-compartment Ca²⁺ dynamics
    * Three-pool vesicle kinetics
    * Ca²⁺-triggered exocytosis

    Parameters
    ----------
    membrane : Membrane, optional
        Membrane model.  Defaults to a standard HH membrane with L-type Ca²⁺.
    calcium : CalciumDynamics, optional
        Calcium dynamics model.  Defaults to :class:`~electroexo.models.calcium.CalciumDynamics`.
    vesicle_pool : VesiclePool, optional
        Vesicle pool model.  Defaults to :class:`~electroexo.models.vesicle.VesiclePool`.
    exocytosis_model : CooperativeExocytosis, optional
        Exocytosis rate model.  Defaults to
        :class:`~electroexo.models.exocytosis.CooperativeExocytosis`.
    dt : float
        Integration time step (ms).  Default 0.01.
    record_every : int
        Record state every *N* steps.  Default 10 (records every 0.1 ms).
    """

    def __init__(
        self,
        membrane: Optional[Membrane] = None,
        calcium: Optional[CalciumDynamics] = None,
        vesicle_pool: Optional[VesiclePool] = None,
        exocytosis_model: Optional[CooperativeExocytosis] = None,
        dt: float = 0.01,
        record_every: int = 10,
    ) -> None:
        self.membrane = membrane if membrane is not None else Membrane()
        self.calcium = calcium if calcium is not None else CalciumDynamics()
        self.vesicle_pool = vesicle_pool if vesicle_pool is not None else VesiclePool()
        self.exocytosis_model = (
            exocytosis_model if exocytosis_model is not None else CooperativeExocytosis()
        )
        self.dt = dt
        self.record_every = record_every

    # ------------------------------------------------------------------
    # Internal ODE system
    # ------------------------------------------------------------------

    def _build_state_vector(self) -> np.ndarray:
        """Assemble the initial state vector."""
        mem0 = self.membrane.initial_state()
        ca0 = self.calcium.initial_state()
        ves0 = self.vesicle_pool.initial_state()
        y0 = np.zeros(10)
        y0[_IDX_V] = mem0.v
        y0[_IDX_M] = mem0.m
        y0[_IDX_H] = mem0.h
        y0[_IDX_N] = mem0.n
        y0[_IDX_D] = mem0.d
        y0[_IDX_F] = mem0.f
        y0[_IDX_CA] = ca0
        y0[_IDX_RESERVE] = ves0.reserve
        y0[_IDX_DOCKED] = ves0.docked
        y0[_IDX_RELEASED] = ves0.released
        return y0

    def _make_derivatives(
        self, stimulus: BaseStimulus
    ):
        """Return the derivative function for the full state vector."""
        membrane = self.membrane
        calcium = self.calcium
        vesicle_pool = self.vesicle_pool
        exocytosis_model = self.exocytosis_model

        def derivatives(t: float, y: np.ndarray) -> np.ndarray:
            v = y[_IDX_V]
            m = y[_IDX_M]
            h = y[_IDX_H]
            n = y[_IDX_N]
            d = y[_IDX_D]
            f = y[_IDX_F]
            ca = y[_IDX_CA]
            reserve = y[_IDX_RESERVE]
            docked = y[_IDX_DOCKED]

            i_ext = stimulus.amplitude_at(t)

            # Membrane
            from electroexo.models.membrane import MembraneState  # avoid circular
            mem_state = MembraneState(v=v, m=m, h=h, n=n, d=d, f=f)
            dmem = membrane.derivatives(mem_state, i_ext)

            # Calcium
            i_ca = membrane.calcium_current(v, d, f)
            dca = calcium.dca_dt(ca, i_ca)

            # Exocytosis rate
            exo_rate = exocytosis_model.exocytosis_rate(ca)

            # Vesicle pools
            from electroexo.models.vesicle import VesicleState  # avoid circular
            ves_state = VesicleState(reserve=reserve, docked=docked)
            d_reserve, d_docked = vesicle_pool.derivatives(ves_state, exo_rate)

            # Release rate
            d_released = exo_rate * max(docked, 0.0)

            dy = np.zeros(10)
            dy[_IDX_V] = dmem.v
            dy[_IDX_M] = dmem.m
            dy[_IDX_H] = dmem.h
            dy[_IDX_N] = dmem.n
            dy[_IDX_D] = dmem.d
            dy[_IDX_F] = dmem.f
            dy[_IDX_CA] = dca
            dy[_IDX_RESERVE] = d_reserve
            dy[_IDX_DOCKED] = d_docked
            dy[_IDX_RELEASED] = d_released

            return dy

        return derivatives

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        t_end: float = 200.0,
        stimulus: Optional[BaseStimulus] = None,
    ) -> SimulationResults:
        """Run the simulation and return a :class:`~electroexo.results.SimulationResults`.

        Parameters
        ----------
        t_end :
            Simulation end time (ms).  Default 200.
        stimulus :
            Electrical stimulus protocol.  Defaults to zero current.

        Returns
        -------
        SimulationResults
        """
        if stimulus is None:
            stimulus = ConstantStimulus(amplitude=0.0)

        y0 = self._build_state_vector()
        f = self._make_derivatives(stimulus)
        solver = RK4Solver(dt=self.dt)

        n_steps = int(t_end / self.dt)
        t = 0.0
        y = y0.copy()

        # Record lists (pre-allocate estimate)
        n_records = n_steps // self.record_every + 2
        t_arr = np.empty(n_records)
        y_arr = np.empty((n_records, 10))
        i_stim_arr = np.empty(n_records)
        exo_rate_arr = np.empty(n_records)

        rec_idx = 0
        t_arr[rec_idx] = t
        y_arr[rec_idx] = y
        i_stim_arr[rec_idx] = stimulus.amplitude_at(t)
        exo_rate_arr[rec_idx] = self.exocytosis_model.exocytosis_rate(y[_IDX_CA])

        for step_i in range(1, n_steps + 1):
            t, y = solver.step(f, t, y)
            # Clamp gate variables to [0, 1] and concentrations to non-negative
            y[_IDX_M] = np.clip(y[_IDX_M], 0.0, 1.0)
            y[_IDX_H] = np.clip(y[_IDX_H], 0.0, 1.0)
            y[_IDX_N] = np.clip(y[_IDX_N], 0.0, 1.0)
            y[_IDX_D] = np.clip(y[_IDX_D], 0.0, 1.0)
            y[_IDX_F] = np.clip(y[_IDX_F], 0.0, 1.0)
            y[_IDX_CA] = max(y[_IDX_CA], 0.0)
            y[_IDX_RESERVE] = max(y[_IDX_RESERVE], 0.0)
            y[_IDX_DOCKED] = max(y[_IDX_DOCKED], 0.0)
            y[_IDX_RELEASED] = max(y[_IDX_RELEASED], 0.0)

            if step_i % self.record_every == 0:
                rec_idx += 1
                if rec_idx >= n_records:
                    # Expand arrays
                    t_arr = np.concatenate([t_arr, np.empty(n_records)])
                    y_arr = np.concatenate([y_arr, np.empty((n_records, 10))], axis=0)
                    i_stim_arr = np.concatenate([i_stim_arr, np.empty(n_records)])
                    exo_rate_arr = np.concatenate([exo_rate_arr, np.empty(n_records)])
                    n_records *= 2

                t_arr[rec_idx] = t
                y_arr[rec_idx] = y
                i_stim_arr[rec_idx] = stimulus.amplitude_at(t)
                exo_rate_arr[rec_idx] = self.exocytosis_model.exocytosis_rate(y[_IDX_CA])

        # Trim to actual recorded length
        n_out = rec_idx + 1
        t_arr = t_arr[:n_out]
        y_arr = y_arr[:n_out]
        i_stim_arr = i_stim_arr[:n_out]
        exo_rate_arr = exo_rate_arr[:n_out]

        return SimulationResults(
            t=t_arr,
            v=y_arr[:, _IDX_V],
            ca=y_arr[:, _IDX_CA],
            docked=y_arr[:, _IDX_DOCKED],
            reserve=y_arr[:, _IDX_RESERVE],
            released=y_arr[:, _IDX_RELEASED],
            m=y_arr[:, _IDX_M],
            h=y_arr[:, _IDX_H],
            n=y_arr[:, _IDX_N],
            d=y_arr[:, _IDX_D],
            f=y_arr[:, _IDX_F],
            i_stim=i_stim_arr,
            exo_rate=exo_rate_arr,
            dt=self.dt * self.record_every,
        )
