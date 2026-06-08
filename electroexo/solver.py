"""
electroexo.solver
==================
Numerical ODE solvers for time-integration of the simulation state.

Two fixed-step solvers are provided:

* :class:`EulerSolver`  — explicit Euler (first-order; fast)
* :class:`RK4Solver`    — classical fourth-order Runge-Kutta (preferred)

Both solvers accept a ``derivatives`` callable with signature::

    derivs = f(t, y)   →  same type as y

where ``y`` is a NumPy array of state variables.

A convenience :func:`solve_ivp_fixed` function runs a solver over a fixed
time grid and returns the trajectory as a NumPy array.

Units
-----
Time : ms  (the user is responsible for consistent unit choice)
"""

from __future__ import annotations

from typing import Callable

import numpy as np

__all__ = ["EulerSolver", "RK4Solver", "solve_ivp_fixed"]

DerivsFn = Callable[[float, np.ndarray], np.ndarray]


class EulerSolver:
    """Explicit (forward) Euler method.

    Parameters
    ----------
    dt : float
        Fixed time step (ms).
    """

    def __init__(self, dt: float = 0.01) -> None:
        self.dt = dt

    def step(
        self, f: DerivsFn, t: float, y: np.ndarray
    ) -> tuple[float, np.ndarray]:
        """Advance the state by one time step.

        Parameters
        ----------
        f : callable
            Derivative function ``f(t, y) → dy/dt``.
        t : float
            Current time (ms).
        y : ndarray
            Current state vector.

        Returns
        -------
        tuple[float, ndarray]
            ``(t_new, y_new)``
        """
        dy = f(t, y)
        return t + self.dt, y + self.dt * dy


class RK4Solver:
    """Classical fourth-order Runge-Kutta method.

    Parameters
    ----------
    dt : float
        Fixed time step (ms).
    """

    def __init__(self, dt: float = 0.01) -> None:
        self.dt = dt

    def step(
        self, f: DerivsFn, t: float, y: np.ndarray
    ) -> tuple[float, np.ndarray]:
        """Advance the state by one time step using the RK4 scheme.

        Parameters
        ----------
        f : callable
            Derivative function ``f(t, y) → dy/dt``.
        t : float
            Current time (ms).
        y : ndarray
            Current state vector.

        Returns
        -------
        tuple[float, ndarray]
            ``(t_new, y_new)``
        """
        dt = self.dt
        k1 = f(t, y)
        k2 = f(t + dt / 2, y + dt / 2 * k1)
        k3 = f(t + dt / 2, y + dt / 2 * k2)
        k4 = f(t + dt, y + dt * k3)
        return t + dt, y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def solve_ivp_fixed(
    f: DerivsFn,
    y0: np.ndarray,
    t_end: float,
    dt: float = 0.01,
    solver: str = "rk4",
    record_every: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate an ODE system over ``[0, t_end]`` with a fixed step.

    Parameters
    ----------
    f :
        Derivative function ``f(t, y) → dy/dt``.
    y0 :
        Initial state vector.
    t_end :
        End time (ms).
    dt :
        Time step (ms).  Default 0.01.
    solver :
        ``"rk4"`` (default) or ``"euler"``.
    record_every :
        Record state every *N* steps to reduce memory use.  Default 1.

    Returns
    -------
    tuple[ndarray, ndarray]
        ``(t_arr, y_arr)`` where *t_arr* has shape ``(n_steps,)`` and
        *y_arr* has shape ``(n_steps, n_vars)``.
    """
    if solver == "rk4":
        _solver = RK4Solver(dt=dt)
    elif solver == "euler":
        _solver = EulerSolver(dt=dt)
    else:
        raise ValueError(f"Unknown solver '{solver}'; choose 'rk4' or 'euler'.")

    n_steps = int(t_end / dt)
    t = 0.0
    y = y0.copy()

    t_list: list[float] = [t]
    y_list: list[np.ndarray] = [y.copy()]

    for step_i in range(1, n_steps + 1):
        t, y = _solver.step(f, t, y)
        if step_i % record_every == 0:
            t_list.append(t)
            y_list.append(y.copy())

    return np.array(t_list), np.array(y_list)
