"""Tests for electroexo.solver."""
import numpy as np
import pytest

from electroexo.solver import EulerSolver, RK4Solver, solve_ivp_fixed


def _exponential_decay(t, y, lam=0.1):
    """dy/dt = -lambda * y  →  exact: y(t) = y0 * exp(-lambda * t)"""
    return -lam * y


def _make_decay_fn(lam=0.1):
    return lambda t, y: -lam * y


class TestEulerSolver:
    def test_single_step_positive(self):
        solver = EulerSolver(dt=0.1)
        f = _make_decay_fn()
        y0 = np.array([1.0])
        t1, y1 = solver.step(f, 0.0, y0)
        assert t1 == pytest.approx(0.1)
        assert y1[0] < y0[0]

    def test_time_advances(self):
        solver = EulerSolver(dt=1.0)
        f = _make_decay_fn()
        t, y = 5.0, np.array([2.0])
        t_new, _ = solver.step(f, t, y)
        assert t_new == pytest.approx(6.0)

    def test_stable_zero_derivative(self):
        solver = EulerSolver(dt=0.5)
        f = lambda t, y: np.zeros_like(y)
        y0 = np.array([3.0, -1.0])
        _, y1 = solver.step(f, 0.0, y0)
        np.testing.assert_allclose(y1, y0)


class TestRK4Solver:
    def test_accuracy_versus_euler(self):
        """RK4 should be more accurate than Euler for the same step size."""
        dt = 1.0
        lam = 0.1
        y0 = np.array([1.0])
        t_end = 10.0
        exact = y0 * np.exp(-lam * t_end)

        euler = EulerSolver(dt=dt)
        rk4 = RK4Solver(dt=dt)
        f = _make_decay_fn(lam)

        t_e, y_e = 0.0, y0.copy()
        t_r, y_r = 0.0, y0.copy()
        for _ in range(int(t_end / dt)):
            t_e, y_e = euler.step(f, t_e, y_e)
            t_r, y_r = rk4.step(f, t_r, y_r)

        err_euler = abs(y_e[0] - exact[0])
        err_rk4 = abs(y_r[0] - exact[0])
        assert err_rk4 < err_euler

    def test_time_advances(self):
        solver = RK4Solver(dt=0.5)
        f = _make_decay_fn()
        t, y = 0.0, np.array([1.0])
        t_new, _ = solver.step(f, t, y)
        assert t_new == pytest.approx(0.5)

    def test_multidimensional_state(self):
        """Check that RK4 handles vector states correctly."""
        solver = RK4Solver(dt=0.01)
        # Two independent decays
        f = lambda t, y: np.array([-0.1 * y[0], -0.5 * y[1]])
        y0 = np.array([1.0, 2.0])
        t, y = 0.0, y0.copy()
        for _ in range(100):
            t, y = solver.step(f, t, y)
        assert y[0] < y0[0]
        assert y[1] < y0[1]


class TestSolveIvpFixed:
    def test_output_shapes(self):
        f = _make_decay_fn()
        y0 = np.array([1.0])
        t_arr, y_arr = solve_ivp_fixed(f, y0, t_end=10.0, dt=0.1, record_every=1)
        assert t_arr.shape[0] == y_arr.shape[0]
        assert y_arr.shape[1] == 1

    def test_first_point_is_initial_condition(self):
        f = _make_decay_fn()
        y0 = np.array([2.5])
        t_arr, y_arr = solve_ivp_fixed(f, y0, t_end=5.0, dt=0.1)
        assert t_arr[0] == pytest.approx(0.0)
        assert y_arr[0, 0] == pytest.approx(2.5)

    def test_record_every_reduces_output(self):
        f = _make_decay_fn()
        y0 = np.array([1.0])
        _, y_full = solve_ivp_fixed(f, y0, t_end=10.0, dt=0.1, record_every=1)
        _, y_thin = solve_ivp_fixed(f, y0, t_end=10.0, dt=0.1, record_every=10)
        assert y_thin.shape[0] < y_full.shape[0]

    def test_euler_solver(self):
        f = _make_decay_fn()
        y0 = np.array([1.0])
        t_arr, y_arr = solve_ivp_fixed(f, y0, t_end=5.0, dt=0.1, solver="euler")
        assert t_arr[-1] == pytest.approx(5.0)

    def test_unknown_solver_raises(self):
        f = _make_decay_fn()
        y0 = np.array([1.0])
        with pytest.raises(ValueError, match="Unknown solver"):
            solve_ivp_fixed(f, y0, t_end=5.0, solver="bogus")
