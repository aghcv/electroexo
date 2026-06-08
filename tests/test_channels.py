"""Tests for electroexo.models.channels."""
import numpy as np
import pytest

from electroexo.models.channels import CalciumChannel, PotassiumChannel, SodiumChannel


class TestSodiumChannel:
    def setup_method(self):
        self.ch = SodiumChannel()

    def test_defaults(self):
        assert self.ch.g_max == 120.0
        assert self.ch.e_rev == 50.0

    def test_m_inf_range(self):
        for v in np.linspace(-100, 50, 30):
            m = self.ch.m_inf(v)
            assert 0.0 <= m <= 1.0, f"m_inf out of range at V={v}"

    def test_h_inf_range(self):
        for v in np.linspace(-100, 50, 30):
            h = self.ch.h_inf(v)
            assert 0.0 <= h <= 1.0, f"h_inf out of range at V={v}"

    def test_m_inf_increases_with_depolarisation(self):
        m_rest = self.ch.m_inf(-65.0)
        m_depol = self.ch.m_inf(0.0)
        assert m_depol > m_rest

    def test_h_inf_decreases_with_depolarisation(self):
        h_rest = self.ch.h_inf(-65.0)
        h_depol = self.ch.h_inf(0.0)
        assert h_depol < h_rest

    def test_tau_m_positive(self):
        for v in np.linspace(-100, 50, 20):
            assert self.ch.tau_m(v) > 0

    def test_tau_h_positive(self):
        for v in np.linspace(-100, 50, 20):
            assert self.ch.tau_h(v) > 0

    def test_current_sign_at_reversal(self):
        m, h = 0.5, 0.8
        i = self.ch.current(self.ch.e_rev, m, h)
        assert abs(i) < 1e-10

    def test_current_inward_below_reversal(self):
        m, h = 0.5, 0.8
        i = self.ch.current(-65.0, m, h)
        assert i < 0  # inward (below reversal)

    def test_initial_state_keys(self):
        state = self.ch.initial_state(-65.0)
        assert "m" in state and "h" in state

    def test_dm_dt_drives_to_steady_state(self):
        v = -20.0
        m_inf = self.ch.m_inf(v)
        # If m < m_inf, dm/dt should be positive
        dm = self.ch.dm_dt(v, m_inf - 0.1)
        assert dm > 0
        # If m > m_inf, dm/dt should be negative
        dm = self.ch.dm_dt(v, m_inf + 0.1)
        assert dm < 0


class TestPotassiumChannel:
    def setup_method(self):
        self.ch = PotassiumChannel()

    def test_defaults(self):
        assert self.ch.g_max == 36.0
        assert self.ch.e_rev == -77.0

    def test_n_inf_range(self):
        for v in np.linspace(-100, 50, 30):
            n = self.ch.n_inf(v)
            assert 0.0 <= n <= 1.0

    def test_n_inf_increases_with_depolarisation(self):
        assert self.ch.n_inf(0.0) > self.ch.n_inf(-65.0)

    def test_tau_n_positive(self):
        for v in np.linspace(-100, 50, 20):
            assert self.ch.tau_n(v) > 0

    def test_current_sign(self):
        n = 0.4
        # At rest (~-65 mV), below reversal (-77 mV) — current should be outward
        i = self.ch.current(-65.0, n)
        assert i > 0  # outward

    def test_initial_state_key(self):
        state = self.ch.initial_state(-65.0)
        assert "n" in state

    def test_dn_dt_drives_to_steady_state(self):
        v = -20.0
        n_inf = self.ch.n_inf(v)
        assert self.ch.dn_dt(v, n_inf - 0.1) > 0
        assert self.ch.dn_dt(v, n_inf + 0.1) < 0


class TestCalciumChannel:
    def setup_method(self):
        self.ch = CalciumChannel()

    def test_defaults(self):
        assert self.ch.g_max == 2.0
        assert self.ch.e_rev == 120.0

    def test_d_inf_range(self):
        for v in np.linspace(-100, 50, 30):
            d = self.ch.d_inf(v)
            assert 0.0 <= d <= 1.0

    def test_f_inf_range(self):
        for v in np.linspace(-100, 50, 30):
            f = self.ch.f_inf(v)
            assert 0.0 <= f <= 1.0

    def test_d_inf_increases_with_depolarisation(self):
        assert self.ch.d_inf(0.0) > self.ch.d_inf(-65.0)

    def test_f_inf_decreases_with_depolarisation(self):
        assert self.ch.f_inf(0.0) < self.ch.f_inf(-65.0)

    def test_tau_positive(self):
        for v in np.linspace(-100, 50, 20):
            assert self.ch.tau_d(v) > 0
            assert self.ch.tau_f(v) > 0

    def test_current_inward_at_rest(self):
        d, f = 0.5, 0.8
        i = self.ch.current(-65.0, d, f)
        assert i < 0  # -65 < 120 mV → inward

    def test_initial_state_keys(self):
        state = self.ch.initial_state(-65.0)
        assert "d" in state and "f" in state
