"""Tests for electroexo.models.membrane."""
import numpy as np
import pytest

from electroexo.models.membrane import Membrane, MembraneState


class TestMembraneState:
    def test_round_trip(self):
        s = MembraneState(v=-65.0, m=0.05, h=0.6, n=0.32, d=0.01, f=0.99)
        arr = s.as_array()
        s2 = MembraneState.from_array(arr)
        assert s.v == s2.v
        assert s.m == s2.m

    def test_array_length(self):
        s = MembraneState(v=-65.0, m=0.05, h=0.6, n=0.32, d=0.01, f=0.99)
        assert len(s.as_array()) == 6


class TestMembrane:
    def setup_method(self):
        self.mem = Membrane()

    def test_initial_state_near_rest(self):
        s = self.mem.initial_state()
        assert abs(s.v - (-65.0)) < 1e-10

    def test_gate_variables_in_range(self):
        s = self.mem.initial_state()
        for val in [s.m, s.h, s.n, s.d, s.f]:
            assert 0.0 <= val <= 1.0

    def test_derivatives_return_correct_type(self):
        s = self.mem.initial_state()
        ds = self.mem.derivatives(s, i_ext=0.0)
        assert isinstance(ds, MembraneState)

    def test_leak_current_at_rest_small(self):
        # At rest, net current should be near zero
        s = self.mem.initial_state()
        i_leak = self.mem.leak_current(s.v)
        assert abs(i_leak) < 5.0  # small but non-zero due to model parameters

    def test_positive_current_depolarises(self):
        s = self.mem.initial_state()
        ds = self.mem.derivatives(s, i_ext=20.0)
        assert ds.v > 0  # positive external current should raise dV/dt

    def test_calcium_current_at_rest_is_inward(self):
        s = self.mem.initial_state()
        i_ca = self.mem.calcium_current(s.v, s.d, s.f)
        assert i_ca < 0  # inward at rest

    def test_custom_channels(self):
        from electroexo.models.channels import CalciumChannel
        ca_ch = CalciumChannel(g_max=5.0)
        mem = Membrane(ca_channel=ca_ch)
        assert mem.ca.g_max == 5.0

    def test_no_current_stability(self):
        """With no external current, dV/dt should be small at rest."""
        s = self.mem.initial_state()
        ds = self.mem.derivatives(s, i_ext=0.0)
        assert abs(ds.v) < 1.0  # nearly at steady state
