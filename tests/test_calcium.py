"""Tests for electroexo.models.calcium."""
import pytest

from electroexo.models.calcium import CalciumDynamics


class TestCalciumDynamics:
    def setup_method(self):
        self.ca = CalciumDynamics()

    def test_defaults(self):
        assert self.ca.ca_rest == 0.1
        assert self.ca.tau_ca == 100.0
        assert self.ca.alpha == 1e-3

    def test_initial_state(self):
        assert self.ca.initial_state() == self.ca.ca_rest

    def test_inward_current_raises_calcium(self):
        # Inward Ca2+ current is negative → should increase [Ca]
        dca = self.ca.dca_dt(ca=0.1, i_ca=-1.0)
        assert dca > 0

    def test_outward_current_lowers_calcium(self):
        # Outward Ca2+ current is positive → should decrease [Ca] (only clearance)
        # With elevated [Ca] and outward current
        dca = self.ca.dca_dt(ca=5.0, i_ca=1.0)
        assert dca < 0

    def test_at_rest_near_zero_derivative(self):
        # No current and [Ca] at rest → dCa/dt ≈ 0
        dca = self.ca.dca_dt(ca=self.ca.ca_rest, i_ca=0.0)
        assert abs(dca) < 1e-12

    def test_clearance_drives_to_rest(self):
        # Elevated [Ca] with no current → should decrease
        dca = self.ca.dca_dt(ca=2.0, i_ca=0.0)
        assert dca < 0

    def test_custom_parameters(self):
        ca = CalciumDynamics(ca_rest=0.05, tau_ca=50.0, alpha=5e-4)
        assert ca.ca_rest == 0.05
        assert ca.tau_ca == 50.0

    def test_large_inward_current(self):
        # Very large inward current → strongly positive dCa/dt
        dca = self.ca.dca_dt(ca=0.1, i_ca=-100.0)
        assert dca > 0.05
