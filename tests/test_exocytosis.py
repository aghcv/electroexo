"""Tests for electroexo.models.exocytosis."""
import pytest

from electroexo.models.exocytosis import AllostericExocytosis, CooperativeExocytosis


class TestCooperativeExocytosis:
    def setup_method(self):
        self.model = CooperativeExocytosis()

    def test_rate_at_zero_calcium_equals_basal(self):
        assert self.model.exocytosis_rate(0.0) == self.model.k_basal

    def test_rate_increases_with_calcium(self):
        r_low = self.model.exocytosis_rate(0.5)
        r_high = self.model.exocytosis_rate(5.0)
        assert r_high > r_low

    def test_rate_saturates_at_k_max(self):
        rate = self.model.exocytosis_rate(1e6)
        assert abs(rate - (self.model.k_basal + self.model.k_max)) < 1e-8

    def test_half_maximal_at_k_d(self):
        half_rate = self.model.k_basal + self.model.k_max / 2.0
        rate_at_kd = self.model.exocytosis_rate(self.model.k_d)
        assert abs(rate_at_kd - half_rate) < 1e-8

    def test_negative_calcium_clamped(self):
        rate = self.model.exocytosis_rate(-1.0)
        assert rate == self.model.k_basal

    def test_cooperativity(self):
        """Higher n_hill → steeper activation curve."""
        model_low = CooperativeExocytosis(n_hill=1.0)
        model_high = CooperativeExocytosis(n_hill=8.0)
        # At subthreshold Ca, high cooperativity → lower rate
        ca_sub = model_low.k_d * 0.5
        assert model_high.exocytosis_rate(ca_sub) < model_low.exocytosis_rate(ca_sub)

    def test_with_basal_rate(self):
        model = CooperativeExocytosis(k_basal=1e-5)
        assert model.exocytosis_rate(0.0) == pytest.approx(1e-5)


class TestAllostericExocytosis:
    def setup_method(self):
        self.model = AllostericExocytosis()

    def test_rate_at_zero_calcium_equals_l0(self):
        assert self.model.exocytosis_rate(0.0) == pytest.approx(self.model.l0)

    def test_rate_increases_with_calcium(self):
        r_low = self.model.exocytosis_rate(0.1)
        r_high = self.model.exocytosis_rate(10.0)
        assert r_high > r_low

    def test_k_max_clamp(self):
        model = AllostericExocytosis(k_max=1e-2)
        rate = model.exocytosis_rate(1e6)
        assert rate <= 1e-2 + 1e-12

    def test_negative_calcium_clamped(self):
        rate = self.model.exocytosis_rate(-5.0)
        assert rate == pytest.approx(self.model.l0)

    def test_n_sites_effect(self):
        model1 = AllostericExocytosis(n_sites=1)
        model5 = AllostericExocytosis(n_sites=5)
        ca = 3.0
        assert model5.exocytosis_rate(ca) > model1.exocytosis_rate(ca)
