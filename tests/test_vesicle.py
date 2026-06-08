"""Tests for electroexo.models.vesicle."""
import pytest

from electroexo.models.vesicle import VesiclePool, VesicleState


class TestVesicleState:
    def test_defaults(self):
        s = VesicleState(reserve=100.0, docked=10.0)
        assert s.released == 0.0


class TestVesiclePool:
    def setup_method(self):
        self.pool = VesiclePool()

    def test_defaults(self):
        assert self.pool.n_reserve == 200.0
        assert self.pool.n_docked == 20.0
        assert self.pool.k_dock == 5e-4

    def test_initial_state(self):
        s = self.pool.initial_state()
        assert s.reserve == self.pool.n_reserve
        assert s.docked == self.pool.n_docked
        assert s.released == 0.0

    def test_no_depletion_reserve_unchanged(self):
        """By default reserve pool is not depleted."""
        s = self.pool.initial_state()
        d_res, _ = self.pool.derivatives(s, exo_rate=0.1)
        assert d_res == 0.0

    def test_high_exo_rate_drains_docked_pool(self):
        s = self.pool.initial_state()
        _, d_docked = self.pool.derivatives(s, exo_rate=1.0)
        # High exo_rate should make d_docked negative (drain faster than replenish)
        assert d_docked < 0

    def test_zero_exo_rate_replenishes_docked(self):
        """With no exocytosis the docked pool is replenished from reserve."""
        pool = VesiclePool(n_docked=0.0)
        s = pool.initial_state()
        _, d_docked = pool.derivatives(s, exo_rate=0.0)
        assert d_docked > 0  # replenishment flux

    def test_deplete_reserve_mode(self):
        pool = VesiclePool(deplete_reserve=True)
        s = pool.initial_state()
        d_res, _ = pool.derivatives(s, exo_rate=0.0)
        assert d_res < 0  # reserve is drained

    def test_custom_parameters(self):
        pool = VesiclePool(n_reserve=500, n_docked=50, k_dock=1e-3)
        s = pool.initial_state()
        assert s.reserve == 500
        assert s.docked == 50
