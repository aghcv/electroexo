"""
electroexo.models.vesicle
===========================
Vesicle pool kinetic model.

Three-pool scheme
-----------------
Vesicles reside in one of three pools:

``R`` — **Reserve pool** (large; slowly replenishes docked pool)
``D`` — **Docked pool** (primed/docked; immediately available for fusion)
``F`` — **Fused** (released; counted as cumulative release events)

The rate equations are::

    dR/dt = −k_dock · R
    dD/dt =  k_dock · R − k_exo(Ca) · D
    cumulative_releases += k_exo(Ca) · D · dt

For simplicity the reserve pool is treated as inexhaustible (constant
replenishment flux keeps ``R`` at its resting level) unless
``deplete_reserve=True``.

Units
-----
Pools : vesicle count (dimensionless)
Rates : ms⁻¹
Time  : ms
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["VesiclePool", "VesicleState"]


@dataclass
class VesicleState:
    """Snapshot of vesicle pool occupancies.

    Attributes
    ----------
    reserve : float
        Number of vesicles in the reserve pool.
    docked : float
        Number of primed/docked vesicles ready for fusion.
    released : float
        Cumulative number of release events since simulation start.
    """

    reserve: float
    docked: float
    released: float = 0.0


class VesiclePool:
    """Three-pool vesicle kinetic model.

    Parameters
    ----------
    n_reserve : float
        Initial reserve-pool size (vesicle count).  Default 200.
    n_docked : float
        Initial docked-pool size (readily-releasable pool, RRP).  Default 20.
    k_dock : float
        Reserve → docked docking rate (ms⁻¹).  Default 5e-4.
    deplete_reserve : bool
        If ``True``, track reserve-pool depletion.  Default ``False``.
    """

    def __init__(
        self,
        n_reserve: float = 200.0,
        n_docked: float = 20.0,
        k_dock: float = 5e-4,
        deplete_reserve: bool = False,
    ) -> None:
        self.n_reserve = n_reserve
        self.n_docked = n_docked
        self.k_dock = k_dock
        self.deplete_reserve = deplete_reserve

    # ------------------------------------------------------------------
    # ODE derivatives
    # ------------------------------------------------------------------

    def derivatives(
        self,
        state: VesicleState,
        exo_rate: float,
    ) -> tuple[float, float]:
        """Compute d(reserve)/dt and d(docked)/dt.

        Parameters
        ----------
        state :
            Current vesicle pool state.
        exo_rate :
            Exocytosis rate (ms⁻¹) applied to the docked pool.

        Returns
        -------
        tuple[float, float]
            ``(d_reserve_dt, d_docked_dt)``
        """
        if self.deplete_reserve:
            d_reserve = -self.k_dock * state.reserve
            # Replenishment of docked pool from reserve
            replenishment = self.k_dock * state.reserve
        else:
            d_reserve = 0.0
            replenishment = self.k_dock * self.n_reserve

        d_docked = replenishment - exo_rate * state.docked
        return d_reserve, d_docked

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def initial_state(self) -> VesicleState:
        """Return resting vesicle pool state."""
        return VesicleState(
            reserve=self.n_reserve,
            docked=self.n_docked,
            released=0.0,
        )
