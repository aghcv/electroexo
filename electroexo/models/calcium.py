"""
electroexo.models.calcium
===========================
Intracellular calcium dynamics model.

The free cytosolic Ca²⁺ concentration [Ca²⁺]_i evolves according to

    d[Ca]/dt = −α · I_Ca / (2 · F · vol) − ([Ca] − [Ca]_rest) / τ_Ca

where:

* *I_Ca*  — Ca²⁺ channel current density (µA/cm²)
* *α*     — conversion factor (µM · cm² / µA / ms)
* *F*     — Faraday's constant (conversion embedded in α)
* *vol*   — effective submembrane shell volume (cm)
* *τ_Ca*  — calcium clearance/buffering time constant (ms)
* *[Ca]_rest* — resting free Ca²⁺ concentration (µM)

The default parameter values correspond to a small spherical neuronal
compartment with fast buffering.

Units
-----
Concentration : µM
Time          : ms
"""

from __future__ import annotations

__all__ = ["CalciumDynamics"]


class CalciumDynamics:
    """Single-compartment intracellular Ca²⁺ dynamics.

    Parameters
    ----------
    ca_rest : float
        Resting free Ca²⁺ concentration (µM).  Default 0.1.
    tau_ca : float
        Calcium clearance time constant (ms).  Default 100.
    alpha : float
        Current-to-concentration conversion factor (µM·cm²/µA/ms).
        Encodes ``1 / (2 · F · depth)``.  Default 1e-3.
    """

    def __init__(
        self,
        ca_rest: float = 0.1,
        tau_ca: float = 100.0,
        alpha: float = 1e-3,
    ) -> None:
        self.ca_rest = ca_rest
        self.tau_ca = tau_ca
        self.alpha = alpha

    # ------------------------------------------------------------------
    # ODE derivative
    # ------------------------------------------------------------------

    def dca_dt(self, ca: float, i_ca: float) -> float:
        """Time derivative of cytosolic Ca²⁺ concentration.

        Parameters
        ----------
        ca :
            Current free Ca²⁺ concentration (µM).
        i_ca :
            Ca²⁺ channel current density (µA/cm²).  Positive = outward
            by convention; inward current is negative and thus *increases*
            [Ca²⁺]_i.

        Returns
        -------
        float
            d[Ca]/dt (µM/ms).
        """
        # Inward Ca²⁺ current is negative in I_Ca convention → sign inversion
        influx = -self.alpha * i_ca
        clearance = -(ca - self.ca_rest) / self.tau_ca
        return influx + clearance

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def initial_state(self) -> float:
        """Return resting calcium concentration (µM)."""
        return self.ca_rest
