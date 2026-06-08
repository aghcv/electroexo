"""
electroexo.models.channels
===========================
Voltage-gated ion channel models used in the electro-exocytosis framework.

Implemented channels
--------------------
* :class:`SodiumChannel`  — fast inactivating Na⁺ channel (Hodgkin-Huxley)
* :class:`PotassiumChannel` — delayed-rectifier K⁺ channel (Hodgkin-Huxley)
* :class:`CalciumChannel` — L-type voltage-gated Ca²⁺ channel

All models follow the Hodgkin-Huxley formalism.  Gate variables obey

    dx/dt = alpha_x(V) * (1 - x) - beta_x(V) * x

which can be recast as

    dx/dt = (x_inf(V) - x) / tau_x(V)

Units
-----
Voltage : mV
Time    : ms
Current density : µA/cm²
Conductance density : mS/cm²
"""

from __future__ import annotations

import numpy as np

__all__ = ["SodiumChannel", "PotassiumChannel", "CalciumChannel"]


# ---------------------------------------------------------------------------
# Sodium channel (Hodgkin-Huxley, 1952)
# ---------------------------------------------------------------------------

class SodiumChannel:
    """Fast, transient Na⁺ channel (Hodgkin-Huxley formalism).

    Parameters
    ----------
    g_max : float
        Maximum conductance density (mS/cm²).  Default 120.
    e_rev : float
        Reversal (Nernst) potential (mV).  Default +50.
    """

    def __init__(self, g_max: float = 120.0, e_rev: float = 50.0) -> None:
        self.g_max = g_max
        self.e_rev = e_rev

    # ------------------------------------------------------------------
    # Steady-state values and time constants
    # ------------------------------------------------------------------

    @staticmethod
    def _alpha_m(v: float) -> float:
        """Opening rate for activation gate *m* (ms⁻¹)."""
        dv = v + 40.0
        if abs(dv) < 1e-7:
            return 1.0
        return 0.1 * dv / (1.0 - np.exp(-dv / 10.0))

    @staticmethod
    def _beta_m(v: float) -> float:
        """Closing rate for activation gate *m* (ms⁻¹)."""
        return 4.0 * np.exp(-(v + 65.0) / 18.0)

    @staticmethod
    def _alpha_h(v: float) -> float:
        """Opening rate for inactivation gate *h* (ms⁻¹)."""
        return 0.07 * np.exp(-(v + 65.0) / 20.0)

    @staticmethod
    def _beta_h(v: float) -> float:
        """Closing rate for inactivation gate *h* (ms⁻¹)."""
        return 1.0 / (1.0 + np.exp(-(v + 35.0) / 10.0))

    def m_inf(self, v: float) -> float:
        """Steady-state activation."""
        a = self._alpha_m(v)
        return a / (a + self._beta_m(v))

    def tau_m(self, v: float) -> float:
        """Activation time constant (ms)."""
        return 1.0 / (self._alpha_m(v) + self._beta_m(v))

    def h_inf(self, v: float) -> float:
        """Steady-state inactivation."""
        a = self._alpha_h(v)
        return a / (a + self._beta_h(v))

    def tau_h(self, v: float) -> float:
        """Inactivation time constant (ms)."""
        return 1.0 / (self._alpha_h(v) + self._beta_h(v))

    # ------------------------------------------------------------------
    # Gate derivatives (for ODE integration)
    # ------------------------------------------------------------------

    def dm_dt(self, v: float, m: float) -> float:
        """Time derivative of *m* gate."""
        return self._alpha_m(v) * (1.0 - m) - self._beta_m(v) * m

    def dh_dt(self, v: float, h: float) -> float:
        """Time derivative of *h* gate."""
        return self._alpha_h(v) * (1.0 - h) - self._beta_h(v) * h

    # ------------------------------------------------------------------
    # Current
    # ------------------------------------------------------------------

    def current(self, v: float, m: float, h: float) -> float:
        """Ionic current density (µA/cm²)."""
        return self.g_max * m ** 3 * h * (v - self.e_rev)

    def conductance(self, m: float, h: float) -> float:
        """Instantaneous conductance density (mS/cm²)."""
        return self.g_max * m ** 3 * h

    def initial_state(self, v_rest: float = -65.0) -> dict:
        """Return gate variables at resting potential."""
        return {"m": self.m_inf(v_rest), "h": self.h_inf(v_rest)}


# ---------------------------------------------------------------------------
# Potassium channel (Hodgkin-Huxley, 1952)
# ---------------------------------------------------------------------------

class PotassiumChannel:
    """Delayed-rectifier K⁺ channel (Hodgkin-Huxley formalism).

    Parameters
    ----------
    g_max : float
        Maximum conductance density (mS/cm²).  Default 36.
    e_rev : float
        Reversal potential (mV).  Default −77.
    """

    def __init__(self, g_max: float = 36.0, e_rev: float = -77.0) -> None:
        self.g_max = g_max
        self.e_rev = e_rev

    # ------------------------------------------------------------------
    # Rate functions
    # ------------------------------------------------------------------

    @staticmethod
    def _alpha_n(v: float) -> float:
        """Opening rate for *n* gate (ms⁻¹)."""
        dv = v + 55.0
        if abs(dv) < 1e-7:
            return 0.1
        return 0.01 * dv / (1.0 - np.exp(-dv / 10.0))

    @staticmethod
    def _beta_n(v: float) -> float:
        """Closing rate for *n* gate (ms⁻¹)."""
        return 0.125 * np.exp(-(v + 65.0) / 80.0)

    def n_inf(self, v: float) -> float:
        """Steady-state activation."""
        a = self._alpha_n(v)
        return a / (a + self._beta_n(v))

    def tau_n(self, v: float) -> float:
        """Activation time constant (ms)."""
        return 1.0 / (self._alpha_n(v) + self._beta_n(v))

    def dn_dt(self, v: float, n: float) -> float:
        """Time derivative of *n* gate."""
        return self._alpha_n(v) * (1.0 - n) - self._beta_n(v) * n

    def current(self, v: float, n: float) -> float:
        """Ionic current density (µA/cm²)."""
        return self.g_max * n ** 4 * (v - self.e_rev)

    def conductance(self, n: float) -> float:
        """Instantaneous conductance density (mS/cm²)."""
        return self.g_max * n ** 4

    def initial_state(self, v_rest: float = -65.0) -> dict:
        """Return gate variable at resting potential."""
        return {"n": self.n_inf(v_rest)}


# ---------------------------------------------------------------------------
# L-type calcium channel
# ---------------------------------------------------------------------------

class CalciumChannel:
    """L-type voltage-gated Ca²⁺ channel (simplified two-gate model).

    The activation gate *d* follows Boltzmann steady-state kinetics and the
    inactivation gate *f* provides voltage-dependent inactivation.

    Parameters
    ----------
    g_max : float
        Maximum conductance density (mS/cm²).  Default 2.
    e_rev : float
        Reversal potential (mV).  Default +120.  Note: the Nernst potential
        for Ca²⁺ is large and positive under typical conditions.
    """

    def __init__(self, g_max: float = 2.0, e_rev: float = 120.0) -> None:
        self.g_max = g_max
        self.e_rev = e_rev

    # ------------------------------------------------------------------
    # Activation gate *d*
    # ------------------------------------------------------------------

    def d_inf(self, v: float) -> float:
        """Steady-state activation."""
        return 1.0 / (1.0 + np.exp(-(v + 10.0) / 8.0))

    def tau_d(self, v: float) -> float:  # noqa: D401
        """Activation time constant (ms)."""
        return 2.0 / (np.exp(-(v + 40.0) / 15.0) + np.exp((v + 40.0) / 15.0)) + 0.5

    def dd_dt(self, v: float, d: float) -> float:
        """Time derivative of *d* gate."""
        return (self.d_inf(v) - d) / self.tau_d(v)

    # ------------------------------------------------------------------
    # Inactivation gate *f*
    # ------------------------------------------------------------------

    def f_inf(self, v: float) -> float:
        """Steady-state inactivation."""
        return 1.0 / (1.0 + np.exp((v + 40.0) / 6.0))

    def tau_f(self, v: float) -> float:
        """Inactivation time constant (ms)."""
        return 80.0 / (0.02 * np.exp(-(v + 65.0) / 20.0) + 0.5) + 2.0

    def df_dt(self, v: float, f: float) -> float:
        """Time derivative of *f* gate."""
        return (self.f_inf(v) - f) / self.tau_f(v)

    # ------------------------------------------------------------------
    # Current
    # ------------------------------------------------------------------

    def current(self, v: float, d: float, f: float) -> float:
        """Ionic current density (µA/cm²)."""
        return self.g_max * d ** 2 * f * (v - self.e_rev)

    def conductance(self, d: float, f: float) -> float:
        """Instantaneous conductance density (mS/cm²)."""
        return self.g_max * d ** 2 * f

    def initial_state(self, v_rest: float = -65.0) -> dict:
        """Return gate variables at resting potential."""
        return {"d": self.d_inf(v_rest), "f": self.f_inf(v_rest)}
