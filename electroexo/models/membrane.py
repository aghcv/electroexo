"""
electroexo.models.membrane
============================
Hodgkin-Huxley-based cell membrane model extended with an L-type Ca²⁺ channel.

The membrane integrates

    C_m dV/dt = I_ext(t) − I_Na − I_K − I_Ca − I_L

where:

* *I_Na* — fast Na⁺ channel current (m³h kinetics)
* *I_K*  — delayed-rectifier K⁺ channel current (n⁴ kinetics)
* *I_Ca* — L-type Ca²⁺ channel current (d²f kinetics)
* *I_L*  — passive leak current

The channel gate variables (*m, h, n, d, f*) are evolved alongside *V* as
part of the full state vector returned by :meth:`Membrane.derivatives`.

Units
-----
Voltage : mV
Time    : ms
Capacitance : µF/cm²
Current density : µA/cm²
Conductance density : mS/cm²
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np

from electroexo.models.channels import CalciumChannel, PotassiumChannel, SodiumChannel

__all__ = ["Membrane", "MembraneState"]


@dataclass
class MembraneState:
    """Snapshot of membrane electrical state variables.

    Attributes
    ----------
    v : float
        Membrane potential (mV).
    m : float
        Na⁺ channel activation gate.
    h : float
        Na⁺ channel inactivation gate.
    n : float
        K⁺ channel activation gate.
    d : float
        Ca²⁺ channel activation gate.
    f : float
        Ca²⁺ channel inactivation gate.
    """

    v: float
    m: float
    h: float
    n: float
    d: float
    f: float

    def as_array(self) -> np.ndarray:
        """Return state as a 1-D NumPy array ``[v, m, h, n, d, f]``."""
        return np.array([self.v, self.m, self.h, self.n, self.d, self.f])

    @classmethod
    def from_array(cls, arr: np.ndarray) -> "MembraneState":
        """Construct from a 1-D array ``[v, m, h, n, d, f]``."""
        return cls(v=arr[0], m=arr[1], h=arr[2], n=arr[3], d=arr[4], f=arr[5])


class Membrane:
    """Conductance-based membrane model with Na⁺, K⁺, Ca²⁺, and leak currents.

    Parameters
    ----------
    c_m : float
        Membrane capacitance density (µF/cm²).  Default 1.0.
    g_leak : float
        Leak conductance density (mS/cm²).  Default 0.3.
    e_leak : float
        Leak reversal potential (mV).  Default −54.4.
    v_rest : float
        Resting membrane potential used to initialise gate variables (mV).
        Default −65.0.
    na_channel : SodiumChannel, optional
        Na⁺ channel instance.  Defaults to :class:`~electroexo.models.channels.SodiumChannel`.
    k_channel : PotassiumChannel, optional
        K⁺ channel instance.  Defaults to :class:`~electroexo.models.channels.PotassiumChannel`.
    ca_channel : CalciumChannel, optional
        Ca²⁺ channel instance.  Defaults to :class:`~electroexo.models.channels.CalciumChannel`.
    """

    def __init__(
        self,
        c_m: float = 1.0,
        g_leak: float = 0.3,
        e_leak: float = -54.4,
        v_rest: float = -65.0,
        na_channel: SodiumChannel | None = None,
        k_channel: PotassiumChannel | None = None,
        ca_channel: CalciumChannel | None = None,
    ) -> None:
        self.c_m = c_m
        self.g_leak = g_leak
        self.e_leak = e_leak
        self.v_rest = v_rest

        self.na = na_channel if na_channel is not None else SodiumChannel()
        self.k = k_channel if k_channel is not None else PotassiumChannel()
        self.ca = ca_channel if ca_channel is not None else CalciumChannel()

    # ------------------------------------------------------------------
    # Currents
    # ------------------------------------------------------------------

    def leak_current(self, v: float) -> float:
        """Passive leak current density (µA/cm²)."""
        return self.g_leak * (v - self.e_leak)

    def total_ionic_current(
        self,
        v: float,
        m: float,
        h: float,
        n: float,
        d: float,
        f: float,
    ) -> float:
        """Sum of all ionic current densities (µA/cm²)."""
        i_na = self.na.current(v, m, h)
        i_k = self.k.current(v, n)
        i_ca = self.ca.current(v, d, f)
        i_l = self.leak_current(v)
        return i_na + i_k + i_ca + i_l

    def calcium_current(self, v: float, d: float, f: float) -> float:
        """Ca²⁺ channel current density (µA/cm²); negative = inward."""
        return self.ca.current(v, d, f)

    # ------------------------------------------------------------------
    # ODE derivatives
    # ------------------------------------------------------------------

    def derivatives(
        self,
        state: MembraneState,
        i_ext: float,
    ) -> MembraneState:
        """Compute dState/dt for all membrane variables.

        Parameters
        ----------
        state :
            Current membrane state.
        i_ext :
            Applied external current density (µA/cm²).

        Returns
        -------
        MembraneState
            Rate of change of each state variable.
        """
        v, m, h, n, d, f = state.v, state.m, state.h, state.n, state.d, state.f

        i_ion = self.total_ionic_current(v, m, h, n, d, f)
        dv = (i_ext - i_ion) / self.c_m

        dm = self.na.dm_dt(v, m)
        dh = self.na.dh_dt(v, h)
        dn = self.k.dn_dt(v, n)
        dd = self.ca.dd_dt(v, d)
        df = self.ca.df_dt(v, f)

        return MembraneState(v=dv, m=dm, h=dh, n=dn, d=dd, f=df)

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def initial_state(self) -> MembraneState:
        """Return the membrane state at rest."""
        v = self.v_rest
        na0 = self.na.initial_state(v)
        k0 = self.k.initial_state(v)
        ca0 = self.ca.initial_state(v)
        return MembraneState(
            v=v,
            m=na0["m"],
            h=na0["h"],
            n=k0["n"],
            d=ca0["d"],
            f=ca0["f"],
        )
