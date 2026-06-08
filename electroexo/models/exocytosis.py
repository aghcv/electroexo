"""
electroexo.models.exocytosis
==============================
Calcium-triggered exocytosis rate models.

Two models are provided:

* :class:`CooperativeExocytosis` — sigmoidal Hill-function model with
  cooperative Ca²⁺ binding (suitable for fast, synchronous release).
* :class:`AllostericExocytosis` — allosteric/primed-release model that
  includes a basal spontaneous release rate (suitable for slow, asynchronous
  release and tonic exocytosis).

Both classes expose a uniform interface:

    rate = model.exocytosis_rate(ca)   →  ms⁻¹

Units
-----
Concentration : µM
Rate          : ms⁻¹
"""

from __future__ import annotations

import numpy as np

__all__ = ["CooperativeExocytosis", "AllostericExocytosis"]


class CooperativeExocytosis:
    """Hill-function model of cooperative Ca²⁺-triggered exocytosis.

    The exocytosis rate is

        k_exo([Ca]) = k_max · [Ca]^n_hill / (K_d^n_hill + [Ca]^n_hill)

    Parameters
    ----------
    k_max : float
        Maximal exocytosis rate (ms⁻¹).  Default 1e-2.
    k_d : float
        Half-maximal Ca²⁺ concentration (µM).  Default 2.0.
    n_hill : float
        Hill coefficient (cooperativity).  Default 4.0.
    k_basal : float
        Basal (spontaneous) exocytosis rate (ms⁻¹).  Default 0.
    """

    def __init__(
        self,
        k_max: float = 1e-2,
        k_d: float = 2.0,
        n_hill: float = 4.0,
        k_basal: float = 0.0,
    ) -> None:
        self.k_max = k_max
        self.k_d = k_d
        self.n_hill = n_hill
        self.k_basal = k_basal

    def exocytosis_rate(self, ca: float) -> float:
        """Return the Ca²⁺-dependent exocytosis rate (ms⁻¹).

        Parameters
        ----------
        ca :
            Free cytosolic Ca²⁺ concentration (µM).
        """
        ca = max(ca, 0.0)
        ca_n = ca ** self.n_hill
        kd_n = self.k_d ** self.n_hill
        return self.k_basal + self.k_max * ca_n / (kd_n + ca_n)


class AllostericExocytosis:
    """Allosteric release model with spontaneous and Ca²⁺-driven rates.

    Based on the allosteric model of Kochubey & Bhatt (2011):

        k_exo([Ca]) = L₀ · (1 + [Ca] / K_Ca)^n_sites

    where *L₀* is the basal opening rate of the release machinery and
    *K_Ca* is the Ca²⁺-binding affinity of each SNARE sensor site.

    Parameters
    ----------
    l0 : float
        Basal exocytosis rate (ms⁻¹).  Default 1e-4.
    k_ca : float
        Ca²⁺ affinity of each sensor site (µM).  Default 3.0.
    n_sites : int
        Number of cooperative Ca²⁺-binding sites.  Default 5.
    k_max : float
        Optional saturation clamp (ms⁻¹); ``None`` means no clamp.
    """

    def __init__(
        self,
        l0: float = 1e-4,
        k_ca: float = 3.0,
        n_sites: int = 5,
        k_max: float | None = None,
    ) -> None:
        self.l0 = l0
        self.k_ca = k_ca
        self.n_sites = n_sites
        self.k_max = k_max

    def exocytosis_rate(self, ca: float) -> float:
        """Return the Ca²⁺-dependent exocytosis rate (ms⁻¹).

        Parameters
        ----------
        ca :
            Free cytosolic Ca²⁺ concentration (µM).
        """
        ca = max(ca, 0.0)
        rate = self.l0 * (1.0 + ca / self.k_ca) ** self.n_sites
        if self.k_max is not None:
            rate = min(rate, self.k_max)
        return rate
