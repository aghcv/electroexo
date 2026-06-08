"""
electroexo.stimulus
====================
Electrical stimulus protocols applied to the membrane.

Each protocol exposes a single method:

    current = stimulus.amplitude_at(t)   →  µA/cm²

Available protocols
-------------------
* :class:`ConstantStimulus`  — constant DC injection
* :class:`StepStimulus`      — rectangular step (on/off)
* :class:`PulseTrainStimulus` — repetitive rectangular pulses
* :class:`SinusoidalStimulus` — sinusoidal AC current
* :class:`RampStimulus`       — linearly ramped current

Units
-----
Current density : µA/cm²
Time            : ms
"""

from __future__ import annotations

import math
from typing import List, Optional

__all__ = [
    "BaseStimulus",
    "ConstantStimulus",
    "StepStimulus",
    "PulseTrainStimulus",
    "SinusoidalStimulus",
    "RampStimulus",
]


class BaseStimulus:
    """Abstract base class for stimulus protocols."""

    def amplitude_at(self, t: float) -> float:
        """Return the stimulus current density at time *t* (µA/cm²).

        Parameters
        ----------
        t :
            Simulation time (ms).
        """
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Concrete stimulus classes
# ---------------------------------------------------------------------------

class ConstantStimulus(BaseStimulus):
    """Constant (DC) current injection.

    Parameters
    ----------
    amplitude : float
        Current density (µA/cm²).  Default 0.
    """

    def __init__(self, amplitude: float = 0.0) -> None:
        self.amplitude = amplitude

    def amplitude_at(self, t: float) -> float:
        return self.amplitude


class StepStimulus(BaseStimulus):
    """Rectangular step stimulus.

    The current is *amplitude* during ``[t_start, t_end)`` and zero otherwise.

    Parameters
    ----------
    amplitude : float
        Current density during the step (µA/cm²).  Default 10.0.
    t_start : float
        Step onset time (ms).  Default 50.0.
    t_end : float
        Step offset time (ms).  Default 150.0.
    """

    def __init__(
        self,
        amplitude: float = 10.0,
        t_start: float = 50.0,
        t_end: float = 150.0,
    ) -> None:
        self.amplitude = amplitude
        self.t_start = t_start
        self.t_end = t_end

    def amplitude_at(self, t: float) -> float:
        return self.amplitude if self.t_start <= t < self.t_end else 0.0


class PulseTrainStimulus(BaseStimulus):
    """Periodic train of rectangular current pulses.

    Parameters
    ----------
    amplitude : float
        Pulse current density (µA/cm²).  Default 10.0.
    pulse_duration : float
        Duration of each pulse (ms).  Default 1.0.
    period : float
        Period between pulse onsets (ms).  Default 20.0.
    t_start : float
        Time of first pulse onset (ms).  Default 0.0.
    t_end : float
        Time after which no further pulses are delivered (ms).
        If ``None``, the train continues indefinitely.
    n_pulses : int | None
        Maximum number of pulses to deliver.  If ``None``, unlimited.
    """

    def __init__(
        self,
        amplitude: float = 10.0,
        pulse_duration: float = 1.0,
        period: float = 20.0,
        t_start: float = 0.0,
        t_end: Optional[float] = None,
        n_pulses: Optional[int] = None,
    ) -> None:
        self.amplitude = amplitude
        self.pulse_duration = pulse_duration
        self.period = period
        self.t_start = t_start
        self.t_end = t_end
        self.n_pulses = n_pulses

    def amplitude_at(self, t: float) -> float:
        if t < self.t_start:
            return 0.0
        if self.t_end is not None and t >= self.t_end:
            return 0.0

        elapsed = t - self.t_start
        pulse_index = int(elapsed / self.period)

        if self.n_pulses is not None and pulse_index >= self.n_pulses:
            return 0.0

        phase = elapsed - pulse_index * self.period
        return self.amplitude if phase < self.pulse_duration else 0.0


class SinusoidalStimulus(BaseStimulus):
    """Sinusoidal alternating current stimulus.

    Parameters
    ----------
    amplitude : float
        Peak current density (µA/cm²).  Default 5.0.
    frequency : float
        Frequency (kHz; 1 kHz = 1 cycle/ms).  Default 0.05 (50 Hz).
    dc_offset : float
        DC offset added to the sinusoid (µA/cm²).  Default 0.0.
    t_start : float
        Stimulus onset (ms).  Default 0.0.
    t_end : float | None
        Stimulus offset (ms); ``None`` means run indefinitely.
    """

    def __init__(
        self,
        amplitude: float = 5.0,
        frequency: float = 0.05,
        dc_offset: float = 0.0,
        t_start: float = 0.0,
        t_end: Optional[float] = None,
    ) -> None:
        self.amplitude = amplitude
        self.frequency = frequency
        self.dc_offset = dc_offset
        self.t_start = t_start
        self.t_end = t_end

    def amplitude_at(self, t: float) -> float:
        if t < self.t_start:
            return 0.0
        if self.t_end is not None and t >= self.t_end:
            return 0.0
        return self.dc_offset + self.amplitude * math.sin(
            2.0 * math.pi * self.frequency * (t - self.t_start)
        )


class RampStimulus(BaseStimulus):
    """Linearly ramped current stimulus.

    The current increases linearly from *amplitude_start* to *amplitude_end*
    over ``[t_start, t_end]`` and is zero outside that window.

    Parameters
    ----------
    amplitude_start : float
        Current density at ramp onset (µA/cm²).  Default 0.0.
    amplitude_end : float
        Current density at ramp offset (µA/cm²).  Default 20.0.
    t_start : float
        Ramp onset (ms).  Default 50.0.
    t_end : float
        Ramp offset (ms).  Default 200.0.
    """

    def __init__(
        self,
        amplitude_start: float = 0.0,
        amplitude_end: float = 20.0,
        t_start: float = 50.0,
        t_end: float = 200.0,
    ) -> None:
        self.amplitude_start = amplitude_start
        self.amplitude_end = amplitude_end
        self.t_start = t_start
        self.t_end = t_end

    def amplitude_at(self, t: float) -> float:
        if t < self.t_start or t > self.t_end:
            return 0.0
        frac = (t - self.t_start) / (self.t_end - self.t_start)
        return self.amplitude_start + frac * (self.amplitude_end - self.amplitude_start)
