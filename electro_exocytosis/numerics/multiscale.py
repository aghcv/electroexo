from __future__ import annotations

from electro_exocytosis.models.pulse import PulseDescriptors
from electro_exocytosis.numerics.solvers import build_multiscale_time_array


class MultiscaleScheduler:
    """
    Manages the multi-rate temporal structure:
    - Fast phase: pulse events (nanosecond, resolved as descriptors)
    - Medium phase: membrane/Ca2+ transients (milliseconds to seconds)
    - Slow phase: EV release and recovery (seconds to hours)
    """

    def __init__(self, pulse_descriptors: PulseDescriptors, t_start: float, t_end: float):
        self.pulse_descriptors = pulse_descriptors
        self.t_start = t_start
        self.t_end = t_end

    def get_time_array(self, dt: float):
        """Return the output time array for simulation."""
        return build_multiscale_time_array(self.t_start, self.t_end, dt)

    def get_pulse_event_times(self) -> list[float]:
        """Return simplified pulse event times beginning at t=0."""
        if self.pulse_descriptors.pulse_number <= 1:
            return [0.0]
        interval = 1.0 / self.pulse_descriptors.repetition_rate_Hz
        return [i * interval for i in range(self.pulse_descriptors.pulse_number)]
