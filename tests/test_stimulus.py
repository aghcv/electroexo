"""Tests for electroexo.stimulus."""
import pytest

from electroexo.stimulus import (
    ConstantStimulus,
    PulseTrainStimulus,
    RampStimulus,
    SinusoidalStimulus,
    StepStimulus,
)


class TestConstantStimulus:
    def test_always_returns_amplitude(self):
        s = ConstantStimulus(amplitude=5.0)
        for t in [-10, 0, 50, 1000]:
            assert s.amplitude_at(t) == 5.0

    def test_default_zero(self):
        s = ConstantStimulus()
        assert s.amplitude_at(100) == 0.0


class TestStepStimulus:
    def setup_method(self):
        self.s = StepStimulus(amplitude=10.0, t_start=50.0, t_end=150.0)

    def test_zero_before_onset(self):
        assert self.s.amplitude_at(49.9) == 0.0

    def test_amplitude_during_step(self):
        assert self.s.amplitude_at(50.0) == 10.0
        assert self.s.amplitude_at(100.0) == 10.0

    def test_zero_after_offset(self):
        assert self.s.amplitude_at(150.0) == 0.0
        assert self.s.amplitude_at(200.0) == 0.0


class TestPulseTrainStimulus:
    def setup_method(self):
        self.s = PulseTrainStimulus(
            amplitude=10.0,
            pulse_duration=1.0,
            period=10.0,
            t_start=0.0,
        )

    def test_first_pulse(self):
        assert self.s.amplitude_at(0.0) == 10.0
        assert self.s.amplitude_at(0.5) == 10.0

    def test_between_pulses(self):
        assert self.s.amplitude_at(1.5) == 0.0
        assert self.s.amplitude_at(9.9) == 0.0

    def test_second_pulse(self):
        assert self.s.amplitude_at(10.0) == 10.0

    def test_n_pulses_limit(self):
        s = PulseTrainStimulus(
            amplitude=10.0,
            pulse_duration=1.0,
            period=10.0,
            t_start=0.0,
            n_pulses=2,
        )
        assert s.amplitude_at(20.0) == 0.0  # 3rd pulse should be blocked

    def test_t_end_cutoff(self):
        s = PulseTrainStimulus(
            amplitude=10.0,
            pulse_duration=1.0,
            period=10.0,
            t_start=0.0,
            t_end=15.0,
        )
        assert s.amplitude_at(20.0) == 0.0

    def test_before_t_start(self):
        s = PulseTrainStimulus(t_start=100.0)
        assert s.amplitude_at(50.0) == 0.0


class TestSinusoidalStimulus:
    def setup_method(self):
        self.s = SinusoidalStimulus(
            amplitude=5.0, frequency=0.05, t_start=0.0, t_end=200.0
        )

    def test_starts_near_zero(self):
        assert abs(self.s.amplitude_at(0.0)) < 1e-10

    def test_before_t_start(self):
        s = SinusoidalStimulus(t_start=50.0)
        assert s.amplitude_at(30.0) == 0.0

    def test_after_t_end(self):
        assert self.s.amplitude_at(201.0) == 0.0

    def test_amplitude_bounded(self):
        import math
        for t in range(0, 200):
            v = self.s.amplitude_at(float(t))
            assert abs(v) <= self.s.amplitude + 1e-10

    def test_dc_offset(self):
        s = SinusoidalStimulus(amplitude=0.0, dc_offset=3.0)
        assert s.amplitude_at(50.0) == pytest.approx(3.0)


class TestRampStimulus:
    def setup_method(self):
        self.s = RampStimulus(
            amplitude_start=0.0, amplitude_end=20.0, t_start=50.0, t_end=200.0
        )

    def test_zero_before_ramp(self):
        assert self.s.amplitude_at(30.0) == 0.0

    def test_start_value(self):
        assert self.s.amplitude_at(50.0) == pytest.approx(0.0)

    def test_end_value(self):
        assert self.s.amplitude_at(200.0) == pytest.approx(20.0)

    def test_midpoint(self):
        assert self.s.amplitude_at(125.0) == pytest.approx(10.0)

    def test_zero_after_ramp(self):
        assert self.s.amplitude_at(250.0) == 0.0
