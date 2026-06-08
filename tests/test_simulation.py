"""
Integration tests for the full SimulationEngine.

These tests exercise complete end-to-end runs and verify biological
plausibility of results (e.g., action potentials fire, calcium rises,
exocytosis increases under stimulation).
"""
import numpy as np
import pytest

from electroexo import SimulationEngine
from electroexo.models.calcium import CalciumDynamics
from electroexo.models.exocytosis import AllostericExocytosis, CooperativeExocytosis
from electroexo.models.membrane import Membrane
from electroexo.models.vesicle import VesiclePool
from electroexo.results import SimulationResults
from electroexo.stimulus import (
    ConstantStimulus,
    PulseTrainStimulus,
    RampStimulus,
    StepStimulus,
)


@pytest.fixture
def engine():
    return SimulationEngine(dt=0.05, record_every=10)


class TestSimulationEngineBasic:
    def test_returns_simulation_results(self, engine):
        results = engine.run(t_end=10.0)
        assert isinstance(results, SimulationResults)

    def test_time_array_starts_at_zero(self, engine):
        results = engine.run(t_end=10.0)
        assert results.t[0] == pytest.approx(0.0)

    def test_time_array_ends_near_t_end(self, engine):
        results = engine.run(t_end=10.0)
        assert results.t[-1] == pytest.approx(10.0, rel=0.05)

    def test_arrays_same_length(self, engine):
        results = engine.run(t_end=10.0)
        n = len(results.t)
        for arr in [results.v, results.ca, results.docked, results.released]:
            assert len(arr) == n

    def test_resting_voltage_without_stimulus(self, engine):
        results = engine.run(t_end=50.0, stimulus=ConstantStimulus(0.0))
        # Voltage should stay near resting potential
        assert abs(results.v[-1] - (-65.0)) < 10.0

    def test_gate_variables_in_range(self, engine):
        results = engine.run(t_end=20.0)
        for arr in [results.m, results.h, results.n, results.d, results.f]:
            assert np.all(arr >= 0.0)
            assert np.all(arr <= 1.0)

    def test_calcium_non_negative(self, engine):
        results = engine.run(t_end=20.0)
        assert np.all(results.ca >= 0.0)

    def test_released_non_decreasing(self, engine):
        results = engine.run(t_end=50.0, stimulus=StepStimulus(10.0, 10.0, 40.0))
        diffs = np.diff(results.released)
        assert np.all(diffs >= -1e-10)

    def test_vesicle_pools_non_negative(self, engine):
        results = engine.run(t_end=20.0)
        assert np.all(results.docked >= 0.0)
        assert np.all(results.reserve >= 0.0)


class TestSimulationEngineStimulation:
    def test_step_stimulus_elicits_spike(self):
        """A suprathreshold step current should elicit at least one AP."""
        engine = SimulationEngine(dt=0.025, record_every=4)
        results = engine.run(
            t_end=100.0, stimulus=StepStimulus(amplitude=10.0, t_start=10.0, t_end=90.0)
        )
        assert results.n_spikes >= 1

    def test_peak_voltage_increases_with_stimulus(self):
        engine = SimulationEngine(dt=0.025, record_every=4)
        r_no = engine.run(t_end=50.0, stimulus=ConstantStimulus(0.0))
        r_stim = engine.run(t_end=50.0, stimulus=StepStimulus(10.0, 5.0, 45.0))
        assert r_stim.peak_voltage > r_no.peak_voltage

    def test_calcium_rises_during_stimulation(self):
        engine = SimulationEngine(dt=0.025, record_every=4)
        r_no = engine.run(t_end=50.0, stimulus=ConstantStimulus(0.0))
        r_stim = engine.run(t_end=50.0, stimulus=StepStimulus(10.0, 5.0, 45.0))
        assert r_stim.peak_calcium > r_no.peak_calcium

    def test_release_increases_during_stimulation(self):
        engine = SimulationEngine(dt=0.025, record_every=4)
        r_stim = engine.run(t_end=100.0, stimulus=StepStimulus(10.0, 10.0, 90.0))
        assert r_stim.total_release_events > 0.0

    def test_pulse_train_stimulus(self):
        engine = SimulationEngine(dt=0.025, record_every=4)
        stim = PulseTrainStimulus(
            amplitude=12.0, pulse_duration=1.0, period=20.0, t_start=10.0, n_pulses=4
        )
        results = engine.run(t_end=100.0, stimulus=stim)
        assert results.n_spikes >= 1

    def test_ramp_stimulus(self):
        engine = SimulationEngine(dt=0.025, record_every=4)
        stim = RampStimulus(amplitude_start=0.0, amplitude_end=15.0, t_start=10.0, t_end=80.0)
        results = engine.run(t_end=100.0, stimulus=stim)
        assert results.peak_voltage > -60.0


class TestSimulationEngineCustomisation:
    def test_custom_exocytosis_model(self):
        engine = SimulationEngine(
            exocytosis_model=AllostericExocytosis(l0=1e-3),
            dt=0.05,
            record_every=10,
        )
        results = engine.run(t_end=10.0)
        assert isinstance(results, SimulationResults)

    def test_deplete_reserve_mode(self):
        engine = SimulationEngine(
            vesicle_pool=VesiclePool(n_reserve=50, deplete_reserve=True),
            dt=0.05,
            record_every=10,
        )
        results = engine.run(
            t_end=100.0, stimulus=StepStimulus(12.0, 10.0, 90.0)
        )
        # Reserve should decrease during heavy stimulation
        assert results.reserve[-1] <= results.reserve[0]

    def test_release_rate_length(self):
        engine = SimulationEngine(dt=0.05, record_every=10)
        results = engine.run(t_end=20.0)
        rate = results.release_rate()
        assert len(rate) == len(results.t) - 1


class TestSimulationResultsProperties:
    @pytest.fixture
    def results(self):
        engine = SimulationEngine(dt=0.05, record_every=10)
        return engine.run(
            t_end=100.0, stimulus=StepStimulus(10.0, 10.0, 90.0)
        )

    def test_total_release_events_non_negative(self, results):
        assert results.total_release_events >= 0.0

    def test_peak_voltage_is_max(self, results):
        assert results.peak_voltage == pytest.approx(np.max(results.v))

    def test_peak_calcium_is_max(self, results):
        assert results.peak_calcium == pytest.approx(np.max(results.ca))

    def test_n_spikes_non_negative(self, results):
        assert results.n_spikes >= 0
