"""
Pulse-train stimulation example.

Delivers a train of short high-frequency current pulses and measures
cumulative vesicle release as a function of pulse count.
"""

from electroexo import SimulationEngine
from electroexo.models.exocytosis import CooperativeExocytosis
from electroexo.models.vesicle import VesiclePool
from electroexo.stimulus import PulseTrainStimulus


def main():
    pool = VesiclePool(n_reserve=200, n_docked=30, deplete_reserve=True)
    exo = CooperativeExocytosis(k_max=5e-2, k_d=1.5, n_hill=4)
    engine = SimulationEngine(
        vesicle_pool=pool,
        exocytosis_model=exo,
        dt=0.025,
        record_every=4,
    )

    # 10 pulses, 2 ms each, at 50 Hz (20 ms inter-pulse interval)
    stimulus = PulseTrainStimulus(
        amplitude=12.0,
        pulse_duration=2.0,
        period=20.0,
        t_start=20.0,
        n_pulses=10,
    )

    results = engine.run(t_end=300.0, stimulus=stimulus)

    print("=== Pulse Train Stimulation ===")
    print(f"  Action potentials   : {results.n_spikes}")
    print(f"  Peak [Ca²⁺]         : {results.peak_calcium:.4f} µM")
    print(f"  Cumulative release  : {results.total_release_events:.2f} vesicles")
    print(f"  Remaining in RRP    : {results.docked[-1]:.2f} vesicles")
    print(f"  Reserve depletion   : {results.reserve[0] - results.reserve[-1]:.2f} vesicles")


if __name__ == "__main__":
    main()
