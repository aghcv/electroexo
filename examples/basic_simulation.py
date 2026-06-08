"""
Basic electro-exocytosis simulation example.

Simulates a neuron receiving a suprathreshold depolarising step, fires
several action potentials, triggers Ca²⁺ influx, and drives exocytosis
from the readily-releasable vesicle pool.
"""

from electroexo import SimulationEngine
from electroexo.stimulus import StepStimulus


def main():
    # ---- Build the simulation engine with default parameters ---------
    engine = SimulationEngine(dt=0.025, record_every=4)

    # ---- Define a 100 ms depolarising step ---------------------------
    stimulus = StepStimulus(amplitude=10.0, t_start=50.0, t_end=150.0)

    # ---- Run for 200 ms ----------------------------------------------
    results = engine.run(t_end=200.0, stimulus=stimulus)

    # ---- Report summary ----------------------------------------------
    print("=== Simulation Results ===")
    print(f"  Duration          : {results.t[-1]:.1f} ms")
    print(f"  Action potentials : {results.n_spikes}")
    print(f"  Peak voltage      : {results.peak_voltage:.1f} mV")
    print(f"  Peak [Ca²⁺]       : {results.peak_calcium:.4f} µM")
    print(f"  Release events    : {results.total_release_events:.2f} vesicles")


if __name__ == "__main__":
    main()
