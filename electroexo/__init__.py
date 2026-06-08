"""
electroexo — A general computational framework for electro-exocytosis simulation.

Provides biophysical models of neuronal membrane dynamics, ion channels,
intracellular calcium handling, vesicle pool kinetics, and calcium-triggered
exocytosis, together with configurable electrical-stimulus protocols and
numerical ODE solvers.

Typical usage
-------------
>>> from electroexo import SimulationEngine
>>> from electroexo.stimulus import StepStimulus
>>> engine = SimulationEngine()
>>> results = engine.run(t_end=200.0, stimulus=StepStimulus(amplitude=10.0,
...                                                         t_start=50.0,
...                                                         t_end=150.0))
>>> print(results.total_release_events)
"""

from electroexo.simulation import SimulationEngine
from electroexo.results import SimulationResults

__version__ = "0.1.0"
__all__ = ["SimulationEngine", "SimulationResults"]
