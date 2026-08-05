# Solution Space Analysis

This example bundle explores the current reduced `electro_exocytosis` model over
a bounded, reasonable nsPEF design space:

- field amplitude: 5-30 kV/cm
- pulse width: 50-200 ns
- pulse count: 5-60 pulses
- repetition rate: 1-20 Hz
- medium conductivity: 0.5-2.0 S/m
- waveform: square, bipolar, exponential
- dosimetry model: legacy, Joule adiabatic, Joule lumped thermal
- cell-state modifiers: 0.7x, 1.0x, and 1.3x around the nominal state

Regenerate the scenario YAML files and aggregate result figures from the repo
root:

```bash
python examples/solution_space_analysis/generate_solution_space.py \
  --out results/solution_space_analysis
```

The script writes `scenario_*.yml` files in this folder and writes aggregate
tables, rankings, and figures under `results/solution_space_analysis`.

The analysis remains exploratory. The model uses reduced, uncalibrated
mechanistic equations, so these figures are for understanding the behavior of
the current simulator rather than for selecting real experimental conditions.
