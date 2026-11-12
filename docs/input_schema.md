# Scenario and parameter input

Use [the short demo](../examples/scenario_demo.yml) as a starting point.
Scenario YAML is loaded with `load_scenario` and validated by
`electro_exocytosis.config.SimulationScenario`. The installed Python package
can emit the complete structural schema.

```bash
python -c 'import json; from electro_exocytosis.config import SimulationScenario; print(json.dumps(SimulationScenario.model_json_schema(), indent=2))'
```

| Section | Fields and units |
| --- | --- |
| `scenario` | Required `name`. Optional `mode` is `cell_based_electro_exocytosis` or `direct_EV_engineering`. |
| `pulse` | Required peak `amplitude_kV_cm` in kV/cm, `pulse_width_ns` in ns, integer `pulse_number`, and `repetition_rate_Hz` in Hz. Optional `waveform` is `square`, `bipolar`, or `exponential`. |
| `exposure` | Required section. Defaults are available for geometry, conductivity in S/m, temperature in °C, cell density per mL, and dosimetry. |
| `cell_state` | Optional cell label and dimensionless modifiers for membrane, calcium handling, baseline release, and stress. |
| `simulation` | Optional `t_start_s`, `t_end_s`, and `output_dt_s` in seconds. End time must follow start time and the output interval must be positive. `numerical_method` is `solve_ivp` or `euler`. |
| `extracellular_medium` | Optional initial volume in mL, viable-fraction setting, and sampling events with time in seconds and volumes in mL. |

The schema describes types and field-level constraints. Pydantic also checks
the time-window relationship when it builds the scenario. See
[extracellular kinetics](extracellular_ev_kinetics.md) for event handling.

Pass sparse parameter overrides with `--parameters-file PATH`. The runtime
deep-merges them with the packaged defaults. See
[parameter defaults and snapshots](parameter_defaults_and_snapshots.md) for
resolution order, parameter status, and output snapshots.
