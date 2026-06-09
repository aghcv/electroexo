# Input schema

## `scenario`
- `name` (`str`, required): scenario label used in outputs.
- `mode` (`cell_based_electro_exocytosis | direct_EV_engineering`, default `cell_based_electro_exocytosis`).

## `pulse`
- `amplitude_kV_cm` (`float`, required, `> 0`): peak electric field.
- `pulse_width_ns` (`float`, required, `> 0`): pulse width in nanoseconds.
- `pulse_number` (`int`, required, `>= 1`): number of pulses.
- `repetition_rate_Hz` (`float`, required, `> 0`): pulse repetition rate.
- `waveform` (`square | bipolar | exponential`, default `square`).

## `exposure`
- `geometry` (`cuvette | dish | flow`, default `cuvette`).
- `medium_conductivity_S_m` (`float`, default `1.5`, `> 0`).
- `temperature_C` (`float`, default `37.0`).
- `cell_density_per_ml` (`float`, default `1e6`, `> 0`).

## `cell_state`
- `cell_type` (`str`, default `generic`).
- `membrane_modifier` (`float`, default `1.0`).
- `calcium_handling_modifier` (`float`, default `1.0`).
- `baseline_EV_release_modifier` (`float`, default `1.0`).
- `stress_sensitivity_modifier` (`float`, default `1.0`).

## `simulation`
- `t_start_s` (`float`, default `0.0`).
- `t_end_s` (`float`, default `7200.0`).
- `output_dt_s` (`float`, default `1.0`).
- `numerical_method` (`solve_ivp | euler`, default `solve_ivp`).

## Notes
- Missing `cell_state` and `simulation` sections are filled with defaults.
- Validation is performed with Pydantic v2 during scenario loading.
- Parameter overrides are handled separately from the scenario YAML.
