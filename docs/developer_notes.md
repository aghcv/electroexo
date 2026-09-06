# Developer notes

## Extending the code
- Replace placeholder equations inside `electro_exocytosis/models/` module by module.
- Keep parameter names synchronized between dataclasses, `default_parameters.yaml`, and `simulation.py`.
- Use `apply_cell_state_modifiers` for new cross-cutting state scalings.
- Keep the packaged defaults authoritative and use sparse run overrides. Update
  `electro_exocytosis.parameter_registry` when a model/submodel or known runtime
  status changes; see `docs/parameter_defaults_and_snapshots.md`.

## Evidence integration
- `EvidenceLoader` currently reads workbook sheets and exposes them as pandas DataFrames.
- Future development can translate workbook rows into parameter overrides or model-selection logic.

## Testing
Run the package in editable mode and execute:

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## Numerical stability
- The current solver uses `solve_ivp` with RK45.
- The pulse is compressed into descriptor-driven initial forcing rather than explicit nanosecond time stepping.
- Keep new couplings bounded to preserve smoke-test stability over multi-hour simulations.
