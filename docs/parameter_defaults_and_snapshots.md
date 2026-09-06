# Parameter defaults, sparse overrides, and run snapshots

The framework already uses a canonical-default plus sparse-override pattern.
The packaged baseline is
`electro_exocytosis/parameters/default_parameters.yaml`. A normal simulation
loads that complete mapping, recursively overlays only values supplied in a
run-specific parameter file, applies configured cell-state modifiers, and
writes the effective run parameters to `parameters_used.yaml`.

For example, a run file can change one extracellular loss term without copying
the other framework parameters:

```yaml
extracellular_kinetics:
  effective_loss_rate_s: 0.00077
```

```bash
electro-exocytosis run examples/basic_scenario.yml \
  --parameters-file my_parameter_changes.yml \
  --out results/my_run
```

This keeps the comprehensive standard in one version-controlled package file
and makes each job input a readable record of intentional differences. The
resolution order is:

```text
packaged defaults -> sparse run override -> cell-state modifier -> derived outputs
```

## Auditable registries

`electro_exocytosis.parameter_registry` provides two read-only table builders:

- `build_parameter_snapshot()` flattens all 241 scalar runtime defaults and can
  join effective overrides and a selected calibration's initial values, search
  bounds, and final values;
- `build_model_registry()` lists the implemented layers, modules, submodules,
  reduced formulations, runtime roles, and their role in a particular fit.

The size-resolved fitting pipeline writes these as
`framework_parameter_snapshot.csv` and `framework_model_registry.csv`. The
parameter snapshot separates mechanistic runtime defaults from the fitted
size-observation layer. For fixed runtime parameters, `effective_value` is the
value used by the simulation; `fit_initial`, `fit_lower_bound`,
`fit_upper_bound`, and `fit_final` are populated only when a parameter was
actually optimized. Blank fit bounds mean “not curated/not fitted,” not an
unbounded biological quantity.

The snapshot also flags known parameters that are declared but not consumed by
the current equations. This is intentional audit information: a value's
presence in the default YAML does not by itself prove that the running model
uses it.

## Current boundary between simulation and observation fitting

The general simulator and the size-resolved calibration currently have
separate parameter namespaces. The latter fits pathway-to-particle source
conversion, lognormal size-kernel descriptors, size-dependent loss, and
optional state/dose adapters. Its `fitted_parameters.yml` is a calibration
artifact; it is not a generic simulator override file. Keeping that boundary
explicit avoids collisions and avoids implying that an observation scale is
an intracellular constitutive rate.

A future consolidated resolver should use qualified paths, for example
`ev_release.baseline_sEV_rate` and
`ev_size_observation.sEV_median_diameter_nm`, and should reject unknown paths.
It should also record the requested value, effective post-modifier value,
scientific validity range, calibration search range, unit, provenance, and
evidence status separately.

## Scope of the remaining work

The software mechanism for canonical defaults, sparse overrides, and complete
snapshots is now present and is not a major redesign. Two larger scientific and
validation tasks remain:

1. curate defensible units, validity domains, and literature provenance for
   every parameter rather than inventing broad numerical ranges; and
2. add strict qualified-path/type/range validation after resolving existing
   aliases and declared-but-unused values, so older parameter files can be
   migrated deliberately instead of silently changing behavior.

Optimizer search bounds in the current fit are engineering constraints. They
must not be reported as biological reference ranges, and multistart endpoint
spread must not be reported as parameter uncertainty.
