# Experimental-to-Computational Observation Bridge

## Purpose

The mechanistic electro-exocytosis simulator represents the response and EV
output of one representative producer cell. Particle instruments commonly
report a concentration in conditioned medium. These quantities must not be
compared directly.

`electro_exocytosis.experimental_bridge.ExperimentalObservationBridge`
provides the standard boundary between experimental particle measurements and
the single-cell model. The bridge is deliberately separate from the
intracellular equations: cell number, culture volume, viability normalization,
assay dilution, recovery, and background are properties of the experiment and
observation process, not constitutive intracellular biology.

## Forward and Inverse Transformations

For measured concentration \(C_{meas}\), background \(C_{bg}\), dilution
back-calculation factor \(D\), recovery fraction \(\eta\), conditioned-medium
volume \(V\), and cell denominator \(N_{basis}\), the single-cell-equivalent
observation is

```text
C_corrected = max(C_meas - C_bg, 0) * D / eta
particles_per_cell = C_corrected * V / N_basis
```

The inverse transformation is

```text
C_meas = C_bg + particles_per_cell * N_basis / V * eta / D
```

The cell denominator is selected explicitly:

- `cell_basis: initial` uses the number of cells entering the experiment. This
  is the default for the current representative-cell model and avoids
  artificially increasing output per cell when treatment reduces viability.
- `cell_basis: viable` uses `initial_cell_count * viability_fraction`. Use this
  only when matched viability or viable-cell counts are measured and the
  scientific endpoint is explicitly per viable cell.

For a time-resolved extracellular stock with uptake or loss, the more complete
observation model is

```text
dC_ext/dt = viable_cell_density(t) * release_rate_per_cell(t)
            - k_loss * C_ext(t)
C_measured = recovery * C_ext / dilution_factor + background
```

The bridge still performs the reversible static unit transformation used by
the legacy FFRCI fitting tools. The version 1.1 simulator now implements the
time-resolved stock equation separately in `models/extracellular_kinetics.py`.
Keeping the two interfaces distinct avoids applying the same cell-count or
volume conversion twice. See `docs/extracellular_ev_kinetics.md`.

## Standard YAML Configuration

```yaml
experimental_bridge:
  initial_cell_count: 5000000
  medium_volume_ml: 5.0
  cell_basis: initial
  viability_fraction: 1.0
  recovery_fraction: 1.0
  dilution_factor: 1.0
  background_concentration_particles_per_ml: 0.0
```

Definitions:

| Field | Required interpretation |
|---|---|
| `initial_cell_count` | Total starting producer cells in the culture represented by the measured sample |
| `medium_volume_ml` | Conditioned-medium volume corresponding to the reported concentration |
| `cell_basis` | `initial` or `viable`; must match the stated endpoint |
| `viability_fraction` | Matched fraction in (0, 1]; only changes conversion for the viable-cell basis |
| `recovery_fraction` | Fraction retained by collection, preparation, and measurement |
| `dilution_factor` | Factor multiplying the instrument result to recover pre-dilution concentration |
| `background_concentration_particles_per_ml` | Blank or medium background to subtract before scaling |

Unknown quantities must remain visible assumptions. Do not silently set them
through biological release parameters.

## Python Interface

```python
from electro_exocytosis.experimental_bridge import ExperimentalObservationBridge

bridge = ExperimentalObservationBridge.from_yaml(
    "examples/ffrci_experimental_bridge.yml"
)

particles_per_cell = bridge.concentration_to_particles_per_cell(6.2e9)
concentration = bridge.particles_per_cell_to_concentration(particles_per_cell)

single_cell_table = bridge.transform_frame(
    experimental_table,
    concentration_column="concentration_particles_per_ml",
)
```

The forward and inverse functions accept scalars or arrays and validate all
scale factors. `to_metadata()` returns the supplied assumptions and derived
cell density for inclusion in result manifests.

## Calibration Command

The FFRCI fit driver accepts the bridge directly:

```bash
MPLCONFIGDIR=/tmp/electroexo-mpl PYTHONPATH=. python \
  tools/fit_ffrci_ev_kinetics.py \
  --experimental-bridge-config examples/ffrci_experimental_bridge.yml \
  --output-dir results/ffrci_ev_kinetics_fit_per_cell
```

Add `--max-particle-diameter-nm 200` for the small-particle sensitivity
analysis. The normalized fit accepts the same bridge option.

Every resulting observation table retains both
`measured_concentration_particles_per_ml` and the transformed calibration
value, and every summary JSON records the complete bridge configuration.

## Normalized Comparisons

When stimulated and sham samples have identical cell count, medium volume,
viability basis, recovery, dilution, and background handling, the shared
population scale cancels:

```text
(particles_per_cell_stim / particles_per_cell_sham)
    = (concentration_stim / concentration_sham)
```

This cancellation is invalid when treatment changes viable cell number, sample
volume, recovery, dilution, or background. Sample-specific bridge settings are
then required. A sham without a collection time cannot establish a
time-matched denominator even when its cell count is known.

## Interpretation Rules

1. Report the original concentration and the transformed single-cell value.
2. Treat cell count, volume, dilution, recovery, and background as measured
   experimental inputs, not fitted intracellular parameters.
3. State whether normalization uses starting cells, viable cells at harvest,
   or viable-cell-hours.
4. Do not interpret particle-diameter bins as biogenesis-defined EV subtypes.
5. Do not infer an extracellular clearance rate from a static unit conversion;
   use the dynamic extracellular stock and fit an unresolved effective loss.
6. Propagate uncertainty in cell count, volume, recovery, and viability into
   parameter intervals once those measurements are available.
7. Reject or sensitivity-test any fit whose constitutive rates change strongly
   when plausible bridge values are varied.

## Current Limitations

- The bridge currently uses one static metadata set per loaded dataset. The
  class can be instantiated per sample, but the FFRCI command-line driver does
  not yet read a row-specific sample manifest.
- The intracellular EV states remain cumulative by design, but the version 1.1
  extracellular stock can rise or fall and is the appropriate supernatant-fit
  target.
- The current model uses reduced particle-equivalent states rather than a
  traceable absolute particle-count calibration.
- Viability, recovery, dilution, background, and culture volume are unconfirmed
  for the supplied FFRCI file.
- The standard scenario interface now supports model-derived time-varying
  viability, extracellular loss, sampling, and medium replacement. Directly
  measured viable-cell trajectories and observation uncertainty are not yet
  accepted as row-specific inputs.
- The legacy FFRCI fitting commands have not yet been migrated to optimize the
  new extracellular stock; their published output folders remain v1.0 analyses.
