# Extracellular EV Kinetics and Observation Layer

## Purpose

The intracellular EV-release model produces nonnegative release fluxes per
representative cell. Its cumulative outputs can only increase. Experimental
particle instruments instead observe a stock in the extracellular medium, and
that stock can increase or decrease as secretion competes with disappearance,
sampling, dilution, and changes in the viable producer population.

Version 1.1 therefore retains the legacy cumulative outputs and adds a distinct
extracellular-medium model in
`electro_exocytosis/models/extracellular_kinetics.py`. No extracellular loss is
enabled by default, so a default run remains compatible with version 1.0.

## Model

For modeled class (m\in\{sEV, mlEV, AB\}), the state (N_m) is the number of
particle equivalents in the complete culture volume:

```text
dN_m/dt = N_initial_cells * viable_fraction(t) * source_scale * release_rate_m(t)
          - k_m * N_m(t)
```

The class loss rate is

```text
k_m = class_multiplier_m * (k_effective + k_uptake + k_degradation + k_adsorption)
```

An optional reduced aggregation term removes small-particle equivalents and
adds fewer medium/large-particle equivalents according to an aggregation yield.
The true and assay-facing concentrations are

```text
C_m(t) = N_m(t) / medium_volume(t)
C_measured(t) = recovery * sum_m C_m(t) / dilution_factor + background
```

The three model classes are mechanistic pathway classes. They must not be
identified with Exoid diameter bins without orthogonal evidence of biogenesis.

## Choosing Loss Parameters

Use `effective_loss_rate_s` for the first fit to a concentration trajectory.
The EV data alone generally constrain the sum of losses, not whether particles
were taken up by cells, degraded, adsorbed to labware, or lost during handling.
Keeping this parameter explicitly unresolved prevents an optimizer from
assigning an unsupported biological interpretation.

Only separate the following terms after adding discriminating measurements:

| Parameter | Interpretation | Useful discriminating experiment |
|---|---|---|
| `uptake_rate_s` | cell-associated removal | labeled EV uptake with and without cells or uptake inhibition |
| `degradation_rate_s` | loss of intact detectable particles | cell-free conditioned-medium stability time course |
| `adsorption_rate_s` | surface and handling loss | spike-in recovery across tubes, plates, and transfers |
| `sEV_to_mlEV_aggregation_rate_s` | small-particle disappearance with size redistribution | matched time-resolved size distributions and orthogonal particle characterization |

`sEV_loss_multiplier`, `mlEV_loss_multiplier`, and `AB_loss_multiplier` permit
class-specific loss. Leave them at one until the data can identify differential
loss. A rate (k) corresponds to half-life `log(2) / k`.

## Sampling and Medium Replacement

Collection events belong in the scenario file:

```yaml
extracellular_medium:
  initial_volume_ml: 5.0
  use_time_varying_viability: true
  sampling_events:
    - time_s: 1800
      sampled_volume_ml: 0.5
      replacement_volume_ml: 0.5
```

The state reported at an event time represents the withdrawn, well-mixed
sample. Removal and particle-free replacement are then applied to subsequent
states. Set both event volumes to zero for destructive, independent cultures.
Do not encode independent biological replicates as serial samples from one
culture.

When `use_time_varying_viability` is true, the source is multiplied by the
model's time-resolved viable fraction. A measured viable-cell trajectory is
preferable for calibration and should eventually replace this model-derived
proxy when available.

## Running an Experiment-Specific Parameter Set

The parameter file is deep-merged over package defaults:

```bash
python -m electro_exocytosis run \
  examples/scenario_rubens_experiment.yml \
  --parameters-file examples/parameters_extracellular_decline.yml \
  --out results/extracellular_decline_example
```

`examples/parameters_extracellular_decline.yml` uses an illustrative 15-minute
effective half-life solely to demonstrate the new behavior. It is not an FFRCI
fit.

## Outputs and Compatibility

`ev_outputs.csv` retains `sEV_cumulative`, `mlEV_cumulative`, and
`AB_cumulative`. It also contains true subtype stocks and concentrations,
total true concentration, assay-facing concentration, culture volume, viable
producer fraction, source rate, and mechanism-specific sink fluxes. The
standard plot `extracellular_ev_kinetics.png` is generated for every run.

Cargo, quality, and manufacturing endpoints still consume the legacy
cumulative quantities in version 1.1. The extracellular stock is the correct
target for fitting supernatant kinetics; changing those downstream interfaces
should be a separate, validated change.

## Calibration Guardrails

1. Treat cell density, volume, sampling, dilution, recovery, and background as
   measured experimental inputs, not constitutive rates.
2. Begin with a source scale and one effective loss rate. Add more degrees of
   freedom only when profile likelihood or held-out prediction supports them.
3. Include an initial extracellular concentration when the first time point is
   not a true pre-stimulation zero.
4. Fit independent cultures jointly but do not connect them through sampling
   events.
5. Report the effective half-life and its uncertainty; do not label it uptake
   or degradation without a mechanism-specific assay.
6. Preserve both latent true concentration and the assay-transformed output in
   every comparison.
