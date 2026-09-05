# FFRCI EV Kinetics: Single-Cell Reanalysis

## Executive Finding

The original FFRCI kinetic fits compared a representative-cell simulation
directly with Exoid particles/mL. That comparison had a unit mismatch. The
baseline release-rate parameters absorbed the missing cell-density factor and
could not be interpreted as per-cell constitutive rates.

The three analyses have now been repeated after converting experimental
concentrations to particle equivalents per initial cell through the standard
experimental observation bridge. This corrects the population scale but does
not correct the separate structural mismatch: simulated release is cumulative,
whereas the experimental supernatant concentration rises and falls.

> Historical-version note: the results in this document were produced with the
> frozen version 1.0 cumulative-output workflow. Version 1.1 subsequently added
> the extracellular stock recommended below, but these folders have not been
> refitted with that state and their reported metrics are unchanged.

## Current Framework Audit

### Representative-cell mechanistic core

The ion, remodeling, and EV-release ODEs describe one representative cell. The
EV output states integrate subtype release fluxes into cumulative small-EV,
medium/large-EV, and apoptotic-body particle equivalents.

### Previously disconnected experimental fields

`ExposureConfig.cell_density_per_ml` existed but was not used to scale the EV
time series. `manufacturing_qc.cell_count` was used only to divide an
already-computed endpoint into a cell-normalized manufacturing metric; it did
not convert representative-cell release into particles/mL. The kinetic fitter
therefore optimized single-cell output directly against concentration.

### New experimental observation bridge

The new bridge performs a reversible and recorded mapping between particles/mL
and particles/cell using cell count, volume, viability basis, recovery,
dilution, and background. It is integrated into both FFRCI fitting tools and is
documented in `docs/experimental_observation_bridge.md`.

## Provisional FFRCI Scaling

The following assumptions were used:

| Quantity | Value | Status |
|---|---:|---|
| Initial cell count | 5,000,000 | supplied by the PI in project discussion |
| Medium volume | 5.0 mL | provisional inference from the repository's previous 1,000,000 cells/mL setting |
| Cell basis | initial cells | chosen to match the representative-cell model |
| Viability fraction | 1.0 | placeholder; no time-matched measurement supplied |
| Recovery fraction | 1.0 | placeholder; no workflow recovery supplied |
| Dilution factor | 1.0 | assumes Exoid values are already back-calculated |
| Background | 0 particles/mL | placeholder; no validated blank subtraction supplied |

The resulting density is 1,000,000 initial cells/mL, so a concentration of
\(10^9\) particles/mL becomes 1,000 particle equivalents per initial cell.

The absolute scaling changes linearly with the metadata:

```text
new_particles_per_cell / reported_particles_per_cell
    = (new_volume / 5 mL) * (5e6 / new_initial_cell_count)
      * (new_dilution / new_recovery)
```

The YAML configuration is `examples/ffrci_experimental_bridge.yml`.

## Repeated Analyses

### 1. All reported particle sizes, absolute per-cell fit

Output: `results/ffrci_ev_kinetics_fit_per_cell`

- Fitted log10 RMSE: 0.3276.
- Median absolute percentage error: 46.5%.
- Observed range after conversion: approximately 2,700 to 13,700 particle
  equivalents per initial cell.
- Dominant fitted baseline medium/large-EV rate: approximately 60.2 in current
  per-cell model units, versus approximately \(6.02\times10^7\) in the invalid
  concentration-scale fit.

### 2. Particles below 200 nm, absolute per-cell fit

Output: `results/ffrci_ev_kinetics_fit_lt200nm_per_cell`

- Fitted log10 RMSE: 0.3292.
- Median absolute percentage error: 46.9%.
- Dominant fitted baseline medium/large-EV rate: approximately 54.9 in current
  per-cell model units.
- Filtering below 200 nm does not materially improve agreement.

### 3. Below-200-nm stimulated/sham ratio fits

Output: `results/ffrci_ev_kinetics_normalized_lt200nm_per_cell`

- `Sham2` reference: log10 fold RMSE 0.1910; median absolute percentage error
  14.9%.
- `sham media` reference: log10 fold RMSE 0.2457; median absolute percentage
  error 34.3%.
- The common cell-density conversion cancels from the ratio because the same
  bridge was provisionally applied to stimulated and sham samples.
- The controls still cannot be treated as time-matched because neither has a
  harvest-time label.

## Interpretation

The unchanged fit errors are expected. Multiplying every observation by a
constant translates log observations without changing trajectory shape; an
adjustable release scale can compensate exactly. The value of the reanalysis
is therefore not a better curve but a corrected constitutive parameter scale.

The optimizer still places several parameters at their allowed bounds and
routes most fitted output through the medium/large-EV pathway—even for the
below-200-nm dataset. This should not be interpreted as evidence that the
measured particles are medium/large EVs. Diameter alone does not identify
biogenesis, and the current observation function sums all modeled subtypes.

For the historical version 1.0 fits, the remaining failure was structural:

- cumulative model outputs cannot decrease;
- no extracellular uptake, degradation, aggregation, adsorption, or sampling
  loss is modeled;
- no time-varying viable-cell population is applied;
- controls are not time matched;
- pulse metadata and the meaning of the two sham labels remain provisional.

## Required Metadata Before Biological Parameter Interpretation

1. Confirm whether five million is the starting, treated, recovered, or viable
   cell count for every sample.
2. Confirm conditioned-medium volume at each collection time.
3. Provide viable-cell counts or viability at 0.5, 1, and 3 hours by condition.
4. Confirm Exoid dilution back-calculation, particle recovery, and background
   subtraction.
5. Identify which sham maps to each treatment and harvest time.
6. Confirm whether time points are destructive cultures, cumulative medium, or
   serial sampling with replacement.

## Implemented Next Step and Remaining Calibration Work

Version 1.1 keeps the bridge and intracellular pool definitions intact and adds
the extracellular particle-stock observation model:

```text
dC_ext,m/dt = viable_cell_density(t) * release_rate_per_cell,m(t)
              - k_loss,m * C_ext,m(t)
```

The implementation supports a shared effective loss, separate named mechanisms,
class multipliers, sampling impulses, medium replacement, assay transformation,
and model-derived viability. The next step is to connect the FFRCI optimizer to
this output, beginning with one shared effective loss and a measured or
interpolated viable-cell trajectory. Add subtype-specific loss or additional
mechanisms only if the confirmed experimental design can identify them.
