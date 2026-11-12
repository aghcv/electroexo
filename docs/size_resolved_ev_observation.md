# Size-Resolved EV Observation and Batch-Mean Fit

## Purpose and scope

The size-resolved observation update connects the existing intracellular EV
release model to experimentally measured particle-diameter distributions.
It does not replace or reinterpret the established `sEV`, `mlEV`, and `AB`
trajectories. Those remain reduced mechanistic release pathways. The update
adds a standalone forward model that:

1. distributes each pathway's release over a latent diameter domain;
2. propagates each latent size bin through extracellular loss;
3. maps the latent distribution into instrument-facing size bins; and
4. fits total concentration and conditional size composition as distinct
   sources of information.

The forward model is implemented in
`electro_exocytosis/models/ev_size_observation.py`. Its inputs are size-grid
edges, cumulative pathway release, cell density, and observation parameters.

## Forward model

### Pathway-to-size lognormal kernel

For pathway \(p\), the release increment over a simulation interval is spread
across latent diameter bins using a lognormal kernel. Its fitted baseline
parameters are a median diameter \(d_{50,p}\), a geometric standard deviation
\(\sigma_{g,p}>1\), and a pathway-specific conversion from one upstream model
unit to particle equivalents.

The probability assigned to a latent bin is the lognormal probability mass
between that bin's edges. The probabilities are normalized over the modeled
latent domain, so they sum to one and conserve the pathway release increment
within that domain. This is a conditional distribution over the represented
range; it is not evidence that physical tail mass outside the range is zero.

### State-conditioned median shifts

The state-conditioned variants allow the pathway median to change with an
upstream treatment-versus-sham state signal:

```text
median_p(t) = baseline_median_p * exp(state_shift_p * state_signal_p(t))
```

The interval midpoint signal is used for the release increment in that
interval. The three signals are constructed from the existing simulator:

- `sEV`: the treatment-minus-sham change in the mean of the ESCRT-dependent
  and ceramide signals;
- `mlEV`: the treatment-minus-sham change in the mean of the budding and
  scission signals; and
- `AB`: the treatment-minus-sham change in apoptotic-blebbing signal.

A positive state-shift coefficient moves the median toward larger diameters
when its state signal is positive; a negative coefficient moves it toward
smaller diameters. These shifts describe a reduced observation kernel. They do
not establish that an observed diameter bin was produced by a particular
biogenesis pathway.

### Broad latent size domain and extracellular loss

A latent grid can extend beyond the measurement window, with finer bins
inside the measured range and wider bins in the tails. This keeps mass outside
the reported window from being assigned to edge measurement bins.

The model applies a smooth first-order loss law at each latent bin's geometric
center \(d_b\):

```text
k_loss(d_b) = log(2) / half_life_150nm * (d_b / 150 nm) ** loss_size_exponent
```

Thus the fitted half-life is defined at 150 nm. A positive exponent makes
larger particles disappear faster, while a negative exponent makes smaller
particles disappear faster. Initial ambient particles and newly released
particles experience the same loss law. Release input and loss are integrated
analytically within each interval, using midpoint viability and midpoint state
signals.

The result retains the ambient concentration, each pathway-attributed latent
concentration, their total, the kernel probabilities, and the bin-specific loss
rates. The present reduced model does not include aggregation, fragmentation,
or other transfers between latent size bins.

### Explicit observation matrix

Latent concentrations are converted to reported bins with an explicit
observed-bin-by-latent-bin matrix. With zero instrument width, a matrix entry is
the fraction of a latent bin's linear-diameter interval that overlaps an
observed bin. With a nonzero log-diameter error, the operator instead assigns
probability using a lognormal measurement-error model centered on each latent
bin.

The forward observation also supports assay recovery and scalar, per-bin, or
time-by-bin background concentration. Columns of the observation matrix may
sum below one when latent particles fall outside the observed diameter window
or measurement error scatters them outside it.

Choose the latent and measured size grids explicitly. A broad latent grid can
preserve tail mass around a narrower measurement window. Estimate instrument
width, recovery, and background from the relevant measurement procedure.

## Calibration workflow

Aggregate repeated measurements at each condition and time. Summarize total
concentration within each histogram before combining repeats so the total
uncertainty retains covariance across size bands. Compare total concentration
and normalized size composition separately. A square-root composition residual
admits zero-valued bins. Use explicit bounds and a recorded seed for multistart
optimization.

Record the input size grids, initial concentration, cell density, recovery and
background, model revision, and all fitted parameter bounds alongside the
outputs. Use the observation and extracellular stock APIs as separate model
components and apply the culture-scale conversion once.
