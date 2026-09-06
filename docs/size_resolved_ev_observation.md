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

The implementation is in
`electro_exocytosis/models/ev_size_observation.py`. The current calibration
driver is `tools/fit_ffrci_size_resolved.py`. Because this layer is standalone,
the legacy cumulative outputs and the version 1.1 three-class extracellular
outputs remain unchanged.

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

The calibration uses a broad 40–2000 nm latent domain. It contains 20 nm bins
through the common measurement window and progressively wider tail bins outside it.
The measured comparison is restricted to the common 80–380 nm range in 20 nm
bands. Modeling a broader latent range prevents probability just outside the
reported window from being forced artificially into the edge measurement
bins.

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

For the current fit, the common measurement bins are exact subsets of the
latent bins, instrument log-diameter error is fixed to zero, recovery is fixed
to one, and background is zero. The resulting response is identity-like within
80–380 nm and zero outside that reported window. Keeping the matrix explicit
makes this provisional assumption visible and provides a direct place to add
instrument pore size, detection efficiency, and calibration data later.

## Experimental alignment

### Provisional population and volume bridge

By default, the fit loads `examples/ffrci_experimental_bridge.yml`. It assumes
5,000,000 initial cells in 5 mL, or 1,000,000 initial cells/mL. The upstream
single-cell release is multiplied by this density, while the simulator's
time-varying viable-producer fraction modifies interval release. The bridge
file records unit recovery, no dilution correction, and no background, but the
current size-resolved driver uses only its initial cell count and volume. Its
observation operator separately fixes recovery to one and background to zero;
bridge-level recovery, dilution, and background are not yet wired into this
fit.

These values are provisional. Absolute pathway source scales should not be
interpreted biologically until cell count, conditioned-medium volume,
viability, dilution, recovery, and background handling are confirmed. The
dynamic concentration prediction is compared directly with particles/mL; it
is not converted to particles per cell and then scaled a second time.

### Batch aggregation and control initialization

Raw measurement histograms are first rebinned to the common 80–380 nm grid.
Within each condition, time, and diameter band, the experimental bridge then
exports the arithmetic mean, sample SD, SE, and non-missing measurement count.
The fit receives one mean distribution per condition and time; acquisition
identifiers are retained only in the raw audit table and do not appear in the
model-facing table or figures.

The three undated cell-containing control histograms are summarized in the
same way. Their all-batch mean distribution provisionally initializes the
extracellular state at time zero. Only common-window bins are populated
initially; the unobserved latent tails begin at zero. This ambient distribution
then experiences the same smooth size-dependent loss as newly released
particles. Because the source file does not establish the experimental-unit
metadata, the records are described conservatively as repeated batch
measurements rather than independent biological replicates.
Treatment and control tables are rebinned and summarized independently; an
acquisition identifier never pairs a treatment histogram to a control
histogram or enters the objective.

The experimental condition labels report 20 kV and 40 kV without confirming
whether these are generator voltages or electric-field strengths. Figures
retain those literal voltage labels. The current simulation provisionally uses
the numerical values as kV/cm; this is an explicit scenario assumption that
must be corrected if the protocol metadata establish another interpretation.

## Calibration objective

The driver deliberately separates concentration scale from distribution
shape:

- The total component contributes one natural-log residual for each
  condition/time batch mean. There are nine such targets: three exposure
  conditions at 0.5, 1, and 3 hours. Total concentration is summed within each
  raw histogram before calculating its mean and SD, which preserves covariance
  across size bands.
- The composition component normalizes each condition/time mean distribution
  to sum to one and compares square-root fractions. This is the geometry
  underlying the Hellinger distance and accepts zero-concentration bins without
  logarithms or pseudocounts.
- Weak penalties regularize pathway medians and widths, the size-loss slope,
  state shifts, and the optional dose correction around their prior centers.

Sample SD is exported and plotted as descriptive uncertainty but does not
inverse-variance weight the objective. With only two or three measurements per
cell, uncertain independence, zero-valued bands, and strong covariance across
diameter bands, such weighting would imply more precision than the data
support. SE is exported for downstream sensitivity analyses but is not the
default figure uncertainty.

The optimizer uses bounded, multistart nonlinear least squares. Model ranking
uses a descriptive BIC-like score calculated from the composite data
residual. Because the total and composition terms do not arise from a complete
sampling likelihood, this score is useful for comparison within this analysis
but is not a formal evidence criterion.

## Fit variants

The default command fits and compares three nested variants:

| Variant | Pathway kernels | Additional term |
|---|---|---|
| `static_kernel` | Fixed lognormal median and width for each pathway | None |
| `state_conditioned` | Median shifts with the three treatment-minus-sham state signals | Three pathway state-shift coefficients |
| `state_conditioned_dose_response` | State-conditioned kernels | Shared linear and quadratic dose-response correction |

The dose correction multiplies all three upstream cumulative-release
trajectories for a condition by

```text
gain(D) = exp(a * (D - 0.72) + b * (D - 0.72) ** 2)
```

where \(D\) is the simulator's dose index. It is included as a diagnostic for
the weak pulse-condition separation in the current intracellular model, not as
a validated constitutive mechanism. A large improvement from this variant is
evidence that the upstream dose-response layer needs revision, not proof of
the correction's biological form.

## Running the fit

From the repository root:

```bash
MPLCONFIGDIR=/tmp/electroexo-mpl PYTHONPATH=. python \
  tools/fit_ffrci_size_resolved.py \
  --experimental-bridge-config examples/ffrci_experimental_bridge.yml \
  --output-dir results/ffrci_size_resolved_fit_v1_1
```

The defaults run 24 starts with at most 100 function evaluations per start
and fit all three variants. One start is the canonical initial vector; the
others use a seeded Latin-hypercube design over the full configured optimizer
bounds. This is enough to make endpoint-stability box plots descriptive; it
does not turn their spread into parameter uncertainty. Use
`--variants` followed by one or more variant
names to run a subset. Pulse width and repetition rate remain explicit command
line assumptions through `--pulse-width-ns` and `--repetition-rate-hz`.

## Outputs

The output directory contains:

- `experimental_batch_summary.csv`: model-facing mean, sample SD, SE, and
  measurement count for every condition, time, and diameter band;
- `observed_vs_predicted_size_bins.csv`: observed and predicted 20 nm bins for
  every fitted variant;
- `total_concentration_fit.csv`: the nine longitudinal total-concentration
  targets and predictions;
- `size_distribution_fit_metrics.csv`: condition/time mean-distribution
  Hellinger distance and mean-diameter error;
- `pathway_size_kernels.csv`: pathway kernel probability, observed-window
  mass, loss rate, and dose gain by condition, time, and latent bin;
- `fitted_parameters.csv`: fitted values, bounds, and bound-contact flags for
  every variant;
- `fitted_parameters.yml`: the selected variant and its decoded parameter
  values;
- `model_comparison.csv`: descriptive composite-residual comparison of the
  three variants;
- `goodness_of_fit_by_condition.csv`: total-concentration and size-distribution
  metrics by condition and for the joint fit;
- `fit_diagnostics_by_target.csv` and `size_bin_fit_diagnostics.csv`:
  target-level and size-bin error components, including Hellinger and
  Wasserstein size distances;
- `optimization_starts.csv` and `optimization_endpoint_parameters.csv`: the
  canonical and full-bound Latin-hypercube starts, all resulting endpoints,
  objective components, bound positions, and the selected endpoint;
- `framework_parameter_snapshot.csv`: all authoritative runtime defaults plus
  the selected size-observation fit's initial values, optimizer bounds, and
  final values;
- `framework_model_registry.csv`: the implemented model/submodel inventory and
  each component's role in this calibration;
- `fit_summary.json`: bridge settings, scenario assumptions, objective
  definition, metrics, selected variant, and interpretation limits;
- `longitudinal_total_fit.png`, `size_profile_fit.png`, and
  `size_time_fit.png`: conventional total, profile, and size–time fit
  diagnostics for the selected variant;
- `size_time_surface_overlay.png`: the observed mean concentration surface
  with the selected fit overlaid as a translucent surface; and
- `size_time_fit_error_contours.png`: signed and absolute bounded normalized
  difference between the fitted model and experimental mean. The signed
  diagnostic is
  `100 * (predicted - observed) / (predicted + observed)`; it is not the
  optimizer's residual or a conventional percent error;
- `fit_error_diagnostics.png` and `goodness_of_fit.png`: complementary
  manuscript-style error, parity, and descriptive error-distribution views;
- `constitutive_and_kinetic_parameter_boxplots.png` and
  `parameter_space_boxplots.png`: multistart endpoint stability in physical and
  normalized optimizer coordinates; and
- `parameter_fit_summary.png`: the selected initial-to-final movement relative
  to each configured optimizer range.

The bounded normalization can make differences in low-concentration tail bins
look visually large. No detection-limit mask is applied because the instrument
limit of detection was not supplied; interpret the contours together with the
absolute concentration profiles and exported tables.

The driver also writes a concise `README.md` into the result directory.
All plots follow `docs/manuscript_figure_policy.md`. Mean line and point plots
show sample SD explicitly; surface, heatmap, and contour views use the mean and
are accompanied by the error-bar profile figure and the exported SD table.

The surface and contour figures are visual connections across the measured
size–time grid, not additional temporal data. Only 0.5, 1, and 3 hours were
observed. Surface faces and contour bands between those rows are rendering
interpolations between grid points; the overlaid dots identify measured mean
grid nodes, not individual observations. One condition/time cell contains two
measurements; the other cells contain three. These figures
must not be read as measured or independently predicted intermediate-time
kinetics. Bins where both observed and predicted concentration are zero are
undefined and masked rather than displayed as perfect agreement.

Parameter endpoint boxes describe repeated numerical optimization from
different initial guesses. They are not confidence intervals, posterior
samples, biological replicate distributions, or evidence that the fitted
parameter is identifiable. Search bounds are likewise engineering choices for
this fit, not validated biological ranges. The point-and-range summary and the
exported endpoint tables should be interpreted before using a fitted value in
a mechanistic argument.

## Interpretation limits

1. The cell-containing controls have no harvest time. Treating their mean as
   the time-zero extracellular distribution is a visible but unverified
   assumption.
2. The provenance and dependence structure of the repeated measurements are
   unknown. They must not be described as independent biological replicates or
   longitudinal cultures without confirmation.
3. The current instrument response is identity-like over the common 80–380 nm
   window. Pore size, detection efficiency, resolution, dilution, recovery,
   and background require instrument and protocol metadata.
4. The optional dose-response correction is a diagnostic adapter. It should
   not be promoted to a biological mechanism without independent validation.
5. Diameter alone cannot identify MVB-associated release, plasma-membrane
   budding, or apoptotic-body biogenesis. The fitted pathway decomposition is
   a non-unique reduced-order explanation that requires orthogonal markers,
   perturbations, imaging, or fractionation for validation.
6. Kernel shape, pathway scale, initial ambient particles, and size-dependent
   loss can compensate for one another. Parameter bounds and weak penalties do
   not remove this practical identifiability problem.
7. One condition/time histogram is absent, leaving 26 treated histograms rather
   than the otherwise expected 27. That mean and SD therefore have n=2; all
   other condition/time summaries have n=3.
8. The upstream equations, viability trajectory, pulse-unit interpretation,
   and assumed 60 ns/1 Hz exposure settings remain provisional. The fitted
   values are exploratory observation-layer parameters, not validated
   constitutive constants.
9. The 3D surfaces and 2D error contours connect observations at 0.5, 1, and
   3 hours. Their continuous appearance does not add observations between
   those harvest times.
