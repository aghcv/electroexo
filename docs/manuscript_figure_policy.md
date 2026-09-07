# Manuscript figure policy

This policy applies to every new or regenerated scientific figure in the
electro-exocytosis pipeline. Existing plotting functions should adopt it when
they are next modified. It separates experimental summarization from plotting
so that internal file names and acquisition identifiers do not leak into a
manuscript figure.

## Experimental observations

- Define the independent experimental unit before aggregation. Repeated
  acquisition channels are not independent biological replicates unless the
  experiment metadata establish that they are.
- Once the experimental unit is established, plot the finite observations as
  the arithmetic mean with the sample standard deviation (SD, `ddof=1`). SD is
  the default because these figures describe variability among the available
  measurement sets. Do not silently substitute standard error (SE).
- Whenever space permits, retain the individual observations as faint
  points behind the mean and SD. This is especially important for fewer than
  five batches.
- A single available observation may be plotted, but it has no SD error bar.
  Record its sample count in the exported summary table and manuscript caption.
- The computational fit is evaluated against the batch mean. The same batch
  aggregation must be used consistently in the objective, metrics, tables, and
  figures.

Use `aggregate_repeated_observations` from `experimental_bridge` for tabular
fitting inputs. Use `summarize_repeated_observations` for arrays and
`plot_observed_mean_sd` for line or point figures:

```python
from electro_exocytosis.visualization.style import (
    FITTED_COLOR,
    FITTED_MODEL_LABEL,
    place_manuscript_legend,
    plot_observed_mean_sd,
)

plot_observed_mean_sd(axis, diameter_nm, batch_concentrations)
axis.plot(
    diameter_nm,
    fitted_concentration,
    color=FITTED_COLOR,
    label=FITTED_MODEL_LABEL,
)
place_manuscript_legend(figure, axis)
```

## Labels and terminology

- Use short scientific labels that describe what is displayed: `Observed mean
  ± SD`, `Fitted model`, and `Fit error` are the shared defaults.
- Do not display source-institution names, spreadsheet names, worksheet names,
  sample codes, or internal acquisition identifiers.
- Experimental conditions may be named by interpretable treatment variables,
  such as pulse count, field strength, and elapsed time, provided those values
  are confirmed by metadata.
- Put physical units in parentheses in axis labels. Avoid redundant panel
  titles and legends.

## Legends

- Legends normally sit outside the plotting region so they cannot obscure data
  or error bars. A multi-panel figure may instead place its one shared legend
  in demonstrably unused panel space when this materially improves the printed
  scale. Such an inset legend must have an opaque white background, a visible
  boundary, and no overlap with observations, error bars, or model curves.
- A single-panel figure has one legend outside the axis, normally on the right.
- A multi-panel figure has exactly one deduplicated figure-level legend shared
  by all panels. Do not repeat the legend in every panel.
- Use `place_manuscript_legend(figure, axes)`. It infers single versus
  multi-panel layout, removes duplicate axis legends, and places the shared
  legend on the right. Use `location="bottom"` only when a right-side legend
  would make a wide panel unreadable.
- Save through `save_manuscript_figure`, which uses a tight bounding box so
  external legends are included in the output.

## Visual encoding

- Observed data use `OBSERVED_COLOR`; fitted trajectories use `FITTED_COLOR`;
  signed or scalar fit errors use `FIT_ERROR_COLOR`.
- Do not use color as the only distinction. The existing `line_style` and
  `line_styles` helpers pair colors with line styles and markers.
- Use linear axes unless multiplicative variation or a broad dynamic range
  makes a logarithmic scale scientifically necessary. State log scales in the
  axis label or caption and never plot nonpositive values on them.
- Apply `manuscript_style_context` around figure construction and
  `style_manuscript_axis` to standardize font sizes, guide lines, tick marks,
  labels, titles, and axis scales.
- Use a white background and the colorblind-safe central palette. New tools
  must not define private palettes without a documented reason.
- Mean curves and point series must include SD bars. Raster, contour, and
  surface views may show the mean without bars when bars would obscure the
  encoded field, but the same figure family must provide variability through
  a companion panel, conventional error-bar figure, or exported SD table.

## Figure notes and captions

- Do not place prose footers, abbreviation lists, sample-count notes, formulas,
  or methodological qualifications inside manuscript image files. Put that
  information in the LaTeX caption and generated data tables, where it remains
  readable and editable.
- Reserve in-figure text for concise panel titles, axis labels, legends, and
  data annotations that directly encode a plotted value.
- `add_figure_note` remains available for exploratory or standalone diagnostic
  exports, but manuscript-bound plots must omit it.

## Physical dimensions and typography

- Generate figures at their intended printed width rather than creating a
  12--16-inch canvas and shrinking it in LaTeX. Use approximately 3.45 inches
  for one-column figures and 7.10 inches for two-column figures.
- Judge type size after manuscript placement. Axis labels and legends should
  remain at least 8 points at final size; tick labels should remain at least
  7.5 points. Higher raster resolution improves line sharpness but does not
  compensate for physically undersized text.
- Reclaim space before reducing type: remove redundant supertitles, share axis
  labels where appropriate, shorten repeated labels, and use otherwise empty
  panel regions for a framed shared legend.

## Fit diagnostics and parameter plots

- Error and goodness-of-fit figures must keep total-concentration and
  size-composition residual families separate unless a statistically explicit
  joint scale is defined. Display every condition–time target when the design
  is small.
- Agreement plots use an identity line and the experimental SD on the observed
  axis. Metrics calculated on fitted data are labeled as in-sample descriptive
  diagnostics, not predictive validation.
- A single fitted parameter is shown as an initial/final point with its search
  bounds, not as a box plot.
- A parameter box plot may summarize repeated optimizer endpoints only when at
  least 12 successful multistart endpoints are available; 20–30 are preferred.
  Show the individual endpoints and the selected solution. State explicitly
  that their spread measures optimization stability, not a confidence
  interval, posterior distribution, or biological variation.
- Parameter values with incompatible units may share a panel only after a
  clearly labeled transformation such as position within the configured
  optimizer range. Search bounds are not biological validity ranges.
- With fewer than 12 successful starts, retain the endpoint points and
  initial/final markers but omit the box glyphs.

## Required generator pattern

Every new or modified plot generator should:

1. aggregate repeated observations in the experimental bridge before fitting;
2. use `manuscript_style_context` and the semantic observed/model/error colors;
3. use generic labels and physical units;
4. place one shared legend with `place_manuscript_legend`, or use a framed
   in-panel legend only when it occupies verified empty space;
5. save through `save_manuscript_figure`; and
6. keep prose qualifications in the caption rather than the image; and
7. include a focused render or smoke test that verifies the output is created.

The current default is 1200 dpi for line art. Color-dense surfaces, heatmaps,
and contours may use `MANUSCRIPT_COLOR_DPI` (600 dpi) to keep file sizes
practical while retaining publication resolution.
