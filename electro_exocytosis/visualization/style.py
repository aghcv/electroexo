from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib")
)

import matplotlib
import numpy as np
from matplotlib.axes import Axes
from matplotlib.container import ErrorbarContainer
from matplotlib.figure import Figure
from matplotlib.legend import Legend

from electro_exocytosis.abbreviations import STANDARD_ABBREVIATIONS

MANUSCRIPT_DPI = 1200
MANUSCRIPT_COLOR_DPI = 600
LANDSCAPE_ASPECT_RATIO = 16 / 9
MANUSCRIPT_LANDSCAPE_FIGSIZE = (8.0, 4.5)
MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE = (12.0, 6.75)
MAX_MONOCHROME_SERIES = 3
LINE_WIDTH = 1.8
MARKER_SIZE = 4.2
MARK_EVERY = 0.12
ERROR_BAR_CAPSIZE = 3.0
ERROR_BAR_LINE_WIDTH = 1.1
LEGEND_FONT_SIZE = 8.0
FIGURE_NOTE_FONT_SIZE = 7.0

# Semantic colors and labels should be stable across analysis pipelines.  The
# labels deliberately describe the plotted quantity rather than the study,
# instrument, file name, or internal acquisition identifiers.
OBSERVED_COLOR = "#4477AA"
FITTED_COLOR = "#EE6677"
FIT_ERROR_COLOR = "#AA3377"
OBSERVED_MEAN_LABEL = "Observed mean"
OBSERVED_MEAN_SD_LABEL = "Observed mean ± SD"
FITTED_MODEL_LABEL = "Fitted model"
FIT_ERROR_LABEL = "Fit error"
LEGEND_LABEL_ORDER = (
    OBSERVED_MEAN_SD_LABEL,
    OBSERVED_MEAN_LABEL,
    FITTED_MODEL_LABEL,
    FIT_ERROR_LABEL,
)

MANUSCRIPT_RCPARAMS: dict[str, Any] = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "font.size": 9.0,
    "axes.labelsize": 9.0,
    "axes.titlesize": 10.0,
    "axes.titleweight": "normal",
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": LEGEND_FONT_SIZE,
    "lines.linewidth": LINE_WIDTH,
    "lines.markersize": MARKER_SIZE,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
}

MONOCHROME_LINE_STYLES = (
    {"color": "#111111", "linestyle": "-", "marker": "o"},
    {"color": "#555555", "linestyle": "--", "marker": "s"},
    {"color": "#888888", "linestyle": "-.", "marker": "^"},
)

COLORBLIND_PALETTE = (
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
    "#000000",
)

LINE_FORMATS = (
    {"linestyle": "-", "marker": "o"},
    {"linestyle": "--", "marker": "s"},
    {"linestyle": "-.", "marker": "^"},
    {"linestyle": ":", "marker": "D"},
    {"linestyle": "-", "marker": "v"},
    {"linestyle": "--", "marker": "P"},
    {"linestyle": "-.", "marker": "X"},
    {"linestyle": ":", "marker": "*"},
)

BAR_HATCHES = ("", "//", "\\\\", "xx", "..", "++", "oo", "**")


@dataclass(frozen=True)
class RepeatedObservationSummary:
    """Pointwise mean, sample SD, and finite observation count.

    ``sd`` uses Bessel's correction (``ddof=1``).  It is ``NaN`` wherever
    fewer than two finite observations are available, because variability
    cannot be estimated from a single observation.
    """

    mean: np.ndarray
    sd: np.ndarray
    n: np.ndarray


def landscape_figsize(
    width_in: float = MANUSCRIPT_LANDSCAPE_FIGSIZE[0],
) -> tuple[float, float]:
    """Return a standard 16:9 landscape figure size."""
    return (width_in, width_in / LANDSCAPE_ASPECT_RATIO)


def manuscript_style_context() -> Any:
    """Return a context manager for the standard manuscript typography."""
    return matplotlib.rc_context(rc=MANUSCRIPT_RCPARAMS)


def summarize_repeated_observations(
    values: Any,
    *,
    axis: int = 0,
) -> RepeatedObservationSummary:
    """Summarize repeated observations as a pointwise mean and sample SD.

    Non-finite values are ignored pointwise.  The returned count therefore
    makes missingness explicit and lets callers avoid presenting an estimated
    error bar when only one observation remains.
    """
    observations = np.asarray(values, dtype=float)
    if observations.ndim == 0:
        raise ValueError("Repeated observations must have at least one dimension.")
    if not -observations.ndim <= axis < observations.ndim:
        raise ValueError(
            f"axis={axis} is invalid for an array with {observations.ndim} dimensions."
        )
    axis = axis % observations.ndim

    finite = np.isfinite(observations)
    count = np.sum(finite, axis=axis)
    total = np.sum(np.where(finite, observations, 0.0), axis=axis)
    mean = np.full(np.shape(total), np.nan, dtype=float)
    np.divide(total, count, out=mean, where=count > 0)

    centered = observations - np.expand_dims(mean, axis=axis)
    squared_deviation = np.sum(np.where(finite, centered**2, 0.0), axis=axis)
    variance = np.full(np.shape(total), np.nan, dtype=float)
    np.divide(squared_deviation, count - 1, out=variance, where=count > 1)
    sd = np.sqrt(variance)
    return RepeatedObservationSummary(mean=mean, sd=sd, n=count.astype(int))


def plot_observed_mean_sd(
    ax: Axes,
    x: Any,
    repeated_values: Any,
    *,
    axis: int = 0,
    label: str = OBSERVED_MEAN_SD_LABEL,
    color: str = OBSERVED_COLOR,
    **errorbar_kwargs: Any,
) -> ErrorbarContainer:
    """Plot the mean and sample SD of repeated observations.

    The default visual is a connected point estimate with capped error bars.
    At positions with only one finite observation, the point is retained but
    no variability bar is drawn.
    """
    summary = summarize_repeated_observations(repeated_values, axis=axis)
    x_values = np.asarray(x)
    if summary.mean.ndim != 1:
        raise ValueError(
            "The summarized observations must be one-dimensional for plotting."
        )
    if x_values.ndim != 1 or x_values.shape != summary.mean.shape:
        raise ValueError(
            "x must be one-dimensional and match the summarized observations."
        )

    plot_kwargs: dict[str, Any] = {
        "color": color,
        "marker": "o",
        "linestyle": "-",
        "linewidth": LINE_WIDTH,
        "markersize": MARKER_SIZE,
        "capsize": ERROR_BAR_CAPSIZE,
        "elinewidth": ERROR_BAR_LINE_WIDTH,
        "capthick": ERROR_BAR_LINE_WIDTH,
        "label": label,
        "zorder": 3,
    }
    plot_kwargs.update(errorbar_kwargs)
    drawable_sd = np.where(summary.n > 1, summary.sd, 0.0)
    return ax.errorbar(x_values, summary.mean, yerr=drawable_sd, **plot_kwargs)


def style_manuscript_axis(
    ax: Axes,
    *,
    x_label: str | None = None,
    y_label: str | None = None,
    title: str | None = None,
    x_scale: str | None = None,
    y_scale: str | None = None,
) -> None:
    """Apply the shared axis typography, scales, and light guide lines."""
    if x_label is not None:
        ax.set_xlabel(x_label)
    if y_label is not None:
        ax.set_ylabel(y_label)
    if title is not None:
        ax.set_title(title)
    if x_scale is not None:
        ax.set_xscale(x_scale)
    if y_scale is not None:
        ax.set_yscale(y_scale)
    ax.tick_params(direction="out", width=0.8, length=3.0)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.65)
    ax.set_axisbelow(True)


def _axes_list(axes: Axes | Iterable[Axes]) -> list[Axes]:
    if isinstance(axes, Axes):
        return [axes]
    return [
        axis
        for axis in np.asarray(list(axes), dtype=object).ravel()
        if isinstance(axis, Axes)
    ]


def _unique_legend_entries(axes: list[Axes]) -> tuple[list[Any], list[str]]:
    handles: list[Any] = []
    labels: list[str] = []
    seen: set[str] = set()
    for axis in axes:
        axis_handles, axis_labels = axis.get_legend_handles_labels()
        for handle, label in zip(axis_handles, axis_labels, strict=True):
            if not label or label.startswith("_") or label in seen:
                continue
            seen.add(label)
            handles.append(handle)
            labels.append(label)
    rank = {label: index for index, label in enumerate(LEGEND_LABEL_ORDER)}
    ordered = sorted(
        zip(handles, labels, strict=True),
        key=lambda entry: rank.get(entry[1], len(rank)),
    )
    return [entry[0] for entry in ordered], [entry[1] for entry in ordered]


def place_manuscript_legend(
    fig: Figure,
    axes: Axes | Iterable[Axes],
    *,
    multi_panel: bool | None = None,
    location: str = "right",
    ncol: int | None = None,
    **legend_kwargs: Any,
) -> Legend | None:
    """Place a single legend outside the data region.

    A single-panel plot receives one axis legend outside its right boundary.
    A multi-panel plot receives one deduplicated figure-level legend shared by
    all panels.  ``location='bottom'`` is available when a right-side legend
    would make a wide panel impractical.
    """
    axis_list = _axes_list(axes)
    if not axis_list:
        raise ValueError("At least one Matplotlib axis is required.")
    if location not in {"right", "bottom"}:
        raise ValueError("location must be either 'right' or 'bottom'.")
    if multi_panel is None:
        multi_panel = len(axis_list) > 1

    handles, labels = _unique_legend_entries(axis_list)
    for axis in axis_list:
        existing = axis.get_legend()
        if existing is not None:
            existing.remove()
    if not handles:
        return None

    common_kwargs: dict[str, Any] = {
        "frameon": False,
        "fontsize": LEGEND_FONT_SIZE,
        "borderaxespad": 0.0,
    }
    common_kwargs.update(legend_kwargs)

    if not multi_panel:
        if location == "right":
            return axis_list[0].legend(
                handles,
                labels,
                loc="upper left",
                bbox_to_anchor=(1.02, 1.0),
                **common_kwargs,
            )
        return axis_list[0].legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.16),
            ncol=ncol or len(labels),
            **common_kwargs,
        )

    if location == "right":
        return fig.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            ncol=ncol or 1,
            **common_kwargs,
        )
    return fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=ncol or len(labels),
        **common_kwargs,
    )


def add_figure_note(
    fig: Figure,
    text: str,
    *,
    x: float = 0.01,
    y: float = 0.01,
    reserve_bottom: float = 0.16,
) -> Any:
    """Add one concise methods-style note below the plotting region."""
    if not text.strip():
        raise ValueError("Figure-note text cannot be empty.")
    fig.subplots_adjust(bottom=max(fig.subplotpars.bottom, reserve_bottom))
    return fig.text(
        x,
        y,
        text.strip(),
        ha="left",
        va="bottom",
        fontsize=FIGURE_NOTE_FONT_SIZE,
        color="#222222",
        wrap=True,
    )


def line_style(
    index: int, series_count: int, include_marker: bool = True
) -> dict[str, Any]:
    """Return a readable manuscript style for one plotted series."""
    if series_count <= MAX_MONOCHROME_SERIES:
        style = dict(MONOCHROME_LINE_STYLES[index % len(MONOCHROME_LINE_STYLES)])
    else:
        style = dict(LINE_FORMATS[index % len(LINE_FORMATS)])
        style["color"] = COLORBLIND_PALETTE[index % len(COLORBLIND_PALETTE)]

    style["linewidth"] = LINE_WIDTH
    if include_marker:
        style["markersize"] = MARKER_SIZE
        style["markevery"] = MARK_EVERY
    else:
        style.pop("marker", None)
    return style


def line_styles(
    series_count: int, include_markers: bool = True
) -> tuple[dict[str, Any], ...]:
    """Return styles for all series in a plot."""
    return tuple(
        line_style(index, series_count, include_marker=include_markers)
        for index in range(series_count)
    )


def bar_colors(series_count: int) -> list[str]:
    """Return bar colors following the same monochrome-then-colorblind rule."""
    return [line_style(index, series_count)["color"] for index in range(series_count)]


def bar_hatch(index: int) -> str:
    """Return a repeatable hatch pattern for bar charts."""
    return BAR_HATCHES[index % len(BAR_HATCHES)]


def save_manuscript_figure(
    fig: Figure,
    path: str | Path,
    *,
    abbreviation_keys: tuple[str, ...] | list[str] | None = None,
    dpi: int | None = None,
    **savefig_kwargs: Any,
) -> None:
    """Save a framework figure using the manuscript-resolution standard."""
    if abbreviation_keys:
        STANDARD_ABBREVIATIONS.add_figure_note(fig, abbreviation_keys)
    selected_dpi = MANUSCRIPT_DPI if dpi is None else int(dpi)
    if selected_dpi <= 0:
        raise ValueError("dpi must be positive")
    savefig_kwargs.setdefault("bbox_inches", "tight")
    savefig_kwargs.setdefault("facecolor", "white")
    fig.savefig(path, dpi=selected_dpi, **savefig_kwargs)
