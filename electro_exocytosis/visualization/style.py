from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib"))

from matplotlib.figure import Figure

MANUSCRIPT_DPI = 1200
LANDSCAPE_ASPECT_RATIO = 16 / 9
MANUSCRIPT_LANDSCAPE_FIGSIZE = (8.0, 4.5)
MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE = (12.0, 6.75)
MAX_MONOCHROME_SERIES = 3
LINE_WIDTH = 1.8
MARKER_SIZE = 4.2
MARK_EVERY = 0.12

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


def landscape_figsize(width_in: float = MANUSCRIPT_LANDSCAPE_FIGSIZE[0]) -> tuple[float, float]:
    """Return a standard 16:9 landscape figure size."""
    return (width_in, width_in / LANDSCAPE_ASPECT_RATIO)


def line_style(index: int, series_count: int, include_marker: bool = True) -> dict[str, Any]:
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


def line_styles(series_count: int, include_markers: bool = True) -> tuple[dict[str, Any], ...]:
    """Return styles for all series in a plot."""
    return tuple(line_style(index, series_count, include_marker=include_markers) for index in range(series_count))


def bar_colors(series_count: int) -> list[str]:
    """Return bar colors following the same monochrome-then-colorblind rule."""
    return [line_style(index, series_count)["color"] for index in range(series_count)]


def bar_hatch(index: int) -> str:
    """Return a repeatable hatch pattern for bar charts."""
    return BAR_HATCHES[index % len(BAR_HATCHES)]


def save_manuscript_figure(fig: Figure, path: str | Path, **savefig_kwargs: Any) -> None:
    """Save a framework figure using the manuscript-resolution standard."""
    fig.savefig(path, dpi=MANUSCRIPT_DPI, **savefig_kwargs)
