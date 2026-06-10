from __future__ import annotations

import pytest

from electro_exocytosis.visualization.style import (
    COLORBLIND_PALETTE,
    LANDSCAPE_ASPECT_RATIO,
    MANUSCRIPT_DPI,
    MANUSCRIPT_LANDSCAPE_FIGSIZE,
    MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE,
    MAX_MONOCHROME_SERIES,
    MONOCHROME_LINE_STYLES,
    bar_colors,
    bar_hatch,
    line_style,
    line_styles,
    landscape_figsize,
)


def test_manuscript_figure_defaults_are_landscape_1200_dpi() -> None:
    assert MANUSCRIPT_DPI == 1200
    assert LANDSCAPE_ASPECT_RATIO == pytest.approx(16 / 9)
    assert MANUSCRIPT_LANDSCAPE_FIGSIZE[0] / MANUSCRIPT_LANDSCAPE_FIGSIZE[1] == pytest.approx(16 / 9)
    assert MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE[0] / MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE[1] == pytest.approx(16 / 9)
    assert landscape_figsize(10.0) == pytest.approx((10.0, 5.625))


def test_line_styles_use_monochrome_for_up_to_three_series() -> None:
    styles = line_styles(MAX_MONOCHROME_SERIES)

    assert [style["color"] for style in styles] == [style["color"] for style in MONOCHROME_LINE_STYLES]
    assert len({style["linestyle"] for style in styles}) == MAX_MONOCHROME_SERIES
    assert len({style["marker"] for style in styles}) == MAX_MONOCHROME_SERIES


def test_line_styles_switch_to_colorblind_palette_beyond_three_series() -> None:
    styles = line_styles(MAX_MONOCHROME_SERIES + 1)

    assert [style["color"] for style in styles] == list(COLORBLIND_PALETTE[: MAX_MONOCHROME_SERIES + 1])
    assert all("marker" in style for style in styles)
    assert "marker" not in line_style(0, 2, include_marker=False)


def test_bar_styles_follow_same_series_rule() -> None:
    assert bar_colors(3) == [style["color"] for style in MONOCHROME_LINE_STYLES]
    assert bar_colors(4) == list(COLORBLIND_PALETTE[:4])
    assert bar_hatch(0) != bar_hatch(1)
