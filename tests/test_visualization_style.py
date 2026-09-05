from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from electro_exocytosis.visualization.style import (
    COLORBLIND_PALETTE,
    FITTED_MODEL_LABEL,
    LANDSCAPE_ASPECT_RATIO,
    MANUSCRIPT_DPI,
    MANUSCRIPT_LANDSCAPE_FIGSIZE,
    MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE,
    MAX_MONOCHROME_SERIES,
    MONOCHROME_LINE_STYLES,
    OBSERVED_MEAN_LABEL,
    OBSERVED_MEAN_SD_LABEL,
    add_figure_note,
    bar_colors,
    bar_hatch,
    line_style,
    line_styles,
    landscape_figsize,
    place_manuscript_legend,
    plot_observed_mean_sd,
    summarize_repeated_observations,
)


def test_manuscript_figure_defaults_are_landscape_1200_dpi() -> None:
    assert MANUSCRIPT_DPI == 1200
    assert LANDSCAPE_ASPECT_RATIO == pytest.approx(16 / 9)
    assert MANUSCRIPT_LANDSCAPE_FIGSIZE[0] / MANUSCRIPT_LANDSCAPE_FIGSIZE[
        1
    ] == pytest.approx(16 / 9)
    assert MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE[0] / MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE[
        1
    ] == pytest.approx(16 / 9)
    assert landscape_figsize(10.0) == pytest.approx((10.0, 5.625))
    assert OBSERVED_MEAN_LABEL == "Observed mean"


def test_line_styles_use_monochrome_for_up_to_three_series() -> None:
    styles = line_styles(MAX_MONOCHROME_SERIES)

    assert [style["color"] for style in styles] == [
        style["color"] for style in MONOCHROME_LINE_STYLES
    ]
    assert len({style["linestyle"] for style in styles}) == MAX_MONOCHROME_SERIES
    assert len({style["marker"] for style in styles}) == MAX_MONOCHROME_SERIES


def test_line_styles_switch_to_colorblind_palette_beyond_three_series() -> None:
    styles = line_styles(MAX_MONOCHROME_SERIES + 1)

    assert [style["color"] for style in styles] == list(
        COLORBLIND_PALETTE[: MAX_MONOCHROME_SERIES + 1]
    )
    assert all("marker" in style for style in styles)
    assert "marker" not in line_style(0, 2, include_marker=False)


def test_bar_styles_follow_same_series_rule() -> None:
    assert bar_colors(3) == [style["color"] for style in MONOCHROME_LINE_STYLES]
    assert bar_colors(4) == list(COLORBLIND_PALETTE[:4])
    assert bar_hatch(0) != bar_hatch(1)


def test_repeated_observations_use_sample_sd_and_track_missing_values() -> None:
    summary = summarize_repeated_observations(
        [[1.0, 4.0, np.nan], [3.0, 8.0, 5.0], [5.0, np.nan, np.nan]]
    )

    np.testing.assert_allclose(summary.mean[:2], [3.0, 6.0])
    np.testing.assert_allclose(summary.sd[:2], [2.0, np.sqrt(8.0)])
    assert np.isnan(summary.sd[2])
    np.testing.assert_array_equal(summary.n, [3, 2, 1])


def test_observed_plot_uses_mean_sd_label_and_retains_single_observation() -> None:
    fig, ax = plt.subplots()
    container = plot_observed_mean_sd(ax, [10.0, 20.0], [[1.0, 2.0], [3.0, np.nan]])

    assert container.get_label() == OBSERVED_MEAN_SD_LABEL
    np.testing.assert_allclose(container.lines[0].get_ydata(), [2.0, 2.0])
    plt.close(fig)


def test_single_panel_legend_is_outside_axis() -> None:
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label=OBSERVED_MEAN_SD_LABEL)
    ax.plot([0, 1], [1, 0], label=FITTED_MODEL_LABEL)

    legend = place_manuscript_legend(fig, ax)

    assert legend is ax.get_legend()
    assert not fig.legends
    assert [text.get_text() for text in legend.get_texts()] == [
        OBSERVED_MEAN_SD_LABEL,
        FITTED_MODEL_LABEL,
    ]
    assert legend.get_bbox_to_anchor()._bbox.x0 > 1.0
    plt.close(fig)


def test_multi_panel_legend_is_shared_and_deduplicated() -> None:
    fig, axes = plt.subplots(1, 2)
    for ax in axes:
        ax.plot([0, 1], [0, 1], label=OBSERVED_MEAN_SD_LABEL)
        ax.plot([0, 1], [1, 0], label=FITTED_MODEL_LABEL)
        ax.legend()

    legend = place_manuscript_legend(fig, axes)

    assert legend is fig.legends[0]
    assert len(fig.legends) == 1
    assert all(ax.get_legend() is None for ax in axes)
    assert [text.get_text() for text in legend.get_texts()] == [
        OBSERVED_MEAN_SD_LABEL,
        FITTED_MODEL_LABEL,
    ]
    plt.close(fig)


def test_semantic_legend_order_is_stable_for_error_bars() -> None:
    fig, ax = plt.subplots()
    ax.errorbar(
        [0, 1],
        [1, 2],
        yerr=[0.1, 0.2],
        label=OBSERVED_MEAN_SD_LABEL,
    )
    ax.plot([0, 1], [2, 1], label=FITTED_MODEL_LABEL)

    legend = place_manuscript_legend(fig, ax)

    assert [text.get_text() for text in legend.get_texts()] == [
        OBSERVED_MEAN_SD_LABEL,
        FITTED_MODEL_LABEL,
    ]
    plt.close(fig)


def test_figure_note_requires_text_and_reserves_space() -> None:
    fig, _ = plt.subplots()
    note = add_figure_note(fig, "Points and bars show mean ± SD.")

    assert note.get_text() == "Points and bars show mean ± SD."
    assert fig.subplotpars.bottom >= 0.16
    with pytest.raises(ValueError, match="cannot be empty"):
        add_figure_note(fig, "  ")
    plt.close(fig)
