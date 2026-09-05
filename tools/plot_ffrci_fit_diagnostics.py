#!/usr/bin/env python3
"""Plot fit-error and parameter-space diagnostics for the FFRi EV calibration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullFormatter


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_DIR = REPO_ROOT / "results" / "ffrci_ev_kinetics_fit"

CONDITION_ORDER = ("1p20kV", "3p40kV", "5p40kV")
CONDITION_COLORS = {
    "1p20kV": "#4477AA",
    "3p40kV": "#EE6677",
    "5p40kV": "#228833",
}

PARAMETER_METADATA = {
    "baseline_sEV_rate": {
        "module": "EV release kinetics",
        "submodule": "Small EV MVB release",
        "parameter_class": "Constitutive release rate",
        "label": "Baseline sEV release",
    },
    "k_ILV_release_s": {
        "module": "EV release kinetics",
        "submodule": "Small EV MVB release",
        "parameter_class": "Kinetic rate constant",
        "label": "ILV release rate",
    },
    "baseline_mlEV_rate": {
        "module": "EV release kinetics",
        "submodule": "Medium and large EV budding",
        "parameter_class": "Constitutive release rate",
        "label": "Baseline mlEV release",
    },
    "k_budding_s": {
        "module": "EV release kinetics",
        "submodule": "Medium and large EV budding",
        "parameter_class": "Kinetic rate constant",
        "label": "Budding rate",
    },
    "baseline_AB_rate": {
        "module": "EV release kinetics",
        "submodule": "Apoptotic body release",
        "parameter_class": "Constitutive release rate",
        "label": "Baseline AB release",
    },
    "k_apoptotic_commitment_s": {
        "module": "EV release kinetics",
        "submodule": "Apoptotic body release",
        "parameter_class": "Kinetic rate constant",
        "label": "Apoptotic commitment rate",
    },
}

SUBMODULE_PARAMETERS = {
    "Small EV MVB release": ("baseline_sEV_rate", "k_ILV_release_s"),
    "Medium and large EV budding": ("baseline_mlEV_rate", "k_budding_s"),
    "Apoptotic body release": ("baseline_AB_rate", "k_apoptotic_commitment_s"),
}


def _require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")


def load_diagnostics(
    results_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    comparison_path = results_dir / "observed_vs_predicted.csv"
    ensemble_path = results_dir / "sample_set_parameter_ensemble.csv"
    best_path = results_dir / "sample_set_best_parameters.csv"
    comparison = pd.read_csv(comparison_path)
    ensemble = pd.read_csv(ensemble_path)
    best = pd.read_csv(best_path)
    _require_columns(
        comparison,
        {
            "condition",
            "time_h",
            "observed_concentration",
            "predicted_fitted",
            "fitted_over_observed",
            "log10_residual",
        },
        comparison_path,
    )
    _require_columns(
        ensemble,
        {
            "condition",
            "replicate",
            "start",
            "parameter",
            "value",
            "multiplier",
            "near_optimal",
        },
        ensemble_path,
    )
    _require_columns(
        best,
        {"condition", "replicate", "parameter", "fitted", "multiplier"},
        best_path,
    )
    return comparison, ensemble, best


def calculate_error_metrics(
    comparison: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observations = comparison.copy()
    observations["fold_error"] = (
        observations["predicted_fitted"] / observations["observed_concentration"]
    )
    observations["absolute_percent_error"] = (
        100.0
        * np.abs(
            observations["predicted_fitted"] - observations["observed_concentration"]
        )
        / observations["observed_concentration"]
    )
    observations["log2_fold_error"] = np.log2(observations["fold_error"])

    rows = []
    scopes = [
        (condition, group) for condition, group in observations.groupby("condition")
    ]
    scopes.append(("Joint", observations))
    for scope, group in scopes:
        residual = group["log10_residual"].to_numpy(dtype=float)
        rows.append(
            {
                "scope": scope,
                "n_observations": len(group),
                "rmse_log10": float(np.sqrt(np.mean(residual**2))),
                "typical_fold_error": float(10.0 ** np.sqrt(np.mean(residual**2))),
                "mean_log10_bias": float(np.mean(residual)),
                "median_absolute_percent_error": float(
                    np.median(group["absolute_percent_error"])
                ),
                "maximum_absolute_fold_error": float(
                    np.max(np.maximum(group["fold_error"], 1.0 / group["fold_error"]))
                ),
            }
        )
    return observations, pd.DataFrame(rows)


def plot_error_diagnostics(
    observations: pd.DataFrame, metrics: pd.DataFrame, output_path: Path
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2))
    observation_unit = (
        str(observations["observation_unit"].iloc[0])
        if "observation_unit" in observations
        else "particles_per_ml"
    )
    value_label = (
        "particle equivalents per cell"
        if observation_unit == "particle_equivalents_per_cell"
        else "concentration (particles/mL)"
    )

    scatter_axis = axes[0, 0]
    all_values = np.concatenate(
        [
            observations["observed_concentration"].to_numpy(dtype=float),
            observations["predicted_fitted"].to_numpy(dtype=float),
        ]
    )
    lower, upper = float(all_values.min() / 1.25), float(all_values.max() * 1.25)
    for condition in CONDITION_ORDER:
        group = observations[observations["condition"] == condition]
        scatter_axis.scatter(
            group["observed_concentration"],
            group["predicted_fitted"],
            s=62,
            color=CONDITION_COLORS[condition],
            label=condition,
            zorder=3,
        )
        for row in group.itertuples():
            annotation_offsets = {
                "1p20kV": (6, 7),
                "3p40kV": (6, 7),
                "5p40kV": (6, -13),
            }
            scatter_axis.annotate(
                f"{row.time_h:g} h",
                (row.observed_concentration, row.predicted_fitted),
                xytext=annotation_offsets[condition],
                textcoords="offset points",
                fontsize=8,
            )
    scatter_axis.plot([lower, upper], [lower, upper], "--", color="#555555", lw=1.4)
    scatter_axis.set(
        xscale="log", yscale="log", xlim=(lower, upper), ylim=(lower, upper)
    )
    scatter_axis.set_xlabel(f"Observed {value_label}")
    scatter_axis.set_ylabel(f"Fitted {value_label}")
    pearson = float(
        observations["observed_concentration"].corr(observations["predicted_fitted"])
    )
    scatter_axis.set_title(f"A  Absolute agreement  Pearson r = {pearson:.2f}")
    scatter_axis.legend(frameon=False, fontsize=8)
    scatter_axis.grid(alpha=0.22)

    residual_axis = axes[0, 1]
    for condition in CONDITION_ORDER:
        group = observations[observations["condition"] == condition].sort_values(
            "time_h"
        )
        residual_axis.plot(
            group["time_h"],
            group["log10_residual"],
            "o-",
            lw=2,
            color=CONDITION_COLORS[condition],
            label=condition,
        )
    residual_axis.axhline(0.0, color="#555555", ls="--", lw=1.2)
    residual_axis.set_xlabel("Time after nsPEF (h)")
    residual_axis.set_ylabel("log10 fitted / observed")
    residual_axis.set_title("B  Signed residual over time")
    residual_axis.grid(alpha=0.22)

    heatmap_axis = axes[1, 0]
    heatmap = (
        observations.pivot(
            index="condition", columns="time_h", values="log2_fold_error"
        )
        .reindex(index=CONDITION_ORDER, columns=(0.5, 1.0, 3.0))
        .to_numpy(dtype=float)
    )
    limit = max(1.0, float(np.nanmax(np.abs(heatmap))))
    image = heatmap_axis.imshow(
        heatmap, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto"
    )
    heatmap_axis.set_xticks(range(3), ["0.5 h", "1 h", "3 h"])
    heatmap_axis.set_yticks(range(3), CONDITION_ORDER)
    heatmap_axis.set_title("C  Fold error by condition and time")
    for row_index, condition in enumerate(CONDITION_ORDER):
        for column_index, time_h in enumerate((0.5, 1.0, 3.0)):
            value = observations.loc[
                (observations["condition"] == condition)
                & (observations["time_h"] == time_h),
                "fold_error",
            ].iloc[0]
            heatmap_axis.text(
                column_index,
                row_index,
                f"{value:.2f}×",
                ha="center",
                va="center",
                color=(
                    "white"
                    if abs(heatmap[row_index, column_index]) > 0.65 * limit
                    else "black"
                ),
                fontsize=10,
                fontweight="bold",
            )
    colorbar = fig.colorbar(image, ax=heatmap_axis, fraction=0.046, pad=0.04)
    colorbar.set_label("log2 fitted / observed")

    metric_axis = axes[1, 1]
    order = [*CONDITION_ORDER, "Joint"]
    ordered_metrics = metrics.set_index("scope").loc[order].reset_index()
    bar_colors = [CONDITION_COLORS[name] for name in CONDITION_ORDER] + ["#666666"]
    bars = metric_axis.bar(
        ordered_metrics["scope"],
        ordered_metrics["typical_fold_error"],
        color=bar_colors,
        alpha=0.9,
    )
    metric_axis.axhline(1.0, color="#555555", ls="--", lw=1.2)
    metric_axis.set_ylabel("Typical fold error (10^RMSE)")
    metric_axis.set_title("D  Error magnitude")
    metric_axis.grid(axis="y", alpha=0.22)
    for bar, value in zip(bars, ordered_metrics["typical_fold_error"], strict=True):
        metric_axis.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.04,
            f"{value:.2f}×",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    metric_axis.set_ylim(0.0, float(ordered_metrics["typical_fold_error"].max() * 1.20))

    figure_context = (
        "single-cell-equivalent output"
        if observation_unit == "particle_equivalents_per_cell"
        else "concentration"
    )
    fig.suptitle(f"Current EV model {figure_context} fit diagnostics", fontsize=16)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _add_parameter_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    enriched = frame.copy()
    for column in ("module", "submodule", "parameter_class", "label"):
        enriched[column] = enriched["parameter"].map(
            {name: metadata[column] for name, metadata in PARAMETER_METADATA.items()}
        )
    if enriched[["module", "submodule", "parameter_class", "label"]].isna().any().any():
        unknown = sorted(set(enriched.loc[enriched["label"].isna(), "parameter"]))
        raise ValueError(f"Missing parameter metadata for {unknown}")
    return enriched


def _distinct_near_optimal_solutions(ensemble: pd.DataFrame) -> pd.DataFrame:
    near = ensemble[ensemble["near_optimal"].astype(bool)].copy()
    near = _add_parameter_metadata(near)
    # p1/p2/p3 concentrations are numerically indistinguishable within each
    # condition. Collapse exact repeated solutions so the boxes show distinct
    # optimizer-compatible trajectories rather than triplicate weighting.
    return near.drop_duplicates(["condition", "start", "parameter", "value"])


def _condition_legend() -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=CONDITION_COLORS[condition],
            markeredgecolor="white",
            markersize=8,
            label=condition,
        )
        for condition in CONDITION_ORDER
    ]


def _box_and_points(
    axis: plt.Axes,
    frame: pd.DataFrame,
    parameter_order: tuple[str, ...],
    *,
    value_column: str,
    ylabel: str,
) -> None:
    grouped_values = [
        frame.loc[frame["parameter"] == parameter, value_column].to_numpy(dtype=float)
        for parameter in parameter_order
    ]
    boxes = axis.boxplot(
        grouped_values,
        positions=np.arange(len(parameter_order)),
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.5},
        whiskerprops={"color": "#555555"},
        capprops={"color": "#555555"},
    )
    for patch in boxes["boxes"]:
        patch.set_facecolor("#D9E6F2")
        patch.set_edgecolor("#555555")
    rng = np.random.default_rng(20260904)
    for parameter_index, parameter in enumerate(parameter_order):
        points = frame[frame["parameter"] == parameter]
        jitter = rng.uniform(-0.13, 0.13, len(points))
        axis.scatter(
            parameter_index + jitter,
            points[value_column],
            c=points["condition"].map(CONDITION_COLORS),
            s=36,
            edgecolor="white",
            linewidth=0.45,
            alpha=0.86,
            zorder=3,
        )
    axis.set_xticks(
        np.arange(len(parameter_order)),
        [PARAMETER_METADATA[parameter]["label"] for parameter in parameter_order],
        rotation=18,
        ha="right",
    )
    axis.set_yscale("log")
    axis.set_ylabel(ylabel)
    axis.grid(axis="y", alpha=0.22)


def plot_parameter_boxplots(ensemble: pd.DataFrame, output_path: Path) -> None:
    near = _distinct_near_optimal_solutions(ensemble)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharey=True)
    for axis, (submodule, parameter_order) in zip(
        axes, SUBMODULE_PARAMETERS.items(), strict=True
    ):
        subset = near[near["submodule"] == submodule]
        _box_and_points(
            axis,
            subset,
            parameter_order,
            value_column="multiplier",
            ylabel="Fitted / default parameter value",
        )
        axis.axhline(1.0, color="#555555", ls="--", lw=1.2)
        axis.set_title(submodule)
    axes[0].legend(handles=_condition_legend(), frameon=False, fontsize=8, loc="best")
    fig.suptitle("Near-optimal parameter space by EV-release submodule")
    fig.text(
        0.5,
        -0.01,
        "Boxes and points use distinct multistart solutions; vertical scale is logarithmic.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_constitutive_and_kinetic_parameters(
    ensemble: pd.DataFrame, output_path: Path
) -> None:
    near = _distinct_near_optimal_solutions(ensemble)
    constitutive = tuple(
        name
        for name, metadata in PARAMETER_METADATA.items()
        if metadata["parameter_class"] == "Constitutive release rate"
    )
    kinetic = tuple(
        name
        for name, metadata in PARAMETER_METADATA.items()
        if metadata["parameter_class"] == "Kinetic rate constant"
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    _box_and_points(
        axes[0],
        near,
        constitutive,
        value_column="value",
        ylabel="Fitted value (current model units)",
    )
    axes[0].set_title("Constitutive release parameters")
    _box_and_points(
        axes[1],
        near,
        kinetic,
        value_column="value",
        ylabel="Fitted kinetic rate (s⁻¹)",
    )
    axes[1].set_title("Kinetic parameters")
    axes[0].legend(handles=_condition_legend(), frameon=False, fontsize=8, loc="best")
    fig.suptitle("Near-optimal constitutive and kinetic parameter values")
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_parameter_tradeoffs(ensemble: pd.DataFrame, output_path: Path) -> None:
    near = _distinct_near_optimal_solutions(ensemble)
    index_columns = ["condition", "replicate", "start"]
    wide = near.pivot_table(
        index=index_columns, columns="parameter", values="value"
    ).reset_index()
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5))
    for axis, (submodule, parameters) in zip(
        axes, SUBMODULE_PARAMETERS.items(), strict=True
    ):
        baseline_parameter, kinetic_parameter = parameters
        subset = wide.dropna(
            subset=[baseline_parameter, kinetic_parameter]
        ).drop_duplicates(["condition", baseline_parameter, kinetic_parameter])
        for condition in CONDITION_ORDER:
            points = subset[subset["condition"] == condition]
            axis.scatter(
                points[baseline_parameter],
                points[kinetic_parameter],
                s=52,
                color=CONDITION_COLORS[condition],
                edgecolor="white",
                linewidth=0.55,
                alpha=0.9,
                label=condition,
            )
        kinetic_values = subset[kinetic_parameter].to_numpy(dtype=float)
        if kinetic_values.max() / kinetic_values.min() > 10.0:
            axis.set_yscale("log")
            axis.yaxis.set_minor_formatter(NullFormatter())
        else:
            center = float(np.median(kinetic_values))
            lower = min(float(kinetic_values.min()), 0.95 * center)
            upper = max(float(kinetic_values.max()), 1.05 * center)
            axis.set_ylim(lower, upper)
            axis.yaxis.set_major_formatter(
                FuncFormatter(lambda value, _: f"{value:.2e}")
            )
        axis.set_xlabel(PARAMETER_METADATA[baseline_parameter]["label"])
        axis.set_ylabel(f"{PARAMETER_METADATA[kinetic_parameter]['label']} (s⁻¹)")
        axis.set_title(submodule)
        axis.grid(alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Constitutive and kinetic parameter tradeoffs")
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_diagnostic_summary(
    results_dir: Path,
    observations: pd.DataFrame,
    metrics: pd.DataFrame,
    ensemble: pd.DataFrame,
) -> None:
    observations.to_csv(
        results_dir / "fit_error_metrics_by_observation.csv", index=False
    )
    metrics.to_csv(results_dir / "fit_error_metrics_by_condition.csv", index=False)
    condition_separation = []
    for time_h, group in observations.groupby("time_h"):
        condition_separation.append(
            {
                "time_h": float(time_h),
                "observed_max_to_min": float(
                    group["observed_concentration"].max()
                    / group["observed_concentration"].min()
                ),
                "predicted_max_to_min": float(
                    group["predicted_fitted"].max() / group["predicted_fitted"].min()
                ),
            }
        )
    near = _distinct_near_optimal_solutions(ensemble)
    summary = {
        "global_range_comparison": {
            "observed_max_to_min": float(
                observations["observed_concentration"].max()
                / observations["observed_concentration"].min()
            ),
            "predicted_max_to_min": float(
                observations["predicted_fitted"].max()
                / observations["predicted_fitted"].min()
            ),
            "pearson_correlation": float(
                observations["observed_concentration"].corr(
                    observations["predicted_fitted"]
                )
            ),
            "spearman_correlation": float(
                observations["observed_concentration"].corr(
                    observations["predicted_fitted"], method="spearman"
                )
            ),
        },
        "condition_separation_by_time": condition_separation,
        "largest_absolute_fold_error": float(
            np.max(
                np.maximum(observations["fold_error"], 1.0 / observations["fold_error"])
            )
        ),
        "largest_error_observation": observations.iloc[
            np.argmax(
                np.maximum(observations["fold_error"], 1.0 / observations["fold_error"])
            )
        ][["condition", "time_h", "fold_error"]].to_dict(),
        "near_optimal_parameter_points_after_duplicate_collapse": int(len(near)),
        "interpretation": [
            "The optimizer aligns the overall concentration scale but cannot match decreases because the current model output is cumulative.",
            "The fitted trajectories separate pulse conditions much less than the observations do.",
            "Parameters at bounds or spanning most of their permitted range are weakly identified by total EV concentration alone.",
        ],
    }
    (results_dir / "fit_diagnostic_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    sample_summary_path = results_dir / "sample_set_fit_summary.json"
    replicate_note = (
        "The parameter boxes summarize optimizer-compatible solutions. Their "
        "interpretation as experimental population distributions depends on the "
        "unconfirmed p1/p2/p3 replicate structure."
    )
    if sample_summary_path.exists():
        sample_summary = json.loads(sample_summary_path.read_text(encoding="utf-8"))
        indistinguishable = int(
            sample_summary.get(
                "condition_time_cells_with_numerically_indistinguishable_replicate_totals",
                0,
            )
        )
        total_cells = int(sample_summary.get("condition_time_cells_total", 0))
        if total_cells and indistinguishable == total_cells:
            replicate_note = (
                "The parameter boxes summarize optimizer-compatible solutions. They are not "
                "experimental population distributions because the available p1/p2/p3 total "
                "concentrations are numerically indistinguishable within each condition/time."
            )
        elif total_cells:
            maximum_spread = 100.0 * float(
                sample_summary.get(
                    "maximum_relative_spread_across_replicate_totals", float("nan")
                )
            )
            replicate_note = (
                "The parameter boxes summarize optimizer-compatible solutions. The filtered "
                f"p1/p2/p3 concentrations differ by up to {maximum_spread:.1f}% within a "
                "condition/time cell because their size distributions differ, but the parameter "
                "ranges should not be treated as biological population distributions until the "
                "replicate structure is confirmed."
            )

    readme = f"""# EV calibration diagnostic figures

These figures use the completed calibration tables and do not rerun or modify
the EV model.

- `fit_error_diagnostics.png` compares absolute agreement, signed residuals,
  condition/time fold errors, and typical error magnitude.
- `parameter_space_boxplots.png` shows fitted/default ratios for each
  constitutive and kinetic parameter within its EV-release submodule.
- `constitutive_and_kinetic_parameter_boxplots.png` shows fitted values in the
  current model units.
- `parameter_tradeoffs_by_submodule.png` shows the constitutive/kinetic
  combinations that produced near-optimal fits.

{replicate_note}
"""
    (results_dir / "fit_diagnostics_README.md").write_text(readme, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    comparison, ensemble, _ = load_diagnostics(args.results_dir)
    observations, metrics = calculate_error_metrics(comparison)
    plot_error_diagnostics(
        observations, metrics, args.results_dir / "fit_error_diagnostics.png"
    )
    plot_parameter_boxplots(ensemble, args.results_dir / "parameter_space_boxplots.png")
    plot_constitutive_and_kinetic_parameters(
        ensemble,
        args.results_dir / "constitutive_and_kinetic_parameter_boxplots.png",
    )
    plot_parameter_tradeoffs(
        ensemble, args.results_dir / "parameter_tradeoffs_by_submodule.png"
    )
    write_diagnostic_summary(args.results_dir, observations, metrics, ensemble)
    print(f"Diagnostic figures written to {args.results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
