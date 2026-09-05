from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = REPO_ROOT / "data" / "experimental" / "ffrci_data_sharing"
DEFAULT_OUT_DIR = REPO_ROOT / "results" / "ffrci_data_audit"

EXOID_FILE = "Exoid_particle number and size in CD4 sham and different fields.csv"
ALL_RNASEQ_FILE = "RNAseq all samples.csv"
DGE_FILE = "RNAseq_CD4 sham vs CD4 3p40kV_DGE (236 genes).csv"
PATHWAY_FILE = "RNAseq_CD4 sham vs CD4 3p40kV_Enriched Pathways.csv"

TREATMENT_PATTERN = re.compile(
    r"^(?P<pulse_count>\d+)p(?P<amplitude_label>\d+)kv\s+"
    r"(?P<harvest_value>\d+)(?P<harvest_unit>min|h)\s+CD4\+\s+p(?P<replicate>\d+)$",
    re.IGNORECASE,
)
CONTROL_PATTERN = re.compile(
    r"^(?P<control_label>Sham2|sham media)\s+CD4\+\s+p(?P<replicate>\d+)$",
    re.IGNORECASE,
)

MECHANISTIC_GENE_SETS = {
    "calcium_and_ion_handling": {
        "ATP2A2",
        "ATP2B1",
        "ATP2C1",
        "CAPN1",
        "CAPN2",
        "KCNN4",
        "MCU",
        "ORAI1",
        "PPP3R1",
        "PRKCA",
        "SLC8A1",
        "SLC8A2",
        "STIM1",
    },
    "membrane_repair_and_cytoskeleton": {
        "ACTN4",
        "ANXA1",
        "ANXA2",
        "ANXA5",
        "CAPN2",
        "CD9",
        "MYH9",
        "RAC1",
        "RHOA",
        "TMEM16F",
        "VIM",
    },
    "ev_biogenesis_and_trafficking": {
        "ALIX",
        "CD9",
        "PDCD6IP",
        "RAB11A",
        "RAB27A",
        "RAB27B",
        "RAB35",
        "SDCBP",
        "SMPD3",
        "TSG101",
    },
    "oxidative_and_metabolic_stress": {
        "GCLC",
        "GPX4",
        "HIF1A",
        "HMOX1",
        "LDHA",
        "MT2A",
        "NQO1",
        "SOD2",
        "SQSTM1",
    },
    "apoptosis_and_injury": {
        "BAX",
        "BCL2",
        "CASP3",
        "CASP8",
        "CDKN1A",
        "FAS",
        "FASLG",
        "GADD45A",
        "GADD45B",
        "PMAIP1",
    },
}


def parse_exoid_file(path: Path) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    dataset_number: int | None = None
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if not row or not any(cell.strip() for cell in row):
                continue
            first = row[0].strip()
            if first.lower().startswith("data set"):
                dataset_number = int(first.split()[-1])
                continue
            if first == "Particle Diameter (nm)":
                continue
            if dataset_number is None or len(row) < 3:
                raise ValueError(f"Unexpected Exoid row before dataset header: {row!r}")
            records.append(
                {
                    "dataset_number": dataset_number,
                    "particle_diameter_nm": float(row[0]),
                    "concentration_particles_per_ml": float(row[1]),
                    "label": row[2].strip(),
                }
            )
    frame = pd.DataFrame.from_records(records)
    if frame.empty:
        raise ValueError(f"No particle records found in {path}")
    label_counts = frame.groupby("dataset_number")["label"].nunique()
    if (label_counts != 1).any():
        raise ValueError("Each Exoid dataset must have exactly one sample label")
    return frame


def parse_sample_label(label: str) -> dict[str, object]:
    treatment = TREATMENT_PATTERN.match(label)
    if treatment:
        values = treatment.groupdict()
        harvest_value = int(values["harvest_value"])
        harvest_hours = harvest_value / 60.0 if values["harvest_unit"].lower() == "min" else float(harvest_value)
        return {
            "sample_type": "treatment",
            "condition": f"{int(values['pulse_count'])}p{int(values['amplitude_label'])}kV",
            "control_label": "",
            "pulse_count": int(values["pulse_count"]),
            "amplitude_label_kV": int(values["amplitude_label"]),
            "harvest_time_h": harvest_hours,
            "replicate": int(values["replicate"]),
        }
    control = CONTROL_PATTERN.match(label)
    if control:
        values = control.groupdict()
        control_label = values["control_label"].lower().replace(" ", "_")
        return {
            "sample_type": "control",
            "condition": control_label,
            "control_label": control_label,
            "pulse_count": 0,
            "amplitude_label_kV": 0,
            "harvest_time_h": np.nan,
            "replicate": int(values["replicate"]),
        }
    raise ValueError(f"Unrecognized Exoid sample label: {label}")


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    if cumulative[-1] <= 0:
        return float("nan")
    target = quantile * cumulative[-1]
    return float(sorted_values[min(np.searchsorted(cumulative, target, side="left"), len(sorted_values) - 1)])


def summarize_exoid(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_number, group in frame.groupby("dataset_number", sort=True):
        diameters = group["particle_diameter_nm"].to_numpy(dtype=float)
        concentrations = group["concentration_particles_per_ml"].to_numpy(dtype=float)
        total = float(concentrations.sum())
        row = {
            "dataset_number": int(dataset_number),
            "label": str(group["label"].iloc[0]),
            **parse_sample_label(str(group["label"].iloc[0])),
            "n_reported_size_bins": int(len(group)),
            "diameter_min_nm": float(diameters.min()),
            "diameter_max_nm": float(diameters.max()),
            "summed_bin_concentration_particles_per_ml": total,
            "modal_diameter_nm": float(diameters[np.argmax(concentrations)]),
            "concentration_weighted_mean_diameter_nm": float(np.average(diameters, weights=concentrations)),
            "weighted_d10_nm": weighted_quantile(diameters, concentrations, 0.10),
            "weighted_d50_nm": weighted_quantile(diameters, concentrations, 0.50),
            "weighted_d90_nm": weighted_quantile(diameters, concentrations, 0.90),
            "fraction_below_100_nm": float(concentrations[diameters < 100].sum() / total),
            "fraction_100_to_200_nm": float(concentrations[(diameters >= 100) & (diameters < 200)].sum() / total),
            "fraction_at_or_above_200_nm": float(concentrations[diameters >= 200].sum() / total),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_conditions(sample_summary: pd.DataFrame) -> pd.DataFrame:
    group_columns = [
        "sample_type",
        "condition",
        "control_label",
        "pulse_count",
        "amplitude_label_kV",
        "harvest_time_h",
    ]
    metric_columns = [
        "summed_bin_concentration_particles_per_ml",
        "concentration_weighted_mean_diameter_nm",
        "weighted_d50_nm",
        "fraction_below_100_nm",
        "fraction_100_to_200_nm",
        "fraction_at_or_above_200_nm",
    ]
    rows: list[dict[str, object]] = []
    normalized = sample_summary.copy()
    normalized["harvest_time_group"] = normalized["harvest_time_h"].fillna(-1.0)
    actual_group_columns = [column if column != "harvest_time_h" else "harvest_time_group" for column in group_columns]
    for keys, group in normalized.groupby(actual_group_columns, dropna=False, sort=True):
        key_map = dict(zip(group_columns, keys, strict=True))
        key_map["harvest_time_h"] = np.nan if key_map["harvest_time_h"] == -1 else key_map["harvest_time_h"]
        row: dict[str, object] = {**key_map, "n_samples": int(len(group))}
        for metric in metric_columns:
            values = group[metric].to_numpy(dtype=float)
            mean = float(np.mean(values))
            sample_sd = float(np.std(values, ddof=1)) if len(values) > 1 else float("nan")
            row[f"{metric}_mean"] = mean
            row[f"{metric}_sample_sd"] = sample_sd
            row[f"{metric}_cv"] = sample_sd / mean if mean != 0 and np.isfinite(sample_sd) else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def audit_rnaseq(data_dir: Path, out_dir: Path) -> dict[str, object]:
    all_samples = pd.read_csv(data_dir / ALL_RNASEQ_FILE, low_memory=False)
    dge_raw = pd.read_csv(data_dir / DGE_FILE)
    dge = dge_raw[dge_raw["EnsemblID"].astype(str).str.startswith("ENSG")].copy()
    pathways = pd.read_csv(data_dir / PATHWAY_FILE)

    cpm_columns = [column for column in all_samples.columns if column.lower().endswith("_cpm")]
    count_columns = [column for column in all_samples.columns if column.lower().endswith("_count")]
    sample_names = [column[: -len("_cpm")] for column in cpm_columns]
    paired_counts = {column[: -len("_count")] for column in count_columns}
    unpaired_cpm = sorted(set(sample_names) - paired_counts)
    unpaired_count = sorted(paired_counts - set(sample_names))

    library_rows = []
    for sample_name in sample_names:
        count = pd.to_numeric(all_samples[f"{sample_name}_count"], errors="coerce")
        cpm = pd.to_numeric(all_samples[f"{sample_name}_cpm"], errors="coerce")
        library_rows.append(
            {
                "sample_name": sample_name,
                "sum_reported_counts": float(count.sum()),
                "sum_cpm": float(cpm.sum()),
                "genes_with_positive_count": int((count > 0).sum()),
                "missing_count_values": int(count.isna().sum()),
                "missing_cpm_values": int(cpm.isna().sum()),
            }
        )
    pd.DataFrame(library_rows).to_csv(out_dir / "rnaseq_sample_library_summary.csv", index=False)

    dge["gene_name_upper"] = dge["gene_name"].astype(str).str.upper()
    mapped_rows = []
    for module, genes in MECHANISTIC_GENE_SETS.items():
        matches = dge[dge["gene_name_upper"].isin(genes)].copy()
        for _, record in matches.iterrows():
            mapped_rows.append(
                {
                    "model_module": module,
                    "gene_name": record["gene_name"],
                    "EnsemblID": record["EnsemblID"],
                    "logFC": record["logFC"],
                    "FDR": record["FDR"],
                    "GroupA_NormCPM": record["GroupA_NormCPM"],
                    "GroupB_NormCPM": record["GroupB_NormCPM"],
                }
            )
    mapped = pd.DataFrame(mapped_rows)
    if mapped.empty:
        mapped = pd.DataFrame(
            columns=["model_module", "gene_name", "EnsemblID", "logFC", "FDR", "GroupA_NormCPM", "GroupB_NormCPM"]
        )
    else:
        mapped = mapped.sort_values(["model_module", "FDR", "gene_name"])
    mapped.to_csv(out_dir / "rnaseq_model_module_gene_mapping.csv", index=False)

    pathway_summary = pathways[
        ["TERM", "NES", "NOM_P_VAL", "FDR_Q_VAL", "FWER_P_VAL", "TAG_%", "GENE_%"]
    ].sort_values("FDR_Q_VAL")
    pathway_summary.to_csv(out_dir / "rnaseq_enriched_pathways_summary.csv", index=False)

    return {
        "all_samples_gene_rows": int(len(all_samples)),
        "all_samples_columns": int(len(all_samples.columns)),
        "rnaseq_sample_count": int(len(sample_names)),
        "cpm_count_pairs_complete": not unpaired_cpm and not unpaired_count,
        "unpaired_cpm_samples": unpaired_cpm,
        "unpaired_count_samples": unpaired_count,
        "dge_rows": int(len(dge)),
        "dge_non_data_footer_rows": int(len(dge_raw) - len(dge)),
        "dge_fdr_below_0_05": int((pd.to_numeric(dge["FDR"], errors="coerce") < 0.05).sum()),
        "dge_positive_logfc": int((pd.to_numeric(dge["logFC"], errors="coerce") > 0).sum()),
        "dge_negative_logfc": int((pd.to_numeric(dge["logFC"], errors="coerce") < 0).sum()),
        "enriched_pathway_rows": int(len(pathways)),
        "pathways_fdr_below_0_05": int((pd.to_numeric(pathways["FDR_Q_VAL"], errors="coerce") < 0.05).sum()),
        "mapped_model_gene_rows": int(len(mapped)),
    }


def make_plots(raw_exoid: pd.DataFrame, sample_summary: pd.DataFrame, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))

    treated = sample_summary[sample_summary["sample_type"] == "treatment"].copy()
    for condition, group in treated.groupby("condition", sort=True):
        by_time = group.groupby("harvest_time_h")["summed_bin_concentration_particles_per_ml"]
        times = np.array(sorted(group["harvest_time_h"].unique()), dtype=float)
        means = np.array([by_time.get_group(time).mean() for time in times], dtype=float)
        stds = np.array(
            [by_time.get_group(time).std(ddof=1) if len(by_time.get_group(time)) > 1 else np.nan for time in times],
            dtype=float,
        )
        axes[0].errorbar(times, means, yerr=stds, marker="o", capsize=3, label=condition)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Harvest time after exposure (h)")
    axes[0].set_ylabel("Summed reported bin concentration (particles/mL)")
    axes[0].set_title("Particle concentration by labeled condition")
    axes[0].legend(frameon=False)

    representative_labels = []
    for condition in ["sham2", "sham_media", "1p20kV", "3p40kV", "5p40kV"]:
        matches = sample_summary[sample_summary["condition"].str.lower() == condition.lower()]
        if matches.empty:
            continue
        if condition in {"sham2", "sham_media"}:
            selected = matches.sort_values("replicate").iloc[0]
        else:
            one_hour = matches[np.isclose(matches["harvest_time_h"], 1.0)]
            selected = (one_hour if not one_hour.empty else matches).sort_values("replicate").iloc[0]
        representative_labels.append((condition, selected["label"]))
    for condition, label in representative_labels:
        group = raw_exoid[raw_exoid["label"] == label]
        axes[1].plot(
            group["particle_diameter_nm"],
            group["concentration_particles_per_ml"],
            label=f"{condition}: {label}",
        )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Particle diameter (nm)")
    axes[1].set_ylabel("Reported concentration (particles/mL)")
    axes[1].set_title("Representative particle size distributions")
    axes[1].legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_audit_summary(
    out_path: Path,
    sample_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    rnaseq_audit: dict[str, object],
) -> None:
    treatment = sample_summary[sample_summary["sample_type"] == "treatment"]
    controls = sample_summary[sample_summary["sample_type"] == "control"]
    missing_replicates = condition_summary[
        (condition_summary["sample_type"] == "treatment") & (condition_summary["n_samples"] < 3)
    ]
    concentration_groups = sample_summary.groupby(["condition", "harvest_time_h"], dropna=False)[
        "summed_bin_concentration_particles_per_ml"
    ]
    all_condition_totals_repeated = all(
        np.allclose(values.to_numpy(dtype=float), float(values.iloc[0]), rtol=1e-9, atol=1.0)
        for _, values in concentration_groups
    )
    lines = [
        "# FFRCI experimental data audit",
        "",
        "## Machine readable inventory",
        "",
        f"- Exoid particle distributions: {len(sample_summary)} sample-level datasets, "
        f"including {len(treatment)} treatment datasets and {len(controls)} control datasets.",
        f"- Labeled treatment combinations: {treatment.groupby(['condition', 'harvest_time_h']).ngroups}.",
        f"- RNA-seq all-sample matrix: {rnaseq_audit['all_samples_gene_rows']} genes across "
        f"{rnaseq_audit['rnaseq_sample_count']} labeled samples.",
        f"- Differential-expression table: {rnaseq_audit['dge_rows']} rows; "
        f"{rnaseq_audit['dge_fdr_below_0_05']} have FDR below 0.05.",
        f"- Enrichment table: {rnaseq_audit['enriched_pathway_rows']} pathways; "
        f"{rnaseq_audit['pathways_fdr_below_0_05']} have FDR q-value below 0.05.",
        "",
        "## Structural findings",
        "",
        "- The Exoid CSV is a concatenation of separately headed datasets, not one ordinary rectangular table. "
        "The parser preserves dataset numbers and sample labels.",
        "- The labels encode pulse-count/amplitude combinations and 0.5, 1, and 3 hour harvest times, "
        "but they do not encode pulse width, repetition rate, exposure geometry, medium, or whether the amplitude is generator voltage or electric-field strength.",
        "- Sham2 and sham media are distinct control labels without harvest times. They must not be pooled until their protocols and intended comparisons are confirmed.",
        "- Each particle distribution is one experimental unit. Size bins within a distribution are correlated and must not be treated as independent replicates.",
        "- RNA-seq contains both CD4 and HUVEC samples and CD4 EGTA series. The supplied DGE table covers only the named CD4 sham versus 3p40kV contrast.",
    ]
    if not missing_replicates.empty:
        descriptions = [
            f"{row.condition} at {row.harvest_time_h:g} h has n={int(row.n_samples)}"
            for row in missing_replicates.itertuples()
        ]
        lines.extend(["- Incomplete labeled treatment replication: " + "; ".join(descriptions) + "."])
    if all_condition_totals_repeated:
        lines.extend(
            [
                "- Within every labeled condition, all replicates have exactly the same summed bin concentration while their size distributions differ. "
                "This is consistent with condition-level normalization or a shared concentration value and must be clarified before replicate-level concentration fitting.",
            ]
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The generated concentration totals are sums of the reported per-size-bin concentrations. "
            "Confirm with the Exoid export settings that summing bins yields the intended total concentration before using this value as a calibration target. "
            "Particle diameter alone does not uniquely distinguish small EVs, medium or large EVs, apoptotic bodies, aggregates, lipoproteins, or debris.",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit the FFRCI Exoid and RNA-seq data shared for model calibration.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    raw_exoid = parse_exoid_file(args.data_dir / EXOID_FILE)
    sample_summary = summarize_exoid(raw_exoid)
    condition_summary = summarize_conditions(sample_summary)

    raw_exoid.to_csv(args.out / "exoid_particle_bins_tidy.csv", index=False)
    sample_summary.to_csv(args.out / "exoid_sample_summary.csv", index=False)
    condition_summary.to_csv(args.out / "exoid_condition_summary.csv", index=False)
    rnaseq_audit = audit_rnaseq(args.data_dir, args.out)
    make_plots(raw_exoid, sample_summary, args.out / "exoid_preliminary_overview.png")
    write_audit_summary(args.out / "audit_summary.md", sample_summary, condition_summary, rnaseq_audit)

    print(f"Wrote FFRCI data audit to {args.out}")


if __name__ == "__main__":
    main()
