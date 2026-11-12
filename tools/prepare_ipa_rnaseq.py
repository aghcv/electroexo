"""Prepare FFRCI RNA-seq tables for paired DE followed by IPA interpretation.

This is a lossless format conversion and descriptive audit, not a DE analysis.
The within-group sham/exposed pair order is user-confirmed; donor identity,
p-values, count rounding, and sample exclusions are not inferred here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/experimental/ffrci_data_sharing"
OUT = SOURCE / "ipa_preparation_20260907"
MATRIX = "RNAseq all samples.csv"
DGE = "RNAseq_CD4 sham vs CD4 3p40kV_DGE (236 genes).csv"
IPA_COLUMNS = ["Ensembl_ID", "Gene_Symbol", "Log Ratio", "p-value", "FDR"]
IMPORT_METADATA_COLUMNS = ["Sample", "cell_type", "treatment", "egta_label", "group", "pair_id"]
PAIRING_BASIS = (
    "user-confirmed within-group positional pairing; donor identity and RNA-to-Exoid batch "
    "crosswalk remain unverified"
)


def positional_pair_id(number: int, cell: str, level: str) -> str:
    """Return the user-confirmed block ID for the known sample layout.

    These are analysis block IDs, not donor IDs. They are deliberately scoped
    to each cell/EGTA comparison so they cannot imply cross-EGTA donor matching.
    """
    if cell == "HUVEC" and 1 <= number <= 9:
        position = ((number - 1) % 3) + 1
        return f"HUVEC_EGTA0_pair{position}"
    starts = {"0": 10, "0.125": 16, "0.25": 22, "0.5": 28}
    if cell == "CD4" and level in starts and starts[level] <= number < starts[level] + 6:
        position = ((number - starts[level]) % 3) + 1
        slug = level.replace(".", "p")
        return f"CD4_EGTA{slug}_pair{position}"
    return ""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sample_metadata(name: str) -> dict[str, str]:
    number = int(name.split("_")[0])
    cell = "HUVEC" if "HUVEC" in name else "CD4"
    treatment = "Sham" if "control" in name.lower() else (
        "Treatment_2" if "Treatment_2" in name else "Treatment_1"
    )
    level = next((x for x in ["0.125", "0.25", "0.5"] if x in name), "0")
    group = f"{cell}_{treatment}_EGTA_{level}"
    pair_id = positional_pair_id(number, cell, level)
    return {
        "Sample": name, "sample_number": str(number), "cell_type": cell,
        "treatment": treatment, "egta_label": level,
        "group": group, "pair_id": pair_id, "donor_id": "", "batch": "",
        "pairing_basis": PAIRING_BASIS if pair_id else "",
        "biological_replicate_confirmed": "yes" if pair_id else "", "rna_source": "",
        "rna_collection_time_h": "", "exposure_protocol": "", "egta_units": "",
        "egta_timing": "", "use_for_de": "pending", "exclusion_reason": "",
        "counts_provenance_confirmed": "", "qc_review_notes": "",
    }


def primary_contrasts() -> list[dict[str, str]]:
    rows = []
    for cell, treatment, level in [
        ("CD4", "Treatment_1", "0"), ("CD4", "Treatment_1", "0.125"),
        ("CD4", "Treatment_1", "0.25"), ("CD4", "Treatment_1", "0.5"),
        ("HUVEC", "Treatment_1", "0"), ("HUVEC", "Treatment_2", "0"),
    ]:
        slug = level.replace(".", "p")
        short = treatment.replace("Treatment_", "T")
        rows.append({
            "contrast_id": f"{cell}_{short}_vs_Sham_EGTA{slug}",
            "analysis_role": "primary_treatment_effect" if cell == "CD4" else "secondary_cell_type",
            "cell_type": cell, "numerator": f"{cell}_{treatment}_EGTA_{level}",
            "denominator": f"{cell}_Sham_EGTA_{level}",
            "positive_log2FC": "higher in exposed cells than their matched sham",
            "model": "paired RNA-seq model with confirmed donor or pair blocking",
            "status": "requires pairing provenance and QC before DE",
        })
    return rows


def write_rows(path: Path, columns: list[str], rows: list[dict], delimiter: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(source: Path, output: Path) -> dict:
    before = {p.name: sha256(p) for p in sorted(source.glob("*.csv"))}
    # Strings preserve the precision and spelling of every source expression value.
    raw = pd.read_csv(source / MATRIX, dtype=str, keep_default_na=False)
    ids = raw.gene_id
    human = ids.str.fullmatch(r"ENSG\d+(?:\.\d+)?")
    spikein = ids.str.fullmatch(r"ERCC-\d+")
    if not (human | spikein).all() or ids.duplicated().any():
        raise ValueError("Expected unique Ensembl or ERCC feature IDs; do not silently collapse rows.")
    sample_names = sorted([c[:-6] for c in raw if c.endswith("_count")], key=lambda s: int(s.split("_")[0]))
    if set(sample_names) != {c[:-4] for c in raw if c.endswith("_cpm")}:
        raise ValueError("CPM and count sample sets do not match.")
    meta = [sample_metadata(s) for s in sample_names]
    contrasts = primary_contrasts()
    output.mkdir(parents=True, exist_ok=True)

    manifest = output / "design/sample_manifest.csv"
    # Re-running preparation must never erase the user's completed design fields.
    if not manifest.exists():
        write_rows(manifest, list(meta[0]), meta, ",")
    saved_meta = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    if set(saved_meta.Sample) != set(sample_names) or saved_meta.Sample.duplicated().any():
        raise ValueError("Existing sample manifest does not match source samples.")
    # Migrate manifests created by earlier package versions, then populate the
    # newly confirmed default block IDs without touching any completed custom
    # pair/donor assignments.
    if "pairing_basis" not in saved_meta:
        saved_meta["pairing_basis"] = ""
    for default in meta:
        selected_row = saved_meta.Sample == default["Sample"]
        if not selected_row.any():
            continue
        if not saved_meta.loc[selected_row, "pair_id"].iloc[0] and default["pair_id"]:
            saved_meta.loc[selected_row, "pair_id"] = default["pair_id"]
        if not saved_meta.loc[selected_row, "pairing_basis"].iloc[0] and default["pair_id"]:
            saved_meta.loc[selected_row, "pairing_basis"] = PAIRING_BASIS
        if "biological_replicate_confirmed" in saved_meta and not saved_meta.loc[selected_row, "biological_replicate_confirmed"].iloc[0] and default["pair_id"]:
            saved_meta.loc[selected_row, "biological_replicate_confirmed"] = "yes"
    saved_meta.to_csv(manifest, index=False, lineterminator="\n")
    # The positional pairing was confirmed by the user. Preserve any later
    # laboratory-record assignments, and allow accepted_pairing.csv to replace
    # only these default block IDs with a more authoritative crosswalk.
    accepted_file = output / "design/accepted_pairing.csv"
    accepted = pd.read_csv(accepted_file, dtype=str, keep_default_na=False) if accepted_file.exists() else pd.DataFrame()
    if not accepted.empty:
        assignments = {}
        definitions = {c["contrast_id"]: c for c in contrasts}
        if accepted.duplicated(["contrast_id", "sham_sample"]).any() or accepted.duplicated(["contrast_id", "exposed_sample"]).any():
            raise ValueError("Duplicate samples within an accepted pairing contrast.")
        source_groups = {r["Sample"]: r["group"] for r in meta}
        for row in accepted.to_dict(orient="records"):
            contrast = definitions[row["contrast_id"]]
            if not row["pair_id"] or not row["pairing_basis"]:
                raise ValueError("Accepted pairing requires a block ID and recorded provenance.")
            for field, group_key in [("sham_sample", "denominator"), ("exposed_sample", "numerator")]:
                sample = row[field]
                if source_groups.get(sample) != contrast[group_key]:
                    raise ValueError(f"Accepted pairing sample is outside its contrast: {sample}")
                value = (row["pair_id"], row["pairing_basis"])
                if sample in assignments and assignments[sample] != value:
                    raise ValueError(f"Conflicting accepted pair assignments: {sample}")
                assignments[sample] = value
        if "pairing_basis" not in saved_meta:
            saved_meta["pairing_basis"] = ""
        for sample, (pair_id, basis) in assignments.items():
            selected_row = saved_meta.Sample == sample
            current = saved_meta.loc[selected_row, "pair_id"].iloc[0]
            current_basis = saved_meta.loc[selected_row, "pairing_basis"].iloc[0]
            if current and current != pair_id and current_basis != PAIRING_BASIS:
                raise ValueError(f"Accepted pairing conflicts with completed manifest: {sample}")
            saved_meta.loc[selected_row, ["pair_id", "pairing_basis"]] = [pair_id, basis]
        saved_meta.to_csv(manifest, index=False, lineterminator="\n")
    meta_lookup = saved_meta.set_index("Sample").to_dict(orient="index")

    for kind in ["count", "cpm"]:
        frame = raw[["gene_id"] + [s + "_" + kind for s in sample_names]].copy()
        frame.columns = ["gene_id"] + sample_names
        dest = output / f"expression/{'reported_counts' if kind == 'count' else 'cpm'}_all_samples.tsv"
        dest.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(dest, sep="\t", index=False, lineterminator="\n")
        check = pd.read_csv(dest, sep="\t", dtype=str, keep_default_na=False)
        if not check.equals(frame):
            raise AssertionError(f"Lossless round-trip failed: {dest}")
    raw[["gene_id", "gene_name", "gene_biotype"]].to_csv(
        output / "expression/gene_annotations.tsv", sep="\t", index=False, lineterminator="\n"
    )
    raw.loc[spikein, ["gene_id", "gene_name", "gene_biotype"]].to_csv(
        output / "expression/spike_in_features.tsv", sep="\t", index=False, lineterminator="\n"
    )
    all_input = output / "all_expression_for_de"
    all_input.mkdir(exist_ok=True)
    all_counts = raw.loc[human, ["gene_id"] + [s + "_count" for s in sample_names]].copy()
    all_counts.columns = ["gene_id"] + sample_names
    all_counts.to_csv(all_input / "human_reported_counts.csv", index=False, lineterminator="\n")
    check = pd.read_csv(all_input / "human_reported_counts.csv", dtype=str, keep_default_na=False)
    if not check.equals(all_counts.reset_index(drop=True)):
        raise AssertionError("All-human count matrix failed lossless round-trip.")
    saved_meta.set_index("Sample").loc[sample_names].reset_index()[IMPORT_METADATA_COLUMNS].to_csv(
        all_input / "sample_metadata.csv", index=False, lineterminator="\n"
    )

    for contrast in contrasts:
        selected = [r["Sample"] for r in meta if r["group"] in {contrast["numerator"], contrast["denominator"]}]
        for group in [contrast["numerator"], contrast["denominator"]]:
            if sum(r["group"] == group for r in meta) != 3:
                raise ValueError(f"Unexpected group size: {group}")
        frame = raw.loc[human, ["gene_id"] + [s + "_count" for s in selected]].copy()
        frame.columns = ["gene_id"] + selected
        folder = output / "paired_de_inputs" / contrast["contrast_id"]
        folder.mkdir(parents=True, exist_ok=True)
        frame.to_csv(folder / "reported_counts.csv", index=False, lineterminator="\n")
        restored = pd.read_csv(folder / "reported_counts.csv", dtype=str, keep_default_na=False)
        if not restored.equals(frame.reset_index(drop=True)):
            raise AssertionError(f"Lossless subset round-trip failed: {folder}")
        metadata_rows = [{"Sample": s, **meta_lookup[s]} for s in selected]
        write_rows(folder / "metadata.csv", list(saved_meta.columns), metadata_rows, ",")
        clean_metadata = pd.DataFrame(metadata_rows)[IMPORT_METADATA_COLUMNS]
        clean_metadata.to_csv(folder / "metadata_for_import.csv", index=False, lineterminator="\n")
        if clean_metadata.pair_id.ne("").all():
            blocks = pd.crosstab(clean_metadata.pair_id, clean_metadata.treatment)
            if blocks.shape != (3, 2) or not (blocks == 1).all().all():
                raise ValueError(f"Expected three complete sham/exposed pairs: {contrast['contrast_id']}")
            contrast["status"] = "pair IDs populated; review count provenance and QC before DE"
        contrast["samples_numerator"] = ";".join(r["Sample"] for r in meta if r["group"] == contrast["numerator"])
        contrast["samples_denominator"] = ";".join(r["Sample"] for r in meta if r["group"] == contrast["denominator"])
    write_rows(output / "design/primary_contrasts.csv", list(contrasts[0]), contrasts, ",")
    candidates = []
    confirmed = []
    for contrast in contrasts:
        sham = contrast["samples_denominator"].split(";")
        exposed_samples = contrast["samples_numerator"].split(";")
        for order, (control_sample, exposed_sample) in enumerate(zip(sham, exposed_samples), start=1):
            pair_id = meta_lookup[control_sample]["pair_id"]
            candidates.append({"contrast_id": contrast["contrast_id"], "within_group_order": order,
                               "sham_sample": control_sample, "exposed_sample": exposed_sample,
                               "candidate_exoid_label_if_order_matches": f"p{order}" if contrast["cell_type"] == "CD4" else "",
                               "evidence": PAIRING_BASIS,
                               "confirmed_pair": "yes", "used_in_analysis_metadata": "yes"})
            confirmed.append({"contrast_id": contrast["contrast_id"], "within_group_order": order,
                              "sham_sample": control_sample, "exposed_sample": exposed_sample,
                              "pair_id": pair_id, "pairing_basis": PAIRING_BASIS,
                              "exoid_label_if_order_matches": f"p{order}" if contrast["cell_type"] == "CD4" else "",
                              "exoid_crosswalk_confirmed": "no"})
    # Preserve this historical proposal; accepted_pairing.csv supersedes it.
    candidate_file = output / "design/pairing_candidates_NOT_CONFIRMED.csv"
    if not candidate_file.exists():
        write_rows(candidate_file, list(candidates[0]), candidates, ",")
    write_rows(output / "design/confirmed_pairing.csv", list(confirmed[0]), confirmed, ",")

    extra = []
    for level in ["0.125", "0.25", "0.5"]:
        slug = level.replace(".", "p")
        extra.append({"contrast_id": f"CD4_EGTA{slug}_interaction",
                      "contrast": f"(Treatment_1 - Sham) at EGTA {level} minus (Treatment_1 - Sham) at EGTA 0",
                      "meaning_of_positive_log2FC": "treatment effect shifts upward with EGTA; not necessarily more expression or more calcium",
                      "prerequisite": "joint model with correct donor and EGTA design; cannot subtract IPA z-scores or p-values"})
        extra.append({"contrast_id": f"CD4_Sham_EGTA{slug}_vs_EGTA0",
                      "contrast": f"Sham at EGTA {level} minus Sham at EGTA 0",
                      "meaning_of_positive_log2FC": "higher sham expression with EGTA",
                      "prerequisite": "cross-EGTA donor correspondence must be confirmed"})
    extra.append({"contrast_id": "HUVEC_T2_vs_T1", "contrast": "HUVEC Treatment_2 minus Treatment_1",
                  "meaning_of_positive_log2FC": "higher in Treatment_2", "prerequisite": "confirm protocol and donor correspondence; secondary QC-sensitive comparison"})
    write_rows(output / "design/additional_contrasts.csv", list(extra[0]), extra, ",")

    counts = raw[[s + "_count" for s in sample_names]].apply(pd.to_numeric)
    cpm = raw[[s + "_cpm" for s in sample_names]].apply(pd.to_numeric)
    if not np.isfinite(counts.to_numpy()).all() or not np.isfinite(cpm.to_numpy()).all() or (counts < 0).any().any() or (cpm < 0).any().any():
        raise ValueError("Expression has missing, nonfinite, or negative values.")
    count_totals = counts.sum(axis=0).to_numpy()
    qc = []
    for i, m in enumerate(meta):
        same = [j for j, other in enumerate(meta) if other["group"] == m["group"]]
        group_median = float(np.median(count_totals[same]))
        ratio = float(count_totals[i] / group_median)
        flags = []
        if count_totals[i] < 1e6:
            flags.append("reported_count_total_below_1million")
        if ratio < 0.25:
            flags.append("reported_count_total_below_25pct_group_median")
        qc.append({"Sample": m["Sample"], "group": m["group"],
                   "reported_count_total": float(count_totals[i]),
                   "human_reported_count_total": float(counts.loc[human].iloc[:, i].sum()),
                   "ERCC_reported_count_total": float(counts.loc[spikein].iloc[:, i].sum()),
                   "sum_cpm": float(cpm.iloc[:, i].sum()),
                   "genes_with_positive_count": int((counts.iloc[:, i] > 0).sum()),
                   "fractional_count_entries": int((np.abs(counts.iloc[:, i] - np.rint(counts.iloc[:, i])) > 1e-8).sum()),
                   "total_relative_to_group_median": ratio,
                   "review_flags": ";".join(flags), "excluded_by_preparation": "no"})
    write_rows(output / "qc/sample_qc.csv", list(qc[0]), qc, ",")
    # Modest descriptive gene detection summaries, explicitly not a DE universe.
    detected = []
    for contrast in contrasts:
        nums = [i for i, m in enumerate(meta) if m["group"] == contrast["numerator"]]
        dens = [i for i, m in enumerate(meta) if m["group"] == contrast["denominator"]]
        mask = human & ((cpm.iloc[:, nums].ge(1).sum(axis=1) >= 2) | (cpm.iloc[:, dens].ge(1).sum(axis=1) >= 2))
        detected.append({"contrast_id": contrast["contrast_id"], "genes_CPM_ge1_in_2_of_3_in_either_group": int(mask.sum()),
                         "interpretation": "descriptive only; use actual fitted-model tested genes as IPA reference"})
    write_rows(output / "qc/detection_summary.csv", list(detected[0]), detected, ",")

    existing = pd.read_csv(source / DGE, dtype=str, keep_default_na=False)
    footer = existing.loc[~existing.EnsemblID.str.fullmatch(r"ENSG\d+"), "EnsemblID"].tolist()
    existing = existing[existing.EnsemblID.str.fullmatch(r"ENSG\d+")].copy()
    if existing.EnsemblID.duplicated().any():
        raise ValueError("Duplicate IDs in source DE table.")
    for column in ["logFC", "PValue", "FDR"]:
        values = pd.to_numeric(existing[column])
        if not np.isfinite(values).all() or (column != "logFC" and ((values < 0) | (values > 1)).any()):
            raise ValueError(f"Invalid source DE values in {column}")
    ipa = existing[["EnsemblID", "gene_name", "logFC", "PValue", "FDR"]].copy()
    ipa.columns = IPA_COLUMNS
    dest = output / "ipa_existing/CD4_T1_vs_Sham_236genes_PRELIMINARY.txt"
    dest.parent.mkdir(parents=True, exist_ok=True)
    ipa.to_csv(dest, sep="\t", index=False, lineterminator="\n")
    assert pd.read_csv(dest, sep="\t", dtype=str, keep_default_na=False).equals(ipa)

    joined = existing.merge(raw, left_on="EnsemblID", right_on="gene_id", validate="one_to_one")
    control = joined[[s + "_cpm" for s in sample_names if s in ["10_Control_CD4", "11_Control_CD4", "12_Control_CD4"]]].apply(pd.to_numeric).mean(axis=1)
    exposed = joined[[s + "_cpm" for s in sample_names if s in ["13_Treatment_1_CD4", "14_Treatment_1_CD4", "15_Treatment_1_CD4"]]].apply(pd.to_numeric).mean(axis=1)
    observed_logfc = pd.to_numeric(joined.logFC)
    ratio = np.log2(exposed / control)
    audit = pd.DataFrame({"Ensembl_ID": joined.EnsemblID, "Gene_Symbol": joined.gene_name_x,
                          "source_log2FC": observed_logfc, "sham_mean_source_CPM": control,
                          "treatment1_mean_source_CPM": exposed, "descriptive_log2_mean_CPM_ratio": ratio,
                          "sign_matches": np.sign(observed_logfc) == np.sign(ratio)})
    audit.to_csv(output / "qc/existing_de_direction_check.tsv", sep="\t", index=False, lineterminator="\n")

    summary = {
        "source_feature_rows": len(raw), "human_gene_rows": int(human.sum()),
        "ERCC_feature_rows": int(spikein.sum()),
        "ERCC_handling": "preserved in full expression matrices; excluded explicitly from human comparison inputs",
        "samples": len(sample_names), "primary_contrasts": len(contrasts),
        "existing_de_genes": len(existing), "existing_de_up": int((observed_logfc > 0).sum()),
        "existing_de_down": int((observed_logfc < 0).sum()), "source_de_footer": footer,
        "existing_de_min_absolute_log2FC": float(observed_logfc.abs().min()),
        "existing_de_max_FDR": float(pd.to_numeric(existing.FDR).max()),
        "direction_matches": int(audit.sign_matches.sum()),
        "direction_correlation": float(np.corrcoef(observed_logfc, ratio)[0, 1]),
        "direction_conclusion": "Numerical evidence strongly supports Treatment_1 relative to Sham for the supplied CD4 DE table; original model remains undocumented.",
        "fractional_count_entries": int(sum(row["fractional_count_entries"] for row in qc)),
        "blank_gene_symbols": int((raw.gene_name == "").sum()),
        "duplicate_nonblank_symbol_rows": int(raw.loc[raw.gene_name != "", "gene_name"].duplicated().sum()),
        "sample_qc_flagged": [row["Sample"] for row in qc if row["review_flags"]],
        "new_differential_expression_computed": False, "statistical_pairing_inferred": False,
        "pairing_status": "user_confirmed_by_within_group_order",
        "pairing_basis": PAIRING_BASIS,
        "accepted_pairing_rows": len(accepted),
        "samples_with_pair_id": int(saved_meta.pair_id.ne("").sum()),
        "samples_removed": 0, "counts_rounded": False,
        "source_files_sha256": before,
    }
    if before != {p.name: sha256(p) for p in sorted(source.glob("*.csv"))}:
        raise AssertionError("Source files changed during preparation.")
    (output / "qc/preparation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_rows(output / "templates/IPA_full_DE_table_template.tsv", IPA_COLUMNS, [])
    (output / "ipa_exports").mkdir(exist_ok=True)
    logs = ["analysis_id", "contrast_id", "ipa_version", "knowledge_base_release", "run_date",
            "de_method_version", "model_formula", "donor_blocking", "included_samples", "excluded_samples",
            "counts_provenance", "input_file_sha256", "tested_genes", "mapped_genes", "unmapped_genes",
            "analysis_ready_genes", "reference_set", "gene_FDR_cutoff", "absolute_log2FC_cutoff",
            "relationship_filters", "species_filters", "tissue_filters", "export_filename", "notes"]
    if not (output / "ipa_exports/analysis_log.csv").exists():
        write_rows(output / "ipa_exports/analysis_log.csv", logs, [], ",")
    hashes = [{"file": str(p.relative_to(output)), "sha256": sha256(p), "bytes": p.stat().st_size}
              for p in sorted(output.rglob("*")) if p.is_file() and p.name != "file_manifest.csv"]
    write_rows(output / "file_manifest.csv", ["file", "sha256", "bytes"], hashes, ",")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.out), indent=2))
