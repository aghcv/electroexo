"""Behavior checks for scientific-data preservation and confirmed pair blocks."""
from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd


MODULE = Path(__file__).resolve().parents[1] / "tools/prepare_ipa_rnaseq.py"
SPEC = importlib.util.spec_from_file_location("prepare_ipa", MODULE)
prep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prep)


def make_source(folder: Path) -> None:
    names = []
    for start, group in [(1, "Control_HUVEC"), (4, "treatment_1_HUVEC"),
                         (7, "Treatment_2_HUVEC"), (10, "Control_CD4"),
                         (13, "Treatment_1_CD4"), (16, "Control_CD4_EGTA_0.125"),
                         (19, "Treatment_1_CD4_EGTA_0.125"), (22, "Control_CD4_0.25_EGTA"),
                         (25, "Treatment_1_CD4_0.25_EGTA"), (28, "Control_CD4_0.5_EGTA"),
                         (31, "Treatment_1_CD4_0.5_EGTA")]:
        names.extend(f"{n}_{group}" for n in range(start, start + 3))
    rows = []
    for gene_id, gene_name, kind, values in [
        ("ENSG00000000001", "TEST", "protein_coding", ["1.25", "5.00"]),
        ("ENSG00000000002", "TEST", "protein_coding", ["8.0", "4"]),
        ("ERCC-00002", "ERCC-00002", "transgene", ["2.0", "2.0"]),
    ]:
        row = {"gene_id": gene_id, "gene_name": gene_name, "gene_biotype": kind}
        for name in sorted(names):  # Mimic lexical sample ordering in the source.
            value = values["treatment" in name.lower()]
            row[name + "_cpm"] = value
            row[name + "_count"] = value
        rows.append(row)
    pd.DataFrame(rows).to_csv(folder / prep.MATRIX, index=False)
    pd.DataFrame([
        {"EnsemblID": "ENSG00000000001", "gene_name": "TEST", "logFC": "2.00", "PValue": "0.0010", "FDR": "0.010"},
        {"EnsemblID": "ENSG00000000002", "gene_name": "TEST", "logFC": "-1.0", "PValue": "0.005", "FDR": "0.020"},
        {"EnsemblID": "# Note: omitted 123 gene(s) outside of your thresholds"},
    ]).to_csv(folder / prep.DGE, index=False)


class PreparationTests(unittest.TestCase):
    def test_lossless_conversion_spikein_separation_and_confirmed_pairs(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp)
            make_source(source)
            original = {p.name: p.read_bytes() for p in source.glob("*.csv")}
            out = source / "prepared"
            summary = prep.prepare(source, out)
            self.assertEqual(summary["human_gene_rows"], 2)
            self.assertEqual(summary["ERCC_feature_rows"], 1)
            self.assertEqual(summary["existing_de_up"], 1)
            self.assertEqual(summary["existing_de_down"], 1)
            self.assertEqual(summary["direction_matches"], 2)
            self.assertFalse(summary["new_differential_expression_computed"])
            full = pd.read_csv(out / "expression/reported_counts_all_samples.tsv", sep="\t", dtype=str)
            self.assertEqual(full.columns[1], "1_Control_HUVEC")
            self.assertEqual(full.loc[0, "10_Control_CD4"], "1.25")
            self.assertEqual(full.loc[0, "13_Treatment_1_CD4"], "5.00")
            self.assertIn("ERCC-00002", full.gene_id.tolist())
            all_human = pd.read_csv(out / "all_expression_for_de/human_reported_counts.csv", dtype=str)
            self.assertEqual(all_human.shape, (2, 34))
            self.assertEqual(all_human.loc[0, "13_Treatment_1_CD4"], "5.00")
            subset = pd.read_csv(out / "paired_de_inputs/CD4_T1_vs_Sham_EGTA0/reported_counts.csv")
            self.assertEqual(subset.shape, (2, 7))
            self.assertTrue(subset.gene_id.str.startswith("ENSG").all())
            manifest = pd.read_csv(out / "design/sample_manifest.csv", keep_default_na=False)
            self.assertEqual(manifest.loc[manifest.Sample == "10_Control_CD4", "pair_id"].iloc[0], "CD4_EGTA0_pair1")
            self.assertEqual(manifest.loc[manifest.Sample == "13_Treatment_1_CD4", "pair_id"].iloc[0], "CD4_EGTA0_pair1")
            self.assertTrue((manifest.pairing_basis != "").all())
            self.assertTrue((manifest.donor_id == "").all())
            ipa = pd.read_csv(out / "ipa_existing/CD4_T1_vs_Sham_236genes_PRELIMINARY.txt", sep="\t", dtype=str)
            self.assertEqual(ipa["Log Ratio"].tolist(), ["2.00", "-1.0"])
            self.assertEqual(ipa["FDR"].tolist(), ["0.010", "0.020"])
            self.assertEqual(len(ipa), 2)  # Footer removed, duplicate symbols not collapsed.
            for name, content in original.items():
                self.assertEqual((source / name).read_bytes(), content)

    def test_rerun_preserves_completed_pairing_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp)
            make_source(source)
            out = source / "prepared"
            prep.prepare(source, out)
            path = out / "design/sample_manifest.csv"
            manifest = pd.read_csv(path, dtype=str, keep_default_na=False)
            manifest.loc[manifest.Sample.isin(["10_Control_CD4", "13_Treatment_1_CD4"]), "pair_id"] = "confirmed_A"
            manifest.to_csv(path, index=False)
            prep.prepare(source, out)
            retained = pd.read_csv(path, dtype=str, keep_default_na=False)
            self.assertEqual(retained.pair_id.tolist().count("confirmed_A"), 2)
            subset = pd.read_csv(out / "paired_de_inputs/CD4_T1_vs_Sham_EGTA0/metadata.csv", keep_default_na=False)
            self.assertEqual(subset.pair_id.tolist().count("confirmed_A"), 2)

    def test_source_metadata_decoding_does_not_assign_batch(self):
        row = prep.sample_metadata("31_Treatment_1_CD4_0.5_EGTA")
        self.assertEqual(row["group"], "CD4_Treatment_1_EGTA_0.5")
        self.assertEqual(row["batch"], "")
        self.assertEqual(row["pair_id"], "CD4_EGTA0p5_pair1")
        self.assertEqual(prep.sample_metadata("4_treatment_1_HUVEC")["treatment"], "Treatment_1")

    def test_explicit_accepted_pairs_propagate_and_conflicts_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp)
            make_source(source)
            out = source / "prepared"
            prep.prepare(source, out)
            rows = [{"contrast_id": "CD4_T1_vs_Sham_EGTA0",
                     "sham_sample": f"{10+i}_Control_CD4",
                     "exposed_sample": f"{13+i}_Treatment_1_CD4",
                     "pair_id": f"approved_pair{i+1}", "pairing_basis": "user approved test fixture"}
                    for i in range(3)]
            pd.DataFrame(rows).to_csv(out / "design/accepted_pairing.csv", index=False)
            result = prep.prepare(source, out)
            self.assertEqual(result["samples_with_pair_id"], 33)
            self.assertEqual(result["accepted_pairing_rows"], 3)
            clean = pd.read_csv(out / "paired_de_inputs/CD4_T1_vs_Sham_EGTA0/metadata_for_import.csv")
            self.assertEqual(clean.pair_id.tolist(), [f"approved_pair{i+1}" for i in range(3)] * 2)
            manifest = pd.read_csv(out / "design/sample_manifest.csv", keep_default_na=False)
            self.assertTrue((manifest.donor_id == "").all())
            rows[0]["pair_id"] = "different_pair"
            pd.DataFrame(rows).to_csv(out / "design/accepted_pairing.csv", index=False)
            with self.assertRaisesRegex(ValueError, "conflicts with completed manifest"):
                prep.prepare(source, out)


if __name__ == "__main__":
    unittest.main()
