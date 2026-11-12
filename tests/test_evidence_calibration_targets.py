from __future__ import annotations

from pathlib import Path

from electro_exocytosis.evidence.evidence_loader import EvidenceLoader


def test_full_text_calibration_targets_are_available_without_workbook() -> None:
    loader = EvidenceLoader(Path("missing_workbook.xlsx"))
    assert not loader.load()

    targets = loader.get_calibration_targets()
    assert len(targets) >= 10
    assert {"target_id", "module", "source_key", "model_targets", "priority"}.issubset(targets.columns)
    assert "bhattacharya2022calcium" in set(targets["source_key"])


def test_full_text_calibration_targets_can_be_filtered_by_module() -> None:
    loader = EvidenceLoader(Path("missing_workbook.xlsx"))
    loader.load()

    repair_targets = loader.get_calibration_targets("A7")
    assert not repair_targets.empty
    assert set(repair_targets["module"]) == {"A7"}
    assert "tau_repair_s" in " ".join(repair_targets["model_targets"].astype(str))
