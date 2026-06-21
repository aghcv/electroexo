from __future__ import annotations

from pathlib import Path

import pandas as pd


CALIBRATION_TARGETS_FILE = Path(__file__).with_name("calibration_targets.csv")


class EvidenceLoader:
    """
    Loads and indexes the nsPEF EV literature evidence workbook.
    Currently reads sheets into DataFrames and maps rows to model layers.
    Future versions will use this to override placeholder parameters.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._loaded = False
        self._module_map = pd.DataFrame()
        self._literature_tracker = pd.DataFrame()
        self._calibration_targets = pd.DataFrame()

    def load(self) -> bool:
        """Load workbook sheets, returning False if the file is unavailable."""
        self._calibration_targets = self._load_calibration_targets()
        if not self.path.exists():
            self._loaded = False
            self._module_map = pd.DataFrame()
            self._literature_tracker = pd.DataFrame()
            return False
        self._module_map = pd.read_excel(self.path, sheet_name="Module Map")
        self._literature_tracker = pd.read_excel(self.path, sheet_name="Literature Tracker")
        self._loaded = True
        return True

    @staticmethod
    def _load_calibration_targets() -> pd.DataFrame:
        """Load curated full-text calibration targets bundled with the package."""
        if not CALIBRATION_TARGETS_FILE.exists():
            return pd.DataFrame()
        return pd.read_csv(CALIBRATION_TARGETS_FILE)

    def get_module_map(self) -> pd.DataFrame:
        """Return the module map sheet."""
        return self._module_map.copy()

    def get_literature_tracker(self) -> pd.DataFrame:
        """Return the literature tracker sheet."""
        return self._literature_tracker.copy()

    def get_calibration_targets(self, module: str | None = None) -> pd.DataFrame:
        """Return full-text-derived calibration targets, optionally filtered by module."""
        targets = self._calibration_targets.copy()
        if module and not targets.empty and "module" in targets.columns:
            targets = targets[targets["module"].astype(str).str.casefold() == module.casefold()]
        return targets

    def get_placeholder_report(self) -> dict[str, object]:
        """Report module placeholder coverage."""
        layers = []
        if self._loaded and not self._module_map.empty and "Layer" in self._module_map.columns:
            layers = [str(value) for value in self._module_map["Layer"].dropna().tolist()]
        return {
            "loaded": self._loaded,
            "placeholder_fraction": 1.0,
            "placeholder_modules": [
                "pulse",
                "dosimetry",
                "electrodynamics",
                "ion_transport",
                "remodeling_repair",
                "ev_release",
                "cargo_potency",
                "injury_quality",
                "manufacturing_qc",
                "cell_state",
            ],
            "layers": layers,
            "calibration_target_count": int(len(self._calibration_targets)),
        }

    def summarize(self) -> str:
        """Return a short text summary of workbook status."""
        if not self._loaded:
            return f"Evidence workbook not loaded: {self.path}"
        return (
            f"Evidence workbook loaded with {len(self._literature_tracker)} literature rows "
            f"and {len(self._module_map)} module-map rows. "
            f"Full-text calibration targets: {len(self._calibration_targets)}."
        )
