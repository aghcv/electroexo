"""Reversible scaling between experimental concentration and single-cell output.

The mechanistic simulator represents one producer cell. Particle instruments
usually report concentration in conditioned medium. This module keeps the
population, volume, viability, processing-recovery, dilution, and background
assumptions outside the intracellular model and makes the conversion explicit.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml


CellBasis = Literal["initial", "viable"]


@dataclass(frozen=True, slots=True)
class ExperimentalObservationBridge:
    """Map particles/mL measurements to particle equivalents per model cell.

    ``dilution_factor`` is the factor by which an instrument result must be
    multiplied to recover the concentration before assay dilution. A ten-fold
    diluted sample therefore uses ``dilution_factor=10``.

    ``recovery_fraction`` is the fraction of particles retained through sample
    preparation and measurement. It is divided out when inferring released
    particles per cell and applied when mapping model output back to a measured
    concentration.
    """

    initial_cell_count: float
    medium_volume_ml: float
    viability_fraction: float = 1.0
    recovery_fraction: float = 1.0
    dilution_factor: float = 1.0
    background_concentration_particles_per_ml: float = 0.0
    cell_basis: CellBasis = "initial"

    def __post_init__(self) -> None:
        if not np.isfinite(self.initial_cell_count) or self.initial_cell_count <= 0:
            raise ValueError("initial_cell_count must be finite and positive")
        if not np.isfinite(self.medium_volume_ml) or self.medium_volume_ml <= 0:
            raise ValueError("medium_volume_ml must be finite and positive")
        if not np.isfinite(self.viability_fraction) or not 0 < self.viability_fraction <= 1:
            raise ValueError("viability_fraction must be in (0, 1]")
        if not np.isfinite(self.recovery_fraction) or not 0 < self.recovery_fraction <= 1:
            raise ValueError("recovery_fraction must be in (0, 1]")
        if not np.isfinite(self.dilution_factor) or self.dilution_factor <= 0:
            raise ValueError("dilution_factor must be finite and positive")
        if (
            not np.isfinite(self.background_concentration_particles_per_ml)
            or self.background_concentration_particles_per_ml < 0
        ):
            raise ValueError(
                "background_concentration_particles_per_ml must be finite and nonnegative"
            )
        if self.cell_basis not in {"initial", "viable"}:
            raise ValueError("cell_basis must be 'initial' or 'viable'")

    @property
    def effective_cell_count(self) -> float:
        """Cell denominator selected for the experimental normalization."""

        if self.cell_basis == "viable":
            return self.initial_cell_count * self.viability_fraction
        return self.initial_cell_count

    @property
    def effective_cell_density_per_ml(self) -> float:
        """Selected cell denominator divided by conditioned-medium volume."""

        return self.effective_cell_count / self.medium_volume_ml

    def concentration_to_particles_per_cell(
        self, concentration_particles_per_ml: float | Sequence[float] | np.ndarray
    ) -> float | np.ndarray:
        """Convert measured particles/mL to released particles per selected cell."""

        concentration = np.asarray(concentration_particles_per_ml, dtype=float)
        if np.any(~np.isfinite(concentration)) or np.any(concentration < 0):
            raise ValueError("Measured concentration must be finite and nonnegative")
        background_corrected = np.clip(
            concentration - self.background_concentration_particles_per_ml,
            0.0,
            None,
        )
        corrected_concentration = (
            background_corrected * self.dilution_factor / self.recovery_fraction
        )
        particles_per_cell = (
            corrected_concentration
            * self.medium_volume_ml
            / self.effective_cell_count
        )
        if particles_per_cell.ndim == 0:
            return float(particles_per_cell)
        return particles_per_cell

    def particles_per_cell_to_concentration(
        self, particles_per_cell: float | Sequence[float] | np.ndarray
    ) -> float | np.ndarray:
        """Map model particles/cell to the corresponding measured particles/mL."""

        per_cell = np.asarray(particles_per_cell, dtype=float)
        if np.any(~np.isfinite(per_cell)) or np.any(per_cell < 0):
            raise ValueError("Particles per cell must be finite and nonnegative")
        measured = (
            per_cell
            * self.effective_cell_count
            / self.medium_volume_ml
            * self.recovery_fraction
            / self.dilution_factor
            + self.background_concentration_particles_per_ml
        )
        if measured.ndim == 0:
            return float(measured)
        return measured

    def transform_frame(
        self,
        frame: pd.DataFrame,
        *,
        concentration_column: str,
        output_column: str = "observed_particles_per_cell",
    ) -> pd.DataFrame:
        """Return a copy with a single-cell observation column added."""

        if concentration_column not in frame:
            raise ValueError(f"Missing concentration column: {concentration_column}")
        transformed = frame.copy()
        transformed[output_column] = self.concentration_to_particles_per_cell(
            transformed[concentration_column].to_numpy(dtype=float)
        )
        return transformed

    def to_metadata(self) -> dict[str, float | str]:
        """Return JSON/YAML-safe settings plus derived density."""

        return {
            **asdict(self),
            "effective_cell_count": self.effective_cell_count,
            "effective_cell_density_per_ml": self.effective_cell_density_per_ml,
            "single_cell_output_unit": "particle_equivalents_per_cell",
        }

    @classmethod
    def from_mapping(
        cls, values: Mapping[str, float | str]
    ) -> "ExperimentalObservationBridge":
        """Build from either a direct mapping or an ``experimental_bridge`` block."""

        nested = values.get("experimental_bridge")
        config = nested if isinstance(nested, Mapping) else values
        allowed = {
            "initial_cell_count",
            "medium_volume_ml",
            "viability_fraction",
            "recovery_fraction",
            "dilution_factor",
            "background_concentration_particles_per_ml",
            "cell_basis",
        }
        unknown = set(config) - allowed
        if unknown:
            raise ValueError(f"Unknown experimental bridge fields: {sorted(unknown)}")
        return cls(**dict(config))

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ExperimentalObservationBridge":
        """Load bridge settings from a YAML file."""

        with Path(path).open("r", encoding="utf-8") as handle:
            values = yaml.safe_load(handle) or {}
        if not isinstance(values, Mapping):
            raise ValueError("Experimental bridge YAML must contain a mapping")
        return cls.from_mapping(values)

