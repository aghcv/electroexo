from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True, slots=True)
class FieldGeometry2D:
    """Geometry and conductivity inputs for a 2D quasi-static field solve.

    Coordinates are expressed in millimetres and conductivity in S/m. The
    potential solution is normalized first (needle electrode = 1, outer return
    boundary = 0), then scaled so the median field in a reference annulus equals
    ``reference_field_kV_cm``.
    """

    x_limits_mm: tuple[float, float] = (-4.0, 4.0)
    y_limits_mm: tuple[float, float] = (-3.0, 3.0)
    nx: int = 201
    ny: int = 151
    tissue_conductivity_S_m: float = 0.20
    vessel_conductivity_S_m: float = 1.00
    electrode_center_mm: tuple[float, float] = (-1.25, 0.0)
    electrode_radius_mm: float = 0.16
    vessel_center_mm: tuple[float, float] = (0.80, 0.35)
    vessel_semi_major_mm: float = 2.15
    vessel_semi_minor_mm: float = 0.28
    vessel_angle_deg: float = 23.0
    reference_annulus_mm: tuple[float, float] = (0.32, 0.50)
    reference_field_kV_cm: float = 40.0

    def validate(self) -> None:
        if self.nx < 21 or self.ny < 21:
            raise ValueError("Field grid must contain at least 21 nodes per axis.")
        if self.nx % 2 == 0 or self.ny % 2 == 0:
            raise ValueError("Use odd grid dimensions so geometric centers are represented.")
        if self.tissue_conductivity_S_m <= 0 or self.vessel_conductivity_S_m <= 0:
            raise ValueError("Conductivities must be positive.")
        if self.electrode_radius_mm <= 0:
            raise ValueError("Electrode radius must be positive.")
        if self.reference_field_kV_cm <= 0:
            raise ValueError("Reference field must be positive.")
        if not 0 < self.reference_annulus_mm[0] < self.reference_annulus_mm[1]:
            raise ValueError("Reference annulus radii must be ordered and positive.")


@dataclass(frozen=True, slots=True)
class FieldSolution2D:
    x_mm: np.ndarray
    y_mm: np.ndarray
    conductivity_S_m: np.ndarray
    potential_normalized: np.ndarray
    potential_V: np.ndarray
    field_x_kV_cm: np.ndarray
    field_y_kV_cm: np.ndarray
    field_magnitude_kV_cm: np.ndarray
    electrode_mask: np.ndarray
    vessel_mask: np.ndarray
    applied_voltage_V: float
    geometry: FieldGeometry2D


def build_conductivity_map(
    geometry: FieldGeometry2D,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build the tissue/vessel conductivity map and geometry masks."""

    geometry.validate()
    x_mm = np.linspace(*geometry.x_limits_mm, geometry.nx)
    y_mm = np.linspace(*geometry.y_limits_mm, geometry.ny)
    xx_mm, yy_mm = np.meshgrid(x_mm, y_mm)

    electrode_x, electrode_y = geometry.electrode_center_mm
    electrode_mask = (
        (xx_mm - electrode_x) ** 2 + (yy_mm - electrode_y) ** 2
        <= geometry.electrode_radius_mm**2
    )

    vessel_x, vessel_y = geometry.vessel_center_mm
    angle = np.deg2rad(geometry.vessel_angle_deg)
    u_mm = np.cos(angle) * (xx_mm - vessel_x) + np.sin(angle) * (yy_mm - vessel_y)
    v_mm = -np.sin(angle) * (xx_mm - vessel_x) + np.cos(angle) * (yy_mm - vessel_y)
    vessel_mask = (
        (u_mm / geometry.vessel_semi_major_mm) ** 2
        + (v_mm / geometry.vessel_semi_minor_mm) ** 2
        <= 1.0
    )
    vessel_mask &= ~electrode_mask

    conductivity = np.full_like(xx_mm, geometry.tissue_conductivity_S_m, dtype=float)
    conductivity[vessel_mask] = geometry.vessel_conductivity_S_m
    return x_mm, y_mm, conductivity, electrode_mask, vessel_mask


def solve_quasistatic_field(geometry: FieldGeometry2D | None = None) -> FieldSolution2D:
    """Solve ``div(sigma grad(phi)) = 0`` on a rectangular 2D tissue section.

    The outer boundary is the grounded return and the embedded circular needle
    electrode is held at unit potential. Harmonic face conductivities enforce
    current continuity across the tissue-vessel interface. The normalized
    solution is scaled to the requested reference field after the linear solve.
    """

    geometry = geometry or FieldGeometry2D()
    x_mm, y_mm, conductivity, electrode_mask, vessel_mask = build_conductivity_map(geometry)
    ny, nx = conductivity.shape
    dx_mm = float(x_mm[1] - x_mm[0])
    dy_mm = float(y_mm[1] - y_mm[0])

    boundary_mask = np.zeros_like(electrode_mask, dtype=bool)
    boundary_mask[0, :] = True
    boundary_mask[-1, :] = True
    boundary_mask[:, 0] = True
    boundary_mask[:, -1] = True
    dirichlet_mask = boundary_mask | electrode_mask
    dirichlet_values = np.zeros_like(conductivity, dtype=float)
    dirichlet_values[electrode_mask] = 1.0

    unknown_mask = ~dirichlet_mask
    unknown_index = np.full((ny, nx), -1, dtype=int)
    unknown_index[unknown_mask] = np.arange(int(unknown_mask.sum()))

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    rhs = np.zeros(int(unknown_mask.sum()), dtype=float)

    for row, col in np.argwhere(unknown_mask):
        equation = int(unknown_index[row, col])
        sigma_center = float(conductivity[row, col])
        diagonal = 0.0
        for drow, dcol, spacing_mm in (
            (0, -1, dx_mm),
            (0, 1, dx_mm),
            (-1, 0, dy_mm),
            (1, 0, dy_mm),
        ):
            neighbor_row = row + drow
            neighbor_col = col + dcol
            sigma_neighbor = float(conductivity[neighbor_row, neighbor_col])
            sigma_face = _harmonic_mean(sigma_center, sigma_neighbor)
            coefficient = sigma_face / (spacing_mm**2)
            diagonal += coefficient
            if unknown_mask[neighbor_row, neighbor_col]:
                rows.append(equation)
                cols.append(int(unknown_index[neighbor_row, neighbor_col]))
                values.append(-coefficient)
            else:
                rhs[equation] += coefficient * float(
                    dirichlet_values[neighbor_row, neighbor_col]
                )
        rows.append(equation)
        cols.append(equation)
        values.append(diagonal)

    matrix = coo_matrix(
        (np.asarray(values), (np.asarray(rows), np.asarray(cols))),
        shape=(len(rhs), len(rhs)),
    ).tocsr()
    normalized_unknowns = spsolve(matrix, rhs)
    potential_normalized = dirichlet_values.copy()
    potential_normalized[unknown_mask] = normalized_unknowns

    grad_y_per_mm, grad_x_per_mm = np.gradient(
        potential_normalized, dy_mm, dx_mm, edge_order=2
    )
    normalized_field_per_mm = np.hypot(grad_x_per_mm, grad_y_per_mm)

    xx_mm, yy_mm = np.meshgrid(x_mm, y_mm)
    electrode_x, electrode_y = geometry.electrode_center_mm
    radius_mm = np.hypot(xx_mm - electrode_x, yy_mm - electrode_y)
    reference_mask = (
        (radius_mm >= geometry.reference_annulus_mm[0])
        & (radius_mm <= geometry.reference_annulus_mm[1])
        & ~electrode_mask
        & ~vessel_mask
    )
    reference_gradient_per_mm = float(np.median(normalized_field_per_mm[reference_mask]))
    if not np.isfinite(reference_gradient_per_mm) or reference_gradient_per_mm <= 0:
        raise RuntimeError("Could not determine a positive reference field gradient.")

    target_field_V_mm = geometry.reference_field_kV_cm * 100.0
    applied_voltage_V = target_field_V_mm / reference_gradient_per_mm
    field_x_kV_cm = -grad_x_per_mm * applied_voltage_V / 100.0
    field_y_kV_cm = -grad_y_per_mm * applied_voltage_V / 100.0
    field_magnitude_kV_cm = np.hypot(field_x_kV_cm, field_y_kV_cm)

    return FieldSolution2D(
        x_mm=x_mm,
        y_mm=y_mm,
        conductivity_S_m=conductivity,
        potential_normalized=potential_normalized,
        potential_V=potential_normalized * applied_voltage_V,
        field_x_kV_cm=field_x_kV_cm,
        field_y_kV_cm=field_y_kV_cm,
        field_magnitude_kV_cm=field_magnitude_kV_cm,
        electrode_mask=electrode_mask,
        vessel_mask=vessel_mask,
        applied_voltage_V=float(applied_voltage_V),
        geometry=geometry,
    )


def sample_field_kV_cm(
    solution: FieldSolution2D,
    points_mm: np.ndarray | list[tuple[float, float]],
) -> np.ndarray:
    """Bilinearly sample field magnitude at ``(x_mm, y_mm)`` points."""

    points = np.asarray(points_mm, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points_mm must have shape (n, 2) with x/y coordinates.")
    interpolator = RegularGridInterpolator(
        (solution.y_mm, solution.x_mm),
        solution.field_magnitude_kV_cm,
        bounds_error=True,
    )
    return np.asarray(interpolator(points[:, [1, 0]]), dtype=float)


def _harmonic_mean(left: float, right: float) -> float:
    return 2.0 * left * right / max(left + right, np.finfo(float).tiny)
