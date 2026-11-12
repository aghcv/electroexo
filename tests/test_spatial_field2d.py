from __future__ import annotations

import numpy as np

from electro_exocytosis.spatial.field2d import (
    FieldGeometry2D,
    sample_field_kV_cm,
    solve_quasistatic_field,
)


def test_quasistatic_field_solution_is_bounded_and_reference_scaled() -> None:
    geometry = FieldGeometry2D(nx=61, ny=51)
    solution = solve_quasistatic_field(geometry)

    assert solution.potential_normalized.shape == (geometry.ny, geometry.nx)
    assert np.isfinite(solution.potential_normalized).all()
    assert np.isfinite(solution.field_magnitude_kV_cm).all()
    assert float(solution.potential_normalized.min()) >= -1e-10
    assert float(solution.potential_normalized.max()) <= 1.0 + 1e-10
    assert np.allclose(solution.potential_normalized[solution.electrode_mask], 1.0)
    assert solution.applied_voltage_V > 0
    assert np.isclose(
        float(solution.conductivity_S_m[solution.vessel_mask][0]),
        geometry.vessel_conductivity_S_m,
    )

    xx_mm, yy_mm = np.meshgrid(solution.x_mm, solution.y_mm)
    electrode_x, electrode_y = geometry.electrode_center_mm
    radius_mm = np.hypot(xx_mm - electrode_x, yy_mm - electrode_y)
    reference_mask = (
        (radius_mm >= geometry.reference_annulus_mm[0])
        & (radius_mm <= geometry.reference_annulus_mm[1])
        & ~solution.electrode_mask
        & ~solution.vessel_mask
    )
    assert np.isclose(
        float(np.median(solution.field_magnitude_kV_cm[reference_mask])),
        geometry.reference_field_kV_cm,
        rtol=0.02,
    )


def test_field_sampler_returns_location_specific_amplitudes() -> None:
    solution = solve_quasistatic_field(FieldGeometry2D(nx=61, ny=51))
    sampled = sample_field_kV_cm(
        solution,
        [(-0.8, 0.7), (0.6, -1.2), (2.8, -2.0)],
    )

    assert sampled.shape == (3,)
    assert np.isfinite(sampled).all()
    assert np.all(sampled > 0)
    assert float(np.ptp(sampled)) > 1.0


def test_conductive_vessel_distorts_the_field_solution() -> None:
    heterogeneous = solve_quasistatic_field(FieldGeometry2D(nx=61, ny=51))
    uniform = solve_quasistatic_field(
        FieldGeometry2D(
            nx=61,
            ny=51,
            vessel_conductivity_S_m=0.20,
        )
    )

    relative_change = np.mean(
        np.abs(
            heterogeneous.field_magnitude_kV_cm
            - uniform.field_magnitude_kV_cm
        )
    ) / np.mean(uniform.field_magnitude_kV_cm)
    assert float(relative_change) > 0.01
