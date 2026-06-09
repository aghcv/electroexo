from __future__ import annotations

from typing import Any

import numpy as np
from scipy.integrate import solve_ivp



def solve_ode_system(rhs, y0, t_span, t_eval, method: str = "RK45") -> Any:
    """Wrapper around scipy.integrate.solve_ivp with error handling."""
    result = solve_ivp(rhs, t_span=t_span, y0=y0, t_eval=t_eval, method=method, vectorized=False, rtol=1e-6, atol=1e-8)
    if not result.success:
        raise RuntimeError(f"ODE solve failed: {result.message}")
    return result



def build_multiscale_time_array(t_start, t_end, dt_coarse) -> np.ndarray:
    """Build time array suitable for multi-scale integration."""
    if dt_coarse <= 0:
        raise ValueError("dt_coarse must be positive")
    n_steps = int(np.floor((t_end - t_start) / dt_coarse)) + 1
    t_array = t_start + np.arange(n_steps, dtype=float) * dt_coarse
    if t_array[-1] < t_end:
        t_array = np.append(t_array, float(t_end))
    return t_array
