from __future__ import annotations

"""Physical constants and unit conversion helpers."""

c_light: float = 299_792_458.0
"""Speed of light in vacuum, m/s."""

epsilon_0: float = 8.854_187_812_8e-12
"""Vacuum permittivity, F/m."""

k_B: float = 1.380_649e-23
"""Boltzmann constant, J/K."""

F_faraday: float = 96_485.332_12
"""Faraday constant, C/mol."""

R_gas: float = 8.314_462_618
"""Universal gas constant, J/(mol*K)."""

N_A: float = 6.022_140_76e23
"""Avogadro constant, 1/mol."""


def ns_to_s(value_ns: float) -> float:
    """Convert nanoseconds to seconds."""
    return value_ns / 1e9



def kV_cm_to_V_m(value_kV_cm: float) -> float:
    """Convert kilovolts per centimeter to volts per meter."""
    return value_kV_cm * 1e5



def celsius_to_kelvin(value_celsius: float) -> float:
    """Convert temperature in degrees Celsius to Kelvin."""
    return value_celsius + 273.15
