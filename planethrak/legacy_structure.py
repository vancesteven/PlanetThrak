"""Literal legacy radial structure helpers used to establish MATLAB parity.

These functions intentionally preserve the simplifying assumptions in
``plot_PTCracking_planets.m``. They are reference fixtures for the Python port,
not the intended final Mars structure model. Modern calculations should consume
PlanetProfile columns through the neutral adapter and compare against these
fixtures explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .radiogenic import past_radiogenic_heat_uthk


@dataclass(frozen=True)
class RadialColumn:
    depth_m: np.ndarray
    pressure_MPa: np.ndarray
    temperature_C: np.ndarray
    radiogenic_power_W: float
    surface_heat_flux_Wm2: float


# Literal Mars constants in plot_PTCracking_planets.m.
LEGACY_MARS_RADIUS_M = 3397.0e3
LEGACY_MARS_MASS_KG = 6.42e23
LEGACY_MARS_SILICATE_MASS_FRACTION = 1.0
LEGACY_MARS_ROCK_DENSITY_KGM3 = 3500.0
LEGACY_MARS_GRAVITY_MS2 = 6.0
LEGACY_THERMAL_CONDUCTIVITY_WMK = 3.0


def legacy_mars_pt_on_depths(
    depth_m,
    time_before_present_yr: float,
) -> RadialColumn:
    """Evaluate the archived Mars P-T approximation on arbitrary depths.

    This helper is intentionally literal. It exists so a modern PlanetProfile
    radial grid can be compared to the historical structure without conflating
    a grid change with a pressure or thermal-physics change.
    """
    if time_before_present_yr < 0:
        raise ValueError("time_before_present_yr must be nonnegative")
    z = np.asarray(depth_m, dtype=float).reshape(-1)
    if z.size < 2 or not np.all(np.isfinite(z)):
        raise ValueError("depth_m must contain at least two finite values")
    if np.any(z < 0) or np.any(z > LEGACY_MARS_RADIUS_M):
        raise ValueError("depth_m must lie within the legacy Mars radius")
    if np.any(np.diff(z) <= 0):
        raise ValueError("depth_m must increase strictly downward")

    power = float(
        past_radiogenic_heat_uthk(
            LEGACY_MARS_MASS_KG,
            LEGACY_MARS_SILICATE_MASS_FRACTION,
            time_before_present_yr,
        )
    )
    heat_flux = power / (4.0 * np.pi * LEGACY_MARS_RADIUS_M**2)
    temperature = z * heat_flux / LEGACY_THERMAL_CONDUCTIVITY_WMK
    pressure = (
        LEGACY_MARS_ROCK_DENSITY_KGM3
        * LEGACY_MARS_GRAVITY_MS2
        * z
        / 1.0e6
    )
    return RadialColumn(
        depth_m=z,
        pressure_MPa=pressure,
        temperature_C=temperature,
        radiogenic_power_W=power,
        surface_heat_flux_Wm2=heat_flux,
    )


def legacy_mars_column(time_before_present_yr: float, dz_m: float = 1000.0) -> RadialColumn:
    """Reproduce the archived Mars P-T column at one time before present.

    Pressure is the original constant-rho/constant-g approximation,
    ``P = 3500*6*z``. Temperature is the steady conductive profile from the
    instantaneous long-lived radiogenic power, with zero-C surface temperature
    and constant 3 W m-1 K-1 conductivity.
    """
    if dz_m <= 0:
        raise ValueError("dz_m must be positive")
    z = np.arange(0.0, LEGACY_MARS_RADIUS_M + 0.5 * dz_m, dz_m)
    z = z[z <= LEGACY_MARS_RADIUS_M]
    return legacy_mars_pt_on_depths(z, time_before_present_yr)
