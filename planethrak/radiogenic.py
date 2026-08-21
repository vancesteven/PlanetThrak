"""Legacy long-lived U-Th-K radiogenic heating used by PlanetThrak."""

from __future__ import annotations

import numpy as np


# Literal constants from get_pastRadiogenicHeat_UThK.m.
_PRESENT_CONCENTRATION_MASS_FRACTION = np.array([0.012, 0.040, 840.0]) * 1.0e-6
_HEAT_PER_CONCENTRATION_WKG = np.array([9.75, 2.60, 3.52e-4]) * 1.0e-5
_DECAY_CONSTANT_YR_INV = np.array([1.551, 0.495, 5.543]) * 1.0e-10


def past_radiogenic_heat_uthk(
    mass_kg: float,
    silicate_mass_fraction: float,
    time_before_present_yr: float | np.ndarray,
) -> np.ndarray:
    """Return legacy U-Th-K radiogenic power in watts.

    ``time_before_present_yr`` follows the archived MATLAB convention: zero is
    present day and positive values go backward in time, so radionuclide
    abundances scale as ``exp(lambda * t)``.
    """

    if mass_kg <= 0:
        raise ValueError("mass_kg must be positive")
    if not 0.0 <= silicate_mass_fraction <= 1.0:
        raise ValueError("silicate_mass_fraction must be between 0 and 1")

    t = np.asarray(time_before_present_yr, dtype=float)
    if np.any(~np.isfinite(t)) or np.any(t < 0):
        raise ValueError("time_before_present_yr must be finite and nonnegative")

    isotope_power = (
        _PRESENT_CONCENTRATION_MASS_FRACTION
        * _HEAT_PER_CONCENTRATION_WKG
        * np.exp(t[..., None] * _DECAY_CONSTANT_YR_INV)
    )
    return mass_kg * silicate_mass_fraction * np.sum(isotope_power, axis=-1)
