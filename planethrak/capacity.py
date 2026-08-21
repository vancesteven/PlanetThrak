"""Global fractured-reactive-rock volume and water-capacity diagnostics.

These functions implement the Task-1 quantity

    V_accessible,reactive = integral A_f * f_reactive dV

without making any assumption that mechanically accessible rock is actually
hydrated.  Flow and reaction efficiency belong downstream.
"""

from __future__ import annotations

import numpy as np


def spherical_shell_volumes(radius_edges_m: np.ndarray) -> np.ndarray:
    """Return exact spherical-shell cell volumes between successive radii."""
    r = np.asarray(radius_edges_m, dtype=float)
    if r.ndim != 1 or r.size < 2:
        raise ValueError("radius_edges_m must be a 1-D array with at least two entries")
    if not np.all(np.isfinite(r)) or np.any(r <= 0):
        raise ValueError("radius edges must be finite and positive")
    if np.any(np.diff(r) <= 0):
        raise ValueError("radius_edges_m must increase strictly outward")
    return (4.0 * np.pi / 3.0) * (r[1:] ** 3 - r[:-1] ** 3)


def fractured_reactive_volume(
    accessibility: np.ndarray,
    reactive_fraction: np.ndarray,
    radius_edges_m: np.ndarray,
    *,
    area_weights: np.ndarray | None = None,
) -> float:
    """Integrate accessible reactive-rock volume over a spherical shell grid.

    Parameters
    ----------
    accessibility
        Dimensionless fracture accessibility in [0,1]. Shape is ``(..., nz)``.
        A binary legacy field and a future continuous susceptibility/probability
        field use the same interface.
    reactive_fraction
        Fraction of each cell occupied by compositionally reactive lithology,
        broadcastable to ``accessibility`` and constrained to [0,1].
    radius_edges_m
        Radial cell edges, increasing outward, with length ``nz+1``.
    area_weights
        Optional non-negative surface-area fractions for the leading dimensions
        of ``accessibility``. They are normalized internally to sum to one.
        If omitted, a one-dimensional radial column is interpreted as global.

    Notes
    -----
    The result is a geometric/material upper bound. It does not include water
    delivery, reaction kinetics, incomplete alteration, or crack sealing.
    """

    a = np.asarray(accessibility, dtype=float)
    f = np.asarray(reactive_fraction, dtype=float)
    if a.ndim < 1:
        raise ValueError("accessibility must have at least one dimension")
    if not np.all(np.isfinite(a)) or np.any(a < 0) or np.any(a > 1):
        raise ValueError("accessibility must be finite and lie in [0,1]")
    if not np.all(np.isfinite(f)) or np.any(f < 0) or np.any(f > 1):
        raise ValueError("reactive_fraction must be finite and lie in [0,1]")

    try:
        af = a * f
    except ValueError as exc:
        raise ValueError("reactive_fraction must be broadcastable to accessibility") from exc

    shell_v = spherical_shell_volumes(radius_edges_m)
    if af.shape[-1] != shell_v.size:
        raise ValueError("last accessibility dimension must match radial shell count")

    if af.ndim == 1:
        return float(np.sum(af * shell_v))

    surface_shape = af.shape[:-1]
    if area_weights is None:
        w = np.full(surface_shape, 1.0 / np.prod(surface_shape), dtype=float)
    else:
        w = np.asarray(area_weights, dtype=float)
        if w.shape != surface_shape:
            raise ValueError("area_weights must match accessibility leading dimensions")
        if not np.all(np.isfinite(w)) or np.any(w < 0) or not np.any(w > 0):
            raise ValueError("area_weights must be finite, non-negative, and not all zero")
        w = w / np.sum(w)

    column_fraction = np.sum(af * shell_v, axis=-1)
    return float(np.sum(w * column_fraction))


def maximum_bound_water_mass(
    fractured_reactive_volume_m3: float,
    rock_density_kg_m3: float,
    hydration_water_mass_fraction: float,
) -> float:
    """Convert accessible reactive-rock volume to a full-hydration water bound."""
    if fractured_reactive_volume_m3 < 0:
        raise ValueError("fractured_reactive_volume_m3 must be non-negative")
    if rock_density_kg_m3 <= 0:
        raise ValueError("rock_density_kg_m3 must be positive")
    if not 0.0 <= hydration_water_mass_fraction <= 1.0:
        raise ValueError("hydration_water_mass_fraction must lie in [0,1]")
    return float(fractured_reactive_volume_m3 * rock_density_kg_m3 * hydration_water_mass_fraction)
