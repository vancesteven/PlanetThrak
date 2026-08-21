"""Neutral radial-column interface for PlanetThrak fracture calculations.

The legacy MATLAB implementation owns its planetary structure internally.  The
Python redesign does the opposite: pressure, temperature, and depth are supplied
by the caller.  This keeps the fracture kernel independent of PlanetProfile and
lets the same code consume archived legacy columns, PlanetProfile outputs, or
other thermal-evolution models.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fronts import CrackingFront
from .intersection import CrackingIntersection, find_cracking_intersection


@dataclass(frozen=True)
class FractureColumn:
    """One radial material column supplied to the fracture calculation.

    ``depth_m`` is measured positive downward.  Pressure must increase
    monotonically with depth for the current legacy-front adapter.  Additional
    material fields are intentionally optional because Phase 1 only requires
    P-T parity; later constitutive models can consume them without changing the
    geometry API.
    """

    depth_m: np.ndarray
    pressure_MPa: np.ndarray
    temperature_C: np.ndarray
    cooling_rate_K_per_yr: np.ndarray | None = None
    pore_pressure_MPa: np.ndarray | None = None
    reactive_fraction: np.ndarray | None = None

    def __post_init__(self) -> None:
        z = np.asarray(self.depth_m, dtype=float)
        p = np.asarray(self.pressure_MPa, dtype=float)
        t = np.asarray(self.temperature_C, dtype=float)
        if z.ndim != 1 or p.ndim != 1 or t.ndim != 1:
            raise ValueError("depth, pressure, and temperature must be 1-D arrays")
        if not (z.size == p.size == t.size) or z.size < 2:
            raise ValueError("depth, pressure, and temperature must have equal length >= 2")
        if not np.all(np.isfinite(z)) or not np.all(np.isfinite(p)) or not np.all(np.isfinite(t)):
            raise ValueError("depth, pressure, and temperature must be finite")
        if np.any(np.diff(z) <= 0):
            raise ValueError("depth_m must increase strictly downward")
        if np.any(np.diff(p) <= 0):
            raise ValueError("pressure_MPa must increase strictly with depth")

        object.__setattr__(self, "depth_m", z)
        object.__setattr__(self, "pressure_MPa", p)
        object.__setattr__(self, "temperature_C", t)

        for name in ("cooling_rate_K_per_yr", "pore_pressure_MPa", "reactive_fraction"):
            value = getattr(self, name)
            if value is None:
                continue
            arr = np.asarray(value, dtype=float)
            if arr.shape != z.shape or not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} must be finite and have the same shape as depth_m")
            if name == "reactive_fraction" and (np.any(arr < 0) or np.any(arr > 1)):
                raise ValueError("reactive_fraction must lie in [0, 1]")
            object.__setattr__(self, name, arr)


@dataclass(frozen=True)
class LegacyAccessibility:
    """Binary accessibility implied by one archived cracking-front boundary."""

    intersection: CrackingIntersection | None
    accessibility: np.ndarray


def legacy_front_accessibility(
    column: FractureColumn,
    front: CrackingFront,
    *,
    no_intersection: str = "none",
) -> LegacyAccessibility:
    """Map a legacy P-T cracking front onto a supplied column.

    The archived PlanetThrak calculation reports the depth where the planetary
    P-T path intersects the cracking boundary.  Its physical interpretation is
    that material shallower than that front has passed through the thermal-
    cracking regime.  This function makes that implied binary field explicit.

    ``no_intersection`` controls the intentionally ambiguous case where the
    curves do not cross over their common pressure range.  ``"none"`` returns
    an all-zero field.  ``"all"`` returns an all-one field.  The caller must
    choose the physically appropriate interpretation rather than having the
    library silently extrapolate the archived lookup table.
    """

    if no_intersection not in {"none", "all"}:
        raise ValueError("no_intersection must be 'none' or 'all'")

    hit = find_cracking_intersection(
        front,
        column.pressure_MPa,
        column.temperature_C,
        column.depth_m,
    )
    if hit is None:
        value = 1.0 if no_intersection == "all" else 0.0
        return LegacyAccessibility(None, np.full(column.depth_m.shape, value, dtype=float))

    accessibility = (column.depth_m <= hit.depth_m).astype(float)
    return LegacyAccessibility(hit, accessibility)
