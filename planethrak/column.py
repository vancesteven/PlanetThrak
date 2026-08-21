"""Neutral radial-column interface for PlanetThrak fracture calculations.

The legacy MATLAB implementation owns its planetary structure internally. The
Python redesign does the opposite: pressure, temperature, depth, and optional
material properties are supplied by the caller. This keeps the fracture kernel
independent of PlanetProfile and lets the same code consume archived legacy
columns, PlanetProfile outputs, or other thermal-evolution models.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fronts import CrackingFront
from .intersection import CrackingIntersection, find_cracking_intersection


@dataclass(frozen=True)
class FractureColumn:
    """One radial material column supplied to the fracture calculation.

    ``depth_m`` is measured positive downward and pressure must increase with
    depth for the current legacy-front adapter. Optional material arrays are
    carried here so later fracture constitutive laws can use self-consistent
    PlanetProfile/Perple_X properties without changing the geometry API.
    """

    depth_m: np.ndarray
    pressure_MPa: np.ndarray
    temperature_C: np.ndarray
    cooling_rate_K_per_yr: np.ndarray | None = None
    pore_pressure_MPa: np.ndarray | None = None
    reactive_fraction: np.ndarray | None = None
    density_kg_m3: np.ndarray | None = None
    gravity_m_s2: np.ndarray | None = None
    thermal_expansivity_Kinv: np.ndarray | None = None
    thermal_conductivity_W_mK: np.ndarray | None = None
    porosity_fraction: np.ndarray | None = None
    bulk_modulus_Pa: np.ndarray | None = None
    shear_modulus_Pa: np.ndarray | None = None
    vp_m_s: np.ndarray | None = None
    vs_m_s: np.ndarray | None = None
    youngs_modulus_Pa: np.ndarray | None = None
    poisson_ratio: np.ndarray | None = None
    fracture_toughness_Pa_sqrt_m: np.ndarray | None = None
    grain_size_m: np.ndarray | None = None

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

        optional = (
            "cooling_rate_K_per_yr",
            "pore_pressure_MPa",
            "reactive_fraction",
            "density_kg_m3",
            "gravity_m_s2",
            "thermal_expansivity_Kinv",
            "thermal_conductivity_W_mK",
            "porosity_fraction",
            "bulk_modulus_Pa",
            "shear_modulus_Pa",
            "vp_m_s",
            "vs_m_s",
            "youngs_modulus_Pa",
            "poisson_ratio",
            "fracture_toughness_Pa_sqrt_m",
            "grain_size_m",
        )
        for name in optional:
            value = getattr(self, name)
            if value is None:
                continue
            arr = np.asarray(value, dtype=float)
            if arr.shape != z.shape or not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} must be finite and have the same shape as depth_m")
            if name in {"reactive_fraction", "porosity_fraction"} and (
                np.any(arr < 0) or np.any(arr > 1)
            ):
                raise ValueError(f"{name} must lie in [0, 1]")
            if name == "pore_pressure_MPa" and np.any(arr < 0):
                raise ValueError("pore_pressure_MPa must be non-negative")
            if name in {"gravity_m_s2", "shear_modulus_Pa", "vp_m_s", "vs_m_s"} and np.any(arr < 0):
                raise ValueError(f"{name} must be non-negative")
            if name in {
                "density_kg_m3",
                "thermal_conductivity_W_mK",
                "bulk_modulus_Pa",
                "youngs_modulus_Pa",
                "fracture_toughness_Pa_sqrt_m",
                "grain_size_m",
            } and np.any(arr <= 0):
                raise ValueError(f"{name} must be positive")
            if name == "poisson_ratio" and (np.any(arr <= -1) or np.any(arr >= 0.5)):
                raise ValueError("poisson_ratio must lie in (-1, 0.5)")
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

    Material shallower than the archived intersection is marked accessible.
    ``no_intersection`` is deliberately explicit because an absent crossing can
    mean either no accessible rock or an all-accessible column depending on the
    relative positions of the two P-T curves.
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
