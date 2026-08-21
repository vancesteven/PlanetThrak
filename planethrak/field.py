"""Validated 3D fracture-accessibility fields for Task-1 Mars calculations.

The field object is intentionally a mechanical/compositional product. It knows
where rock is accessible and how much locally reactive lithology is present, but
it does not assume that water arrived or that reaction went to completion.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .capacity import (
    equiangular_latlon_area_weights,
    fractured_reactive_volume,
    maximum_bound_water_mass,
)


@dataclass(frozen=True)
class FractureField:
    """Full-sphere latitude-longitude-radius fracture/accessibility field.

    Arrays use shape ``(nlat, nlon, nz)``. ``radius_edges_m`` increases outward,
    so radial index zero is the deepest shell. Latitude and longitude values are
    cell centres in degrees. The current global integration path requires a
    full-sphere equiangular surface grid, matching the geometry expected by the
    downstream spherical-harmonic forward models.
    """

    latitude_deg: np.ndarray
    longitude_deg: np.ndarray
    radius_edges_m: np.ndarray
    accessibility: np.ndarray
    reactive_fraction: np.ndarray

    def __post_init__(self) -> None:
        lat = np.asarray(self.latitude_deg, dtype=float)
        lon = np.asarray(self.longitude_deg, dtype=float)
        r = np.asarray(self.radius_edges_m, dtype=float)
        a = np.asarray(self.accessibility, dtype=float)
        f = np.asarray(self.reactive_fraction, dtype=float)

        if lat.ndim != 1 or lon.ndim != 1 or r.ndim != 1:
            raise ValueError("latitude, longitude, and radius_edges must be 1-D")
        if lat.size < 2 or lon.size < 1 or r.size < 2:
            raise ValueError("field coordinate arrays are too short")
        expected = (lat.size, lon.size, r.size - 1)
        if a.shape != expected:
            raise ValueError(f"accessibility must have shape {expected}")
        try:
            f = np.broadcast_to(f, expected).astype(float, copy=True)
        except ValueError as exc:
            raise ValueError(f"reactive_fraction must broadcast to shape {expected}") from exc

        if not all(np.all(np.isfinite(x)) for x in (lat, lon, r, a, f)):
            raise ValueError("field coordinates and values must be finite")
        if np.any(np.diff(lat) <= 0) or (lon.size > 1 and np.any(np.diff(lon) <= 0)):
            raise ValueError("latitude and longitude cell centres must increase")
        if np.any(np.diff(r) <= 0) or np.any(r <= 0):
            raise ValueError("radius_edges_m must be positive and increase outward")
        if np.any(a < 0) or np.any(a > 1):
            raise ValueError("accessibility must lie in [0,1]")
        if np.any(f < 0) or np.any(f > 1):
            raise ValueError("reactive_fraction must lie in [0,1]")

        # Validate the latitude cells immediately instead of waiting until a
        # volume integral requests area weights.
        equiangular_latlon_area_weights(lat, lon.size)

        if lon.size > 1:
            dlon = np.diff(lon)
            if not np.allclose(dlon, dlon[0], rtol=0, atol=1e-10):
                raise ValueError("longitude_deg must be uniformly spaced")
            if not np.isclose(float(dlon[0]) * lon.size, 360.0, rtol=0, atol=1e-9):
                raise ValueError("longitude cells must span a full 360 degrees")

        object.__setattr__(self, "latitude_deg", lat)
        object.__setattr__(self, "longitude_deg", lon)
        object.__setattr__(self, "radius_edges_m", r)
        object.__setattr__(self, "accessibility", a)
        object.__setattr__(self, "reactive_fraction", f)

    @property
    def area_weights(self) -> np.ndarray:
        """Exact spherical cell-area fractions for the equiangular grid."""
        return equiangular_latlon_area_weights(self.latitude_deg, self.longitude_deg.size)

    @property
    def accessible_reactive_fraction(self) -> np.ndarray:
        """Cell-wise ``A_f * f_reactive`` entering the Task-1 volume integral."""
        return self.accessibility * self.reactive_fraction

    def fractured_reactive_volume_m3(self) -> float:
        """Return ``integral A_f f_reactive dV`` over the whole field."""
        return fractured_reactive_volume(
            self.accessibility,
            self.reactive_fraction,
            self.radius_edges_m,
            area_weights=self.area_weights,
        )

    def maximum_bound_water_mass_kg(
        self,
        rock_density_kg_m3: float,
        hydration_water_mass_fraction: float,
    ) -> float:
        """Return the full-hydration water ceiling for accessible reactive rock."""
        return maximum_bound_water_mass(
            self.fractured_reactive_volume_m3(),
            rock_density_kg_m3,
            hydration_water_mass_fraction,
        )
