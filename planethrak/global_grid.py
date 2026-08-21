"""Build global fracture-accessibility fields from gridded P-T columns.

This is the first direct bridge from a 3D thermal/interior model to the Task-1
mechanical accessibility field ``A_f(theta, phi, r, t)``. The routine remains
agnostic about where the P-T arrays came from; PlanetProfile is one intended
producer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .column import FractureColumn, legacy_front_accessibility
from .field import FractureField
from .fronts import CrackingFront


@dataclass(frozen=True)
class FractureGridResult:
    field: FractureField
    cracking_depth_m: np.ndarray


def fracture_field_from_pt_grids(
    *,
    latitude_deg,
    longitude_deg,
    radius_edges_m,
    pressure_MPa,
    temperature_K,
    reactive_fraction,
    front: CrackingFront,
    no_intersection: str = "none",
) -> FractureGridResult:
    """Evaluate an archived/general P-T cracking front column by column.

    Parameters
    ----------
    pressure_MPa, temperature_K
        Arrays with shape ``(nlat, nlon, nz)``. Radial index increases outward,
        matching ``radius_edges_m``. The function reverses each column before
        passing it to :class:`FractureColumn`, whose public convention is depth
        increasing downward.
    reactive_fraction
        Scalar or array broadcastable to the P-T grid. It is carried into the
        returned :class:`FractureField` but does not affect mechanical cracking.
    no_intersection
        Explicit legacy behavior, ``'none'`` or ``'all'``. No automatic physical
        interpretation is made when the curves do not cross.

    Returns
    -------
    FractureGridResult
        ``field.accessibility`` is in the same deep-to-shallow radial ordering as
        the input arrays. ``cracking_depth_m`` has shape ``(nlat,nlon)`` and is
        NaN when no intersection was found.
    """
    lat = np.asarray(latitude_deg, dtype=float)
    lon = np.asarray(longitude_deg, dtype=float)
    edges = np.asarray(radius_edges_m, dtype=float)
    p = np.asarray(pressure_MPa, dtype=float)
    tk = np.asarray(temperature_K, dtype=float)
    if lat.ndim != 1 or lon.ndim != 1 or edges.ndim != 1:
        raise ValueError("latitude, longitude, and radius_edges_m must be 1-D")
    if edges.size < 3 or np.any(np.diff(edges) <= 0):
        raise ValueError("radius_edges_m must increase outward and define at least two shells")
    expected = (lat.size, lon.size, edges.size - 1)
    if p.shape != expected or tk.shape != expected:
        raise ValueError(f"pressure_MPa and temperature_K must have shape {expected}")
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(tk)):
        raise ValueError("P-T grids must be finite")

    try:
        reactive = np.broadcast_to(np.asarray(reactive_fraction, dtype=float), expected).copy()
    except ValueError as exc:
        raise ValueError(f"reactive_fraction must broadcast to shape {expected}") from exc

    centers = 0.5 * (edges[:-1] + edges[1:])
    surface_radius = float(edges[-1])
    depth_deep_to_shallow = surface_radius - centers
    if np.any(depth_deep_to_shallow <= 0):
        raise ValueError("radial shell centers must lie below the reference surface")

    accessibility = np.zeros(expected, dtype=float)
    depth_map = np.full((lat.size, lon.size), np.nan, dtype=float)

    for i in range(lat.size):
        for j in range(lon.size):
            # FractureColumn expects shallow -> deep, opposite the radius order.
            column = FractureColumn(
                depth_m=depth_deep_to_shallow[::-1],
                pressure_MPa=p[i, j, ::-1],
                temperature_C=tk[i, j, ::-1] - 273.15,
                reactive_fraction=reactive[i, j, ::-1],
            )
            result = legacy_front_accessibility(
                column,
                front,
                no_intersection=no_intersection,
            )
            accessibility[i, j, :] = result.accessibility[::-1]
            if result.intersection is not None:
                depth_map[i, j] = result.intersection.depth_m

    field = FractureField(
        latitude_deg=lat,
        longitude_deg=lon,
        radius_edges_m=edges,
        accessibility=accessibility,
        reactive_fraction=reactive,
    )
    return FractureGridResult(field=field, cracking_depth_m=depth_map)
