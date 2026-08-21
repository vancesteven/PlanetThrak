"""Portable radial-profile artifacts for PlanetProfile/PlanetThrak coupling.

The file format is intentionally a simple NumPy ``.npz`` container with named
SI-unit arrays. It is not a serialization of PlanetProfile internals. This
keeps archived fracture calculations reproducible when PlanetProfile evolves
and gives collaborators a small, inspectable handoff artifact that can also be
consumed by other forward models such as pyLOV3D.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .column import FractureColumn
from .planetprofile import fracture_column_from_arrays, fracture_column_from_planetprofile


SCHEMA_VERSION = 1

_OPTIONAL_FIELDS = (
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


def save_fracture_column_npz(path, column: FractureColumn, **metadata) -> Path:
    """Write a neutral radial column to ``.npz``.

    Required arrays are stored in explicit units. Optional material fields are
    included only when populated. Metadata values should be scalar strings or
    numbers and are stored under the ``meta_`` prefix.
    """
    out = Path(path)
    arrays = {
        "schema_version": np.array(SCHEMA_VERSION, dtype=np.int64),
        "depth_m": np.asarray(column.depth_m, dtype=float),
        "pressure_MPa": np.asarray(column.pressure_MPa, dtype=float),
        "temperature_K": np.asarray(column.temperature_C, dtype=float) + 273.15,
    }
    for name in _OPTIONAL_FIELDS:
        value = getattr(column, name)
        if value is not None:
            arrays[name] = np.asarray(value, dtype=float)
    for key, value in metadata.items():
        if isinstance(value, (str, int, float, np.integer, np.floating, bool)):
            arrays[f"meta_{key}"] = np.asarray(value)
        else:
            raise TypeError(f"metadata {key!r} must be a scalar string/number/bool")
    np.savez_compressed(out, **arrays)
    return out


def save_planetprofile_npz(path, planet, *, mask=None, reactive_fraction=None, **metadata) -> Path:
    """Duck-typed convenience wrapper for a completed PlanetProfile object.

    ``body_radius_m`` and ``body`` metadata are added automatically when the
    corresponding PlanetProfile attributes are available, unless the caller
    supplied explicit values.
    """
    column = fracture_column_from_planetprofile(
        planet,
        mask=mask,
        reactive_fraction=reactive_fraction,
    )
    bulk = getattr(planet, "Bulk", None)
    if "body_radius_m" not in metadata and bulk is not None:
        radius = getattr(bulk, "R_m", None)
        if radius is not None:
            metadata["body_radius_m"] = float(radius)
    if "body" not in metadata:
        name = getattr(planet, "name", None)
        if name is not None:
            metadata["body"] = str(name)
    return save_fracture_column_npz(path, column, **metadata)


def load_fracture_column_npz(path) -> tuple[FractureColumn, dict[str, object]]:
    """Load a neutral radial column and scalar metadata from ``.npz``."""
    with np.load(Path(path), allow_pickle=False) as data:
        version = int(np.asarray(data["schema_version"]).item())
        if version != SCHEMA_VERSION:
            raise ValueError(f"unsupported profile schema_version={version}; expected {SCHEMA_VERSION}")
        kwargs = {}
        for name in _OPTIONAL_FIELDS:
            if name in data:
                kwargs[name] = np.asarray(data[name], dtype=float)
        column = fracture_column_from_arrays(
            depth_m=np.asarray(data["depth_m"], dtype=float),
            pressure_MPa=np.asarray(data["pressure_MPa"], dtype=float),
            temperature_K=np.asarray(data["temperature_K"], dtype=float),
            **kwargs,
        )
        metadata = {}
        for key in data.files:
            if key.startswith("meta_"):
                value = np.asarray(data[key]).item()
                metadata[key[5:]] = value
    return column, metadata
