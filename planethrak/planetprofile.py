"""Adapters from PlanetProfile profile arrays to :class:`FractureColumn`.

This module deliberately uses duck typing and NumPy only. PlanetThrak does not
import PlanetProfile, so the fracture code can remain independently testable.
Current PlanetProfile profile arrays are mapped explicitly, including seismic
moduli needed by both fracture mechanics and the shared pyLOV3D state.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .column import FractureColumn


def _optional_array(value: Any, n: int, name: str) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size != n:
        raise ValueError(f"{name} length {arr.size} does not match required length {n}")
    return arr


def fracture_column_from_arrays(
    *,
    depth_m,
    pressure_MPa,
    temperature_K,
    mask=None,
    cooling_rate_K_per_yr=None,
    pore_pressure_MPa=None,
    reactive_fraction=None,
    density_kg_m3=None,
    gravity_m_s2=None,
    thermal_expansivity_Kinv=None,
    thermal_conductivity_W_mK=None,
    porosity_fraction=None,
    bulk_modulus_Pa=None,
    shear_modulus_Pa=None,
    vp_m_s=None,
    vs_m_s=None,
    youngs_modulus_Pa=None,
    poisson_ratio=None,
    fracture_toughness_Pa_sqrt_m=None,
    grain_size_m=None,
) -> FractureColumn:
    """Build a fracture column from PlanetProfile-like arrays.

    Required arrays may be ordered either shallow-to-deep or deep-to-shallow;
    the output is sorted by increasing depth. ``temperature_K`` is converted to
    degrees C at the package boundary. ``mask`` can select, for example, only
    silicate/crustal nodes while excluding an iron core.
    """
    z = np.asarray(depth_m, dtype=float).reshape(-1)
    p = np.asarray(pressure_MPa, dtype=float).reshape(-1)
    tk = np.asarray(temperature_K, dtype=float).reshape(-1)
    if not (z.size == p.size == tk.size):
        raise ValueError("depth_m, pressure_MPa, and temperature_K must have equal length")
    n = z.size
    if n < 2:
        raise ValueError("at least two radial samples are required")

    optionals = {
        "cooling_rate_K_per_yr": _optional_array(cooling_rate_K_per_yr, n, "cooling_rate_K_per_yr"),
        "pore_pressure_MPa": _optional_array(pore_pressure_MPa, n, "pore_pressure_MPa"),
        "reactive_fraction": _optional_array(reactive_fraction, n, "reactive_fraction"),
        "density_kg_m3": _optional_array(density_kg_m3, n, "density_kg_m3"),
        "gravity_m_s2": _optional_array(gravity_m_s2, n, "gravity_m_s2"),
        "thermal_expansivity_Kinv": _optional_array(thermal_expansivity_Kinv, n, "thermal_expansivity_Kinv"),
        "thermal_conductivity_W_mK": _optional_array(thermal_conductivity_W_mK, n, "thermal_conductivity_W_mK"),
        "porosity_fraction": _optional_array(porosity_fraction, n, "porosity_fraction"),
        "bulk_modulus_Pa": _optional_array(bulk_modulus_Pa, n, "bulk_modulus_Pa"),
        "shear_modulus_Pa": _optional_array(shear_modulus_Pa, n, "shear_modulus_Pa"),
        "vp_m_s": _optional_array(vp_m_s, n, "vp_m_s"),
        "vs_m_s": _optional_array(vs_m_s, n, "vs_m_s"),
        "youngs_modulus_Pa": _optional_array(youngs_modulus_Pa, n, "youngs_modulus_Pa"),
        "poisson_ratio": _optional_array(poisson_ratio, n, "poisson_ratio"),
        "fracture_toughness_Pa_sqrt_m": _optional_array(fracture_toughness_Pa_sqrt_m, n, "fracture_toughness_Pa_sqrt_m"),
        "grain_size_m": _optional_array(grain_size_m, n, "grain_size_m"),
    }

    if mask is None:
        keep = np.ones(n, dtype=bool)
    else:
        keep = np.asarray(mask, dtype=bool).reshape(-1)
        if keep.size != n:
            raise ValueError("mask must have the same length as the profile arrays")
    if np.count_nonzero(keep) < 2:
        raise ValueError("mask must retain at least two radial samples")

    z = z[keep]
    p = p[keep]
    tk = tk[keep]
    optionals = {name: None if arr is None else arr[keep] for name, arr in optionals.items()}

    order = np.argsort(z)
    z = z[order]
    p = p[order]
    tk = tk[order]
    optionals = {name: None if arr is None else arr[order] for name, arr in optionals.items()}

    return FractureColumn(
        depth_m=z,
        pressure_MPa=p,
        temperature_C=tk - 273.15,
        **optionals,
    )


def _planetprofile_seismic_arrays(planet):
    seismic = getattr(planet, "Seismic", None)
    if seismic is None:
        return None, None, None, None
    ks = getattr(seismic, "KS_GPa", None)
    gs = getattr(seismic, "GS_GPa", None)
    vp = getattr(seismic, "VP_kms", None)
    vs = getattr(seismic, "VS_kms", None)
    return (
        None if ks is None else np.asarray(ks, dtype=float) * 1e9,
        None if gs is None else np.asarray(gs, dtype=float) * 1e9,
        None if vp is None else np.asarray(vp, dtype=float) * 1e3,
        None if vs is None else np.asarray(vs, dtype=float) * 1e3,
    )


def fracture_column_from_planetprofile(
    planet,
    *,
    mask=None,
    reactive_fraction=None,
    cooling_rate_K_per_yr=None,
    youngs_modulus_Pa=None,
    poisson_ratio=None,
    fracture_toughness_Pa_sqrt_m=None,
    grain_size_m=None,
) -> FractureColumn:
    """Duck-typed adapter for a completed PlanetProfile ``PlanetStruct``.

    The function intentionally does not infer which phase is reactive. Callers
    should supply ``mask`` and/or ``reactive_fraction`` from the composition
    model. If ``z_m`` is unavailable, depth is reconstructed from ``Bulk.R_m``
    and ``r_m``.

    When PlanetProfile seismic K and G are present they are carried directly in
    SI units. If E and nu are not supplied explicitly, they are derived from
    K and G via the isotropic relations. Callers should mask out fluid/core
    samples when those derived solid-rock quantities are required.
    """
    if getattr(planet, "P_MPa", None) is None or getattr(planet, "T_K", None) is None:
        raise ValueError("PlanetProfile object must contain populated P_MPa and T_K arrays")

    z = getattr(planet, "z_m", None)
    if z is None:
        r = getattr(planet, "r_m", None)
        bulk = getattr(planet, "Bulk", None)
        radius = None if bulk is None else getattr(bulk, "R_m", None)
        if r is None or radius is None:
            raise ValueError("PlanetProfile object must provide z_m or both r_m and Bulk.R_m")
        z = float(radius) - np.asarray(r, dtype=float)

    K, G, vp, vs = _planetprofile_seismic_arrays(planet)
    if youngs_modulus_Pa is None and poisson_ratio is None and K is not None and G is not None:
        denom = 3.0 * K + G
        with np.errstate(divide="ignore", invalid="ignore"):
            youngs_modulus_Pa = 9.0 * K * G / denom
            poisson_ratio = (3.0 * K - 2.0 * G) / (2.0 * denom)

    return fracture_column_from_arrays(
        depth_m=z,
        pressure_MPa=planet.P_MPa,
        temperature_K=planet.T_K,
        mask=mask,
        cooling_rate_K_per_yr=cooling_rate_K_per_yr,
        pore_pressure_MPa=getattr(planet, "Ppore_MPa", None),
        reactive_fraction=reactive_fraction,
        density_kg_m3=getattr(planet, "rho_kgm3", None),
        gravity_m_s2=getattr(planet, "g_ms2", None),
        thermal_expansivity_Kinv=getattr(planet, "alpha_pK", None),
        thermal_conductivity_W_mK=getattr(planet, "kTherm_WmK", None),
        porosity_fraction=getattr(planet, "phi_frac", None),
        bulk_modulus_Pa=K,
        shear_modulus_Pa=G,
        vp_m_s=vp,
        vs_m_s=vs,
        youngs_modulus_Pa=youngs_modulus_Pa,
        poisson_ratio=poisson_ratio,
        fracture_toughness_Pa_sqrt_m=fracture_toughness_Pa_sqrt_m,
        grain_size_m=grain_size_m,
    )
