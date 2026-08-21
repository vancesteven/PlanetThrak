"""Intersection of a planetary P-T column with a thermal-cracking boundary.

This module mirrors the numerical role of MATLAB ``get_Pz_cracking`` while
removing globals and plotting side effects.  It intentionally operates on
supplied P-T-depth arrays so the same kernel can consume either legacy radial
profiles or future PlanetProfile columns.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from .fronts import CrackingFront


@dataclass(frozen=True)
class CrackingIntersection:
    pressure_MPa: float
    temperature_C: float
    depth_m: float


def _ascending_xy(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return strictly increasing x with y reordered to match."""

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size:
        raise ValueError("x and y must be equal-length 1-D arrays")
    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    if np.any(np.diff(xs) <= 0):
        raise ValueError("independent-variable values must be unique")
    return xs, ys


def find_cracking_intersection(
    front: CrackingFront,
    profile_pressure_MPa: np.ndarray,
    profile_temperature_C: np.ndarray,
    profile_depth_m: np.ndarray,
) -> CrackingIntersection | None:
    """Find the shallowest P-T intersection with the cracking boundary.

    MATLAB's legacy implementation constructs cubic splines for T(P) for both
    the cracking front and planetary column, then uses ``fzero``.  Here we use
    SciPy's not-a-knot ``CubicSpline`` (the corresponding cubic-spline default)
    and bracketed Brent roots.  Multiple roots are handled deterministically by
    returning the shallowest intersection in the supplied monotonic pressure
    column.

    Returns ``None`` when the two curves have no intersection over their common
    pressure range.
    """

    p_prof, t_prof = _ascending_xy(profile_pressure_MPa, profile_temperature_C)
    p_depth, depth = _ascending_xy(profile_pressure_MPa, profile_depth_m)
    p_front, t_front = _ascending_xy(front.pressure_MPa, front.temperature_C)

    if not np.all(np.isfinite(t_prof)) or not np.all(np.isfinite(depth)):
        raise ValueError("profile temperature/depth arrays must be finite")

    lo = max(float(p_prof[0]), float(p_front[0]))
    hi = min(float(p_prof[-1]), float(p_front[-1]))
    if not lo < hi:
        return None

    front_spline = CubicSpline(p_front, t_front, extrapolate=False)
    profile_spline = CubicSpline(p_prof, t_prof, extrapolate=False)
    depth_spline = CubicSpline(p_depth, depth, extrapolate=False)

    def delta_t(p: float) -> float:
        return float(front_spline(p) - profile_spline(p))

    # Include all knot locations in the common interval and densify each
    # interval.  This makes root finding robust to an interior crossing that
    # does not coincide with a knot while preserving the legacy cubic curves.
    knots = np.unique(
        np.concatenate(
            [
                p_front[(p_front >= lo) & (p_front <= hi)],
                p_prof[(p_prof >= lo) & (p_prof <= hi)],
                np.asarray([lo, hi]),
            ]
        )
    )
    samples: list[float] = []
    for a, b in zip(knots[:-1], knots[1:]):
        samples.extend(np.linspace(a, b, 17, endpoint=False).tolist())
    samples.append(float(hi))
    p_test = np.asarray(samples, dtype=float)
    f_test = np.asarray([delta_t(p) for p in p_test])

    roots: list[float] = []
    atol = 1e-10
    for p, f in zip(p_test, f_test):
        if abs(f) <= atol:
            roots.append(float(p))
    for a, b, fa, fb in zip(p_test[:-1], p_test[1:], f_test[:-1], f_test[1:]):
        if fa * fb < 0:
            roots.append(float(brentq(delta_t, float(a), float(b))))

    if not roots:
        return None

    # Deduplicate numerically coincident roots and choose the shallowest one.
    roots = sorted(roots)
    unique_roots = [roots[0]]
    for root in roots[1:]:
        if not np.isclose(root, unique_roots[-1], rtol=1e-10, atol=1e-10):
            unique_roots.append(root)

    candidates = [
        CrackingIntersection(
            pressure_MPa=p,
            temperature_C=float(front_spline(p)),
            depth_m=float(depth_spline(p)),
        )
        for p in unique_roots
    ]
    return min(candidates, key=lambda item: item.depth_m)
