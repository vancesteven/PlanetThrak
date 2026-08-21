"""Regression anchor for the archived MATLAB Mars thermal-cracking series.

This test intentionally validates only the historical 1-mm, 1 C/yr legacy
calculation stored in Planets.mat.  It does not validate the physical fidelity
of the archived constant-density/constant-gravity Mars structure.
"""

from pathlib import Path

import numpy as np
from scipy.io import loadmat

from planethrak.fronts import load_front
from planethrak.intersection import find_cracking_intersection
from planethrak.legacy_structure import legacy_mars_column


ROOT = Path(__file__).resolve().parents[1]
TOLERANCE_KM = 1.0e-9  # 1 micrometre; observed 2026-08-21 max = 1.45519152e-14 km


def _records(value):
    if isinstance(value, dict):
        if "name" in value:
            yield value
        else:
            for item in value.values():
                yield from _records(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _records(item)
    elif isinstance(value, np.ndarray):
        for item in value.flat:
            yield from _records(item)


def _find_planet(planets, name: str) -> dict:
    wanted = name.strip().lower()
    for record in _records(planets):
        if str(record.get("name", "")).strip().lower() == wanted:
            return record
    raise KeyError(f"planet {name!r} not found in archived structure")


def _one_dimensional_mars_depth(value, n_time: int) -> np.ndarray:
    z = np.squeeze(np.asarray(value, dtype=float))
    if z.ndim == 1 and z.size == n_time:
        return z
    if z.ndim == 2:
        if z.shape == (1, n_time):
            return z[0].ravel()
        if z.shape == (n_time, 1):
            return z[:, 0].ravel()
        candidates = [np.asarray(row, dtype=float).ravel() for row in z if np.asarray(row).size == n_time]
        if len(candidates) == 1:
            return candidates[0]
    raise ValueError(f"cannot reduce archived Mars cracking-depth shape {z.shape} to {n_time} times")


def test_mars_1mm_1cyr_archive_depth_series_matches_to_micrometre():
    data = loadmat(ROOT / "Planets.mat", simplify_cells=True)
    mars = _find_planet(data["Planets"], "Mars")
    times = np.asarray(mars["t_yr"], dtype=float).ravel()
    archived_m = _one_dimensional_mars_depth(mars["z_cracking_1mm_m"], times.size)
    front = load_front(ROOT / "PcT1mm_1oCyr.ext")

    predicted_m = np.empty_like(archived_m)
    for i, time_yr in enumerate(times):
        column = legacy_mars_column(float(time_yr))
        hit = find_cracking_intersection(
            front,
            column.pressure_MPa,
            column.temperature_C,
            column.depth_m,
        )
        predicted_m[i] = 0.0 if hit is None else hit.depth_m

    error_km = (predicted_m - archived_m) / 1.0e3
    assert times.size == 46
    assert np.all(np.isfinite(error_km))
    assert np.max(np.abs(error_km)) <= TOLERANCE_KM
