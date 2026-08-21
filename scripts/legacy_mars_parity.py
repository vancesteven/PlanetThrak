#!/usr/bin/env python3
"""Compare Python legacy Mars cracking depths directly with ``Planets.mat``.

This diagnostic closes the loop between the new Python kernel and the archived
MATLAB product. It intentionally reproduces the historical 1-mm, 1 C/yr lookup
path before any modern PlanetProfile structure is introduced.

Use ``--assert-km`` only after inspecting the first comparison and choosing a
scientifically justified numerical tolerance.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from planethrak.fronts import load_front
from planethrak.intersection import find_cracking_intersection
from planethrak.legacy_structure import legacy_mars_column


ROOT = Path(__file__).resolve().parents[1]


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
    z = np.asarray(value, dtype=float)
    z = np.squeeze(z)
    if z.ndim == 1 and z.size == n_time:
        return z
    if z.ndim == 2:
        candidates = [row for row in z if np.asarray(row).size == n_time]
        if len(candidates) == 1:
            return np.asarray(candidates[0], dtype=float).ravel()
        if z.shape[0] == 1 and z.shape[1] == n_time:
            return z[0].ravel()
        if z.shape[1] == 1 and z.shape[0] == n_time:
            return z[:, 0].ravel()
    raise ValueError(f"cannot reduce archived Mars cracking-depth shape {z.shape} to {n_time} times")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mat", type=Path, default=ROOT / "Planets.mat")
    parser.add_argument("--table", type=Path, default=ROOT / "PcT1mm_1oCyr.ext")
    parser.add_argument(
        "--assert-km",
        type=float,
        default=None,
        help="fail if maximum absolute depth error exceeds this tolerance",
    )
    parser.add_argument(
        "--show-all",
        action="store_true",
        help="print every archived/Python depth pair rather than summary plus worst rows",
    )
    args = parser.parse_args()

    data = loadmat(args.mat, simplify_cells=True)
    if "Planets" not in data:
        raise KeyError(f"Planets variable not found in {args.mat}")
    mars = _find_planet(data["Planets"], "Mars")
    times = np.asarray(mars["t_yr"], dtype=float).ravel()
    archived = _one_dimensional_mars_depth(mars["z_cracking_1mm_m"], times.size)
    front = load_front(args.table)

    predicted = np.zeros(times.size, dtype=float)
    for i, time_yr in enumerate(times):
        column = legacy_mars_column(float(time_yr))
        hit = find_cracking_intersection(
            front,
            column.pressure_MPa,
            column.temperature_C,
            column.depth_m,
        )
        # MATLAB get_Pz_cracking catches a failed intersection and returns 0.
        predicted[i] = 0.0 if hit is None else hit.depth_m

    error_km = (predicted - archived) / 1.0e3
    abs_error_km = np.abs(error_km)
    finite = np.isfinite(abs_error_km)
    if not np.all(finite):
        raise ValueError("non-finite archived or predicted cracking depths encountered")

    print("PlanetThrak legacy Mars parity: Python vs archived MATLAB")
    print(f"archive: {args.mat}")
    print(f"front:   {args.table}")
    print(f"samples: {times.size}")
    print(f"max |depth error|: {abs_error_km.max():.9g} km")
    print(f"rms depth error:   {np.sqrt(np.mean(error_km**2)):.9g} km")
    print(f"mean depth bias:   {np.mean(error_km):+.9g} km")

    order = np.argsort(abs_error_km)[::-1]
    rows = np.arange(times.size) if args.show_all else order[: min(8, times.size)]
    print("# time_yr, archived_km, python_km, error_km")
    for i in rows:
        print(
            f"{times[i]:.17g}, {archived[i]/1e3:.12g}, "
            f"{predicted[i]/1e3:.12g}, {error_km[i]:+.12g}"
        )

    if args.assert_km is not None:
        if args.assert_km < 0:
            raise ValueError("--assert-km must be non-negative")
        if abs_error_km.max() > args.assert_km:
            raise SystemExit(
                f"FAIL: max depth error {abs_error_km.max():.9g} km "
                f"> tolerance {args.assert_km:.9g} km"
            )
        print(f"PASS: max depth error <= {args.assert_km:g} km")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
