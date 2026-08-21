#!/usr/bin/env python3
"""Extract archived Mars cracking fields from PlanetThrak MAT files.

This is the bridge from the historical MATLAB archive to strict Python parity.
It reads ``Planets.mat`` with SciPy's ``simplify_cells=True`` and prints the
stored Mars time/cracking-depth arrays in a machine-readable way.  No model is
run and no archived data are modified.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]


def _records(value):
    if isinstance(value, dict):
        # A single MATLAB struct or a mapping containing struct elements.
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


def find_planet(planets, name: str) -> dict:
    wanted = name.strip().lower()
    for record in _records(planets):
        if str(record.get("name", "")).strip().lower() == wanted:
            return record
    raise KeyError(f"planet {name!r} not found in archived structure")


def _summary(record: dict, key: str) -> None:
    value = record.get(key)
    if value is None:
        print(f"{key}: MISSING")
        return
    arr = np.asarray(value)
    if arr.ndim == 0:
        print(f"{key}: {arr.item()!r}")
        return
    flat = np.asarray(arr, dtype=float).ravel()
    finite = flat[np.isfinite(flat)]
    if finite.size:
        print(
            f"{key}: shape={arr.shape} n={flat.size} "
            f"first={flat[0]:.12g} last={flat[-1]:.12g} "
            f"min={finite.min():.12g} max={finite.max():.12g}"
        )
    else:
        print(f"{key}: shape={arr.shape} n={flat.size} all_nonfinite")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mat", type=Path, default=ROOT / "Planets.mat")
    parser.add_argument("--planet", default="Mars")
    parser.add_argument(
        "--dump",
        action="store_true",
        help="also print paired t_yr,z_cracking_1mm_m rows for direct regression",
    )
    args = parser.parse_args()

    data = loadmat(args.mat, simplify_cells=True)
    if "Planets" not in data:
        raise KeyError(f"Planets variable not found in {args.mat}")
    record = find_planet(data["Planets"], args.planet)

    print(f"archive: {args.mat}")
    print(f"planet: {record.get('name', args.planet)}")
    for key in (
        "t_yr",
        "z_cracking_1mm_m",
        "z_cracking_10mm_m",
        "P_z_1mm_MPa",
        "P_z_10mm_MPa",
        "P_MPa",
        "T_oC",
    ):
        _summary(record, key)

    if args.dump:
        t = np.asarray(record["t_yr"], dtype=float).ravel()
        z = np.asarray(record["z_cracking_1mm_m"], dtype=float)
        if z.ndim > 1:
            # Legacy files may hold multiple ocean-depth rows. Mars normally
            # has a single row; preserve the first for a deterministic dump.
            z = z[0]
        z = z.ravel()
        if t.size != z.size:
            raise ValueError(f"t_yr length {t.size} != z_cracking_1mm_m length {z.size}")
        print("# time_before_present_yr,z_cracking_1mm_m")
        for ti, zi in zip(t, z):
            print(f"{ti:.17g},{zi:.17g}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
