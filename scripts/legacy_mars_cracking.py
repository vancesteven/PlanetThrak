#!/usr/bin/env python3
"""Run the first end-to-end Python reproduction of the legacy Mars path.

This is a parity diagnostic, not yet a validation anchor.  It combines the
literal archived Mars P-T approximation with one or more shipped cracking-front
tables and prints the resulting intersection depth.  The next gate is to compare
these numbers directly with the archived MATLAB ``Planets.mat`` output or a fresh
MATLAB run.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from planethrak.fronts import load_front
from planethrak.intersection import find_cracking_intersection
from planethrak.legacy_structure import legacy_mars_column


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TABLES = (
    "PcT1mm_1oCyr.ext",
    "PcT10mm_1oCyr.ext",
    "PcTp1mm_1oCyr.ext",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--time-yr",
        type=float,
        default=0.0,
        help="time before present in years, matching the legacy MATLAB convention",
    )
    parser.add_argument(
        "--table",
        action="append",
        default=None,
        help="cracking-front table relative to repo root; repeat for multiple tables",
    )
    args = parser.parse_args()

    tables = args.table or list(DEFAULT_TABLES)
    column = legacy_mars_column(args.time_yr)
    print("PlanetThrak legacy Mars Python parity diagnostic")
    print(f"time before present: {args.time_yr:.6g} yr")
    print(f"radiogenic power: {column.radiogenic_power_W:.12e} W")
    print(f"surface heat flux: {column.surface_heat_flux_Wm2:.12e} W/m^2")

    for relpath in tables:
        path = ROOT / relpath
        front = load_front(path)
        hit = find_cracking_intersection(
            front,
            column.pressure_MPa,
            column.temperature_C,
            column.depth_m,
        )
        if hit is None:
            print(f"{relpath}: no P-T intersection")
        else:
            print(
                f"{relpath}: depth={hit.depth_m/1e3:.9f} km  "
                f"P={hit.pressure_MPa:.9f} MPa  T={hit.temperature_C:.9f} C"
            )

    print("guard rail: compare against archived/fresh MATLAB before calling this parity")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
