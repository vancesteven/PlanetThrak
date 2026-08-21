#!/usr/bin/env python3
"""Compare a neutral radial profile against the archived Mars cracking result.

The input is a ``planethrak.profile_io`` NPZ artifact.  The same archived
PcT1mm_1oCyr front is applied to both the supplied profile and the historical
constant-rho/constant-g Mars column, isolating the effect of replacing the
planetary P-T structure before any new fracture constitutive law is introduced.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from planethrak.column import legacy_front_accessibility
from planethrak.fronts import load_front
from planethrak.legacy_structure import legacy_mars_column
from planethrak.planetprofile import fracture_column_from_arrays
from planethrak.profile_io import load_fracture_column_npz


ROOT = Path(__file__).resolve().parents[1]


def _depth_km(result) -> float:
    return 0.0 if result.intersection is None else result.intersection.depth_m / 1e3


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path, help="neutral radial profile .npz")
    parser.add_argument("--time-yr", type=float, default=4.5e9)
    parser.add_argument("--front", type=Path, default=ROOT / "PcT1mm_1oCyr.ext")
    args = parser.parse_args()

    modern, metadata = load_fracture_column_npz(args.profile)
    front = load_front(args.front)
    modern_hit = legacy_front_accessibility(modern, front)

    legacy = legacy_mars_column(args.time_yr)
    legacy_column = fracture_column_from_arrays(
        depth_m=legacy.depth_m,
        pressure_MPa=legacy.pressure_MPa,
        temperature_K=legacy.temperature_C + 273.15,
    )
    legacy_hit = legacy_front_accessibility(legacy_column, front)

    z_modern = _depth_km(modern_hit)
    z_legacy = _depth_km(legacy_hit)
    print("PlanetThrak radial-structure comparison using the archived cracking front")
    print(f"profile: {args.profile}")
    if metadata:
        print("metadata:")
        for key in sorted(metadata):
            print(f"  {key}: {metadata[key]}")
    print(f"time before present: {args.time_yr:.9g} yr")
    print(f"legacy cracking depth: {z_legacy:.9f} km")
    print(f"supplied-profile depth: {z_modern:.9f} km")
    print(f"delta supplied - legacy: {z_modern - z_legacy:+.9f} km")
    print()
    print("Interpretation guard rail:")
    print("  Both columns use the same archived PcT front. The depth difference")
    print("  therefore measures only the effect of changing the planetary P-T")
    print("  structure. It is not yet a generalized fracture-mechanics result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
