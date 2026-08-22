#!/usr/bin/env python3
"""Compare a neutral PlanetProfile/PlanetThrak radial artifact with legacy Mars.

The comparison is staged so pressure/gravity structure and thermal/geotherm
changes are reported separately before their combined P-T effect.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from planethrak.fronts import load_front
from planethrak.profile_io import load_fracture_column_npz
from planethrak.staged_comparison import format_stage_summary, staged_mars_cracking_depths


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path, help="neutral radial .npz artifact")
    parser.add_argument("--time-yr", type=float, default=0.0, help="time before present in years")
    parser.add_argument("--front", type=Path, default=ROOT / "PcT1mm_1oCyr.ext")
    args = parser.parse_args()

    column, metadata = load_fracture_column_npz(args.profile)
    front = load_front(args.front)
    stages = staged_mars_cracking_depths(
        column,
        front,
        time_before_present_yr=args.time_yr,
    )

    print("PlanetThrak staged Mars cracking comparison")
    print(f"profile: {args.profile}")
    print(f"front:   {args.front}")
    print(f"time:    {args.time_yr:.9g} yr before present")
    if metadata:
        for key in sorted(metadata):
            print(f"meta.{key}: {metadata[key]}")
    print()
    print(format_stage_summary(stages))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
