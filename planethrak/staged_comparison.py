"""Controlled comparisons between legacy and modern Mars radial columns.

The key scientific purpose is to avoid one opaque legacy-vs-modern difference.
Given a candidate PlanetProfile-derived column, this module evaluates four
states on the *same depth grid*:

1. legacy pressure + legacy temperature,
2. candidate pressure + legacy temperature,
3. legacy pressure + candidate temperature,
4. candidate pressure + candidate temperature.

The first two isolate pressure/gravity structure, the first/third isolate the
thermal profile, and the fourth contains their combined effect. Composition-
dependent fracture constitutive properties remain a later, separately testable
axis.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .column import FractureColumn
from .fronts import CrackingFront
from .intersection import find_cracking_intersection
from .legacy_structure import legacy_mars_pt_on_depths


@dataclass(frozen=True)
class CrackingDepthStages:
    legacy_m: float | None
    pressure_only_m: float | None
    thermal_only_m: float | None
    modern_pt_m: float | None

    @staticmethod
    def _shift(value: float | None, reference: float | None) -> float | None:
        if value is None or reference is None:
            return None
        return value - reference

    @property
    def pressure_shift_m(self) -> float | None:
        return self._shift(self.pressure_only_m, self.legacy_m)

    @property
    def thermal_shift_m(self) -> float | None:
        return self._shift(self.thermal_only_m, self.legacy_m)

    @property
    def total_shift_m(self) -> float | None:
        return self._shift(self.modern_pt_m, self.legacy_m)

    @property
    def interaction_m(self) -> float | None:
        """Non-additive P-T interaction in cracking depth.

        Defined as total - pressure-only - thermal-only, all relative to the
        common legacy baseline.
        """
        if (
            self.total_shift_m is None
            or self.pressure_shift_m is None
            or self.thermal_shift_m is None
        ):
            return None
        return self.total_shift_m - self.pressure_shift_m - self.thermal_shift_m


def _depth(front: CrackingFront, column: FractureColumn) -> float | None:
    hit = find_cracking_intersection(
        front,
        column.pressure_MPa,
        column.temperature_C,
        column.depth_m,
    )
    return None if hit is None else float(hit.depth_m)


def staged_mars_cracking_depths(
    candidate: FractureColumn,
    front: CrackingFront,
    *,
    time_before_present_yr: float,
) -> CrackingDepthStages:
    """Decompose legacy-to-candidate cracking-depth changes into P and T axes."""
    legacy = legacy_mars_pt_on_depths(candidate.depth_m, time_before_present_yr)

    common = dict(depth_m=candidate.depth_m)
    legacy_column = FractureColumn(
        pressure_MPa=legacy.pressure_MPa,
        temperature_C=legacy.temperature_C,
        **common,
    )
    pressure_only = FractureColumn(
        pressure_MPa=candidate.pressure_MPa,
        temperature_C=legacy.temperature_C,
        **common,
    )
    thermal_only = FractureColumn(
        pressure_MPa=legacy.pressure_MPa,
        temperature_C=candidate.temperature_C,
        **common,
    )

    return CrackingDepthStages(
        legacy_m=_depth(front, legacy_column),
        pressure_only_m=_depth(front, pressure_only),
        thermal_only_m=_depth(front, thermal_only),
        modern_pt_m=_depth(front, candidate),
    )


def format_stage_summary(stages: CrackingDepthStages) -> str:
    """Return a compact human-readable stage summary in kilometres."""
    def km(value):
        return "none" if value is None else f"{value / 1e3:.6f} km"

    return "\n".join(
        [
            f"legacy:        {km(stages.legacy_m)}",
            f"pressure-only: {km(stages.pressure_only_m)}  shift={km(stages.pressure_shift_m)}",
            f"thermal-only:  {km(stages.thermal_only_m)}  shift={km(stages.thermal_shift_m)}",
            f"modern P-T:    {km(stages.modern_pt_m)}  shift={km(stages.total_shift_m)}",
            f"P-T interaction: {km(stages.interaction_m)}",
        ]
    )
