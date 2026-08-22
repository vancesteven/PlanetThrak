import numpy as np

from planethrak.column import FractureColumn
from planethrak.fronts import load_front
from planethrak.legacy_structure import legacy_mars_pt_on_depths
from planethrak.staged_comparison import staged_mars_cracking_depths


def test_identical_candidate_has_zero_stage_shifts():
    z = np.arange(0.0, 60_001.0, 1000.0)
    legacy = legacy_mars_pt_on_depths(z, 0.0)
    candidate = FractureColumn(
        depth_m=z,
        pressure_MPa=legacy.pressure_MPa,
        temperature_C=legacy.temperature_C,
    )
    front = load_front("PcT1mm_1oCyr.ext")
    out = staged_mars_cracking_depths(candidate, front, time_before_present_yr=0.0)
    assert out.legacy_m is not None
    assert out.pressure_shift_m == 0.0
    assert out.thermal_shift_m == 0.0
    assert out.total_shift_m == 0.0
    assert out.interaction_m == 0.0


def test_pressure_and_temperature_axes_are_separated():
    z = np.arange(0.0, 60_001.0, 1000.0)
    legacy = legacy_mars_pt_on_depths(z, 0.0)
    candidate = FractureColumn(
        depth_m=z,
        pressure_MPa=legacy.pressure_MPa * 0.9,
        temperature_C=legacy.temperature_C * 1.1,
    )
    front = load_front("PcT1mm_1oCyr.ext")
    out = staged_mars_cracking_depths(candidate, front, time_before_present_yr=0.0)
    assert out.legacy_m is not None
    assert out.pressure_only_m is not None
    assert out.thermal_only_m is not None
    assert out.modern_pt_m is not None
    assert out.pressure_shift_m != 0.0
    assert out.thermal_shift_m != 0.0
    assert np.isfinite(out.interaction_m)
