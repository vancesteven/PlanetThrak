import numpy as np
import pytest

from planethrak.fronts import CrackingFront
from planethrak.intersection import find_cracking_intersection


def test_linear_equivalent_pt_curves_intersect_at_expected_depth():
    # Front: T = 300 - P over the interval, with P in MPa.
    front = CrackingFront(
        pressure_MPa=np.array([100.0, 200.0, 300.0]),
        temperature_C=np.array([200.0, 100.0, 0.0]),
    )
    # Planetary column: T = P - 100, depth = 100 m/MPa * P.
    p = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
    t = p - 100.0
    z = 100.0 * p
    hit = find_cracking_intersection(front, p, t, z)
    assert hit is not None
    np.testing.assert_allclose(hit.pressure_MPa, 200.0, rtol=0, atol=1e-9)
    np.testing.assert_allclose(hit.temperature_C, 100.0, rtol=0, atol=1e-9)
    np.testing.assert_allclose(hit.depth_m, 20_000.0, rtol=0, atol=1e-6)


def test_no_overlap_returns_none():
    front = CrackingFront(
        pressure_MPa=np.array([100.0, 200.0]),
        temperature_C=np.array([0.0, 100.0]),
    )
    p = np.array([0.0, 50.0])
    t = np.array([0.0, 10.0])
    z = np.array([0.0, 1000.0])
    assert find_cracking_intersection(front, p, t, z) is None


def test_decreasing_front_pressure_is_supported_like_legacy_tables():
    front = CrackingFront(
        pressure_MPa=np.array([300.0, 200.0, 100.0]),
        temperature_C=np.array([0.0, 100.0, 200.0]),
    )
    p = np.array([0.0, 100.0, 200.0, 300.0])
    t = p - 100.0
    z = 50.0 * p
    hit = find_cracking_intersection(front, p, t, z)
    assert hit is not None
    np.testing.assert_allclose(hit.pressure_MPa, 200.0, rtol=0, atol=1e-9)


def test_singleton_archive_fixture_is_not_silently_interpolated():
    front = CrackingFront(
        pressure_MPa=np.array([3.5723]),
        temperature_C=np.array([-90.0]),
        source="PcTp1mm_1oCGyr.ext",
    )
    p = np.array([0.0, 10.0])
    t = np.array([0.0, 100.0])
    z = np.array([0.0, 1000.0])
    with pytest.raises(ValueError, match="at least two P-T points"):
        find_cracking_intersection(front, p, t, z)
