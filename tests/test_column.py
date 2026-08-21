import numpy as np
import pytest

from planethrak.column import FractureColumn, legacy_front_accessibility
from planethrak.fronts import CrackingFront


def _example_column():
    p = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
    return FractureColumn(
        depth_m=100.0 * p,
        pressure_MPa=p,
        temperature_C=p - 100.0,
        reactive_fraction=np.linspace(0.1, 0.5, p.size),
    )


def test_legacy_accessibility_is_shallow_side_of_intersection():
    front = CrackingFront(
        pressure_MPa=np.array([100.0, 200.0, 300.0]),
        temperature_C=np.array([200.0, 100.0, 0.0]),
    )
    result = legacy_front_accessibility(_example_column(), front)
    assert result.intersection is not None
    np.testing.assert_allclose(result.intersection.depth_m, 20_000.0, atol=1e-6)
    np.testing.assert_array_equal(result.accessibility, [1.0, 1.0, 1.0, 0.0, 0.0])


def test_no_intersection_requires_explicit_interpretation():
    column = _example_column()
    front = CrackingFront(
        pressure_MPa=np.array([500.0, 600.0]),
        temperature_C=np.array([0.0, 100.0]),
    )
    np.testing.assert_array_equal(
        legacy_front_accessibility(column, front, no_intersection="none").accessibility,
        np.zeros(column.depth_m.size),
    )
    np.testing.assert_array_equal(
        legacy_front_accessibility(column, front, no_intersection="all").accessibility,
        np.ones(column.depth_m.size),
    )


def test_column_rejects_nonmonotonic_pressure_and_invalid_reactive_fraction():
    with pytest.raises(ValueError, match="pressure_MPa"):
        FractureColumn(
            depth_m=np.array([0.0, 1.0, 2.0]),
            pressure_MPa=np.array([0.0, 2.0, 1.0]),
            temperature_C=np.array([0.0, 1.0, 2.0]),
        )
    with pytest.raises(ValueError, match="reactive_fraction"):
        FractureColumn(
            depth_m=np.array([0.0, 1.0, 2.0]),
            pressure_MPa=np.array([0.0, 1.0, 2.0]),
            temperature_C=np.array([0.0, 1.0, 2.0]),
            reactive_fraction=np.array([0.0, 1.1, 0.0]),
        )
