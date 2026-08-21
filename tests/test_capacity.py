import numpy as np
import pytest

from planethrak.capacity import (
    equiangular_latlon_area_weights,
    fractured_reactive_volume,
    maximum_bound_water_mass,
    spherical_shell_volumes,
)


def test_full_accessibility_recovers_exact_shell_volume():
    r = np.array([3.0, 4.0, 5.0])
    shell = spherical_shell_volumes(r)
    got = fractured_reactive_volume(np.ones(2), np.ones(2), r)
    np.testing.assert_allclose(got, np.sum(shell), rtol=0, atol=1e-12)


def test_surface_weights_and_reactive_fraction_are_integrated_separately():
    r = np.array([10.0, 11.0])
    # Two surface columns. First is fully accessible/reactive, second is not.
    a = np.array([[1.0], [0.0]])
    f = np.ones_like(a)
    got = fractured_reactive_volume(a, f, r, area_weights=np.array([0.25, 0.75]))
    expected = 0.25 * spherical_shell_volumes(r)[0]
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)


def test_equiangular_weights_cover_sphere_and_downweight_poles():
    # Native lmax=2 latitude centres: -67.5, -22.5, 22.5, 67.5 deg.
    lat = np.array([-67.5, -22.5, 22.5, 67.5])
    w = equiangular_latlon_area_weights(lat, 8)
    assert w.shape == (4, 8)
    np.testing.assert_allclose(np.sum(w), 1.0, rtol=0, atol=1e-15)
    assert w[0, 0] < w[1, 0]
    assert w[-1, 0] < w[-2, 0]
    np.testing.assert_allclose(w[0], w[-1], rtol=0, atol=1e-15)


def test_latlon_weighted_global_integral_recovers_full_shell_volume():
    lat = np.array([-67.5, -22.5, 22.5, 67.5])
    w = equiangular_latlon_area_weights(lat, 8)
    r = np.array([100.0, 110.0])
    a = np.ones((4, 8, 1))
    f = np.ones_like(a)
    got = fractured_reactive_volume(a, f, r, area_weights=w)
    expected = spherical_shell_volumes(r)[0]
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)


def test_multidimensional_integral_requires_explicit_area_weights():
    r = np.array([100.0, 110.0])
    with pytest.raises(ValueError, match="area_weights"):
        fractured_reactive_volume(np.ones((2, 1)), np.ones((2, 1)), r)


def test_reactive_fraction_scales_capacity_without_implying_hydration():
    r = np.array([100.0, 110.0])
    full = fractured_reactive_volume(np.array([1.0]), np.array([1.0]), r)
    half = fractured_reactive_volume(np.array([1.0]), np.array([0.5]), r)
    np.testing.assert_allclose(half, 0.5 * full, rtol=0, atol=1e-12)


def test_water_mass_bound_is_simple_full_hydration_ceiling():
    got = maximum_bound_water_mass(2.0, 3000.0, 0.13)
    assert got == 780.0


def test_invalid_fraction_fields_are_rejected():
    r = np.array([1.0, 2.0])
    with pytest.raises(ValueError, match="accessibility"):
        fractured_reactive_volume(np.array([1.1]), np.array([1.0]), r)
    with pytest.raises(ValueError, match="reactive_fraction"):
        fractured_reactive_volume(np.array([1.0]), np.array([-0.1]), r)
