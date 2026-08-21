import numpy as np
import pytest

from planethrak.capacity import spherical_shell_volumes
from planethrak.field import FractureField


def _grid():
    lat = np.array([-67.5, -22.5, 22.5, 67.5])
    lon = np.arange(-157.5, 180.0, 45.0)
    radius_edges = np.array([3_300_000.0, 3_320_000.0, 3_340_000.0])
    return lat, lon, radius_edges


def test_full_accessible_reactive_field_recovers_exact_shell_volume():
    lat, lon, r = _grid()
    shape = (lat.size, lon.size, r.size - 1)
    field = FractureField(
        latitude_deg=lat,
        longitude_deg=lon,
        radius_edges_m=r,
        accessibility=np.ones(shape),
        reactive_fraction=np.ones(shape),
    )
    expected = np.sum(spherical_shell_volumes(r))
    np.testing.assert_allclose(field.fractured_reactive_volume_m3(), expected, rtol=2e-15)
    np.testing.assert_allclose(np.sum(field.area_weights), 1.0, rtol=0, atol=1e-15)


def test_reactive_fraction_broadcasts_over_full_field():
    lat, lon, r = _grid()
    shape = (lat.size, lon.size, r.size - 1)
    field = FractureField(
        latitude_deg=lat,
        longitude_deg=lon,
        radius_edges_m=r,
        accessibility=np.ones(shape),
        reactive_fraction=np.array([0.25, 0.75]),
    )
    shell = spherical_shell_volumes(r)
    expected = 0.25 * shell[0] + 0.75 * shell[1]
    np.testing.assert_allclose(field.fractured_reactive_volume_m3(), expected, rtol=2e-15)
    assert field.reactive_fraction.shape == shape


def test_water_mass_method_is_full_hydration_ceiling_only():
    lat, lon, r = _grid()
    shape = (lat.size, lon.size, r.size - 1)
    field = FractureField(lat, lon, r, np.ones(shape), np.full(shape, 0.5))
    volume = field.fractured_reactive_volume_m3()
    got = field.maximum_bound_water_mass_kg(3000.0, 0.13)
    assert got == pytest.approx(volume * 3000.0 * 0.13, rel=1e-15)


def test_field_rejects_invalid_shape_and_fraction():
    lat, lon, r = _grid()
    with pytest.raises(ValueError, match="accessibility must have shape"):
        FractureField(lat, lon, r, np.ones((lat.size, lon.size, 1)), 1.0)
    shape = (lat.size, lon.size, r.size - 1)
    with pytest.raises(ValueError, match="accessibility"):
        FractureField(lat, lon, r, np.full(shape, 1.01), 1.0)
    with pytest.raises(ValueError, match="reactive_fraction"):
        FractureField(lat, lon, r, np.ones(shape), -0.1)
