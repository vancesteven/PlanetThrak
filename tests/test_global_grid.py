import numpy as np

from planethrak.fronts import CrackingFront
from planethrak.global_grid import fracture_field_from_pt_grids


def _surface_grid():
    lat = np.array([-45.0, 45.0])
    lon = np.array([-135.0, -45.0, 45.0, 135.0])
    R = 3_390_000.0
    edges = np.array([R - 3000.0, R - 2000.0, R - 1000.0, R])
    return lat, lon, edges


def test_two_pt_columns_produce_different_cracking_depths_and_accessibility():
    lat, lon, edges = _surface_grid()
    # Shell-center depths are [2500, 1500, 500] m in deep->shallow order.
    depth = np.array([2500.0, 1500.0, 500.0])
    p1 = depth / 100.0  # [25, 15, 5] MPa
    # Constant front T=15 C. Column 0 has T=P, crossing at depth 1500 m;
    # column 1 has T=P+10, crossing at depth 500 m.
    front = CrackingFront(
        pressure_MPa=np.array([0.0, 30.0]),
        temperature_C=np.array([15.0, 15.0]),
    )

    p = np.broadcast_to(p1, (lat.size, lon.size, depth.size)).copy()
    tc = np.empty_like(p)
    tc[:, 0::2, :] = p[:, 0::2, :]
    tc[:, 1::2, :] = p[:, 1::2, :] + 10.0

    result = fracture_field_from_pt_grids(
        latitude_deg=lat,
        longitude_deg=lon,
        radius_edges_m=edges,
        pressure_MPa=p,
        temperature_K=tc + 273.15,
        reactive_fraction=0.5,
        front=front,
    )

    np.testing.assert_allclose(result.cracking_depth_m[:, 0::2], 1500.0, atol=1e-8)
    np.testing.assert_allclose(result.cracking_depth_m[:, 1::2], 500.0, atol=1e-8)
    np.testing.assert_array_equal(result.field.accessibility[:, 0, :], [[0, 1, 1], [0, 1, 1]])
    np.testing.assert_array_equal(result.field.accessibility[:, 1, :], [[0, 0, 1], [0, 0, 1]])
    np.testing.assert_allclose(result.field.reactive_fraction, 0.5)


def test_global_grid_integral_uses_same_accessibility_field():
    lat, lon, edges = _surface_grid()
    depth = np.array([2500.0, 1500.0, 500.0])
    p1 = depth / 100.0
    p = np.broadcast_to(p1, (lat.size, lon.size, depth.size)).copy()
    tc = p.copy()
    front = CrackingFront(np.array([0.0, 30.0]), np.array([15.0, 15.0]))

    result = fracture_field_from_pt_grids(
        latitude_deg=lat,
        longitude_deg=lon,
        radius_edges_m=edges,
        pressure_MPa=p,
        temperature_K=tc + 273.15,
        reactive_fraction=1.0,
        front=front,
    )
    # Every column has identical [0,1,1] accessibility, so the global volume
    # must equal the exact volumes of the two shallow shells.
    r = edges
    expected = (4.0 * np.pi / 3.0) * (r[-1] ** 3 - r[1] ** 3)
    np.testing.assert_allclose(result.field.fractured_reactive_volume_m3(), expected, rtol=2e-15)
