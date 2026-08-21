import numpy as np

from planethrak.legacy_structure import legacy_mars_column


def test_legacy_mars_pressure_is_21_mpa_per_km():
    col = legacy_mars_column(0.0)
    np.testing.assert_allclose(col.pressure_MPa[1], 21.0, rtol=0, atol=1e-14)
    np.testing.assert_allclose(col.pressure_MPa[10], 210.0, rtol=0, atol=1e-12)


def test_legacy_mars_present_radiogenic_power_and_flux():
    col = legacy_mars_column(0.0)
    np.testing.assert_allclose(col.radiogenic_power_W, 3.3170856e12, rtol=5e-15)
    np.testing.assert_allclose(col.surface_heat_flux_Wm2, 0.022874716547976, rtol=5e-15)


def test_legacy_mars_column_heats_more_in_the_past():
    present = legacy_mars_column(0.0)
    ancient = legacy_mars_column(4.5e9)
    assert ancient.surface_heat_flux_Wm2 > present.surface_heat_flux_Wm2
    assert ancient.temperature_C[10] > present.temperature_C[10]
    # The archived pressure approximation is time independent.
    np.testing.assert_allclose(ancient.pressure_MPa, present.pressure_MPa)
