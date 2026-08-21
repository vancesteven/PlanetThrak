from types import SimpleNamespace

import numpy as np

from planethrak.legacy_structure import legacy_mars_column
from planethrak.planetprofile import (
    fracture_column_from_arrays,
    fracture_column_from_planetprofile,
)


def test_array_adapter_reproduces_legacy_pt_column_exactly():
    legacy = legacy_mars_column(0.0, dz_m=50_000.0)
    column = fracture_column_from_arrays(
        depth_m=legacy.depth_m,
        pressure_MPa=legacy.pressure_MPa,
        temperature_K=legacy.temperature_C + 273.15,
    )
    np.testing.assert_allclose(column.depth_m, legacy.depth_m, rtol=0, atol=0)
    np.testing.assert_allclose(column.pressure_MPa, legacy.pressure_MPa, rtol=0, atol=0)
    np.testing.assert_allclose(column.temperature_C, legacy.temperature_C, rtol=0, atol=1e-12)


def test_planetprofile_duck_adapter_maps_current_planetstruct_names_and_elasticity():
    K_GPa = np.array([70.0, 72.0, 75.0])
    G_GPa = np.array([30.0, 31.0, 32.0])
    vp_kms = np.sqrt((K_GPa + 4.0 * G_GPa / 3.0) * 1e9 / np.array([2900.0, 2920.0, 2950.0])) / 1e3
    vs_kms = np.sqrt(G_GPa * 1e9 / np.array([2900.0, 2920.0, 2950.0])) / 1e3
    planet = SimpleNamespace(
        z_m=np.array([0.0, 1000.0, 2000.0]),
        P_MPa=np.array([0.1, 5.0, 10.0]),
        T_K=np.array([220.0, 240.0, 260.0]),
        rho_kgm3=np.array([2900.0, 2920.0, 2950.0]),
        g_ms2=np.array([3.71, 3.70, 3.69]),
        alpha_pK=np.array([2e-5, 2e-5, 2e-5]),
        kTherm_WmK=np.array([3.0, 3.1, 3.2]),
        phi_frac=np.array([0.1, 0.05, 0.0]),
        Ppore_MPa=np.array([0.0, 2.0, 4.0]),
        Seismic=SimpleNamespace(KS_GPa=K_GPa, GS_GPa=G_GPa, VP_kms=vp_kms, VS_kms=vs_kms),
    )
    column = fracture_column_from_planetprofile(
        planet,
        reactive_fraction=np.array([0.2, 0.4, 0.6]),
    )
    np.testing.assert_allclose(column.temperature_C, planet.T_K - 273.15)
    np.testing.assert_allclose(column.density_kg_m3, planet.rho_kgm3)
    np.testing.assert_allclose(column.gravity_m_s2, planet.g_ms2)
    np.testing.assert_allclose(column.thermal_expansivity_Kinv, planet.alpha_pK)
    np.testing.assert_allclose(column.thermal_conductivity_W_mK, planet.kTherm_WmK)
    np.testing.assert_allclose(column.porosity_fraction, planet.phi_frac)
    np.testing.assert_allclose(column.pore_pressure_MPa, planet.Ppore_MPa)
    np.testing.assert_allclose(column.bulk_modulus_Pa, K_GPa * 1e9)
    np.testing.assert_allclose(column.shear_modulus_Pa, G_GPa * 1e9)
    np.testing.assert_allclose(column.vp_m_s, vp_kms * 1e3)
    np.testing.assert_allclose(column.vs_m_s, vs_kms * 1e3)
    expected_E = 9.0 * K_GPa * G_GPa / (3.0 * K_GPa + G_GPa) * 1e9
    expected_nu = (3.0 * K_GPa - 2.0 * G_GPa) / (2.0 * (3.0 * K_GPa + G_GPa))
    np.testing.assert_allclose(column.youngs_modulus_Pa, expected_E)
    np.testing.assert_allclose(column.poisson_ratio, expected_nu)


def test_adapter_sorts_deep_to_shallow_input_and_applies_mask():
    column = fracture_column_from_arrays(
        depth_m=np.array([3000.0, 2000.0, 1000.0, 0.0]),
        pressure_MPa=np.array([30.0, 20.0, 10.0, 0.1]),
        temperature_K=np.array([300.0, 280.0, 260.0, 240.0]),
        mask=np.array([False, True, True, True]),
        reactive_fraction=np.array([0.0, 0.6, 0.4, 0.2]),
    )
    np.testing.assert_array_equal(column.depth_m, [0.0, 1000.0, 2000.0])
    np.testing.assert_allclose(column.reactive_fraction, [0.2, 0.4, 0.6])


def test_adapter_can_reconstruct_depth_from_radius():
    planet = SimpleNamespace(
        z_m=None,
        r_m=np.array([3.39e6, 3.389e6, 3.388e6]),
        P_MPa=np.array([0.1, 5.0, 10.0]),
        T_K=np.array([220.0, 240.0, 260.0]),
        Bulk=SimpleNamespace(R_m=3.39e6),
        rho_kgm3=None,
        g_ms2=None,
        alpha_pK=None,
        kTherm_WmK=None,
        phi_frac=None,
        Ppore_MPa=None,
        Seismic=None,
    )
    column = fracture_column_from_planetprofile(planet)
    np.testing.assert_allclose(column.depth_m, [0.0, 1000.0, 2000.0])
