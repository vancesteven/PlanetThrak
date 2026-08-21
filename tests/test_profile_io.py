from pathlib import Path

import numpy as np

from planethrak.planetprofile import fracture_column_from_arrays
from planethrak.profile_io import load_fracture_column_npz, save_fracture_column_npz


def test_profile_npz_roundtrip(tmp_path: Path):
    column = fracture_column_from_arrays(
        depth_m=[0.0, 1000.0, 2000.0],
        pressure_MPa=[0.1, 5.0, 10.0],
        temperature_K=[220.0, 240.0, 260.0],
        density_kg_m3=[2900.0, 2920.0, 2950.0],
        gravity_m_s2=[3.71, 3.70, 3.69],
        reactive_fraction=[0.2, 0.4, 0.6],
        pore_pressure_MPa=[0.0, 2.0, 4.0],
        bulk_modulus_Pa=[70e9, 72e9, 75e9],
        shear_modulus_Pa=[30e9, 31e9, 32e9],
        vp_m_s=[6000.0, 6100.0, 6200.0],
        vs_m_s=[3200.0, 3250.0, 3300.0],
    )
    path = save_fracture_column_npz(
        tmp_path / "mars_column.npz",
        column,
        body="Mars",
        source="unit-test",
        time_yr=4.5e9,
    )
    loaded, meta = load_fracture_column_npz(path)
    np.testing.assert_allclose(loaded.depth_m, column.depth_m)
    np.testing.assert_allclose(loaded.pressure_MPa, column.pressure_MPa)
    np.testing.assert_allclose(loaded.temperature_C, column.temperature_C, atol=1e-12)
    np.testing.assert_allclose(loaded.density_kg_m3, column.density_kg_m3)
    np.testing.assert_allclose(loaded.gravity_m_s2, column.gravity_m_s2)
    np.testing.assert_allclose(loaded.reactive_fraction, column.reactive_fraction)
    np.testing.assert_allclose(loaded.pore_pressure_MPa, column.pore_pressure_MPa)
    np.testing.assert_allclose(loaded.bulk_modulus_Pa, column.bulk_modulus_Pa)
    np.testing.assert_allclose(loaded.shear_modulus_Pa, column.shear_modulus_Pa)
    np.testing.assert_allclose(loaded.vp_m_s, column.vp_m_s)
    np.testing.assert_allclose(loaded.vs_m_s, column.vs_m_s)
    assert meta["body"] == "Mars"
    assert meta["source"] == "unit-test"
    assert float(meta["time_yr"]) == 4.5e9
