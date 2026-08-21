from pathlib import Path

import numpy as np

from planethrak.matlab_fronts import load_matlab_cracking_front


ROOT = Path(__file__).resolve().parents[1]


def test_kic_06_mat_front_exact_archive_values():
    front = load_matlab_cracking_front(ROOT / "CrackingFrontKicp6_1mm_1oCyr.mat")
    np.testing.assert_allclose(
        front.temperature_C,
        [-90, 10, 110, 210, 310, 410, 510],
        rtol=0,
        atol=0,
    )
    np.testing.assert_allclose(front.pressure_MPa[0], 287.73794424, rtol=0, atol=5e-9)
    np.testing.assert_allclose(front.pressure_MPa[-1], 139.03068172, rtol=0, atol=5e-9)


def test_higher_fracture_toughness_shifts_archived_front_to_lower_pressure():
    k06 = load_matlab_cracking_front(ROOT / "CrackingFrontKicp6_1mm_1oCyr.mat")
    k09 = load_matlab_cracking_front(ROOT / "CrackingFrontKicp9_1mm_1oCyr.mat")
    np.testing.assert_array_equal(k06.temperature_C, k09.temperature_C)
    assert np.all(k09.pressure_MPa < k06.pressure_MPa)
    np.testing.assert_allclose(k09.pressure_MPa[0], 250.61855646, rtol=0, atol=5e-9)
    np.testing.assert_allclose(k09.pressure_MPa[-1], 102.72179226, rtol=0, atol=5e-9)
