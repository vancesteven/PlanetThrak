from pathlib import Path

import numpy as np

from planethrak.fronts import load_front


ROOT = Path(__file__).resolve().parents[1]


def test_load_legacy_1mm_1cyr_table_exact_values():
    front = load_front(ROOT / "PcT1mm_1oCyr.ext")
    np.testing.assert_allclose(
        front.pressure_MPa,
        [287.78, 262.95, 238.14, 213.34, 188.56, 163.80, 139.08],
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        front.temperature_C,
        [-90, 10, 110, 210, 310, 410, 510],
        rtol=0,
        atol=0,
    )


def test_all_shipped_pct_tables_are_two_column_finite_fronts():
    # PcT*.ext includes the PcTp* effective-pressure variants as well.
    paths = sorted(ROOT.glob("PcT*.ext"))
    # Avoid accidental test success if the legacy fixtures are moved/renamed.
    assert len(paths) >= 9
    for path in paths:
        front = load_front(path)
        assert front.pressure_MPa.size >= 2
        assert np.all(np.isfinite(front.pressure_MPa))
        assert np.all(np.isfinite(front.temperature_C))
        assert np.unique(front.pressure_MPa).size == front.pressure_MPa.size
