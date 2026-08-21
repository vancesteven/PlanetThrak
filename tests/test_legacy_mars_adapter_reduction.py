"""Reduction test from the neutral PlanetProfile-style adapter to legacy Mars.

This is the bridge between Phase 1 archive reproduction and Phase 2 modern
PlanetProfile structure.  It proves that replacing the hard-coded structure
container with the neutral adapter changes no cracking result when the same
legacy P-T path is supplied.
"""

from pathlib import Path

import numpy as np
from scipy.io import loadmat

from planethrak.fronts import load_front
from planethrak.legacy_structure import legacy_mars_column
from planethrak.planetprofile import fracture_column_from_arrays
from planethrak.column import legacy_front_accessibility


ROOT = Path(__file__).resolve().parents[1]


def _records(value):
    if isinstance(value, dict):
        if "name" in value:
            yield value
        else:
            for item in value.values():
                yield from _records(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _records(item)
    elif isinstance(value, np.ndarray):
        for item in value.flat:
            yield from _records(item)


def _mars_record():
    data = loadmat(ROOT / "Planets.mat", simplify_cells=True)
    for record in _records(data["Planets"]):
        if str(record.get("name", "")).strip().lower() == "mars":
            return record
    raise AssertionError("Mars record not found in Planets.mat")


def test_neutral_adapter_reproduces_all_archived_mars_depths():
    mars = _mars_record()
    times = np.asarray(mars["t_yr"], dtype=float).ravel()
    archived = np.asarray(mars["z_cracking_1mm_m"], dtype=float).squeeze().ravel()
    assert archived.size == times.size
    front = load_front(ROOT / "PcT1mm_1oCyr.ext")

    predicted = np.empty_like(times)
    for i, time_yr in enumerate(times):
        legacy = legacy_mars_column(float(time_yr))
        column = fracture_column_from_arrays(
            depth_m=legacy.depth_m,
            pressure_MPa=legacy.pressure_MPa,
            temperature_K=legacy.temperature_C + 273.15,
        )
        result = legacy_front_accessibility(column, front)
        predicted[i] = 0.0 if result.intersection is None else result.intersection.depth_m

    # The adapter path should add no material numerical error.  A 1-micron
    # tolerance is deliberately far above the observed ~1e-11 m residual yet
    # negligible relative to any physical cracking-depth uncertainty.
    np.testing.assert_allclose(predicted, archived, rtol=0.0, atol=1.0e-6)
