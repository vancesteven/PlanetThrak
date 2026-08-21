import numpy as np
import pytest

from planethrak.radiogenic import past_radiogenic_heat_uthk


def test_mars_legacy_radiogenic_power_reference_values():
    # Literal Mars mass/f_Si values from plot_PTCracking_planets.m.
    times = np.array([0.0, 1.0e9, 4.5e9])
    power = past_radiogenic_heat_uthk(6.42e23, 1.0, times)
    np.testing.assert_allclose(
        power,
        [3.3170856e12, 4.883076921677151e12, 2.5339102075894004e13],
        rtol=5e-15,
        atol=0,
    )


def test_radiogenic_power_increases_backward_in_time():
    power = past_radiogenic_heat_uthk(6.42e23, 1.0, np.array([0.0, 1e9, 4e9]))
    assert np.all(np.diff(power) > 0)


def test_radiogenic_inputs_are_guarded():
    with pytest.raises(ValueError):
        past_radiogenic_heat_uthk(-1.0, 1.0, 0.0)
    with pytest.raises(ValueError):
        past_radiogenic_heat_uthk(1.0, 1.1, 0.0)
    with pytest.raises(ValueError):
        past_radiogenic_heat_uthk(1.0, 1.0, -1.0)
