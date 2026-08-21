"""Legacy PlanetThrak thermal-cracking P-T front tables.

The archived ``PcT*.ext`` files are treated as scientific reference fixtures.
Column 1 is pressure in Pa and column 2 is temperature in degrees Celsius, as
used by ``plot_PTCracking_planets.m`` before pressure is divided by 1e6.

Most files contain an interpolable P-T cracking curve, but the archive also
contains at least one legitimate single-point effective-pressure artifact.  The
loader preserves such files exactly.  Code that needs an interpolated cracking
boundary must explicitly require ``is_interpolable`` rather than silently
reinterpreting an incomplete lookup as a physical curve.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class CrackingFront:
    """Archived pressure-temperature cracking boundary or boundary fragment."""

    pressure_MPa: np.ndarray
    temperature_C: np.ndarray
    source: str | None = None

    def __post_init__(self) -> None:
        p = np.asarray(self.pressure_MPa, dtype=float)
        t = np.asarray(self.temperature_C, dtype=float)
        if p.ndim != 1 or t.ndim != 1 or p.size != t.size:
            raise ValueError("pressure_MPa and temperature_C must be equal-length 1-D arrays")
        if p.size < 1:
            raise ValueError("a cracking-front archive fixture requires at least one P-T point")
        if not np.all(np.isfinite(p)) or not np.all(np.isfinite(t)):
            raise ValueError("cracking-front arrays must be finite")
        if np.any(p < 0):
            raise ValueError("cracking-front pressure cannot be negative")
        if np.unique(p).size != p.size:
            raise ValueError("cracking-front pressure values must be unique")
        object.__setattr__(self, "pressure_MPa", p)
        object.__setattr__(self, "temperature_C", t)

    @property
    def is_interpolable(self) -> bool:
        """Whether the archived fixture contains enough points for T(P)."""
        return self.pressure_MPa.size >= 2


def load_front(path: str | Path) -> CrackingFront:
    """Load one archived two-column ``PcT*.ext`` table.

    The MATLAB files store pressure in Pa; the Python public interface uses MPa.
    No smoothing or regridding is performed here.  Single-row archive files are
    preserved as non-interpolable ``CrackingFront`` objects.
    """

    path = Path(path)
    data = np.loadtxt(path, dtype=float)
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] != 2:
        raise ValueError(f"expected two columns in {path}, found {data.shape[1]}")
    return CrackingFront(
        pressure_MPa=data[:, 0] / 1.0e6,
        temperature_C=data[:, 1],
        source=str(path),
    )
