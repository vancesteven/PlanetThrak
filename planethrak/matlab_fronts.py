"""Load archived MATLAB cracking-front structures as validation fixtures.

The repository includes small ``CrackingFrontKic*.mat`` products that preserve
fracture-toughness sensitivity from the historical DeMartin-style generator.
The generator source itself is not present, so these products are valuable
independent anchors for a future reconstructed fracture-mechanics kernel.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import loadmat

from .fronts import CrackingFront


def load_matlab_cracking_front(
    path: str | Path,
    *,
    grain_label: str = "1mm",
) -> CrackingFront:
    """Load ``CrackingFront.T_<grain>_oC`` and ``Pc_<grain>_MPa`` from MAT.

    No interpolation, smoothing, or unit conversion is applied because the MAT
    archive already stores pressure in MPa and temperature in degrees C.
    """
    path = Path(path)
    data = loadmat(path, simplify_cells=True)
    if "CrackingFront" not in data:
        raise KeyError(f"CrackingFront variable not found in {path}")
    record = data["CrackingFront"]
    if not isinstance(record, dict):
        # ``simplify_cells`` normally returns a dict, but keep a clear error for
        # archive variants rather than guessing scipy's internal struct layout.
        raise TypeError(f"unexpected CrackingFront representation in {path}: {type(record)!r}")

    t_key = f"T_{grain_label}_oC"
    p_key = f"Pc_{grain_label}_MPa"
    if t_key not in record or p_key not in record:
        raise KeyError(f"{path} does not contain {t_key!r} and {p_key!r}")

    return CrackingFront(
        pressure_MPa=np.asarray(record[p_key], dtype=float).reshape(-1),
        temperature_C=np.asarray(record[t_key], dtype=float).reshape(-1),
        source=str(path),
    )
