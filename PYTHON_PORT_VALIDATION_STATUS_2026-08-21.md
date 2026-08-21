# PlanetThrak Python validation status

Date: 2026-08-21
Branch: `python-port`

## Purpose

The Python port is being developed in two deliberately separated stages. First, it must reproduce the archived MATLAB PlanetThrak calculation. Only after that parity gate is passed will the code replace the historical radial approximations with PlanetProfile/Perple_X structure and generalize thermal cracking to a local material-state calculation.

The archived Mars 1-mm, 1 C/yr path has now passed its direct MATLAB-product parity gate at essentially floating-point roundoff. This validates the Python implementation as a reproduction of that specific archived calculation. It does not validate the physical assumptions of the historical constant-density, constant-gravity, constant-conductivity Mars model, which remain deliberately isolated as legacy fixtures.

## Implemented Phase-1 components

- `planethrak.fronts`: parses archived two-column `PcT*.ext` thermal-cracking boundaries without smoothing or regridding.
- `planethrak.intersection`: reproduces the numerical role of MATLAB `get_Pz_cracking` using cubic T(P) splines and bracketed roots.
- `planethrak.radiogenic`: literal U-Th-K radiogenic-heating constants and time convention from the MATLAB archive.
- `planethrak.legacy_structure`: literal historical Mars constant-density, constant-gravity pressure and conductive-temperature column, retained only as a parity fixture.
- `scripts/legacy_mat_anchor.py`: extracts archived Mars cracking fields from `Planets.mat`.
- `scripts/legacy_mars_cracking.py`: runs the Python legacy path at a requested time.
- `scripts/legacy_mars_parity.py`: compares Python cracking depths directly with every archived Mars 1-mm time sample and reports maximum, RMS, and mean depth error.

## Implemented generalized interfaces

- `FractureColumn`: neutral P-T-depth/material-column API independent of PlanetProfile internals.
- `legacy_front_accessibility`: makes the mechanically accessible side of a legacy cracking front explicit rather than returning only a depth.
- `FractureField`: validated full-sphere `(lat, lon, radial shell)` accessibility/reactive-lithology container.
- exact latitude-longitude surface-area weights for equiangular full-sphere grids.
- exact spherical-shell integration of

  `V_fr,react = integral A_f * f_reactive dV`.

- a maximum structurally bound-water mass ceiling that is explicitly labeled a full-hydration upper bound, not a reaction prediction.

## Test order

Run the tests in this order so failures remain interpretable:

```bash
python -m pip install -e '.[dev]'
pytest -q tests/test_fronts.py tests/test_intersection.py tests/test_radiogenic.py
pytest -q tests/test_legacy_structure.py
pytest -q tests/test_column.py tests/test_capacity.py tests/test_field.py
pytest -q
```

The first group tests immutable archive conventions and the root-finding kernel. The second tests the literal Mars parity fixture. The third tests generalized APIs and global integration independently of the archive. The final command detects unintended interactions.

## Archive parity gate: PASSED

The direct comparison was run locally on 2026-08-21 with

```bash
python scripts/legacy_mars_parity.py --show-all
```

using:

- archive: `/Users/svance/Library/CloudStorage/Dropbox/PlanetThrak/Planets.mat`
- cracking front: `/Users/svance/Library/CloudStorage/Dropbox/PlanetThrak/PcT1mm_1oCyr.ext`
- samples: 46 archived Mars time points from 0 to 4.5 Gyr

Reported errors were:

- maximum absolute depth error: `1.45519152e-14 km`
- RMS depth error: `5.10100234e-15 km`
- mean depth bias: `+1.75967453e-15 km`

The largest discrepancy is about `1.5e-11 m`, so the Python and archived MATLAB depth series are identical for practical purposes and differ only at floating-point roundoff. There is no time-dependent drift or systematic bias in the 46-point series.

A regression tolerance of `1e-9 km` is adopted for this legacy anchor. This corresponds to one micrometre in depth, roughly five orders of magnitude looser than the observed maximum numerical residual while remaining negligible relative to every physical uncertainty in the fracture model. The margin is intentional so harmless platform/library-level floating-point differences do not break the regression while any scientifically meaningful implementation change still does.

The hard-gate command is therefore:

```bash
python scripts/legacy_mars_parity.py --assert-km 1e-9
```

The exact Python, NumPy, SciPy, and platform versions used for the first reported run were not captured in the console output and should be recorded on the next full validation run.

## Interpretation of the passed gate

This result establishes a narrow but important claim: for the archived Mars `PcT1mm_1oCyr.ext` path, the Python spline/intersection/legacy-structure workflow reproduces the values stored in `Planets.mat` to roundoff. It is therefore safe to use the Python implementation as the reference implementation for controlled changes to the planetary structure model.

It does **not** establish that the archived cracking depths are physically accurate for Mars. In particular, the historical `rho=3500 kg/m3`, `g=6 m/s2`, constant-thermal-conductivity structure remains an approximation to be replaced and compared explicitly against self-consistent PlanetProfile/Perple_X columns.

## Next gate after MATLAB parity

1. Feed the *legacy* Mars column through the neutral PlanetProfile-style radial adapter and require reduction to the passed Phase-1 result.
2. Feed modern PlanetProfile Mars P-T-density/composition columns through the same interface and quantify how the inferred cracking front changes relative to the historical `rho=3500 kg/m3`, `g=6 m/s2`, constant-k model.
3. Expand to longitude-latitude columns and compute `A_f(theta,phi,z,t)` and fractured-reactive-rock volume with uncertainty.
4. Couple accessibility to spatially variable reactive ultramafic abundance, but keep hydration efficiency downstream of mechanical accessibility until fluid supply and reaction kinetics are applied.
5. Use the resulting common physical alteration state as the shared input to gravity, seismic, magnetic/EM, and pyLOV3D tidal forward models.
