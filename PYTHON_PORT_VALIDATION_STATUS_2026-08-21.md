# PlanetThrak Python validation status

Date: 2026-08-21
Branch: `python-port`

## Purpose

The Python port is being developed in two deliberately separated stages. First, it must reproduce the archived MATLAB PlanetThrak calculation. Only after that parity gate is passed will the code replace the historical radial approximations with PlanetProfile/Perple_X structure and generalize thermal cracking to a local material-state calculation.

No current Python result should be described as a validated replacement for the MATLAB calculation until the archive parity diagnostic below has passed with a documented tolerance.

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

## Archive parity gate

After the unit suite:

```bash
python scripts/legacy_mat_anchor.py --dump
python scripts/legacy_mars_parity.py --show-all
```

Inspect the error distribution before selecting a hard tolerance. Only then run, for example,

```bash
python scripts/legacy_mars_parity.py --assert-km <JUSTIFIED_TOLERANCE_KM>
```

and record the exact command, environment, maximum error, RMS error, and bias here. Do not choose the tolerance merely to make the test pass.

## Current result status

The code and tests above are committed, but the full suite and direct `Planets.mat` parity diagnostic have not yet been executed in the development environment after the latest changes. No pass count or parity tolerance is therefore claimed here yet.

## Next gate after MATLAB parity

1. Define the PlanetProfile radial-column adapter without importing PlanetProfile deeply into PlanetThrak.
2. Feed the *legacy* Mars column through that adapter and require reduction to the Phase-1 result.
3. Feed modern PlanetProfile Mars P-T-density/composition columns through the same interface and quantify how the inferred cracking front changes relative to the historical `rho=3500 kg/m3`, `g=6 m/s2`, constant-k model.
4. Expand to longitude-latitude columns and compute `A_f(theta,phi,z,t)` and fractured-reactive-rock volume with uncertainty.
5. Keep actual hydration efficiency downstream of mechanical accessibility until fluid supply and reaction kinetics are applied.
