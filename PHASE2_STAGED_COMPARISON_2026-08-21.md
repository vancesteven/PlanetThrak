# PlanetThrak Phase-2 staged Mars structure comparison

Date: 2026-08-21
Branch: `python-port`

Phase 1 is closed for the archived Mars 1-mm, 1 C/yr path. The Python implementation reproduces all 46 archived `Planets.mat` cracking depths to floating-point roundoff and the hard `1e-9 km` regression gate passes.

Phase 2 must therefore change the *physics* without changing the validated numerical kernel invisibly.

## Controlled decomposition

For any PlanetProfile-derived candidate radial column, evaluate four cracking depths on the same depth grid:

1. legacy pressure + legacy temperature,
2. candidate pressure + legacy temperature,
3. legacy pressure + candidate temperature,
4. candidate pressure + candidate temperature.

The pressure-only shift isolates the effect of replacing the historical `rho=3500 kg/m3`, `g=6 m/s2` pressure approximation. The thermal-only shift isolates the modern geotherm/thermal structure. The full P-T shift contains both. The residual

`interaction = total_shift - pressure_shift - thermal_shift`

reports their non-additive coupling.

This decomposition intentionally precedes composition-dependent fracture toughness, elastic moduli, thermal expansion, grain size, pore pressure, and damage. Those are subsequent axes rather than being folded into one opaque modern-vs-legacy change.

## Implemented API

- `legacy_mars_pt_on_depths`: evaluates the archived Mars P-T assumptions on any strictly increasing depth grid.
- `staged_mars_cracking_depths`: evaluates the four states above using one archived cracking front.
- `scripts/compare_profile_to_legacy.py`: loads a neutral `.npz` radial artifact and prints the staged comparison.

## Validation order

```bash
pytest -q tests/test_staged_comparison.py
pytest -q
```

Then, after producing the first Mars PlanetProfile artifact:

```bash
python scripts/compare_profile_to_legacy.py path/to/mars_profile.npz --time-yr 0
```

Repeat at selected lookback times only if the candidate thermal structure is physically defined at those epochs. A present-day PlanetProfile geotherm must not be relabeled as a 4.5-Gyr thermal history.

## Next scientific gate

The first real PlanetProfile comparison should report:

- legacy and modern cracking depth,
- pressure-only contribution,
- thermal-only contribution,
- P-T interaction,
- the exact PlanetProfile/Perple_X provenance of the radial artifact,
- mass and C/MR^2 consistency of the exported structure.

Only after this gate should the fracture boundary itself be generalized away from the archived lookup tables.
