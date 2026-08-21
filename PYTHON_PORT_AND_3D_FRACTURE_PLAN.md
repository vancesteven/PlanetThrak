# PlanetThrak Python port and 3D fracture-generalization plan

Date: 2026-08-20

## Scientific objective

Port the Vance et al. PlanetThrak thermal-cracking calculation from MATLAB to a tested Python library that can consume PlanetProfile/Perple_X interior states and, after exact legacy reproduction, generalize the cracking calculation from one-dimensional hard-coded planetary columns to spatially variable Mars models. The end product should predict where rock is mechanically accessible to water as a function of position, depth, time, composition, thermal history, and effective stress. It should not duplicate the groundwater-flow or reactive-transport calculations that consume this accessibility field.

## What the current MATLAB code actually does

The current repository is a compact research archive rather than a reusable library. `plot_PTCracking_planets.m` hard-codes planetary cases, constructs simplified pressure and conductive-temperature profiles, loads precomputed DeMartin-style cracking-front curves for 0.1, 1, and 10 mm grains and three cooling-rate labels, and finds the intersection of each planetary P-T path with the selected cracking-front P-T curve. `SerpHeatOceanPlanets.m` then converts the time-dependent cracking depth into bounding serpentinization heat/H2 production estimates. For Mars specifically, the present pressure approximation is hard-coded as `3500*6*z/1e6` MPa, and the thermal profile uses a constant conductivity and globally averaged heat flow. Those approximations were appropriate for the original comparative-planet calculation but should not be inherited into the new Mars application when PlanetProfile already supplies self-consistent P-T-density structure.

The `.ext` cracking-front tables are therefore part of the legacy scientific reference and should be preserved as validation fixtures. The generator referred to in comments as `Calc_z_CrackingDeMartin2004.m` is not present in the current repository, so Phase 1 should reproduce the shipped lookup-table behavior exactly before attempting to reconstruct or replace that upstream fracture calculation.

## Phase 1: exact Python reproduction

Create a small package, provisionally `planethrak`, with no plotting side effects and no global variables. The first implementation should intentionally preserve the legacy physics and data conventions.

Suggested modules:

- `planethrak.fronts`: read and validate the shipped P-T cracking-front tables; represent grain size and cooling-rate metadata explicitly.
- `planethrak.structure`: reproduce the legacy two-layer pressure and conductive-temperature helpers for reference tests only.
- `planethrak.intersection`: robustly solve for the P-T intersection and corresponding cracking depth, replacing MATLAB `spline`/`fzero` with explicit monotonic interpolation and bracketed root finding.
- `planethrak.radiogenic`: port `get_pastRadiogenicHeat_UThK.m` with isotope constants and units documented.
- `planethrak.serpentinization_bounds`: port the legacy accessible-volume-to-heat/H2 bookkeeping as an explicitly labeled upper/bounding calculation.
- `planethrak.legacy_cases`: Earth, Mars, Europa, Enceladus, and other archived body cases needed only to reproduce published/archived results.

Acceptance gates:

1. Parse every shipped `.ext` cracking-front table without modification.
2. Reproduce archived Mars cracking depths for each available grain-size/cooling-rate combination to a stated numerical tolerance.
3. Reproduce a representative icy-body case to ensure the two-layer gravity/pressure branch is also preserved.
4. Reproduce selected `Planets.mat`/`ThermalCracking_Planets.mat` fields from the same inputs. Where MAT-file contents are too coarse or provenance is ambiguous, record that limitation rather than weakening the test silently.
5. Unit-test all units at API boundaries: Pa/MPa, m/km, K/degC, s/yr/Gyr.

The Python port should be considered scientifically equivalent to legacy PlanetThrak only after these gates pass.

## Phase 2: separate the fracture criterion from the planetary structure model

The key redesign is to make cracking a local material-state calculation instead of a property of a hard-coded planet. Define a fracture-state input such as

`FractureState(P, T, Tdot, rho, g, grain_size, E, nu, alpha_tensor_or_mismatch, K_IC, pore_pressure, preexisting_damage, phase_fractions)`.

The initial generalized criterion can still use the legacy cracking-front tables, but the API should expose a dimensionless fracture susceptibility or margin, for example

`F = K_I / K_IC`

or an equivalent signed distance from the calibrated legacy P-T cracking boundary. `F >= 1` then denotes material expected to be thermally microfractured under the selected constitutive model. This makes the later replacement of the lookup table by a reconstructed DeMartin/Vance fracture-mechanics kernel possible without changing downstream code.

Important uncertainty axes to expose rather than hard-code are grain size, cooling rate, thermal-expansion anisotropy/mismatch, Young's modulus, Poisson ratio, fracture toughness, pore pressure/effective stress, and pre-existing damage. Composition-dependent values should be supplied by the caller, not hidden inside PlanetThrak.

## Phase 3: PlanetProfile/Perple_X adapter

PlanetProfile should own hydrostatic/thermal structure and composition. Perple_X-derived phase equilibria should provide density and elastic/mineralogical properties where appropriate. PlanetThrak should consume a neutral radial or 3D state rather than importing PlanetProfile internals deeply.

Minimum fields for a radial column:

- radius/depth
- pressure
- temperature and, for evolutionary calculations, cooling rate or T(t)
- density and gravity
- phase/mineral fractions
- elastic moduli or E/nu
- thermal-expansion properties
- fracture toughness or an externally supplied lithology-dependent prior
- pore pressure/effective stress when available

The adapter returns fracture susceptibility versus depth and a mechanically accessible depth/front. A reduction test must demonstrate that feeding the old simplified Mars P-T state through the adapter reproduces Phase-1 legacy results before using modern PlanetProfile states.

## Phase 4: 3D Mars fracture-accessibility field

For each sample from the global composition/thermal/gravity model, evaluate the local fracture criterion column-by-column to obtain

- `A_f(lon, lat, z, t)`: fracture accessibility/susceptibility,
- `b_f(lon, lat, t)`: depth of the mechanically accessible/fractured domain,
- uncertainty envelopes over material properties and thermal histories.

This is the natural interface to the groundwater/reactive-transport work: Task 1 supplies `A_f` and `b_f`; the flow/reaction model determines whether water actually reaches those fractures and how much reaction occurs. Mechanical accessibility must therefore not be equated with hydration fraction.

## Phase 5: combine fracture accessibility with 3D composition and gravity

The new 3D composition model makes PlanetThrak scientifically more powerful than the original depth-only calculation. For each posterior sample, define a dry/reference mineralogical density and elastic structure, then propagate candidate alteration through density, rigidity, seismic, magnetic, and electrical-property contrasts. The composition model identifies where reactive ultramafic/olivine-rich material exists, while the fracture model identifies where that material is mechanically accessible.

Useful Task-1 products include:

- global volume of fractured reactive lithology versus time,
- an upper bound on alterable rock and structurally bound water if accessible rock were fully hydrated,
- spatial priors on hydration fraction for the joint geophysical inversion,
- predicted density-deficit fields relative to a compositionally consistent dry reference,
- predicted rigidity and seismic-property fields for pyLOV3D and seismic forward models,
- predicted magnetite-bearing alteration source geometry for magnetic forward models.

For a simple two-endmember diagnostic, an inferred alteration fraction can be written schematically as

`f_h ~ (rho_dry - rho_obs) / (rho_dry - rho_hyd)`

but production inference must include porosity, crustal thickness/compensation, temperature, and compositional covariance. The recent GMM-3 tests show why: formal coefficient uncertainties can be far below plausible hydration signals, so the limiting uncertainty is geological attribution, not coefficient precision.

## Phase 6: feedbacks, only after the one-way model is validated

Do not begin with a fully coupled damage-reaction model. First validate one-way coupling:

composition + thermal history -> fracture accessibility -> candidate alteration -> geophysical observables.

Only then add optional feedbacks such as hydration-induced volume change, reaction heat, changes in thermal conductivity/elasticity, crack sealing, and reactive cracking. These feedbacks should be separate constitutive modules so they can be enabled individually and tested against the one-way baseline.

## Validation strategy

Validation should have four levels:

1. MATLAB legacy parity for archived PlanetThrak cases.
2. Constitutive tests for monotonic trends with grain size, cooling rate, toughness, pore pressure, and thermal-expansion mismatch.
3. Earth analog validation using the Mid-Continent Rift, where gravity, magnetic/EM, and independent seismic constraints can test whether the same inferred alteration geometry is compatible with all observables.
4. Mars closure tests against mass/MoI, InSight seismic structure, GMM-3 gravity, crustal magnetic maps, and pyLOV3D tidal response.

## Immediate implementation sequence

1. Freeze checksums and metadata for every legacy `.ext` and selected `.mat` validation artifact.
2. Write table loaders and the P-T intersection solver.
3. Build archived Mars parity tests.
4. Port radiogenic/thermal helpers only as needed to reproduce the archive.
5. Define the neutral `FractureState`/result API.
6. Add a PlanetProfile radial adapter.
7. Replace the hard-coded Mars P-T approximation with PlanetProfile columns and quantify the difference.
8. Generalize to a longitude-latitude set of columns from the Task-1 3D posterior.
9. Couple fracture accessibility to compositionally available reactive-rock fraction and generate global alteration-capacity maps.
10. Only after these gates, add reaction/damage feedbacks.
