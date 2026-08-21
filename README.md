# PlanetThrak: Planetary Thermal Cracking

Thermal fracturing and serpentinization calculations used in Vance et al. (2007, *Astrobiology*), Vance et al. (2016, *Geophysical Research Letters*), and Vance, Melwani Daswani (2020, *Philosophical Transactions of the Royal Society*).

The original MATLAB functions and lookup tables generate simplified radial structures for planets and spherical moons. They are used to estimate the extent of microfracturing and bounding hydrogen production through retrogressive hydration of olivine- and pyroxene-bearing rock. The archived MATLAB calculation is run with `plot_PTCracking_planets.m`.

## Python port

The `python-port` branch is a staged port, not a rewrite of the physics in one jump. Phase 1 preserves the archived `PcT*.ext` cracking-front tables as scientific fixtures and reproduces the MATLAB pressure-temperature/intersection workflow with tests. Only after legacy parity is established will the fracture calculation be generalized for PlanetProfile and spatially variable Mars models.

Current Python pieces include:

- `planethrak.fronts`: validated loading of archived cracking-front P-T tables;
- `planethrak.intersection`: side-effect-free cubic-spline P-T intersection kernel;
- `planethrak.radiogenic`: literal U-Th-K radiogenic-heating port;
- `planethrak.legacy_structure`: literal legacy Mars P-T column used only for MATLAB parity.

Install and run the tests from the repository root with:

```bash
python -m pip install -e '.[dev]'
pytest -q
```

The scientific roadmap for the port and the 3D fracture-accessibility extension is in `PYTHON_PORT_AND_3D_FRACTURE_PLAN.md`.

The long-term architecture separates planetary structure from fracture physics: PlanetProfile/Perple_X supplies self-consistent P-T-density-mineralogical columns; PlanetThrak evaluates local fracture susceptibility and the mechanically accessible depth. A later 3D application will combine fracture accessibility with the spatial abundance of reactive lithologies to estimate the time-dependent volume of fractured reactive rock. That quantity is an upper bound/prior on alteration capacity, not an assumption that every accessible fracture hydrates.
