# API overview

This is the auto-generated API reference for CALLIOPE's public interface. For the underlying physics, see the [model overview](../../Explanations/model.md); for the package layout in detail, see [code architecture](../../Explanations/code_architecture.md).

## Module overview

```
src/calliope/
├── __init__.py          # version string only
├── constants.py         # molar masses, element list, plotting colours
├── oxygen_fugacity.py   # OxygenFugacity (IW buffer parameterisations)
├── chemistry.py         # ModifiedKeq (six redox couples)
├── solubility.py        # Solubility{H2O,CO2,CO,CH4,N2,S2}
├── solve.py             # equilibrium_atmosphere() and supporting functions
└── structure.py         # calculate_mantle_mass() simple core/mantle split
```

## API reference

### Solver

- [`calliope.solve`](calliope.solve.md) - `equilibrium_atmosphere()` (the public entry point), the residual function `func`, the L2-norm objective `obj`, the speciation tree `get_partial_pressures`, the column-mass aggregation `atmosphere_mass`, the dissolved-mass aggregation `dissolved_mass`, and the elemental-target builders `get_target_from_params` / `get_target_from_pressures`.

### Chemistry

- [`calliope.chemistry`](calliope.chemistry.md) - `ModifiedKeq` class with eight model methods (`schaefer_C`, `schaefer_H`, `schaefer_CH4`, `janaf_CO`, `janaf_H2`, `janaf_SO2`, `janaf_H2S`, `janaf_NH3`) returning the modified equilibrium constant $G_\mathrm{eq}(T, f_{\mathrm{O}_2})$.

### Solubility

- [`calliope.solubility`](calliope.solubility.md) - `Solubility` base class plus six concrete subclasses with multiple alternative compositions per species. Returns dissolved-volatile concentration in ppmw given partial pressure (and, where relevant, $T_\mathrm{magma}$, $p_\mathrm{tot}$, $\Delta\mathrm{IW}$).

### Oxygen fugacity

- [`calliope.oxygen_fugacity`](calliope.oxygen_fugacity.md) - `OxygenFugacity` class with two model methods (`oneill` default, `fischer`). Returns $\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T) + \Delta\mathrm{IW}$.

### Structure

- [`calliope.structure`](calliope.structure.md) - `calculate_mantle_mass(radius, mass, core_frac)` simple Earth-density core/mantle split. Provided as a fallback; PROTEUS-coupled runs source `M_mantle` from Zalmoxis or SPIDER instead.

### Constants

- [`calliope.constants`](calliope.constants.md) - Static reference values: planetary masses and radii, ocean-of-H scale, molar masses, the eleven-species `volatile_species` list, the five-element `element_list`, and the PROTEUS-shared plotting colour scheme `dict_colors`.

## What is not in the public API

CALLIOPE intentionally has a flat module structure with no underscore-prefixed sub-packages or private submodules. Every function or class declared in the modules above is part of the public API. The exception is `solve.func` and `solve.obj` (the residual and objective functions passed to `scipy.optimize`); these are public for testing and inspection purposes but should not be called directly outside of an `equilibrium_atmosphere` flow because they assume the input dict is in a particular state (e.g. zero-clipped, IW-shifted, included-flagged).
