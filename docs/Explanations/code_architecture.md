# Code architecture

CALLIOPE is a small package: seven Python modules totaling under 1000 raw lines, no compiled extensions, no thread or process pools. This page maps each module to its responsibility and shows how they call each other.

## Package layout

```
src/calliope/
├── __init__.py              # version string only
├── constants.py             # molar masses, element list, plotting colours
├── oxygen_fugacity.py       # IW buffer parameterisations
├── chemistry.py             # ModifiedKeq class, six redox couples
├── solubility.py            # five solubility classes
├── solve.py                 # speciation, mass balance, outer solver
└── structure.py             # simple core/mantle mass split
```

Plus the test suite:

```
tests/
├── test_init.py             # version smoke test
├── test_core.py             # per-module unit tests
└── test_stoichiometry.py    # cross-module mass-balance tests
```

And the development tooling under `tools/` (Jupyter notebooks for fitting JANAF tables to Schaefer-style functional forms; the JANAF source data files for H$_2$S, NH$_3$, S$_2$).

## Module responsibilities

### `constants.py`

Exports the static data CALLIOPE needs:

- `M_earth`, `R_earth`, `R_core_earth`, `M_core_earth`, `ocean_moles`: planetary reference values for unit conversions.
- `const_G`, `mol`, `R_gas`: physical constants.
- `molar_mass`: dictionary mapping species names to kg/mol.
- `volatile_species`: ordered list of the eleven species CALLIOPE handles.
- `element_list`: ordered list of the five elements (H, O, C, N, S).
- `dict_colors`: colour scheme used in PROTEUS plots (so figures across the ecosystem stay consistent).

Pure data; no side effects, no logic.

### `oxygen_fugacity.py`

Single class `OxygenFugacity` with two model methods (`fischer` default, `oneill` legacy). Stateless: instantiate once, call repeatedly with `(T, fO2_shift)`. See [Oxygen fugacity](oxygen_fugacity.md) for the equations.

### `chemistry.py`

Single class `ModifiedKeq` with eight model methods (`schaefer_C`, `schaefer_H`, `schaefer_CH4`, `janaf_CO`, `janaf_H2`, `janaf_SO2`, `janaf_H2S`, `janaf_NH3`). Composes with `OxygenFugacity` to produce the modified equilibrium constant $G_\mathrm{eq}(T, f_{\mathrm{O}_2})$. See [Equilibrium chemistry](equilibrium_chemistry.md).

### `solubility.py`

Six concrete subclasses of an abstract `Solubility` base: `SolubilityH2O`, `SolubilityCO2`, `SolubilityCO`, `SolubilityCH4`, `SolubilityN2`, `SolubilityS2`. Each carries multiple alternative composition methods selectable through the constructor. The `__call__` interface is uniform (`solubility(p, *args)` returns ppmw) but the `*args` differ per species: H$_2$O takes only $p$, CO$_2$ takes $(p, T)$, S$_2$ takes $(p, T, \Delta\mathrm{IW})$, etc. This irregularity is by design; flat-out wrapping into a uniform signature would hide the actual physical dependences.

See [Solubility laws](solubility.md) for the equations.

### `solve.py`

The orchestration layer. Two solver entry points, the buffered mode and the authoritative-O mode, share the speciation tree and aggregation functions:

- `equilibrium_atmosphere(target, ddict, ...)`: buffered-mode outer driver. Takes a four-key `target` and the IW-buffer offset in `ddict`, solves the $4\times 4$ H/C/N/S mass-balance system.
- `equilibrium_atmosphere_authoritative_O(target_d, ddict, ...)`: authoritative-O outer driver. Takes a five-key `target_d` (H, C, N, S, O) and solves the $5\times 5$ system for the four pressures plus $\Delta\mathrm{IW}$. See [Authoritative-oxygen mode](authoritative_oxygen.md) for the augmented mass balance.

Plus seven shared helpers:

- `get_partial_pressures(pin, ddict)`: walks the eleven-species speciation tree from the four primary pressures.
- `atmosphere_mass(pin, ddict)`: applies Bower et al. (2019)[^cite-bower2019] Eq. 2 to every species and aggregates atomic-mass tallies per element.
- `dissolved_mass(pin, ddict)`: applies the chosen solubility law for each soluble species and aggregates atomic-mass tallies per element.
- `get_target_from_params(ddict)`: translates `hydrogen_earth_oceans`, `CH_ratio`, `nitrogen_ppmw`, `sulfur_ppmw` into kg-per-element targets (four-key, for the buffered mode).
- `get_target_from_pressures(ddict)`: back-computes kg-per-element targets from prescribed initial atmospheric pressures.
- `get_initial_pressures(target_d)`: cold-start guess generator for the four primary pressures.
- `get_initial_pressures_with_fO2(target_d, fO2_hint, ...)`: cold-start guess generator for the five-unknown authoritative-O solver, with a separate restart redraw for $\Delta\mathrm{IW}$.

Plus the inner-loop residual functions `func` (four-vector for buffered) and `func_authoritative_O` (five-vector), the L2-norm objectives `obj` and `obj_authoritative_O`, and three private versions (`_get_partial_pressures`, `_atmosphere_mass`, `_dissolved_mass`) that take `fO2_shift` as a positional argument so both solver modes can share the physics path.

`get_partial_pressures` is *the* call that defines what species CALLIOPE knows about. To add a new species, you add an entry in `volatile_species` (constants.py), an entry in `molar_mass` (constants.py), a `ModifiedKeq` model method (chemistry.py) if it is a derived species, a `Solubility` subclass (solubility.py) if it has dissolved-melt physics, and the speciation step in `get_partial_pressures` plus the atomic-mass tally in `atmosphere_mass`. Tests need the analytical-vs-code consistency check in `tests/test_stoichiometry.py`.

### `structure.py`

Single function `calculate_mantle_mass(radius, mass, core_frac)` that returns `M_mantle = M_planet - M_core` assuming Earth-like core density. This is a fallback for use cases where a full structure solve (Zalmoxis, SPIDER) is not warranted. The PROTEUS coupling does *not* go through this function; it gets `M_mantle` directly from Zalmoxis or SPIDER via `hf_row`.

## Call graph

```
            equilibrium_atmosphere              equilibrium_atmosphere_authoritative_O
                       │                                          │
                       ▼                                          ▼
              get_initial_pressures              get_initial_pressures_with_fO2
                       │                                          │
                       └────────────┬─────────────────────────────┘
                                    ▼
                    scipy.optimize.fsolve / minimize  ◄─── ModifiedKeq.__call__
                                    │                              ▲
                                    ▼                              │
                       func / func_authoritative_O                 │
                                    │                              │
                  ┌─────────────────┴─────────────────┐            │
                  ▼                                   ▼            │
        _atmosphere_mass                     _dissolved_mass       │
                  │                                   │            │
                  │     ┌─────────────────────────────┤            │
                  │     ▼                             ▼            │
                  │     _get_partial_pressures  Solubility{...}.__call__
                  │             │                                  │
                  ▼             ▼                                  │
        atmosphere_mean_molar_mass     ModifiedKeq.__call__  ◄─────┘
                                              │
                                              ▼
                                       OxygenFugacity.__call__
```

The graph is a DAG with a single feedback loop (the outer solver → residual → speciation → modified-equilibrium-constant cycle) and two parallel entry points that diverge at the cold-start generator and re-merge at the residual function. The private versions (`_get_partial_pressures`, `_atmosphere_mass`, `_dissolved_mass`) take `fO2_shift` as a positional argument so both modes share the physics path; the public versions read it from `ddict` for backward compatibility. No module imports any other module's private state; all couplings are through the `pin`/`ddict` dictionaries.

## What is *not* in CALLIOPE

By design:

- **No interior structure** (`structure.py` is a one-line fallback, not a real model).
- **No time stepping** (CALLIOPE is a steady-state equilibrium solver).
- **No I/O** (no NetCDF writes, no file reads, no logging beyond the module-level `log` handle).
- **No threading** (deterministic single-threaded; PROTEUS provides parallelism by sweeping at the outer level).
- **No JAX, no PyTorch** (pure NumPy + SciPy).
- **No state** (every call is functional; no class-level caches).

These omissions are deliberate: CALLIOPE is a calculation kernel, not a framework. Anything that needs framework-level features (logging, status files, retries, parallel sweeps, plotting orchestration) lives in the PROTEUS wrapper.

## Performance

A single `equilibrium_atmosphere()` call takes:

- ~10 ms with a good warm start (`p_guess` from previous iteration);
- ~100 ms-1 s with a cold start (Monte-Carlo guess phase);
- the residual evaluation (`func`) is $\le$1 ms; the bulk of the time is in `fsolve` step iteration.

For batch use cases (sensitivity sweeps, parameter studies), wrap a Python loop around `equilibrium_atmosphere()`; the cost of the function call is small enough that running thousands of cases serially is straightforward. For Monte-Carlo studies with $10^5$+ cases, parallelise at the outer Python level (e.g. `multiprocessing.Pool.map`); CALLIOPE itself is not threaded but is process-safe.

## See also

- [API reference](../Reference/api/index.md) for the auto-generated per-symbol documentation.
- [Source on GitHub](https://github.com/FormingWorlds/CALLIOPE/tree/main/src/calliope) for the actual implementation.

[^cite-bower2019]: D. J. Bower, D. Kitzmann, A. S. Wolf, P. Sanan, C. Dorn, A. V. Oza, *[Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations](https://doi.org/10.1051/0004-6361/201935710)*, Astronomy & Astrophysics, 631, A103, 2019. [SciX](https://scixplorer.org/abs/2019A%26A...631A.103B/abstract).
