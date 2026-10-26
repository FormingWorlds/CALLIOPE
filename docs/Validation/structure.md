# Validation: `src/calliope/structure.py`

## What this page checks

This page checks `calculate_mantle_mass` in `src/calliope/structure.py`.
It splits the mass of a planet into a core and a mantle, from the planet radius, the planet mass and the core radius as a fraction of the planet radius.
Validation here means that we run the code on a case where another source gives the right answer.
The test then checks that the code gives that answer within a stated tolerance.
Here the right answer is the core mass fraction of the Earth from Wang et al. (2018).
The mantle mass sets how much melt can hold dissolved volatiles in a standalone CALLIOPE run.
The case below says why the check exists, what it tests, how we know the right answer and why it matters, and a box "Test details" ends it.

## Mantle mass of the Earth

**Why this check exists.** The core mass is the core volume times a core density. A core volume with $r^2$ in place of $r^3$, a wrong mass unit (g in place of kg) or a mantle that is not the planet mass minus the core gives the wrong mantle mass.

**What it tests.** For the Earth radius, the Earth mass and a core radius of 0.55 times the planet radius, the code must give a mantle mass of $(1 - 0.325)\,M_\oplus$ = 4.03e24 kg, to a relative 1e-6. It must differ by more than $0.1\,M_\oplus$ from the result of a core volume with $r^2$, and it must lie between 1e24 and 1e25 kg.

**How we know the right answer.** Wang et al. (2018)[^cite-wang2018] give a core mass fraction of the Earth of 32.5 ± 0.3 wt%. The code builds its core density from this fraction and a core radius of 0.55 Earth radii, so for Earth inputs the mantle is $0.675\,M_\oplus$ by construction. The test therefore checks the form of the formula and its scale, not an independent measurement.

**Why it matters.** A standalone CALLIOPE run dissolves volatiles in this mantle mass. PROTEUS does not call `calculate_mantle_mass`: it passes the mantle mass of its own interior model, as the [code architecture](../Explanations/code_architecture.md#structurepy) page explains.

![Mantle mass over Earth mass against core radius fraction, with the values at 0.55](../assets/validation/structure_earth_mantle_mass.svg#only-light){ width="80%" }
![Mantle mass over Earth mass against core radius fraction, with the values at 0.55](../assets/validation/structure_earth_mantle_mass_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The line is `calculate_mantle_mass` for the Earth against the core radius fraction and gives context. The test asserts the marks at 0.55: the code value (circle), the value $1 - 0.325$ (square) and the result with an $r^2$ core volume (cross). A pass means that the circle agrees with the square to a relative 1e-6, too small to see. The cross must lie more than 0.1 above it on this axis, that is more than 0.1 Earth masses.

??? info "Test details"
    - Reference: Wang, Lineweaver and Ireland (2018)[^cite-wang2018], arxiv:1708.08718: Earth core mass fraction 32.5 ± 0.3 wt% ([arxiv:1708.08718](https://arxiv.org/abs/1708.08718)).
    - Tolerance: relative 1e-6; more than $0.1\,M_\oplus$ from the $r^2$ result; $10^{24} < M_\mathrm{mantle} < 10^{25}$ kg.
    - Test id: `tests/test_structure.py::test_calculate_mantle_mass_recovers_wang_2018_earth_core_fraction`
    - Command: `pytest tests/test_structure.py -k wang_2018`
    - Companion test: `test_calculate_mantle_mass_closure_holds_for_earth_like` checks that the mantle and the core are both positive. Its sum check holds by definition, because the test takes the core as the planet mass minus the mantle.

??? note "Worked numbers"
    `calculate_mantle_mass` subtracts the core mass from the planet mass. The core mass is the volume of a sphere of radius `core_frac * R_planet` times a core density from the Earth:

    ```
    core_rho = 3 * earth_fm * M_earth / (4 pi * (earth_fr * R_earth)^3)
    ```

    with `earth_fm = 0.325` and `earth_fr = 0.55`. For `R = R_earth`, `M = M_earth`, `core_frac = 0.55` the mantle is

    ```
    (1 - 0.325) * M_earth = 0.675 * 5.972e24 kg = 4.031e24 kg
    ```

    A core volume with `r^2` in place of `r^3` gives a core of 5.6e17 kg, so the mantle would be almost the whole planet mass, 5.972e24 kg. The scale check `1e24 < mantle < 1e25` also fails on a factor of 1000 from a g/kg slip.

## Related code

- `src/calliope/structure.py` (`calculate_mantle_mass`): cites arxiv:1708.08718 for the Earth core mass fraction 0.325.

## References

 [^cite-wang2018]: H. S. Wang, C. H. Lineweaver, T. R. Ireland, *[The elemental abundances (with uncertainties) of the most Earth-like planet](https://doi.org/10.1016/j.icarus.2017.08.024)*, Icarus, 299, 460-474, 2018. [SciX](https://scixplorer.org/abs/2018Icar..299..460W/abstract).
