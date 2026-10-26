# Validation: `src/calliope/solve.py`

## What this page checks

This page checks the equilibrium solve in `src/calliope/solve.py`.
The solve finds the partial pressures of the gas and the dissolved masses for a given inventory of H, C, N, S and O.
Validation here means that we run the code on a case where another source gives the right answer.
The test then checks that the code gives that answer within a stated tolerance.
On this page the right answer comes from a second CALLIOPE solve that must agree with the first, and from a closed-form result for the noble gases.
If the solve is wrong, every outgassing step of a PROTEUS run gets the wrong atmosphere.
Each case below says why the check exists, what it tests, how we know the right answer and why it matters.
The box "Test details" at the end of each case gives the reference, the tolerance, the test id and the command.
The cases of the solve with a metal reservoir are on the [metal reservoir page](metal.md).

**The 2 ways to call the solve.**
`equilibrium_atmosphere(target, ddict)` takes the oxygen fugacity as input: you give the shift $\Delta\mathrm{IW}$ from the iron-wüstite buffer in `ddict['fO2_shift_IW']`, and the solve returns the partial pressures, the mass of each species and the O budget `O_kg_total` that goes with them.
`equilibrium_atmosphere_authoritative_O(target, ddict)` takes the O budget as input: you give `O_kg_total` in `target['O']`, and the solve finds the $\Delta\mathrm{IW}$ that matches it.
It treats the oxygen fugacity as a fifth unknown and solves the 5 by 5 mass balance of 4 partial pressures and the oxygen fugacity.
The [authoritative oxygen](../Explanations/authoritative_oxygen.md) page explains this second way.

## Round trip between the 2 solves

**Why this check exists.** The 2 calls share one chemistry, but each has its own O mass balance. An error in the O mass balance of either call, or in the inverse solve, makes the 2 calls disagree on the same planet.

**What it tests.** An Earth-like planet: $T$ = 1800 K, a fully molten mantle of 4.03e24 kg, budgets of 1.5e20 kg H, 1.5e19 kg C, 8.0e18 kg N and 8.0e20 kg S, every volatile species included, and the default IW buffer. The solve with $\Delta\mathrm{IW}$ = 2 gives an O budget of 1.2254e21 kg; the test requires it between 1e18 and 1e22 kg. The solve with that O budget must then find $\Delta\mathrm{IW}$ = 2 to within 0.05 dex. It finds 2.00004.

**How we know the right answer.** The second solve starts from the O budget of the first, so the exact answer is the $\Delta\mathrm{IW}$ of the first solve (a check of one code against itself).

**Why it matters.** PROTEUS uses the first call with `planet.fO2_source = "user_constant"` and the second with `"from_O_budget"`. The 2 settings must describe the same planet the same way.

??? info "Test details"
    - Reference: cross-check of the 2 CALLIOPE calls: `equilibrium_atmosphere` against `equilibrium_atmosphere_authoritative_O`.
    - Tolerance: 0.05 dex on the recovered $\Delta\mathrm{IW}$; O budget between 1e18 and 1e22 kg.
    - Test id: `tests/test_solve.py::test_round_trip_self_consistency_at_earth_fiducial`
    - Command: `pytest tests/test_solve.py -k round_trip_self_consistency`

## Noble gases between melt and atmosphere

**Why this check exists.** The solve puts each noble gas into the melt by its Henry constant and into the atmosphere by its partial pressure. A wrong formula for the mass of the gas column, a wrong prefactor of the dissolved mass, or the constant of one gas used for another gives the wrong split.

**What it tests.** An Earth-scale magma ocean ($T$ = 1800 K, fully molten mantle of 4.0e24 kg, gravity 9.81 m s$^{-2}$, radius 6.37e6 m) with H, C, N, O and S and the noble gases. For He, Ar and Xe the ratio of the dissolved mass to the atmospheric mass after the solve must equal the closed-form ratio to a relative 1e-6. For He the ratio must lie between 0 and 10.

**How we know the right answer.** For a noble gas the dissolved mass is $10^{-6}\,M_\mathrm{mantle}\,\Phi\,k\,p$ and the atmospheric mass is $10^5\,(A/g)\,(M/\mu)\,p$. Here $p$ \[bar\] is its partial pressure, $k$ \[ppmw bar$^{-1}$\] its Henry constant of Jambon et al. (1986), $\Phi$ the melt fraction, $A = 4\pi R^2$ the surface area, $g$ the gravity, $M$ its molar mass and $\mu$ the mean molar mass of the gas. The ratio of the 2 does not depend on $p$ (an analytical limit). The value of $k$ itself is checked on the [solubility page](solubility.md#noble-gas-henry-constants).

**Why it matters.** This split decides how much of each noble gas the atmosphere holds, and so how much escape can remove.

??? info "Test details"
    - Reference: analytical limit, the ratio of the Henry coefficient to the coefficient of the hydrostatic column.
    - Tolerance: relative 1e-6 for He, Ar and Xe; $0 < M_\mathrm{liquid}/M_\mathrm{atm} < 10$ for He.
    - Test id: `tests/test_noble_gases.py::test_partition_ratio_matches_henry_versus_hydrostatic_column`
    - Command: `pytest tests/test_noble_gases.py -k partition_ratio_matches_henry`

**Known limit.** No case on this page compares CALLIOPE with atmodeller. The comparison at the Earth case is in `scripts/cross_backend/` and on the [cross-code comparison](../Explanations/cross_backend_comparison.md) page, not in the test suite.

## Other tests of the solve

`solve.py` is a large file, so its tests are split by topic:

- `tests/test_authoritative_O.py` and the files `test_authoritative_O_*.py`: the contract of the second call, its monotonicity and its input checks.
- `tests/test_equilibrium_paths.py`: the first call on compositions with many species.
- `tests/test_partial_species.py`: solves with some elements left out.
- `tests/test_stoichiometry.py`: the stoichiometric ratios of the reactions.
- `tests/test_targets.py`: the element budgets.
- `tests/test_noble_gases.py`: the noble gases (He, Ne, Ar, Kr, Xe): their split, mass closure, their effect on the mean molar mass, and their budgets.
- `tests/test_invariants.py`: closure of every element and species.
- `tests/test_invariants_hypothesis.py`: property-based tests at the slow tier.

This split is an exception to the rule of one test file per source file, for sources above 500 lines whose topics are independent.

## Related code and pages

- `src/calliope/solve.py`: the 2 calls of the solve.
- [Authoritative oxygen](../Explanations/authoritative_oxygen.md): the explanation page of the second call.
- [Cross-code comparison](../Explanations/cross_backend_comparison.md): CALLIOPE with the Fischer buffer against atmodeller with the Hirschmann buffer at the Earth case, with a residual of 0.16 dex in $\Delta\mathrm{IW}$.
