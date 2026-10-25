# Validation: `src/calliope/solve.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.solve` against published or self-consistency sources.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_solve.py::test_round_trip_self_consistency_at_earth_fiducial` | Cross-implementation cross-check: CALLIOPE forward `equilibrium_atmosphere` vs CALLIOPE inverse `equilibrium_atmosphere_authoritative_O` | `src/calliope/solve.py` (forward and authoritative-O entry points) | Pins the round-trip property at the Earth-fiducial input: forward solve at `fO2 = IW + 2` produces an O budget; the authoritative-O inverse from that budget recovers `fO2_shift_derived` within 0.05 dex. |
| `tests/test_noble_gases.py::test_partition_ratio_matches_henry_versus_hydrostatic_column` | Analytical limit: closed-form Henry-versus-hydrostatic-column partition ratio | `src/calliope/solve.py` (noble gas unknowns and reservoir tallies) | Pins the dissolved-to-atmospheric mass ratio of He, Ar, and Xe after the joint solve against the independently computed ratio of the Henry coefficient to the hydrostatic column coefficient, anchoring the noble gas partitioning to first principles and the Jambon et al. (1986) constants. |

## Re-derivation note

`calliope.solve` exposes two entry points into the equilibrium-chemistry
solver:

- `equilibrium_atmosphere(target, ddict)` (legacy / forward): user
  supplies `fO2_shift_IW` in `ddict`; the solver returns the equilibrium
  partial pressures and per-species kg.
- `equilibrium_atmosphere_authoritative_O(target, ddict)` (authoritative-O
  inverse): user supplies the target `O_kg_total` in `target['O']`;
  the solver inverts to find the `fO2_shift_IW` that matches it, then
  forwards into the same equilibrium engine.

The round-trip is the contract: the inverse must recover the forward-mode
`fO2_shift_IW` from the forward-mode `O_kg_total`. The test pins this at
the Earth-fiducial input (`T_magma = 1800 K`, `Phi = 1.0`, Earth-like
H/C/N/S budget, all volatile species included, `fO2 = IW + 2`) and asserts
the recovered `fO2_shift_derived` matches the input within 0.05 dex.

Anchor type: forward-inverse closure of one engine. Both entry points
share the same equilibrium chemistry. The legacy mode takes
`fO2_shift_IW` as a control variable and returns the implied `O_kg_total`;
the authoritative-O mode treats `fO2_shift_IW` as a fifth unknown and
solves the 5x5 mass-balance system (four partial pressures plus fO2) for
the `fO2_shift_derived` that reproduces the supplied `O_kg_total`. Their
round-trip agreement is the property the test pins. A regression that
broke either the forward O mass-balance or the inverse solve would lose
the round-trip within 0.1 dex; the 0.05 dex envelope catches a
coefficient-only bug.

## Cross-cutting topical test files

`solve.py` is large (>1200 LOC) and its full test surface is split across
several topical cross-cutting files for readability:

- `tests/test_authoritative_O.py` and siblings: authoritative-O entry
  point contract, monotonicity (`test_authoritative_O_monotonicity.py`),
  input validation (`test_authoritative_O_validation.py`).
- `tests/test_equilibrium_paths.py`: forward solver behaviour on
  multi-species compositions.
- `tests/test_partial_species.py`: partial-species (some elements
  excluded) branches.
- `tests/test_stoichiometry.py`: stoichiometric ratios across the
  published reactions.
- `tests/test_targets.py`: target-element-budget computation.
- `tests/test_noble_gases.py`: noble gas (He, Ne, Ar, Kr, Xe)
  partitioning, mass closure, the mean-molar-mass coupling, and the
  target builders.
- `tests/test_invariants.py`: per-element / per-species closure
  invariants.
- `tests/test_invariants_hypothesis.py`: property-based fuzz tests at
  the slow tier.

This is a documented exception to the strict 1:1 source-to-test
mirroring rule for sources >500 LOC where the test topics are
independent enough that consolidation would hurt readability.

## Metal reservoir

With `metal_mass` above zero the element balance of H, C, N and S gains a core term $D_i M_\mathrm{metal} w_i^\mathrm{melt}$, where $D_i$ is the mass-basis metal over melt partition coefficient from the partition file (`src/calliope/metal.py`). `<e>_kg_total` stays gas plus melt and excludes the core; the bulk inventory of an element is `<e>_kg_total + <e>_kg_core`.

### Gas-free core fraction

**What it tests:** with melt and metal alone the core holds $f = Dx/(1 + Dx)$ of an element, $x = M_\mathrm{metal}/M_\mathrm{melt}$ (0.325/0.675).

- Reference: the mass balance of melt and metal (analytical).
- Tolerance: relative 1e-12; $f = Dx$, the law without the melt back-reaction, lies more than 10 % above it from $D = 0.3$ on.
- Test id: `tests/test_metal.py::test_gas_free_core_fraction_curve`, `tests/test_metal.py::TestMetalReservoir` (unit tier)
- Command: `pytest tests/test_metal.py -k "gas_free"`

![Core fraction of an element against D](../assets/validation/metal_gas_free_core_fraction.svg#only-light){ width="80%" }
![Core fraction of an element against D](../assets/validation/metal_gas_free_core_fraction_dark.svg#only-dark){ width="80%" }

### Equal D is a larger melt

**What it tests:** with the same mass-basis $D$ for H, C, N and S, no O or Si in the metal and $\Delta\mathrm{IW}$ a parameter, the solve with metal equals the solve without metal and a melt mass larger by $D M_\mathrm{metal}$, for $D$ = 1 and 30. The gas holds the rest of each budget, so the core mass is $(M - M_\mathrm{gas})\,Dx/(1 + Dx)$.

- Reference: the identity of the element balance (analytical).
- Tolerance: partial pressures and melt plus core masses relative 1e-6, with both solves from the same random state and the root accepted at relative 1e-10 plus 1 kg; core over gas-free core equals $1 - M_\mathrm{gas}/M$ to 2e-5.
- Test id: `tests/test_metal_solve.py::test_metal_with_equal_d_is_a_larger_melt` (smoke tier)
- Command: `pytest tests/test_metal_solve.py`

### Tagawa et al. (2021) H runs

**What it tests:** the placeholder H law, [Tagawa et al. (2021)](https://doi.org/10.1038/s41467-021-22035-0) Eq. 2, against the 6 laser-heated diamond anvil cell runs of their SI Table 1.

- Reference: Tagawa et al. (2021), SI Table 1 (published table).
- Tolerance: each run inside its 1 sigma (largest 0.61 sigma); their Eq. 10 misses every run by more than 2 sigma.
- Test id: `tests/test_metal.py::test_tagawa_runs_lie_within_one_sigma_of_eq2` (unit tier)
- Command: `pytest tests/test_metal.py -k tagawa`

![Measured against fitted log10 K_D of the Tagawa runs](../assets/validation/metal_tagawa2021_h_runs.svg#only-light){ width="80%" }
![Measured against fitted log10 K_D of the Tagawa runs](../assets/validation/metal_tagawa2021_h_runs_dark.svg#only-dark){ width="80%" }

### Tsutsumi et al. (2025) single-stage cores

**What it tests:** the molar D of [Tsutsumi et al. (2025)](https://doi.org/10.1038/s41467-025-65729-5) on the H- and C-free metal basis, converted with the KLB-1 melt of Takahashi (1986), for their 3 single-stage cases (40 GPa and 3750 K, 50 GPa and 3500 K, 54 GPa and 3350 K; 710 ppm H$_2$O and 120 ppm C in the bulk silicate Earth).

- Reference: Tsutsumi et al. (2025), SI Table 1 (published table).
- Tolerance: core H and C equal to the table at its rounding, 0.005 wt%; D times the melt mass fraction without the basis conversion misses the 40 GPa case.
- Test id: `tests/test_metal.py::TestBasis`, `tests/test_metal.py::test_tsutsumi_core_values_plot` (unit tier)
- Command: `pytest tests/test_metal.py -k "tsutsumi or Basis"`

![Core H and C of the Tsutsumi single-stage cases](../assets/validation/metal_tsutsumi2025_core_h_c.svg#only-light){ width="80%" }
![Core H and C of the Tsutsumi single-stage cases](../assets/validation/metal_tsutsumi2025_core_h_c_dark.svg#only-dark){ width="80%" }

## Anchor types

- Self-consistency cross-implementation (round-trip forward vs inverse).
- Future addition: cross-backend cross-check (CALLIOPE vs atmodeller at
  the Earth-fiducial) at the slow tier, once `scripts/cross_backend/`
  fixtures land in `tests/`.

## Cross-references

- `src/calliope/solve.py`: forward and authoritative-O entry points.
- `docs/Explanations/authoritative_oxygen.md`: user-facing concept page
  on the authoritative-O entry point.
- `docs/Explanations/cross_backend_comparison.md`: empirical comparison
  of CALLIOPE Fischer vs atmodeller Hirschmann at the Earth fiducial
  (ΔIW = 0.16 dex residual after the buffer-default flip).
