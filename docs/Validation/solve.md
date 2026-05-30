# Validation: `src/calliope/solve.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.solve` against published or self-consistency sources.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_solve.py::test_round_trip_self_consistency_at_earth_fiducial` | Cross-implementation cross-check: CALLIOPE forward `equilibrium_atmosphere` vs CALLIOPE inverse `equilibrium_atmosphere_authoritative_O` | `src/calliope/solve.py` (forward and authoritative-O entry points) | Pins the round-trip property at the Earth-fiducial input: forward solve at `fO2 = IW + 2` produces an O budget; the authoritative-O inverse from that budget recovers `fO2_shift_derived` within 0.05 dex. |

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
- `tests/test_invariants.py`: per-element / per-species closure
  invariants.
- `tests/test_invariants_hypothesis.py`: property-based fuzz tests at
  the slow tier.

This is a documented exception to the strict 1:1 source-to-test
mirroring rule for sources >500 LOC where the test topics are
independent enough that consolidation would hurt readability.

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
