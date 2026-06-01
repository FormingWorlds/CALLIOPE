---
description: CALLIOPE-specific code review criteria for the generator-evaluator pattern. Applies domain expertise (chemistry, fO2 buffers, solubility, PROTEUS coupling) to all code review in this repo.
---

# CALLIOPE Code Review Criteria

When reviewing CALLIOPE code (either your own or via code-reviewer agents), apply these domain-specific checks in addition to standard code quality review.

> **Discovery note.** CALLIOPE keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `calliope-tests.md` explicitly. **Before opening any review pass, read both this file and `calliope-tests.md`.**

## Physics plausibility

- Temperature must be positive everywhere (Kelvin). Flag any code path where T could reach zero or go negative.
- Total pressure must be positive; partial pressures must be non-negative and sum to total pressure within solver tolerance.
- Mole fractions must sum to 1.0. Flag any composition-returning function that doesn't enforce or verify normalization.
- Henry's-law solubility outputs must be non-negative. Flag any solubility law that could return a negative or `nan` partial pressure for a physically valid input.
- `log10(fO2)` outputs must be finite (no `nan`, `inf`, or `complex`). Flag any buffer evaluation that does not clip / sanitize before return.
- Equilibrium constants `Keq` must be positive. The modified-Keq formulation in `chemistry.py` returns `log10(Keq)`; verify the exponentiation step has not been accidentally inverted.
- Mantle mass closure: `M_mantle = M_planet - M_core`; flag any path that lets `M_mantle <= 0` reach return.

## Unit convention boundaries

CALLIOPE has a mixed unit convention:

- **Solver inputs**: T in K, P in bar, mass in kg, mole fractions dimensionless.
- **Henry's-law constants**: literature values in `(mol / kg) / bar^n` (variable `n` per fit); the conversion to internal units lives in `solubility.py`.
- **Oxygen fugacity**: `log10(fO2)` (dimensionless) and `fO2_shift_IW` (dex relative to the IW buffer).
- **PROTEUS interop**: the returned dictionary uses kg for species masses, bar for partial pressures, K for temperatures.

When reviewing code that crosses these boundaries (e.g. a new solubility law, a new buffer fit, a new PROTEUS-side caller), verify the unit is correct at each conversion site. The `bar` vs `Pa` boundary is the recurring trap: literature fits are nearly always in bar, internal SI bookkeeping is in Pa.

## Buffer-default flip safety

When the default value of a dispatched-by-name physics path (IW buffer, modified-Keq formulation, solubility law) is changed, the change has fan-out across:

1. Tests pinning a buffer-specific reference value (e.g. `EARTH_VOLATILE_O_REF_KG` was tied to the O'Neill 2002 default and had to update when the default flipped to Fischer 2011).
2. Documentation pages citing the buffer-specific number (the cross-backend comparison, the oxygen-fugacity reference page, the tutorials).
3. PROTEUS-side tests pinning against CALLIOPE outputs.
4. Frozen reference fixtures used by the cross-backend harness.

Required workflow for any buffer / law default flip:

1. `git grep` the old buffer / law name; update every reference value tied to it.
2. Update test discrimination guards (Section 2 rule 4 of `calliope-tests.md`) so the test would fail loudly under the wrong-default regression.
3. Update `docs/Explanations/cross_backend_comparison.md` and any tutorial that quotes a buffer-specific number.
4. Note the flip explicitly in the release notes.

A PR that changes a default value but does not touch the test reference values is a red flag during review.

## Solver intermediate-state types

`solve.py`'s root-finder and equilibrium solvers operate on intermediate vectors that can drift through `complex` / `nan` / `inf` during iteration. The current source has `clip` hardening against these; flag any new solver code path that:

- Does not check `np.isreal` and `np.isfinite` on intermediate state before passing to a downstream step.
- Catches `RuntimeWarning` from numpy without preserving the warning category in a log line, or silences it altogether.
- Coerces a complex intermediate to a real float without first asserting the imaginary part is negligible.

The unit tests of `solve.py` verify the intermediate-state defenses fire (see `calliope-tests.md` Section 16); the code review's job is to make sure the defenses are present before the test asks them to fire.

## PROTEUS coupling patterns

CALLIOPE is called by PROTEUS through the chemistry / outgassing step. Four coupling patterns need explicit care during review.

### 1. Authoritative-O IC reconciliation

Under Path C (`fO2_source = "from_O_budget"`), PROTEUS supplies an authoritative oxygen mass `O_kg_total` and CALLIOPE inverts the IW shift that recovers it. PROTEUS stashes the user-supplied `O_kg_total` as `O_kg_user_ic` in `hf_row` before the first CALLIOPE call; the runtime helper `check_ic_oxygen_budget` compares CALLIOPE's solver-derived `O_kg_total` against the sentinel and hard-fails on >50% divergence.

Any change to `solve.py`'s authoritative-O path (`equilibrium_atmosphere_authoritative_O`) must preserve this contract:

- The return dict MUST contain `O_kg_total` (the solver-recovered value, used by PROTEUS to verify the reconciliation).
- The return dict MUST NOT silently substitute a clipped or fallback value for `O_kg_total` without surfacing the substitution in the return dict (e.g. via a `O_kg_total_clipped` flag).

A regression that breaks reconciliation will surface as a >50% divergence hard-fail in PROTEUS; the test on the CALLIOPE side must catch it before it ships.

### 2. fO2_shift_IW echo-back pattern

PROTEUS sometimes overrides `hf_row['fO2_shift_IW']` (e.g. when iterating Path C) and must restore the original value after the CALLIOPE call. CALLIOPE returns `fO2_shift_IW_derived` as a separate key so the user-supplied and solver-derived values are distinguishable in the helpfile CSV.

Required for any change that touches the chemistry-step return contract:

- The return dict MUST keep `fO2_shift_IW_derived` distinct from `fO2_shift_IW`.
- A new field that overloads the meaning of `fO2_shift_IW` (e.g. silently treating it as a target rather than a user input) is a red flag.

The save/restore pattern lives on the PROTEUS side; CALLIOPE's responsibility is to not silently mutate the input.

### 3. Per-species mass closure across modules

PROTEUS's mass-conservation invariant `M_atm <= M_planet` depends on CALLIOPE summing per-species masses (`H2O_kg_total`, `CO2_kg_total`, etc.) consistently. A regression in CALLIOPE's `_kg_total` aggregation (e.g. forgetting the oxygen contribution from H2O when O is treated as a separate element) would silently violate the PROTEUS-side invariant. The asymmetry between "elements as buffered reservoirs" (chemistry view) and "elements as tracked masses" (PROTEUS view) is the lesson from PROTEUS issue #677.

Flag any chemistry-step code that:

- Adds a new species without updating the corresponding `_kg_total` aggregation.
- Treats oxygen as a separate accounting entity from H2O / CO2 / SO2 in a way that changes the per-species totals.
- Returns a per-species total that is not the sum of `kg_atm + kg_liquid + kg_solid` for that species.

### 4. Buffer-default flip impact on downstream PROTEUS

A buffer-default flip in CALLIOPE (the 2026-05 Fischer-vs-O'Neill flip is the canonical example) changes the IW value silently for any PROTEUS run that uses the CALLIOPE default. PROTEUS-side tests that pin a specific IW value are sensitive to this. The rule:

- Announce buffer-default flips in CALLIOPE's release notes.
- Bump the CALLIOPE version pin in PROTEUS's `pyproject.toml` to the post-flip release.
- Run PROTEUS's unit + smoke suite against the new CALLIOPE before merging the version bump.
- Update any PROTEUS-side hardcoded IW value tied to the old default.

A CALLIOPE PR that flips a default but does not anticipate the PROTEUS-side fallout is a red flag during review.

## Config mutability

`Config` (or any dataclass / attrs object) used to carry user input must not be mutated at runtime after IC. Flag any code that sets `config.X.Y = value` outside of config initialization. Use local variables instead.

## Cross-module constant duplication

Physical constants (`R_gas`, `M_earth`, `R_earth`, `N_avogadro`, `M_O`, `M_H`) are defined in `src/calliope/constants.py`. When reviewing code that uses a physical constant, check that the import is from `calliope.constants` and not re-derived. A new constant introduced as a literal in a body (e.g. `5.4e-26 * T` for a Boltzmann-related coefficient) is a red flag.

## Test marker discipline

Every test file must begin with a module-level `pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]` (unit/30 s, smoke/60 s, integration/300 s, slow/3600 s). Per-function markers are additive but do not replace the module-level marker; CI runs `pytest -m "(unit or smoke) and not skip"` and any file missing the tier marker ships untested.

## Test quality (cross-reference)

Test-content rules (anti-happy-path, discriminating-value guards, physics-invariant tiering, `physics_invariant` / `reference_pinned` certification markers, adversarial-review trigger, mocking discipline, `importorskip` + module-constant-monkeypatch traps, buffer-flip propagation, hypothesis seed stability, solver intermediate-state assertions) live in [`calliope-tests.md`](calliope-tests.md). When reviewing tests, apply both files: this one for marker discipline and review-pass gate, the deep-dive for the content contract.

## Sister rules (cross-link)

- [`.github/copilot-instructions.md`](../../copilot-instructions.md) "Testing Standards" -- high-level rules visible to all readers. Repo-root `CLAUDE.md` is a symlink to this file.
- [`calliope-tests.md`](calliope-tests.md) -- test quality deep-dive; the canonical source for anti-happy-path patterns and the validation certification markers.
