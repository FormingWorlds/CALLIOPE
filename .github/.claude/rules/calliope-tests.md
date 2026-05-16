---
description: CALLIOPE test quality deep-dive. Anti-happy-path patterns, discriminating-value guards, physics-invariant tiering, validation certification markers, adversarial-review trigger. Extends the Testing Standards section in `.github/copilot-instructions.md`.
---

# CALLIOPE Test Quality Rules

This file is the canonical deep-dive on test quality. The high-level summary lives in [`.github/copilot-instructions.md`](../../copilot-instructions.md) under "Testing Standards". The two files MUST stay in sync. If you change one, mirror the change in the other.

> **Discovery note.** CALLIOPE keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `calliope-code-review.md` explicitly so AI tooling and human readers know to load them. **When opening or editing any file under `tests/**` or `src/calliope/**`, read this file first.**

Sister rule files:

- [`.github/copilot-instructions.md`](../../copilot-instructions.md): high-level rules, applied repo-wide.
- [`.github/.claude/rules/calliope-code-review.md`](calliope-code-review.md): review-pass gate, domain-aware code review (buffer-flip propagation, solver intermediate-state types, PROTEUS-coupling patterns). Test-marker discipline lives there too.

CALLIOPE is scientific simulation code and the test suite is held to physics-grade rigor. Tests exist to catch real bugs. A test that asserts the wrong thing, or that passes for the wrong reason, is worse than no test because it generates false confidence. The rules below codify what "real test" means here.

---

## 1. Anti-happy-path rules (every new test)

Every new test function MUST include:

1. **At least one edge case**: a boundary value (Phi = 0 or 1, T = T_solidus, P = 0, fO2_shift_IW = 0), an empty input, or an extreme physical parameter.
2. **At least one path that exercises the error contract**:
   - If the function under test has documented validation (raises on negative T, refuses to dispatch with an unknown buffer name), test that the error fires AND that no side effect ran.
   - If the function has no validation (closed-form mathematics: thermodynamic relations, equilibrium constants), exercise the **limit-input behavior** (single-species composition is a degenerate fixed point of the multi-species solver) and assert the corresponding mathematical invariant.
   - "No validation in source therefore no error test" is not an exemption; the limit-input substitute is.
3. **Assertion values NOT trivially derivable from the implementation**: discriminating numeric pins (see Section 2 below) or property-based assertions (monotonicity, conservation, symmetry, boundedness).

### Forbidden patterns

These are flagged by `tools/check_test_quality.py` and rejected at PR time.

- **Single-assert test functions**. Two or more assertions per test; the second usually pins the invariant the first hand-waves over. Exception: a single assertion of a hard-fail invariant (mass closure within `1e-12`) is acceptable if the test is the only test of that invariant in the file.
- **Weak assertions when they stand alone as the sole meaningful check in the test.** The shapes are:
  - `assert result is not None`
  - `assert result > 0`
  - `assert len(result) > 0`
  - `assert isinstance(result, dict)`
  - `assert result is None` where the function returns `None` implicitly

  Required carve-out: the three-class discrimination guard (Section 2) uses `assert val > 0` as the sign-error guard and `assert lo < val < hi` as the scale-error guard alongside a primary `pytest.approx(...)` pin. Those secondary lines look like weak assertions in isolation; they are NOT flagged when paired with a stronger primary assertion in the same test. The linter applies the carve-out automatically: weak shapes are flagged only when the test has exactly one `assert` statement (`len(asserts) == 1`) and that assertion is itself the weak shape.
- **Tests with no function-level docstring**. The docstring states which physical scenario or contract clause is being verified.
- **`==` adjacent to a float literal**. Use `pytest.approx(val, rel=...)` or `np.testing.assert_allclose(actual, expected, rtol=..., atol=...)`. Comparing two floats with `==` is a known flake source even for "exact" identities like 0.0 (-0.0 vs +0.0, NaN propagation).
- **Tests asserting on a fixture's implicit default**: e.g. `assert fixture_returning_none() is None`. This is trivially true. Delete the test; do not strengthen it by adding more `is None` assertions.

---

## 2. Discriminating test values

The test contract is: a regression that introduces a plausible bug must fail the test. "Plausible bug" means off-by-one exponent, wrong sign, swapped factor of 2, missing factor of pi, dimensionally-wrong unit, **wrong-buffer / wrong-law selection**. Pick input values where the wrong-formula result is far from the correct one.

### Bad / good examples

| Pattern | Bad (any-exponent-passes) | Good (discriminates) |
|---|---|---|
| `log10(fO2) = A/T + B` (IW buffer) | Test at `T = 1500` only (degenerate against any reasonable A) | Test at `T = 1500` AND `T = 2500`; assert the difference matches the Fischer-vs-O'Neill slope ratio |
| Henry's-law solubility | Composition all equal (1/N each, symmetric) | Asymmetric composition (one dominant species + traces) so a swapped Henry constant changes the dominant species more than the test tolerance |
| Equilibrium constant interpolation | Test at grid nodes (interpolation is identity there) | Test at off-grid temperatures where bilinear vs nearest-neighbor differ |
| Stoichiometric closure | One species at unit pressure (the closure is trivial) | Multi-species with non-trivial fO2_shift_IW so each Keq matters |

### Discrimination guard (REQUIRED for pinned-value tests)

When a test pins a numeric value, include explicit assertions that the wrong-formula result would differ from the correct one for **each plausible bug class**. "Each plausible bug class" means at minimum:

1. **Exponent or factor error** (off-by-one exponent, missing factor of 2 / pi). `abs(val - wrong_value)` discriminates.
2. **Sign error** (`-x` vs `+x`). `abs()` hides this; assert the sign explicitly with `val > 0` or `val < 0`.
3. **Unit-conversion error** (Pa vs bar, K vs C, log10 vs ln). Pin the absolute scale with the unit named in the comment.
4. **Wrong-buffer / wrong-law selection** (Fischer vs O'Neill, Dasgupta vs Iacono-Marziano). When the function dispatches by name, the discrimination guard MUST include a value that distinguishes the chosen path from a sibling path.

**Carve-out for conservation-style invariants.** When the primary assertion IS a conservation closure (mass closure, energy balance, sum-equals-total), the equality form `sum(parts) == pytest.approx(total)` already discriminates exponent / factor errors by construction. The exponent guard is satisfied by the conservation equality itself; sign and scale guards remain mandatory.

Canonical pattern:

```python
def test_iw_buffer_fischer_at_3000K_matches_published_value():
    """Pin IW(Fischer 2011) at T = 3000 K against the original Table 1 fit."""
    val = oxygen_fugacity.iw_buffer('fischer', T=3000.0)
    expected = -6.57  # log10(fO2) at the IW buffer, Fischer+2011 Table 1
    assert val == pytest.approx(expected, rel=1e-3)
    # Wrong-buffer discrimination: O'Neill 2002 at the same T gives ~-7.52.
    # A regression that silently dispatches to 'oneill' instead of 'fischer'
    # would land outside the tolerance.
    wrong_oneill = -7.52
    assert abs(val - wrong_oneill) > 0.3
    # Sign guard: log10(fO2) at IW is always negative under standard conditions.
    assert val < 0
    # Scale guard: order of magnitude is -7, not -70 (forgotten log10) or
    # -0.7 (factor-10 unit slip). Pin the magnitude.
    assert -10 < val < -3
```

The guard lines are mandatory whenever the test's primary assertion is a `pytest.approx` against a hand-calculated or published value. Property-based assertions (monotonicity, conservation, symmetry) do not need a separate guard because they are already discriminating across the input space.

---

## 3. Physics-invariant assertions (tiered)

### When required

Every unit test on a **physics module** must assert at least one of the four invariants below. Physics modules are:

```
src/calliope/chemistry.py
src/calliope/oxygen_fugacity.py
src/calliope/solubility.py
src/calliope/solve.py
src/calliope/structure.py
```

Per-source-file granularity: each of the five physics files needs at least one `@pytest.mark.physics_invariant` test and at least one `@pytest.mark.reference_pinned` test in `tests/test_<file>.py`. Granularity is per source file, not per directory.

Utility modules are exempt from the physics-invariant requirement but still subject to all anti-happy-path rules:

```
src/calliope/__init__.py     (re-exports)
src/calliope/_version.py     (auto-generated by setuptools-scm)
src/calliope/constants.py    (pure physical constants, no derivation)
```

### The four invariant families

1. **Conservation**
   - Mass closure: `sum(species_kg_atm + species_kg_liquid + species_kg_solid) ≈ species_kg_total` per species.
   - Element closure: per-element mass balance across the gas / melt / solid partitioning.
   - Stoichiometric closure: `sum(mole_fractions) == 1.0` for any returned composition vector.
2. **Positivity / boundedness**
   - `T > 0` Kelvin everywhere, `P > 0` Pa everywhere.
   - Partial pressures non-negative; mole fractions in `[0, 1]`.
   - Henry's-law solubility non-negative; mantle melt fraction in `[0, 1]`.
   - `log10(fO2)` finite (no `nan` / `inf` / `complex`).
3. **Monotonicity or symmetry**
   - `log10(fO2)` decreasing with `1/T` along an isobaric buffer.
   - CO2 solubility increasing with pressure at fixed T (Henry's law in the linear regime).
   - Doubling pressure at fixed mole fractions doubles partial pressures.
   - Swapping two non-reacting species in the input list leaves all other outputs unchanged.
4. **Pinned numeric value with a discrimination guard**: see Section 2. This is acceptable as the sole invariant when a closed-form result or published table value is the contract.

Property-based assertions (monotonicity, conservation, symmetry, boundedness) are preferred over point-value pins when both are possible. They hold for any valid input and so catch bugs across the entire input space.

### Validation certification markers

Two markers track validation quality independently of line coverage:

- **`@pytest.mark.physics_invariant`** -- this test asserts at least one of the four invariants. Tag every qualifying test in a physics-source test file. `tools/check_test_quality.py` warns when a physics-source test asserts no invariant and is not tagged.
- **`@pytest.mark.reference_pinned`** -- this test pins behavior against a **published benchmark** (paper, figure, table; cite explicitly in the test docstring), an **analytical limit** (Henry's-law linear regime, single-species degenerate solve, IW buffer at a tabulated reference T), or a **cross-implementation cross-check** (CALLIOPE vs atmodeller at the Earth fiducial; see `docs/Explanations/cross_backend_comparison.md`).
  - **Per-source-file**: each of the five physics source files must have at least one `reference_pinned` test in `tests/test_<file>.py`. Anchor type is one of {published benchmark, analytical limit, cross-implementation cross-check}; the specific paper or limit is chosen by the test author and recorded in `docs/Validation/<file>.md`.
  - **Tracking**: each physics source gets a page at `docs/Validation/<file>.md`, created when the first reference_pinned test for that source lands. The page records: the source under test, the reference cited, the test ids carrying the marker, and the date of last comparison against the source.
  - **Audit**: `python tools/check_test_quality.py --reference-pinned-audit` reports the physics source files missing a `reference_pinned` test. This is the punch list for follow-up validation work.

Both markers are registered in `pyproject.toml` under `[tool.pytest.ini_options] markers`. They do not gate CI on their own; their coverage is a separate KPI surfaced in the PR summary comment.

---

## 4. Mocking discipline

- Default to `unittest.mock` for ALL external calls in unit tests: atmodeller cross-backend calls, file I/O, HTTP, subprocess.
- Mock at the **narrowest scope**: patch the specific function (`unittest.mock.patch('calliope.solve.some_helper')`), not the whole module.
- A mocked physics function MUST return **physically plausible** values. A mock that returns `0.0` or `1.0` for everything will mask sign / clamp / fallback bugs.
- NEVER mock the function under test. If you're tempted to, the test is asking the wrong question.
- Smoke tests use the real CALLIOPE solver on minimal compositions; integration and slow tests use the full multi-species CHNS solver. The rules in this file still apply to those tiers, but the mocking discipline is relaxed because the real call is the contract.

---

## 5. Optional-dependency imports

Any test that imports an optional dependency MUST call `pytest.importorskip` at module top so `pip install --no-deps` CI runs do not fail collection:

```python
import pytest

hypothesis = pytest.importorskip('hypothesis')
# ... or for a module-level helper that requires the dep:
pytest.importorskip('atmodeller')
```

Optional deps recognized by the linter (`OPTIONAL_DEPS` constant in `tools/check_test_quality.py`):

- `hypothesis` (used in property-based / fuzz tests; lives in `[develop]` extras).
- `atmodeller` (used in cross-backend comparison work; lives in neither `[dependencies]` nor `[develop]` extras, must always be guarded if imported into a test file).

The lint script enforces this. Rule key `missing_importorskip`: any module-top `import <optional_dep>` or `from <optional_dep> import ...` that is not preceded by a module-scope `pytest.importorskip('<optional_dep>')` is flagged.

---

## 6. Module-level constants and `monkeypatch`

When the source under test reads an env var or a class-level default into a **module-level constant** at import time, `monkeypatch.setenv` alone is not sufficient. The constant is frozen at the original import.

Pattern:

```python
# Source: src/calliope/oxygen_fugacity.py
DEFAULT_BUFFER = 'fischer'  # frozen at import after the 2026-05 default flip
```

```python
# Test (wrong):
monkeypatch.setattr('os.environ', {'CALLIOPE_BUFFER': 'oneill'})   # too late

# Test (right):
monkeypatch.setattr('calliope.oxygen_fugacity.DEFAULT_BUFFER', 'oneill', raising=False)
```

A related pattern is the `_with_calliope_buffer` context manager in `scripts/cross_backend/runners.py`: it mutates `OxygenFugacity.__init__.__defaults__` and restores on exit. The pattern is **single-threaded only**; never apply it inside library code (it has process-wide visibility and races under pytest-xdist).

When in doubt, do both the env-var monkeypatch and the constant monkeypatch. The lint script does NOT currently flag this pattern (it would require source-side analysis to know which constants are env-derived); this is a discipline rule enforced via the >50 LOC review trigger and the recurring-trap table in Section 16.

---

## 7. Marker discipline and timeouts

### Module-level marker is mandatory

Every test file MUST begin with:

```python
import pytest

pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]
```

with budgets:

- `unit` -> `timeout(30)` (target wall-time per test is `< 100 ms`; the 30 s cap is a defensive net).
- `smoke` -> `timeout(60)` (target `< 30 s`).
- `integration` -> `timeout(300)`.
- `slow` -> `timeout(3600)`.

PR CI runs `pytest -m "(unit or smoke) and not skip"`. Tests without the tier marker are invisible to CI and shipped untested. The lint script blocks any file missing the module-level `pytestmark`.

### Per-function markers

Per-function `@pytest.mark.<tier>` markers are **additive**, not a replacement for the module-level marker. They are useful when a file's tests span multiple tiers (rare; prefer separate files).

### Timeout is a safety net, not a target

The `timeout` ceiling exists so a future regression that introduces a hang (real solver call, infinite loop, network retry) surfaces as a specific-test failure rather than a generic job timeout. Current test wall times are 100x below the ceiling; if you find yourself needing the full 30 s for a unit test, something has gone wrong and you should reduce scope or move the test to a slower tier.

---

## 8. Float and numerical comparison

- NEVER use `==` for floats. Use `pytest.approx(val, rel=1e-5)` or `np.testing.assert_allclose(actual, expected, rtol=..., atol=...)`.
- State the tolerance rationale in a comment when the choice is non-obvious. E.g. "`rtol=1e-3` because the Fischer 2011 fit reports four significant figures".
- For pinned numeric values, include a **discrimination guard** (Section 2).
- For property-based assertions, use `pytest.approx` against the exact symbolic relation, with the tightest tolerance the implementation can hit (typically `rel=1e-12` for closed-form algebra; looser for solver outputs).

---

## 9. Voice rule for test artifacts

The repo-wide voice rule (zero AI-process disclosure in any public artifact) applies to test code with the same strictness as to source. The voice rule is **scoped** to public artifacts other contributors and external readers see; it does NOT apply to the rule documents themselves (this file, `calliope-code-review.md`, `copilot-instructions.md`), which legitimately name the procedures they define.

In scope (the voice rule is BANNED here):

- Test-skip reasons (`@pytest.mark.skip(reason='...')`).
- Test-file docstrings.
- Test-function and test-class names.
- Test-function docstrings.
- Parametrize ids (`@pytest.mark.parametrize('name', [...], ids=[...])`).
- Log-capture assertions (regex against `caplog.records`).
- Commit messages on test-touching commits (subject AND body).
- **Pull-request titles and bodies on test-touching PRs**.
- GitHub Actions job names and step names that ship to the PR Checks tab.
- Inline source comments and docstrings on `src/calliope/**`.
- Log strings that ship with the repo.
- **All public-facing documentation** (anything under `docs/`, the repo README, CONTRIBUTING.md, tutorials, wiki pages). Public docs apply the rule silently; they do NOT enumerate the banned phrases, name the voice rule, advertise the existence of `.github/.claude/` rule infrastructure, or cross-reference `.github/.claude/rules/*.md` files. A user docs page that describes the testing infrastructure must do so without naming the AI-process rules that produced it.

Out of scope (these may NAME the procedures they define):

- This file (`calliope-tests.md`).
- `calliope-code-review.md`.
- `copilot-instructions.md`.

Banned phrases inside the in-scope artifacts: "audit", "review pass", "adversarial review", "Phase X" (when "X" is an AI-organized roadmap label, not a real project phase), "T1.x", "Group A/B/C/D" (when AI-organized work groups), `claude-config/...` paths, "Generated with Claude", AI-tool names, em-dashes, en-dashes (except in bibliographic page ranges within citations), process meta-commentary ("after careful analysis").

Write the OUTCOME (what the test verifies; what the PR achieves) never the PROCESS (how the rule was derived; which review caught what). First-person Tim voice. Going-forward only, no history rewrite.

---

## 10. Fixture and parameter conventions

- Use SI units in test parameters unless the function under test explicitly expects config units (bar, K, ppmw).
- Use `@pytest.mark.parametrize` when the same logic spans multiple physical regimes (Earth-like, Mars-like, sub-Neptune, high-fO2, low-fO2). Each parametrize id must read like a physical scenario, not a tuple of numbers.
- Set seeds for any randomness:
  ```python
  np.random.seed(42)
  random.seed(42)
  ```
  Hypothesis tests use `@settings(derandomize=True)` or an explicit `--hypothesis-seed` to keep replays stable across versions (see Section 16 trap).
- Use `tmp_path` (pytest fixture) for temporary files. Do not produce large outputs in the test path.

---

## 11. Documentation per test

- **File-level docstring**: name the source file under test (`Tests for src/calliope/<file>.py`), list the invariants and contract clauses the file exercises, link to `docs/How-to/build_tests.md`. Required.
- **Function-level docstring**: state the physical scenario or contract clause in plain language. Required (lint-enforced).
- **Inline comments**: explain **why** a specific input range was chosen ("T=1500 K and T=2500 K so the Fischer-vs-O'Neill slope difference is resolved well above tolerance").

---

## 12. Naming

- Test names describe behavior, not the called function: `test_iw_buffer_monotonic_with_temperature`, NOT `test_iw_buffer`.
- Test names use snake_case and read as full sentences.
- Group related tests in classes (`class TestIWBuffer:`) when they share setup; use the class to thread a single fixture through several scenarios.
- Test file names mirror source 1:1: `src/calliope/<file>.py` -> `tests/test_<file>.py`. Two documented exceptions to the 1:1 rule:
  - **Cross-cutting fuzz / init tests** (`test_invariants_hypothesis.py`, `test_init.py`): tests that span multiple source files or test package-level concerns.
  - **Topical sub-files of a large physics source**: when a physics source exceeds ~500 LOC and its tests split into independent topics that would not benefit from consolidation, topical sub-files are acceptable alongside the primary `tests/test_<file>.py`. The primary file must still exist and carry at least one `reference_pinned` and one `physics_invariant` test; the topical sub-files cover the remaining surface. The current exemption is `solve.py` (>1200 LOC), whose tests split across `test_authoritative_O.py`, `test_equilibrium_paths.py`, `test_partial_species.py`, `test_stoichiometry.py`, `test_targets.py`, `test_invariants.py`. The primary `tests/test_solve.py` carries the round-trip self-consistency anchor.

---

## 13. Adversarial review trigger

A pull request that adds or substantially modifies **> 50 lines of test code across all its commits** triggers an independent review pass before merge. This is a discipline rule, not CI-automated: the author runs the review pass via a `code-reviewer` agent before pushing the final test-touching commit. The denominator is PR-level, not per-commit: `git diff origin/main...HEAD -- 'tests/**'` is the source of truth. Splitting one large change into 49 + 49 + 49 line commits does NOT dodge the trigger.

The reviewer's mandate:

- Cite the anti-happy-path rule (Section 1) and the discrimination-guard requirement (Section 2).
- Flag single-assert tests, weak `is not None` patterns, missing module-level marker, missing `physics_invariant` tag on a physics-source test, missing `reference_pinned` tag on a per-source benchmark test, dead tests (passes for the wrong reason), tests that mock the function under test.
- Verify discriminating values: re-derive the expected value from a plausible wrong formula and assert the test fails with that wrong formula.
- Verify physics-source coverage of the four invariants: which of the four does this test exercise? If none, why is the test in `tests/test_<physics_file>.py`?

The reviewer is a separate process from the test author. For Claude-Code workflow this means spawning a `proteus-review` skill or a `code-reviewer` agent with the test files in scope; the review must complete and surface findings before the test commit is pushed.

The reviewer's findings are addressed in a follow-up commit (not amended into the test commit). The follow-up subject line is in plain language describing the OUTCOME ("sharpen IW-buffer assertions to distinguish Fischer from O'Neill", NOT "address review findings").

---

## 14. Tooling

The repo provides:

- `bash tools/validate_test_structure.sh` -- structural check (marker presence, file naming).
- `python tools/check_test_quality.py --check` -- CI mode: AST scan for the forbidden patterns in Section 1 and the marker requirement in Section 7. Fails the PR if violations exceed the baseline.
- `python tools/check_test_quality.py --baseline` -- after a deliberate sweep, regenerates `tools/test_quality_baseline.json`. Only run when you have intentionally reduced violations.
- `python tools/check_test_quality.py --reference-pinned-audit` -- prints physics source files missing a `reference_pinned` test.
- `python tools/update_coverage_threshold.py` -- ratchet the fast PR gate upward when measured coverage exceeds the current `fail_under`. Capped at the 90% ecosystem ceiling.
- `ruff check src/ tests/` and `ruff format src/ tests/` -- run before commit.

The lint script is wired into PR CI (`tests.yaml`). The step runs in **blocking** mode: any regression above the baseline fails the PR.

---

## 15. Coverage strategy (operator's view)

CALLIOPE uses two coverage gates with explicit sub-targets. The fast gate is for PR cycle time; the full nightly gate is the long-running KPI.

| Gate | Tests | Target | When |
|---|---|---|---|
| Fast gate (`tool.calliope.coverage_fast`) | unit + smoke | ratcheting toward **90%** (PROTEUS-ecosystem ceiling) | Every PR |
| Full gate (`tool.coverage.report`) | unit + smoke + integration + slow | **90%** | Nightly |

The ratchet is one-way (`tools/update_coverage_threshold.py`), capped at 90%. Never manually decrease the threshold. The CI guard in `.github/workflows/tests.yaml` rejects any PR that lowers `[tool.coverage.report].fail_under` below `min(base_ref, 90.0)`.

What this means for adding tests:

- A new closed-form helper in a utility module: a unit test is sufficient.
- A new function in a physics source: a unit test (counts toward both gates), plus a `physics_invariant` tag if it qualifies. If the function feeds a published benchmark, a `reference_pinned` test goes with it (counts toward the per-source-file inventory in `docs/Validation/<file>.md`).
- A new cross-backend comparison: a slow-tier test that calls atmodeller via `pytest.importorskip` and pins against a `scripts/cross_backend/` fixture.

---

## 16. Failure modes to recognize on review

These are real patterns that have shipped in the past. The lint script catches some of them mechanically; reviewers catch the rest.

| Pattern | Example | Why it slipped | Fix |
|---|---|---|---|
| **Buffer-default flip propagation** | A test pins `EARTH_VOLATILE_O_REF_KG = 1.241e22` against the legacy O'Neill default. The default flips to Fischer in `oxygen_fugacity.py`; the constant becomes 1.260e22. The test still passes against an outdated reference, silently. | The reference value is tied to the buffer choice but the test doesn't name the buffer in its docstring or assert against a buffer-discriminating value. | Discrimination guard (Section 2 rule 4) must include the alternative-buffer value so a regression that silently dispatches to the wrong buffer fails the test. Cite the buffer name AND the reference paper in the test docstring. |
| **Hypothesis seed and version stability** | `@given(...)` test passes on hypothesis 6.0 with the default seed strategy; on hypothesis 6.100 the strategy produces a different sequence and the test surfaces a previously-hidden flake or stops covering the previously-hidden bug class. | Hypothesis seed semantics changed between minor releases; the test author relied on implicit determinism. | Add `pytest.importorskip('hypothesis')` at module top. Use `@settings(derandomize=True)` or pass `--hypothesis-seed=<fixed>` in CI. Document the chosen seed in the test docstring. |
| **`solve.py` intermediate-state types** | Solver loop produces a `complex` or `nan` intermediate; the silently-coerced final output is a real float that looks plausible but is wrong. The unit test only checks the final output. | The output check is too late; the bug lives in an intermediate step. | Tests of `solve.py` must assert that intermediate state is real-valued at each step: `assert np.all(np.isreal(intermediate))` and `assert not np.any(np.isnan(intermediate))`. The solver-side defense (`clip` hardening against complex / NaN / inf) belongs in source, but tests verify the defense actually fires. |
| **Silent skip in helper** | `def _enum_for(field): ...; if actual is None: continue` masks broken introspection | Helper hides a real failure as a no-op | Hard assertion: `assert actual is not None, ...` |
| **Log-line-only assertion** | Test captures a log line and asserts on its text; a regression that changes the code path but keeps the log still passes | Logs are not the contract | Capture the call kwarg and assert on the value passed in |
| **Module-level constant patched only via env var** | `monkeypatch.setenv('CALLIOPE_BUFFER', ...)` on a source that read it at import time | Constants are frozen at import; setenv is too late | `monkeypatch.setattr('mod.CONST', ...)` in addition to setenv |
| **Optional dep imported unconditionally** | `import hypothesis` at module top | `pip install --no-deps` build skips the optional install | `pytest.importorskip('hypothesis')` at module top |
| **Stale marker after refactor** | File renamed without re-applying the module-level `pytestmark` | CI marker filter passed because of per-function markers; coverage tier became invisible | Restore module-level `pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]` |
| **Trivially-true on implicit None** | `def fixture(): pass`; `def test_x(fixture): assert fixture is None` | Fixture returned None implicitly; test passes for the wrong reason | Delete the test |

When you spot a new variant of these, add it here.

---

## 17. Sister rules (cross-link)

- `.github/copilot-instructions.md` "Testing Standards" -- the high-level summary readers without `tests/**` context see first.
- `.github/.claude/rules/calliope-code-review.md` "Test marker discipline" -- the review-pass gate that backs up the rules in this file. Also contains domain-aware physics checks (buffer-flip propagation, solver intermediate-state types, PROTEUS-coupling patterns) that apply when reviewing the **source** code that tests cover.

Any change to the rule set: update both files in the same commit and call out the cross-reference in the commit body.
