# Build a new test

[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/CALLIOPE?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/CALLIOPE)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/nightly.yml)
[![tests](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/FormingWorlds/CALLIOPE/main/.github/badges/tests-total.json)](https://proteus-framework.org/testing)

This page is the practical contributor guide for adding or modifying a test in CALLIOPE.
The conceptual framing (marker scheme, badge system, coverage gates, AST linter) lives in the [testing suite](../Explanations/testing.md) explainer; read it first if you have not yet.

## Run the suite locally

From the repository root, with `pip install -e ".[develop]"` already done:

```console
pytest                                      # everything pytest can collect (318)
pytest -m unit                              # 195 fast unit tests
pytest -m smoke                             # 107 minimal-config solver tests
pytest -m integration                       # 5 full CHNS solves
pytest -m slow                              # 11 sweeps and hypothesis fuzz
pytest -m "(unit or smoke) and not skip"    # PR-gate selection
```

A single test file or a single test:

```console
pytest tests/test_chemistry.py
pytest tests/test_chemistry.py::test_modified_keq_janaf_H2_matches_closed_form_at_2000K_with_oneill
```

Stop at first failure and show local variables:

```console
pytest -x --showlocals
```

## Where the test goes

CALLIOPE follows a 1:1 source-to-test mirroring rule: each file in `src/calliope/` has a same-named companion in `tests/`.

- New unit test for a function in `src/calliope/oxygen_fugacity.py` → `tests/test_oxygen_fugacity.py`.
- New unit test for a function in `src/calliope/solubility.py` → `tests/test_solubility.py`.

The cross-cutting test files are the documented exception:

- `tests/test_invariants.py`, `tests/test_invariants_hypothesis.py`: contract clauses that span multiple sources (mass closure end-to-end, partial-species behaviour, stoichiometry of an arbitrary recipe).
- `tests/test_authoritative_O*.py`, `tests/test_equilibrium_paths.py`, `tests/test_partial_species.py`, `tests/test_stoichiometry.py`, `tests/test_targets.py`: solver-architecture tests that span `solve.py`, `chemistry.py`, and `oxygen_fugacity.py` together.

Place a test in a cross-cutting file only when it genuinely cannot be expressed against a single source.

## Module-level marker and timeout

Every test file begins with:

```python
import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]
```

The timeouts per tier:

| Tier | Timeout |
|---|---|
| `unit` | 30 s |
| `smoke` | 60 s |
| `integration` | 300 s |
| `slow` | 3600 s |

Per-function markers (`@pytest.mark.skip` on a stale parametrisation, for example) are additive and do not replace the module-level marker.

Choose the tier from the actual content of the test:

- `unit`: in-process logic, equilibrium-constant fits, solubility laws, structure formulae. No call into `equilibrium_atmosphere`. Should run in under 100 ms.
- `smoke`: real `equilibrium_atmosphere` call on a minimal configuration. Under 30 s.
- `integration`: full multi-species CHNS solve with mass-conservation invariants.
- `slow`: parameter sweeps and convergence studies that exceed one minute wall.

If a slower test would push the PR gate over its 10-minute budget, mark it `integration` or `slow` and rely on the nightly suite.

## Anti-happy-path rules

Every new test function must contain:

1. **At least one edge case**: a boundary value, an empty input, an extreme physical parameter (`T` at the calibration window edges, `p` near zero, mass fractions near 0 and 1).
2. **At least one path that exercises the error contract**: a documented exception, a guard return, or a graceful clamp. If the function under test has no validation, exercise the limit-input behaviour and assert the mathematical invariant ($e = 0$ for an eccentricity-dependent routine, $T \to 0$ for a Boltzmann factor).
3. **Assertion values that are not trivially derivable from the implementation**: discriminating numeric pins, property-based assertions (monotonicity, conservation, symmetry). Avoid point checks at $T = 1$ where $T^n$ is the same for every $n$.

Forbidden patterns flagged by the AST linter:

- A single-assert test function.
- A standalone weak assertion (`assert result is not None`, `assert result > 0`, `assert len(result) > 0`, `assert isinstance(result, dict)`) as the sole meaningful check. A weak assertion *alongside* a strong primary one (a sign guard next to a `pytest.approx` value pin, for example) is allowed and is not flagged.
- A test with no function-level docstring.
- `==` adjacent to a float literal.
- A test asserting on a fixture's implicit default.

The full ruleset and reasoning live in `.github/.claude/rules/calliope-tests.md`.

## The discrimination guard

A pinned numeric value alone does not discriminate the correct formula from the most plausible wrong one.
The discrimination guard is a follow-up assertion that names the wrong formula and shows the gap is larger than the tolerance.

Example from `tests/test_oxygen_fugacity.py`:

```python
def test_oxygen_fugacity_fischer_value_at_2000K_matches_published_fit():
    """At T=2000 K the Fischer 2011 IW buffer evaluates to -7.14981.

    Discrimination: the O'Neill & Eggins 2002 buffer at the same T
    gives -7.4078, which differs from Fischer by 0.26 dex. The pin
    tolerance is 1e-4, so the wrong buffer would not pass.
    """
    fischer_2000 = log10_fO2_IW_fischer_2011(2000.0)
    assert fischer_2000 == pytest.approx(-7.14981, rel=1e-4)

    # Discrimination guard against the most plausible wrong formula.
    oneill_2000 = log10_fO2_IW_oneill_2002(2000.0)
    assert abs(fischer_2000 - oneill_2000) > 0.2
```

For a solubility law, the guard names the alternative law:

```python
# Discrimination guard against Dixon et al. 1995 (Hawaiian basalt):
dixon = 965.0 * (p_bar ** 0.5)
assert abs(sossi - dixon) > 1000.0   # ppmw at p = 100 bar
```

For a stoichiometric closure, the guard names the most plausible off-by-one:

```python
# Discrimination guard: the wrong stoichiometry (forgetting the 0.5 on O2)
# would give 2.5 instead of 1.0 on the LHS.
assert abs(lhs - 1.0) > 0.5
```

## Physics-invariant marker

Tag any test on a physics source that asserts one of the four invariant families with `@pytest.mark.physics_invariant`.
The marker is per-function, not module-level.

The four families with one-line examples:

| Family | Example |
|---|---|
| Conservation | `assert kg_atm + kg_liquid + kg_solid == pytest.approx(kg_total, rel=1e-6)` |
| Positivity / boundedness | `assert 0.0 <= mole_frac <= 1.0 for mole_frac in result.values()` |
| Monotonicity / symmetry | `assert log10_fO2(T=2000) < log10_fO2(T=1500)` along an isobaric buffer |
| Pinned value with discrimination guard | (see the example above) |

The five physics sources (`chemistry.py`, `oxygen_fugacity.py`, `solubility.py`, `solve.py`, `structure.py`) must each carry at least one such test in the matching `tests/test_<source>.py`.
The utility sources (`__init__.py`, `_version.py`, `constants.py`) are exempt from the invariant requirement but still subject to the anti-happy-path rules above.

Structural tests (ordering, autonomy, mutation-in-place, pass-through assignment) in a physics-source test file should not carry the marker.

## Reference-pinned marker

Tag tests that pin behaviour against an external anchor with `@pytest.mark.reference_pinned`.
The anchor is one of:

- a **published benchmark** (cite paper + figure + table in the docstring),
- an **analytical limit** (the ideal-gas limit at low pressure, the Stefan-Boltzmann black-body limit),
- a **cross-implementation check** (CALLIOPE vs atmodeller at a shared fiducial).

Each of the five physics sources carries at least one reference-pinned test.
The current anchor list is in the [testing suite](../Explanations/testing.md#reference-pinned-validation) explainer.

When the first reference-pinned test for a new source lands, create the matching `docs/Validation/<source>.md` page with:

- the source under test,
- the cited paper or analytical limit,
- the closed-form re-derivation of the pinned value (no skipped algebra),
- the wrong-formula or wrong-buffer discrimination number.

Tests carrying `reference_pinned` typically also carry `physics_invariant` (the published-value pin is itself the invariant).

## Optional-dependency imports

Tests that import an optional dependency call `pytest.importorskip('<name>')` at module top, before the import:

```python
import pytest
hypothesis = pytest.importorskip('hypothesis')
from hypothesis import given, strategies as st
```

The optional dependencies recognised by the linter:

- `hypothesis` (property-based fuzz tests),
- `atmodeller` (cross-backend comparison runners).

The PR Docker image installs with `pip install --no-deps`; without `importorskip`, a top-level `import hypothesis` makes the whole test module fail to collect even though the rest of the suite would run.
This trap has recurred multiple times and is now linter-enforced.

## Float comparisons

Never use `==` for floats.
Use `pytest.approx(val, rel=1e-5)` (or `abs=...`) or `numpy.testing.assert_allclose(actual, expected, rtol=..., atol=...)`.

```python
assert fischer_2000 == pytest.approx(-7.14981, rel=1e-4)
np.testing.assert_allclose(result, expected, rtol=1e-6)
```

State the tolerance rationale in a comment when the choice is non-obvious:

```python
# rtol=1e-3 because the Cp lookup truncates to 4 significant figures.
```

## Mocking discipline

Default to `unittest.mock` for all external calls in unit tests: atmodeller, file I/O, network.
Mock at the narrowest scope: a specific function, not a whole module.

```python
from unittest.mock import patch

@patch('calliope.solve.atmodeller_solve_external')
def test_authoritative_O_round_trip(mock_external):
    mock_external.return_value = _physically_plausible_fixture()
    ...
```

A mocked physics function must return physically plausible values; a mock that returns `0.0` or `1.0` for everything can mask real bugs.
Never mock the function under test.
Smoke, integration, and slow tiers use the real solver.

## Module-level constants and `monkeypatch`

When the source under test reads an environment variable into a module-level constant at import time:

```python
DATA_DIR = Path(os.environ.get('CALLIOPE_DATA', ...))
```

`monkeypatch.setenv` is not sufficient because the constant is frozen at the import that already happened.
Patch the constant directly:

```python
monkeypatch.setattr('calliope.solve.DATA_DIR', tmp_path, raising=False)
```

Patch both the env var (for downstream code that re-reads it) and the constant (for code that reads only the constant).

## Voice rule for test artifacts

The repo-wide voice rule applies to test-skip reasons, test-file and function docstrings, test-function and class names, parametrize ids, log-capture assertions, commit messages on test-touching commits, and pull-request titles and bodies on test-touching PRs.
Banned phrases inside those artifacts:

- "audit", "review pass", AI-roadmap labels (`Phase X`, `Stage X.Y`, `Iteration N`),
- AI-tool names ("Generated with Claude"),
- em-dashes, en-dashes (carve-out: bibliographic page ranges in citations).

Write the outcome, never the process.
A test docstring should state the physical scenario the test verifies; it should not narrate how the test was authored.

The rule files themselves (this page, `.github/.claude/rules/calliope-tests.md`, `.github/.claude/rules/calliope-code-review.md`) are out of scope: they may legitimately name the procedures they define.

## Marker validation

`tools/validate_test_structure.sh` runs in the PR gate ahead of pytest and rejects any test file whose tests have no marker.
Run it locally before pushing:

```console
bash tools/validate_test_structure.sh
```

## Test-quality lint

`tools/check_test_quality.py` is an AST linter that walks `tests/test_*.py` and enforces seven rules (single-assert, weak-assertion, missing docstring, float-eq-literal, missing module-level pytestmark, no assertions, missing importorskip).
It runs in two modes:

```console
python tools/check_test_quality.py --check      # CI mode: compare to baseline
python tools/check_test_quality.py --baseline   # regenerate baseline
```

The baseline (`tools/test_quality_baseline.json`) records the current per-rule violation counts and is the floor.
`--check` exits non-zero if any rule's count exceeds the baseline; the CI workflow blocks the PR on that exit code.
Regenerate the baseline only after a deliberate sweep that reduced violations.
The script refuses to regenerate if the new total exceeds the old; override with `CALLIOPE_TEST_QUALITY_ALLOW_REGRESS=1` only when a new rule was added that surfaces pre-existing violations.

Two advisory modes report gaps without failing CI:

```console
python tools/check_test_quality.py --reference-pinned-audit
python tools/check_test_quality.py --physics-invariant-audit
```

The first lists physics sources whose matching `tests/test_<source>.py` has no `@pytest.mark.reference_pinned` test.
The second lists physics-source tests that assert no invariant and are not tagged `@pytest.mark.physics_invariant`.

## Coverage and the ratchet

The PR gate enforces a fast coverage threshold read live from `[tool.calliope.coverage_fast].fail_under` in `pyproject.toml`:

```console
FAST_FAIL_UNDER=$(python -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['tool']['calliope']['coverage_fast']['fail_under'])")
pytest -m "(unit or smoke) and not skip" --cov=calliope --cov-fail-under=${FAST_FAIL_UNDER}
```

The nightly gate enforces `[tool.coverage.report].fail_under` over the full tier set.
Both thresholds ratchet upward, capped at 90 % (`tools/update_coverage_threshold.py` enforces `ECOSYSTEM_CEILING = 90.0`); neither may be manually decreased.

```console
python tools/update_coverage_threshold.py     # one-way ratchet
```

A pre-flight PR step rejects any change that drops `[tool.coverage.report].fail_under` below `min(base, 90.0)`.

## The 50-line review trigger

A pull request that adds or substantially modifies more than 50 lines of test code across all its commits triggers an independent review pass before merge.
The denominator is PR-level (`git diff origin/main...HEAD -- 'tests/**'`); splitting into many sub-50-line commits does not dodge the trigger.
The reviewer cites the anti-happy-path rule, the discrimination-guard requirement, and the physics-invariant tier; flags single-assert tests, weak `is not None` patterns, missing module-level marker, missing `physics_invariant` tag on a physics-source test, and dead tests (tests that pass for the wrong reason).

## Adding a new physics source

When a new `src/calliope/<file>.py` lands:

1. Create the matching `tests/test_<file>.py` with the module-level `pytestmark` (typically `unit` with a 30 s timeout).
2. Write at least one test that asserts one of the four invariant families and tag it `@pytest.mark.physics_invariant`.
3. Plan a `@pytest.mark.reference_pinned` test against a published benchmark, analytical limit, or cross-implementation check.
4. When that test lands, create `docs/Validation/<file>.md` with the anchor, the re-derivation, and the discrimination number. Link it from `mkdocs.yml` if the docs site should surface it.
5. Update `PHYSICS_SOURCES` in `tools/check_test_quality.py` to include the new file name.
6. Run the full PR checks locally:

```console
ruff check --fix src/ tests/
ruff format src/ tests/
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check
python tools/check_test_quality.py --reference-pinned-audit
pytest -m "(unit or smoke) and not skip" --cov=calliope
```

## Linting

CALLIOPE uses [ruff](https://docs.astral.sh/ruff/):

```console
ruff check src/ tests/
ruff format --check src/ tests/
ruff check --fix src/ tests/        # auto-fix
ruff format src/ tests/             # auto-format
```

Both run on every commit via the pre-commit hook (`pre-commit install -f` after `pip install -e ".[develop]"`) and again in the code-style CI workflow.

## See also

- [Testing suite](../Explanations/testing.md): the conceptual framing of the marker scheme, badges, coverage gates, and AST linter.
- `.github/.claude/rules/calliope-tests.md`: the canonical deep-dive for anti-happy-path patterns, the discrimination guard, physics-invariant tiering, buffer-flip propagation, hypothesis seed stability, and solver intermediate-state assertions.
- `.github/.claude/rules/calliope-code-review.md`: the companion review-pass rules.
- [PROTEUS ecosystem testing standard](https://proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/): the repository-wide rules that every PROTEUS-ecosystem submodule follows.
