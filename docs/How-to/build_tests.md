# Build a new test

[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/CALLIOPE?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/CALLIOPE)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/nightly.yml)
[![tests](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/FormingWorlds/CALLIOPE/main/.github/badges/tests-total.json)](https://proteus-framework.org/testing)

This page is for contributors who are adding a new test to CALLIOPE.
For background on the marker scheme and the badge system, see the [testing suite](../Explanations/testing.md) explanation.

## Run the suite locally

From the repository root, with `pip install -e .[develop]` already done:

```console
pytest                                      # everything pytest can collect
pytest -m unit                              # 96 fast unit tests
pytest -m smoke                             # 25 minimal-config solver tests
pytest -m integration                       # 5 full CHNS solves
pytest -m "(unit or smoke) and not skip"    # PR-gate selection
```

A single test file or a single test:

```console
pytest tests/test_core.py
pytest tests/test_stoichiometry.py::TestEquilibriumChemistry::test_SO2_equilibrium
```

To stop at the first failure and show local variables:

```console
pytest -x --showlocals
```

## Choose a marker

Apply exactly one marker per test, either at module level (`pytestmark = pytest.mark.X`) or per class (`@pytest.mark.X`):

- `unit`: in-process tests of an individual helper, equilibrium-constant fit, solubility law, or structure formula. No call into `equilibrium_atmosphere`. Sub-second.
- `smoke`: real `equilibrium_atmosphere` call on a minimal configuration. Sub-30 s.
- `integration`: full multi-species CHNS solve with mass-conservation invariants.
- `slow`: parameter sweeps and convergence studies that exceed 1 minute wall.

The PR gate runs `(unit or smoke) and not skip`. Mark a slower test as `integration` or `slow` if it would push the gate over the 10-minute budget.

## Test-quality rules

Every new test should:

1. Cover at least one **edge case** (boundary in $T$, $f_{\mathrm{O}_2}$, $\Phi$, or $p$);
2. Cover at least one **physically unreasonable** input that must raise (negative pressure, $T \le 0$, mass fraction above 1);
3. Use **discriminating values**: pick inputs where the correct formula gives a different answer than the most plausible wrong formulas (avoid $T = 1$ where $T^n$ is the same for all $n$);
4. Compare against an **analytical expectation**, not against the model output frozen at some point in time;
5. Use `pytest.approx` (or `np.testing.assert_allclose`) for all float comparisons.

The full rule set lives in `.claude/rules/proteus-tests.md` in the PROTEUS repo and applies to every PROTEUS submodule.

## Marker validation

`tools/validate_test_structure.sh` runs in the PR gate ahead of pytest and rejects any test file whose tests have no marker.
Run it locally before pushing:

```console
bash tools/validate_test_structure.sh
```

## Linting

CALLIOPE uses [ruff](https://docs.astral.sh/ruff/):

```console
ruff check .
ruff format --check .
```

Both run on every commit via the pre-commit hook (`pre-commit install` after `pip install -e .[develop]`) and again in the code-style CI workflow.
