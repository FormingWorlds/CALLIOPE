# Testing

[![Tests](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml/badge.svg)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml)
[![Coverage](https://codecov.io/gh/FormingWorlds/CALLIOPE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/FormingWorlds/CALLIOPE/tree/main)

CALLIOPE uses [pytest](https://docs.pytest.org/en/latest/) for automated tests and [coverage.py](https://coverage.readthedocs.io/) for coverage measurement. The badges above reflect the live state of the `main` branch on GitHub Actions and Codecov.

## Running the test suite

From the repository root, with `pip install -e .[develop]` already done:

```console
pytest
```

A passing run prints something like:

```
=========== 50 passed in 15.29s ============
```

To run a single test file or a single test:

```console
pytest tests/test_core.py
pytest tests/test_stoichiometry.py::TestEquilibriumChemistry::test_SO2_equilibrium
```

To stop at the first failure and show local variables:

```console
pytest -x --showlocals
```

## What the suite covers

The `tests/` directory contains six files:

| File | Coverage |
|---|---|
| `test_init.py` | Smoke check that the package imports and exposes `__version__`. |
| `test_core.py` | Per-component tests of the oxygen-fugacity buffer, modified equilibrium constants, the simple structure model, the H$_2$O solubility laws, and the constants table. |
| `test_stoichiometry.py` | Atom-by-atom mass tallies in `atmosphere_mass()`, equilibrium-constant consistency for SO$_2$ / H$_2$S / NH$_3$, CH$_4$ pressure correction, and end-to-end mass conservation through `equilibrium_atmosphere()`. |
| `test_partial_species.py` | Species-exclusion (`is_included=0`) branches of `get_partial_pressures`, `atmosphere_mass`, and `dissolved_mass`; the `SolubilityN2.libourel` Henry's-law mode. |
| `test_targets.py` | The `get_target_from_params` (oceans + CH ratio + ppmw) and `get_target_from_pressures` (round-trip from initial pressures) entry points, including the no-S / no-C / no-N short-circuit branches. |
| `test_equilibrium_paths.py` | `equilibrium_atmosphere` edge paths: warm-start `p_guess`, `opt_solver=False`, `print_result` log emission, `RuntimeError` on convergence failure, `hide_warnings` toggle. |

Every test follows the project [test-quality rules](https://github.com/FormingWorlds/PROTEUS/blob/main/.claude/rules/proteus-tests.md): one edge case + one physically-unreasonable input per class, discriminating values that distinguish the correct formula from plausible wrong formulas, `pytest.approx` for all float comparisons.

## Coverage

To produce a coverage report:

```console
coverage run -m pytest
coverage report          # text summary in the terminal
coverage html            # HTML report under htmlcov/
```

The floor is **95% combined line + branch coverage**, enforced by `[tool.coverage.report] fail_under = 95` in `pyproject.toml`. The local measurement on this branch is 99%; the floor sits a few points below so unrelated refactors that touch a defensive path don't immediately fail CI. 95% is already a high bar and is not intended to be raised further; the goal is regression protection, not chasing 100%.

CI uploads `coverage.xml` to [Codecov](https://app.codecov.io/gh/FormingWorlds/CALLIOPE/tree/main) on every push to `main` (Linux + Python 3.13 only, to avoid double-counting across the OS-and-version matrix). The Codecov badge at the top of this page is rendered live from that report.

## Linting

CALLIOPE uses [ruff](https://docs.astral.sh/ruff/) for both linting and formatting:

```console
ruff check .
ruff format --check .
```

The pre-commit hook (installed by `pre-commit install` after `pip install -e .[develop]`) runs both on every commit; CI also runs both on every PR.

## Adding a new test

When you add a chemistry path, a solubility law, or any function that takes physical inputs, the new test should:

1. Cover at least one **edge case** (boundary in $T$, $f_{\mathrm{O}_2}$, $\Phi$, or $p$);
2. Cover at least one **physically unreasonable** input that must raise (negative pressure, $T \le 0$, mass fraction above 1);
3. Use **discriminating values**: pick inputs where the correct formula gives a different answer than the most plausible wrong formulas (avoid $T = 1$ where $T^n$ is the same for all $n$);
4. Compare against an **analytical expectation**, not against the model output frozen at some point in time.

These rules are enforced informally during code review; they sit in `.claude/rules/proteus-tests.md` in the PROTEUS repo and apply to every PROTEUS submodule including CALLIOPE.

## Next step

If you have changed any equilibrium constant or solubility coefficient, run [`tests/test_stoichiometry.py`](https://github.com/FormingWorlds/CALLIOPE/blob/main/tests/test_stoichiometry.py) and verify that the analytical-vs-code consistency tests still pass before committing.
