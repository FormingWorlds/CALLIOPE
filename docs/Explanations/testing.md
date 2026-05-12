# Testing suite

[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/CALLIOPE?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/CALLIOPE)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/nightly.yml)
[![tests](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/FormingWorlds/CALLIOPE/main/.github/badges/tests-total.json)](https://proteus-framework.org/testing)

Tests verify that the code does what was written; physical correctness is judged by data, not by tests.
The suite catches regressions in equilibrium chemistry, solubility laws, oxygen-fugacity buffers, and the hybrid solver, but it cannot certify that those formulae match nature.
That judgement belongs to the validation runs and the published comparisons against measured magma-ocean outgassing data.

## Internal 4-marker scheme

Every test in the suite carries exactly one of four markers, applied either at module level (`pytestmark = pytest.mark.X`) or per class (`@pytest.mark.X`).

- `unit` (96 tests): in-process tests of individual chemistry helpers, equilibrium-constant fits, solubility laws, and structure formulae. No real `equilibrium_atmosphere` call. Sub-second per test.
- `smoke` (25 tests): real `equilibrium_atmosphere` invocations on minimal configurations (single composition, default species set, one solve). Sub-30 s per test, exercises the hybrid solver code paths.
- `integration` (5 tests): full multi-species CHNS solves with mass-conservation invariants, all eleven species active.
- `slow`: long parameter sweeps and convergence studies. Currently empty.

## Local commands

```console
pytest -m unit                              # 96 fast unit tests
pytest -m smoke                             # 25 minimal-config solver tests
pytest -m integration                       # 5 full CHNS solves
pytest -m slow                              # empty in CALLIOPE
pytest -m "(unit or smoke) and not skip"    # PR-gate selection
pytest -m "not skip"                        # everything that should ever run
```

Coverage:

```console
pytest --cov=src/calliope --cov-report=term -m "not skip"
pytest --cov=src/calliope --cov-report=html -m "not skip"   # htmlcov/
```

## Public-facing badges versus internal taxonomy

Public-facing badges (README, project website) collapse `smoke + integration + slow` into a single `Integration Tests` category, because a 4-way taxonomy is confusing to non-developer readers.
The 4-marker internal scheme remains for CI infrastructure granularity: the PR gate runs `(unit or smoke)`, the nightly runs everything, and the test-count badge fetches the JSON files written by the publish-test-badges workflow.

## Badge system

Three JSON files at `.github/badges/tests-{total,unit,integration}.json` are rewritten by `.github/workflows/publish-test-badges.yml` on every push to `main` (paths-filtered to source, tests, tools, and pyproject.toml).
Shields.io fetches them live via the endpoint URL embedded in the test-count badge.
The publish workflow auto-commits the badges with `[skip ci]` and retries the push up to three times to absorb concurrent main-branch updates.

## Coverage gate

`[tool.coverage.report] fail_under` in `pyproject.toml` sets the minimum combined line + branch coverage for the nightly run.
The PR gate has a pre-flight step that fetches the base branch's `pyproject.toml` and refuses any PR that drops `fail_under` below the value on `main`; both states tolerate a base branch with no `fail_under` declared.
The 90 % ceiling caps the ratchet so that defensive paths do not become unfailable.

## Coverage union estimation

The PR gate downloads the most recent `nightly-coverage` artifact from `main` (`coverage.xml + coverage.json + nightly-timestamp.txt`, 14 d retention) and line-ORs the unit-tier coverage with the nightly's full coverage to produce a union estimate.
The result is written to `$GITHUB_STEP_SUMMARY` as informational output; it does not gate the PR.
A staleness threshold of 48 h and a grace band of 0.3 % apply to the warn / fail / ok decision.

## Canonical specification

The repository-wide rules that every PROTEUS-ecosystem submodule follows are at [proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/](https://proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/).
