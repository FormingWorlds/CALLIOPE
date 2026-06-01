# Testing suite

[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/CALLIOPE?label=coverage&logo=codecov&color=brightgreen)](https://app.codecov.io/gh/FormingWorlds/CALLIOPE)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/tests.yaml?branch=main&label=Unit%20Tests&color=brightgreen)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/nightly.yml?branch=main&label=Integration%20Tests&color=brightgreen)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/nightly.yml)
[![tests](https://img.shields.io/endpoint?url=https://proteus-framework.org/CALLIOPE/badges/tests-total.json)](https://proteus-framework.org/testing)

Tests verify that the code does what was written; physical correctness is judged by data, not by tests.
The suite catches regressions in equilibrium chemistry, solubility laws, oxygen-fugacity buffers, and the hybrid solver, but it cannot certify that those formulae match nature.
That judgement belongs to the validation runs and the published comparisons against measured magma-ocean outgassing data.

This page describes the suite as a whole.
Contributors writing or modifying tests should read it together with the [Build a new test](../How-to/build_tests.md) how-to.

## Test quality contract

Five layers enforce test rigor across the suite:

1. A four-marker tier scheme (`unit`, `smoke`, `integration`, `slow`) selects what runs in the PR gate versus the nightly.
2. Two validation markers (`physics_invariant`, `reference_pinned`) tag tests that carry physical meaning beyond pure code coverage.
3. A 1:1 mirroring rule pairs every physics source with a same-named test file.
4. An AST linter (`tools/check_test_quality.py`) rejects seven weak-test patterns on every PR.
5. A coverage ratchet capped at 90 % keeps the gate moving upward over time.

Layers 1, 4, and 5 are blocking on PRs.
Layers 2 and 3 are advisory: the linter reports gaps but does not fail the build.

## The four-marker tier scheme

Every test in the suite carries exactly one tier marker, applied either at module level (`pytestmark = pytest.mark.X`) or per class (`@pytest.mark.X`).

| Marker | What it tests | Per-test budget | CI surface |
|---|---|---|---|
| `unit` | Python logic, individual helpers, equilibrium-constant fits, solubility laws, structure formulae. No real `equilibrium_atmosphere` call. | < 100 ms | PR + nightly |
| `smoke` | Real `equilibrium_atmosphere` invocations on minimal configurations (single composition, default species set, one solve). | < 30 s | PR + nightly |
| `integration` | Full multi-species CHNS solves with mass-conservation invariants, all eleven species active. | minutes | Nightly only |
| `slow` | Long parameter sweeps and convergence studies (authoritative-O monotonicity regimes, hypothesis fuzz, cross-buffer property checks). | up to an hour | Nightly only |
| `skip` | Placeholder, deliberately disabled. | n/a | Never |

Live counts per tier are shown by the `tests` badge at the top of this page (the total) and by `pytest -m <tier> --collect-only -q` locally.

Tests without a tier marker are invisible to CI.
The PR gate runs `pytest -m "(unit or smoke) and not skip"`; the nightly runs everything in `(unit or smoke or integration or slow) and not skip`.

## Module-level marker and timeout

Every test file declares its tier and its wall-time ceiling at module top:

```python
pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]
```

| Tier | Timeout |
|---|---|
| `unit` | 30 s |
| `smoke` | 60 s |
| `integration` | 300 s |
| `slow` | 3600 s |

The timeout is a defensive ceiling, not a target.
A unit test that takes 25 s of wall time has either picked the wrong tier or has a leak somewhere; the ceiling catches future regressions that introduce a hang.
Per-function markers (for example `@pytest.mark.skip` on one stale parametrisation) are additive and do not replace the module-level marker.

## Physics-invariant tiering

A unit test on any physics source must assert at least one of four invariant families:

- **Conservation**: per-species mass closure (`kg_atm + kg_liquid + kg_solid ≈ kg_total` within rtol), stoichiometric closure (`sum(mole_fractions) == 1.0`), element closure across the H/C/N/O/S budgets.
- **Positivity or boundedness**: `T > 0`, `P > 0`, mole fractions in `[0, 1]`, `log10(fO2)` finite over the calibration window, partial pressures non-negative.
- **Monotonicity or symmetry**: `log10(fO2)` decreasing with `1/T` along an isobaric buffer; CO$_2$ solubility increasing with $P$ at fixed $T$; swapping two non-reacting species leaving outputs unchanged.
- **Pinned numeric value with a discrimination guard**: a closed-form value or published table entry pinned via `pytest.approx`, plus an explicit follow-up assertion showing that the most plausible wrong formula would differ from the correct one by more than the tolerance.

Tests that meet one or more of these are tagged `@pytest.mark.physics_invariant`.
The marker is per-function, not module-level: structural tests in the same file (for example, an ordering or pass-through-assignment check) should not carry it.

Physics sources in CALLIOPE are the five files `chemistry.py`, `oxygen_fugacity.py`, `solubility.py`, `solve.py`, `structure.py`.
Utility sources (`__init__.py`, `_version.py`, `constants.py`) are exempt from the invariant requirement but remain subject to the anti-happy-path rules below.

## Reference-pinned validation

Tests that pin behaviour against an external anchor are tagged `@pytest.mark.reference_pinned`.
The anchor is one of:

- a **published benchmark** (cite paper + figure + table),
- an **analytical limit** (for example, the ideal-gas limit at low pressure, or the Stefan-Boltzmann black-body limit),
- a **cross-implementation check** (CALLIOPE vs atmodeller at a shared Earth fiducial).

Each of the five physics sources carries at least one reference-pinned test, recorded on the matching `docs/Validation/<file>.md` page.
The current anchors:

| Source | Anchor | Test |
|---|---|---|
| `chemistry.py` | JANAF Thermochemical Tables[^cite-chase1998] (4th ed.), $K_{eq}$ for $\mathrm{H_2O} \to \mathrm{H_2} + 0.5\,\mathrm{O_2}$ at 2000 K with the O'Neill & Eggins 2002[^cite-oneilleggins2002] buffer | `tests/test_chemistry.py::test_modified_keq_janaf_H2_matches_closed_form_at_2000K_with_oneill` |
| `oxygen_fugacity.py` | Fischer et al. 2011[^cite-fischer2011] (EPSL 304, 496) IW buffer at 2000 K | `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_fischer_value_at_2000K_matches_published_fit` |
| `solubility.py` | Sossi et al. 2023[^cite-sossi2023] peridotite H$_2$O fit; Gaillard et al. 2022[^cite-gaillard2022] (EPSL 117255) S$_2$ fit | `tests/test_solubility.py::TestSolubilityH2O::test_peridotite_default_matches_sossi_2023_fit` and `TestSolubilityS2_xFeO::test_default_call_matches_gaillard_2022_earth_mantle_value` |
| `solve.py` | Self-consistency between `equilibrium_atmosphere` (buffered) and `equilibrium_atmosphere_authoritative_O` (the authoritative-O entry point) at the Earth fiducial | `tests/test_solve.py::test_round_trip_self_consistency_at_earth_fiducial` |
| `structure.py` | Wang, Lineweaver & Ireland 2018[^cite-wang2018] (arxiv:1708.08718) Earth core mass fraction 0.325 | `tests/test_structure.py::test_calculate_mantle_mass_recovers_wang_2018_earth_core_fraction` |

The marker is not the same thing as physical correctness: a reference-pinned test certifies that *this implementation* reproduces *that anchor*; it does not certify that the anchor is the right physics for every astrophysical regime.

## One-to-one source-to-test mirroring

Each source file in `src/calliope/` has a same-named companion in `tests/`:

| Source | Test |
|---|---|
| `src/calliope/chemistry.py` | `tests/test_chemistry.py` |
| `src/calliope/oxygen_fugacity.py` | `tests/test_oxygen_fugacity.py` |
| `src/calliope/solubility.py` | `tests/test_solubility.py` |
| `src/calliope/solve.py` | `tests/test_solve.py` |
| `src/calliope/structure.py` | `tests/test_structure.py` |
| `src/calliope/__init__.py` | `tests/test_init.py` |
| `src/calliope/constants.py` | (covered in `tests/test_core.py`; constants are imported across many sources) |

Cross-cutting tests are the documented exception, not the rule:

- `tests/test_invariants.py`, `tests/test_invariants_hypothesis.py`: contract clauses that span multiple sources (mass closure end-to-end, partial-species behaviour, stoichiometry of an arbitrary recipe).
- `tests/test_authoritative_O.py`, `tests/test_authoritative_O_validation.py`, `tests/test_authoritative_O_monotonicity.py`: the authoritative-O entry point, which touches `solve.py`, `chemistry.py`, and `oxygen_fugacity.py` together.
- `tests/test_equilibrium_paths.py`, `tests/test_partial_species.py`, `tests/test_stoichiometry.py`, `tests/test_targets.py`: solver-architecture tests that span the same three files.

When a new physics source is added, its 1:1 test file is created at the same time; the matching `docs/Validation/<file>.md` page is added when the first reference-pinned test for that source lands.

## AST test-quality linter

`tools/check_test_quality.py` walks `tests/test_*.py` as an AST and enforces seven rules:

| Rule | What it flags |
|---|---|
| `missing_module_pytestmark` | Test file with no module-level `pytestmark` (a tier marker is required). |
| `missing_docstring` | Test function with no docstring. |
| `single_assert` | Function with exactly one assertion (anti-happy-path: a single assert is rarely enough to discriminate the correct formula from plausible wrong ones). |
| `no_assertions` | Function with zero assertions (only valid for tests that exercise an exception path with `pytest.raises`). |
| `weak_assert_*` (three sub-rules) | Standalone `assert result is not None`, `assert result > 0`, `assert len(result) > 0`, etc., as the sole meaningful check. A weak assertion *alongside* a strong primary one (the sign guard in a discrimination pattern, for example) is not flagged. |
| `float_eq_literal` | `==` adjacent to a numeric literal in a test body (use `pytest.approx`). |
| `missing_importorskip` | An optional dependency (`hypothesis`, `atmodeller`) imported at module top without a preceding `pytest.importorskip('<name>')`. The PR Docker image uses `pip install --no-deps`; without `importorskip`, collection fails. |

The linter runs in two modes:

- `python tools/check_test_quality.py --baseline` walks the suite and writes the per-rule violation counts to `tools/test_quality_baseline.json`.
  This is the floor.
  Regenerate the baseline only after a deliberate sweep that reduced violations.
- `python tools/check_test_quality.py --check` (CI mode) walks the suite, compares the current counts to the baseline, and exits non-zero if any rule's count exceeds the baseline.
  The CI workflow runs this and blocks the PR on regression.

The baseline ratchets one way: the linter refuses to regenerate it if the new total exceeds the old.
Override with `CALLIOPE_TEST_QUALITY_ALLOW_REGRESS=1` only when a new rule was added that surfaces pre-existing violations.

Two advisory modes report gaps without failing CI:

- `python tools/check_test_quality.py --reference-pinned-status`: lists physics sources whose matching `tests/test_<source>.py` has no `@pytest.mark.reference_pinned` test.
- `python tools/check_test_quality.py --physics-invariant-status`: lists physics-source tests that assert no invariant and are not tagged `@pytest.mark.physics_invariant`.

## Local commands

```console
pytest -m unit                              # fast unit tests
pytest -m smoke                             # minimal-config solver tests
pytest -m integration                       # full multi-species CHNS solves
pytest -m slow                              # sweeps and hypothesis fuzz
pytest -m "(unit or smoke) and not skip"    # PR-gate selection
pytest -m "not skip"                        # everything that should ever run
```

Coverage:

```console
pytest --cov=calliope --cov-report=term -m "not skip"
pytest --cov=calliope --cov-report=html -m "not skip"   # htmlcov/
```

Lint and structure:

```console
bash tools/validate_test_structure.sh         # module-level marker validator
python tools/check_test_quality.py --check    # AST linter against baseline
python tools/check_test_quality.py --reference-pinned-status
python tools/check_test_quality.py --physics-invariant-status
```

## Public-facing badges versus internal taxonomy

Public-facing badges (README, project website) collapse `smoke + integration + slow` into a single `Integration Tests` category, because a four-way taxonomy is confusing to non-developer readers.
The four-marker internal scheme remains for CI infrastructure granularity: the PR gate runs `(unit or smoke)`, the nightly runs everything, and the test-count badge fetches the JSON files written into the documentation site during the docs deploy.

## Badge system

The documentation deploy (`.github/workflows/docs.yaml`) regenerates three JSON files, `tests-{total,unit,integration}.json`, from `pytest --collect-only` and writes them into the published site under `badges/`.
Shields.io fetches them live from the site via the endpoint URL embedded in the test-count badge.
The counts refresh on every documentation build, so they track the suite without running it.

## Coverage gates

Two gates are declared in `pyproject.toml`:

| Gate | Tests included | Threshold key | Where it runs |
|---|---|---|---|
| Fast | `unit + smoke` | `[tool.calliope.coverage_fast].fail_under` | PR `Run unit + smoke tests` step |
| Full | `unit + smoke + integration + slow` | `[tool.coverage.report].fail_under` | Nightly `coverage report` |

Both gates ratchet toward 90 %, capped at 90 % (`tools/update_coverage_threshold.py` enforces `ECOSYSTEM_CEILING = 90.0`); neither may be manually decreased.
The PR gate has a pre-flight step that fetches the base branch's `pyproject.toml` and rejects any PR that drops `[tool.coverage.report].fail_under` below `min(base, 90.0)`.
A one-time ratchet down to the 90 % ceiling is allowed; any drop below 90 % is blocked.

## Coverage union estimation

The PR gate downloads the most recent `nightly-coverage` artifact from `main` (`coverage.xml + coverage.json + nightly-timestamp.txt`, 14-day retention) and line-ORs the unit-tier coverage with the nightly's full coverage to produce a union estimate.
The result is written to `$GITHUB_STEP_SUMMARY` as informational output; it does not gate the PR.
A staleness threshold of 48 hours and a grace band of 0.3 % apply to the warn / fail / ok decision.

## PR validation pipeline

`.github/workflows/tests.yaml` runs on every PR (`if: draft == false`) over a 2-OS by 2-Python matrix (`ubuntu-latest`, `macos-latest` x `3.12`, `3.13`).
The full step sequence:

1. **Validate test markers** (`bash tools/validate_test_structure.sh`): rejects any test file without a module-level `pytestmark`.
2. **Run test-quality lint** (`python tools/check_test_quality.py --check`): blocking; rejects regression against `tools/test_quality_baseline.json`.
3. **Pre-flight fail_under ratchet check**: rejects any PR that drops `[tool.coverage.report].fail_under` below `min(base, 90.0)`.
4. **Run unit + smoke tests**: `pytest -m "(unit or smoke) and not skip" --cov=calliope --cov-fail-under=${FAST_FAIL_UNDER}`, where `FAST_FAIL_UNDER` is read live from `[tool.calliope.coverage_fast].fail_under`.

Steps 1, 2, and 3 are gated to `ubuntu-latest, python 3.12` to avoid four redundant runs; step 4 runs across the full matrix.

Nightly (`.github/workflows/nightly.yml`) runs the full suite, uploads coverage to Codecov, and enforces `--cov-fail-under=90`.

## Canonical specification

The repository-wide rules that every PROTEUS-ecosystem submodule follows are at [proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/](https://proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/).

## References

[^cite-chase1998]: M. W. Chase, *[NIST-JANAF Thermochemical Tables, 4th edition](https://janaf.nist.gov/)*, Journal of Physical and Chemical Reference Data Monograph 9, 1998.
[^cite-fischer2011]: R. A. Fischer, A. J. Campbell, G. A. Shofner, O. T. Lord, P. Dera, V. B. Prakapenka, *[Equation of state and phase diagram of FeO](https://doi.org/10.1016/j.epsl.2011.02.025)*, Earth and Planetary Science Letters, 304, 496-502, 2011. [SciX](https://scixplorer.org/abs/2011E%26PSL.304..496F/abstract).
[^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
[^cite-oneilleggins2002]: H. St. C. O'Neill, S. M. Eggins, *[The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts](https://doi.org/10.1016/S0009-2541(01)00414-4)*, Chemical Geology, 186, 151-181, 2002. [SciX](https://scixplorer.org/abs/2002ChGeo.186..151O/abstract).
[^cite-sossi2023]: P. A. Sossi, P. M. E. Tollan, J. Badro, D. J. Bower, *[Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets](https://doi.org/10.1016/j.epsl.2022.117894)*, Earth and Planetary Science Letters, 601, 117894, 2023. [SciX](https://scixplorer.org/abs/2023E%26PSL.60117894S/abstract).
[^cite-wang2018]: H. S. Wang, C. H. Lineweaver, T. R. Ireland, *[The elemental abundances (with uncertainties) of the most Earth-like planet](https://doi.org/10.1016/j.icarus.2017.08.024)*, Icarus, 299, 460-474, 2018. [SciX](https://scixplorer.org/abs/2018Icar..299..460W/abstract).
