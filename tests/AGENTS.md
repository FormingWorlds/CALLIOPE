# CALLIOPE test instructions

<!-- fwl-tests-core:begin sha256=0d210e346c572a96 -->
## Test rules shared by the PROTEUS ecosystem

Each test file starts with a module-level tier marker and a timeout. CI selects tests by marker, so a file without one runs in no CI job; the timeout stops a hang, it is not a target.

```python
pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]
```

| tier | what it tests | timeout |
|---|---|---|
| `unit` | Python logic with the heavy physics mocked; aim for < 100 ms | 30 s |
| `smoke` | the real binaries or solvers, one step, low resolution | 60 s |
| `integration` | several modules coupled | 300 s |
| `slow` | full physics validation | 3600 s |

Each test carries exactly one tier marker, so that the tier filters select it in the one CI job meant for it; a function marker for a second tier breaks that. Split a file whose tests need different tiers. `skip` excludes a test from every CI job.

Every new test covers an edge case (a boundary value, an empty input, an extreme physical parameter), exercises the error contract (a documented exception with the check that no side effect ran, a guard, a clamp, or the limit input of the formula), and asserts values that do not follow trivially from the implementation. A pinned value of 1 that every exponent reproduces checks nothing.

`python tools/check_test_quality.py --check` fails when the count of any rule rises above `tools/test_quality_baseline.json`; the offenders it prints are the first few of that rule in the tree, not necessarily yours. The rules: a file without a tier marker, a test without a docstring, a test with one assertion or none, a weak assertion as the only one (`is None`, `is not None`, `> 0`, `len(...) > 0`, `isinstance`), `==` next to a non-zero float literal, and an optional dependency imported without `pytest.importorskip`. `bash tools/validate_test_structure.sh` runs the repository's own structure check; the repository part of this file says what it checks.

Physics tests carry markers so their coverage is tracked apart from line coverage:
- `@pytest.mark.physics_invariant` on each test function that asserts a conservation law, a bound (T > 0, fractions in [0, 1]), a monotonicity or symmetry, or a pinned value with a discrimination guard. The marker goes on the function, not the module: structural tests in the same file do not carry it.
- `@pytest.mark.reference_pinned` (together with `physics_invariant`) on a test that pins a published benchmark, an analytical limit or a cross-implementation result; cite the paper, table or figure in the docstring.
- A discrimination guard asserts that the most plausible wrong formula (a missing factor, a swapped exponent, the wrong unit) gives a result outside the tolerance of the pinned value.

Floats: compare with `pytest.approx` or `np.testing.assert_allclose` and a tolerance you can justify from the method, not the one that makes the test pass.

Mocks: mock at the narrowest scope (the one external call), and return physically plausible values, so the code under test runs its real branches. Set random seeds and write files only under `tmp_path`.

A module-level constant read from an environment variable at import time does not change with `monkeypatch.setenv`; patch the constant with `monkeypatch.setattr`.
<!-- fwl-tests-core:end -->

## CALLIOPE specifics

Structure: `src/calliope/<file>.py` is tested in `tests/test_<file>.py`. Three kinds of file do not follow this: cross-cutting tests (`test_invariants_hypothesis.py`, `test_init.py`, `test_element_balance_grid.py`, `test_cross_code_consistency.py`), tests of a repository file or of `tests/conftest.py` (`test_precommit_config.py`, `test_validation_index.py`, `test_conftest_*.py`), and topical files of a physics source above about 500 lines whose topics are independent, for example `test_authoritative_O*.py`, `test_equilibrium_paths.py`, `test_noble_gases.py`, `test_metal_solve.py` and `test_metal_fields.py` for `solve.py`. The primary file `tests/test_<file>.py` still exists and carries at least one `physics_invariant` and one `reference_pinned` test. `bash tools/validate_test_structure.sh` checks that every test carries exactly one of `unit`, `smoke`, `integration`, `slow` and `skip` (module, class or function).

Tiers: `unit` is Python logic with the solver mocked; `smoke` is the real solver on a minimal composition (aim for < 30 s); `integration` is the multi-species CHNS solve; `slow` is full physics validation. A file whose slowest test takes 7.5 s or more in CI uses four times that duration, rounded up to a full minute, as its timeout. pytest-timeout is a required plugin, and a hook in `tests/conftest.py` stops the run when the timeout mark of a test has no positive budget. An autouse fixture seeds the global NumPy generator with 0 before every test, so solver restarts repeat. `docs/How-to/build_tests.md` gives the commands per tier.

### Physics sources and invariants

The physics sources are listed in `AGENTS.md` and in `PHYSICS_SOURCES` of `tools/check_test_quality.py`. Each unit test of a physics source asserts at least one of these; utility sources are exempt from this rule but not from the rest:

- Conservation: per-species closure `kg_atm + kg_liquid == kg_total`, element closure across gas, melt and metal, mole fractions summing to 1.
- Bounds: T > 0, P > 0, partial pressures and solubilities non-negative, mole fractions in [0, 1], `log10(fO2)` finite (no NaN, infinity or complex value).
- Monotonicity or symmetry: `log10(fO2)` along an isobaric buffer decreasing with `1/T`; CO2 solubility rising with P at fixed T; partial pressures doubling with the pressure at fixed mole fractions; outputs unchanged when two non-reacting species swap places.
- A pinned value with a discrimination guard, when a closed-form result or a published table value is the contract.

A property is preferred over a point pin where both are possible. Each case has a section on a page in `docs/Validation/` that names the source, the reference and the test ids; the [Validation overview](../docs/Validation/index.md) lists the cases, and `test_validation_index.py` fails on a pinned test missing from it. `python tools/check_test_quality.py --reference-pinned-status` and `--physics-invariant-status` list the gaps (advisory).

### Edge cases and discriminating values

Edge cases here: `Phi` 0 or 1, T at the solidus, P = 0, `fO2_shift_IW = 0`, `metal_mass = 0`, and a single-species composition, which is a degenerate fixed point of the multi-species solve. A function without input checks (a closed-form relation, an equilibrium constant) is tested at its limit input against the invariant there; the missing check is no exemption.

The plausible bugs here are a wrong exponent or factor, a wrong sign, a unit slip (Pa against bar, K against Celsius, log10 against ln) and a wrong buffer or law (Hirschmann, Fischer or O'Neill, one solubility law or its sibling). A pin with `pytest.approx` adds a guard for each class: a wrong-value distance, an explicit sign (`val < 0` for `log10(fO2)` at IW), a scale window with the unit in the comment (`-10 < val < -3`), and a value that separates the chosen buffer or law from its sibling. When the primary assertion is a conservation closure, the closure covers the exponent class; the sign and scale guards stay.

Discriminating inputs: the IW buffer at 1500 and 2500 K against the Fischer to O'Neill slope ratio; Henry's law with one dominant species and traces, so a swapped constant moves it beyond the tolerance; an equilibrium constant at an off-grid temperature; stoichiometric closure with several species at a non-zero `fO2_shift_IW`.

A single assertion is accepted only for a hard-fail invariant (mass closure at `1e-12`) that no other test in the file checks. Delete a test that asserts a fixture's implicit default.

### Tolerances, mocks, constants

- State the reason for a tolerance in a comment when it is not obvious ("`rtol=1e-3` because the fit gives four significant figures"). Closed-form algebra holds to `rel=1e-12`; solver output needs a looser tolerance.
- Unit tests mock atmodeller, file, HTTP and subprocess access with `unittest.mock`, patching one function (`calliope.solve.<helper>`), not a module. A mock returns physically plausible values, not `0.0` or `1.0` everywhere. Never mock the function under test. Smoke, integration and slow tests use the real solver.
- `hypothesis` and `atmodeller` are the optional dependencies the linter knows (`OPTIONAL_DEPS`, rule `missing_importorskip`); atmodeller is in no extra, so its import is always guarded; a new optional dependency is added there.
- A solve takes its IW buffer from `fO2_model=` or `ddict['fO2_model']`, else `calliope.solve.DEFAULT_FO2_MODEL` read at each call, and passes it to every `OxygenFugacity`, `ModifiedKeq` and `SolubilityS2` it builds; select the buffer of a test with `fO2_model=` or the key. A patch of `calliope.solve.DEFAULT_FO2_MODEL` changes every solve in the process; one of `calliope.oxygen_fugacity.DEFAULT_FO2_MODEL` changes nothing, and one of `__init__.__defaults__` only the direct calls of those 3 classes.
- Test parameters use SI units unless the function takes config units (bar, K, ppmw), named in a comment; parametrize ids name the physical scenario (Earth-like, Mars-like, high fO2). A hypothesis test uses `@settings(derandomize=True)` or a fixed `--hypothesis-seed`, named in its docstring. Write no large outputs.

### Docstrings and names

The file docstring starts with `Tests for src/calliope/<file>.py`, lists the invariants and contract clauses it checks, and links `docs/How-to/build_tests.md`. Each test docstring states the physical scenario or contract clause. A comment explains why an input was chosen ("1500 and 2500 K so the Fischer to O'Neill slope difference is resolved"). Names describe behaviour in snake_case (`test_iw_buffer_monotonic_with_temperature`, not `test_iw_buffer`); a class groups tests that share setup. Test names, skip reasons, parametrize ids and log-capture assertions follow the public-text rule of `AGENTS.md` and do not use the words "audit", "review pass" or "adversarial review"; the pages under `docs/` do not cite the agent instruction files.

### Coverage gates

The fast gate (`[tool.calliope.coverage_fast]`, unit and smoke, every pull request) and the full gate (`[tool.coverage.report]`, all tiers, nightly) move one way toward 90 % (`ECOSYSTEM_CEILING` in `tools/update_coverage_threshold.py`, which raises them: `--target fast` or `--target full`); never lower one by hand. `python tools/check_test_quality.py --baseline` rewrites `tools/test_quality_baseline.json`; run it only after a change that reduced violations.

### Review of test changes

A pull request that adds or changes more than 50 lines under `tests/` (`git diff origin/main...HEAD -- tests/`, over all its commits) gets an independent review of its tests before merge. The reviewer checks the rules above, recomputes a pinned value with the most plausible wrong formula (wrong buffer, unit slip) and confirms the test fails with it, and names the invariant each physics test asserts. Fixes go in a follow-up commit whose subject states the outcome. Recurring faults:

- A reference value tied to a buffer whose test does not name the buffer and has no sibling-buffer guard, so a default change passes silently.
- A solver test that checks only the final output, while a complex or NaN intermediate was coerced to a plausible float; assert `np.isreal` and `np.isfinite` on the intermediates the defence covers.
- A helper that turns a failed lookup into `continue`; assert instead.
- An assertion on a log line where the contract is a value; capture the call argument.
