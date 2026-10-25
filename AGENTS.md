# CALLIOPE agent instructions

CALLIOPE computes the equilibrium between a magma ocean and its atmosphere: the partial pressures of H2O, CO2, O2, H2, CO, CH4, N2, S2, SO2, H2S, NH3 and the noble gases, the volatiles dissolved in the melt, and optionally a metal reservoir in equilibrium with the melt. PROTEUS calls it in the outgassing step (`outgas.module = 'calliope'`); it also runs standalone for equilibrium-chemistry studies. Before a first edit:

- Tests for `src/calliope/<file>.py` go in `tests/test_<file>.py`; the test rules are in `tests/AGENTS.md`.
- Physics sources: `chemistry.py`, `metal.py`, `oxygen_fugacity.py`, `solubility.py`, `solve.py`, `structure.py`. Utility sources: `__init__.py`, `_version.py`, `constants.py`. The [API reference](docs/Reference/api/index.md) describes each module.
- These commands decide whether a change is ready (CI runs the tests, the structure check, the test-quality lint, the agent-file check and pre-commit):

```bash
pytest -m "(unit or smoke) and not skip"
ruff check --fix src/ tests/ tools/ scripts/ && ruff format src/ tests/ tools/ scripts/
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check
python tools/agents/check_agents_md.py
pre-commit run --all-files
zensical build --clean  # when a change touches docs/ or mkdocs.yml
```

Coverage as the nightly measures it, and one test with its locals on failure:

```bash
coverage run -m pytest -m "(unit or smoke or integration or slow) and not skip" && coverage report && coverage html
pytest -x -v --showlocals tests/test_<file>.py::test_function  # add --pdb for the debugger
```

Entry points: the Python API (`from calliope.solve import equilibrium_atmosphere, equilibrium_atmosphere_authoritative_O`). CALLIOPE has no command-line interface.

<!-- fwl-core:begin sha256=1adfa510a4f117ec -->
## PROTEUS ecosystem rules

PROTEUS couples separate module repositories into one model of the evolution of rocky planets and their atmospheres. Its config selects one module per process; the config values are the lowercase names:

- `interior_energetics.module`: `aragog`, `spider` (or the built-in `boundary`); `interior_struct.module`: `zalmoxis`, `spider`
- `atmos_clim.module`: `agni`, `janus`, both with SOCRATES spectral radiative transfer
- `outgas.module`: `calliope`, `atmodeller`; `atmos_chem.module`: `vulcan`
- `star.module`: `mors`; `escape.module`: `zephyrus`, `boreas`; `orbit.module` (tides): `obliqua`, `lovepy`
- fwl-io downloads the reference data the modules use

A module change to anything PROTEUS calls or reads (a function signature, a config key, an output column, a unit) can break the coupled model while the module's own tests pass. Search `src/proteus/` in PROTEUS for the name before you change it, and name the affected PROTEUS call sites in the pull request.

### Physics and numerics

- Units differ between and within modules: MORS and VULCAN use cgs; CALLIOPE and atmodeller work with pressures in bar; AGNI takes surface pressure in bar at setup and holds Pa inside; PROTEUS states the unit of every `hf_row` key in `GetHelpfileKeys` (`src/proteus/utils/coupler.py`: atmospheric pressures in bar, interior pressures in Pa, model time in years, orbital periods in s), and its config reference pages state each key's unit (for example stellar mass in M_sun, initial partial surface pressures in bar, stellar age in Gyr). State the unit of every physical quantity in its docstring and convert explicitly at the boundary: a unit mismatch between two modules passes the tests of both and shows only in the coupled run.
- A physics test must fail for the most plausible wrong formula: check a conservation law, a bound, a monotonicity, or a published or analytical value, and assert that the wrong result falls outside the tolerance.
- Take each physical constant from one source per repository (in Python `scipy.constants` or the module's constants file). Two retyped values of one constant differ at round-off and hide real differences between code paths.
- Do not loosen a solver tolerance, a conservation check or a clamp to make a run or a test pass. Find the cause; a loosened check also hides the next defect.

### Validation

A new or changed physics routine gets the validation cases named in its confirmed plan: an analytical law, a published figure or table, or a cross-code check. Each case is a `reference_pinned` test that runs the real code on real data, with no mocks, asserts the tolerance, and writes its light and dark plot through the `validation_plot` fixture of proteus-mpl. `pytest -m reference_pinned --validation-plots` makes all plots again. Line plots are SVG, maps are AVIF; each file is at most 60 kB, under `docs/assets/validation/`. Each case has a subsection on its page in the Validation tab, with the reference, the tolerance, the test id and the command. Data read from a published figure comes from 2 independent figlab reads that agree within the stated error, as a CSV with its provenance under `tests/validation/data/`.

### Code

- Keep an inline comment block to 2 lines or fewer, and never more than 4, because a long block drifts from the code it describes. An explanation that needs more goes in the docstring, the commit message or the pull-request description.

### Branches and pull requests

- Work on a feature branch `<initials>/<short-description>`; `main` changes only through a reviewed pull request.
- Fill every section of the repository's pull-request template, where it has one.
- Before you push, run the checks listed at the top of this `AGENTS.md` and in `tests/AGENTS.md`.

### Where knowledge goes

Rules for every contributor go in the `AGENTS.md` files, the reason for a change in its commit message and pull-request description, and the scientific validation of a module in its `docs/Validation/` pages, where the repository has them. Do not add memory, notes or decision-log files to the repository: nobody maintains them, and they go stale.
<!-- fwl-core:end -->

<!-- fwl-voice:begin sha256=b4fe24baeef8b734 -->
### Commit messages and public text

Commit messages, pull-request text, code comments, docstrings, test names, test skip reasons, parametrize ids, log strings that ship with the code, and CI job and step names describe the change and the current state of the code. They name no tool used to write the change, carry no tool-attribution trailer, use no internal plan, phase or work-group labels, and use no em or en dashes (a page range in a citation is the exception).
<!-- fwl-voice:end -->

## Environment

Python 3.12 or 3.13 and git. `pip install -e ".[develop]"` and `pre-commit install -f`; there is nothing to compile. Next to PROTEUS, install CALLIOPE editable from the PROTEUS tree (`pip install -e CALLIOPE/.`) so a change takes effect without a reinstall. Serve the docs with `zensical serve` after `pip install -e '.[docs]'`. `scripts/cross_backend/` compares CALLIOPE with atmodeller outside CI.

Style: lines of at most 96 characters (ruff `line-length`); `snake_case` names, constants in `UPPER_CASE`; standard type hints; a docstring states the physical scenario. The pre-commit hook runs `ruff check --fix` on the files a commit touches; no CI job checks `ruff format`, so run it before a commit.

## CI

- `tests.yaml` (every pull request, Ubuntu and macOS, Python 3.12 and 3.13): `pytest -m "(unit or smoke) and not skip"` with coverage against the fast gate, and on one job `validate_test_structure.sh`, `check_test_quality.py --check` (blocking), the agent-file check, and a guard that rejects a pull request that lowers `[tool.coverage.report].fail_under` below `min(base, 90)` or raises it above 90.
- `code-style.yaml` runs pre-commit on the changed files; `docs.yaml` builds and deploys the docs on a push to `main`, so build them locally when a change touches `docs/` or `mkdocs.yml`.
- `nightly.yml`: all tiers, coverage to Codecov, the full gate enforced, and a job that installs atmodeller and runs `tests/test_cross_code_consistency.py`.
- `publish-test-badges.yml` refreshes the test-count badges on the `badges` branch; `publish.yaml` releases to PyPI on a tag.

## Physics and coupling contract

- Units: temperature in K, partial pressures in bar, masses in kg, mole fractions dimensionless; `fO2_shift_IW` in dex relative to the IW buffer; solubility laws take a partial pressure in bar and return ppmw of the melt; the metal reservoir takes `P_eq` in GPa and `T_eq` in K. Literature fits are in bar, SI bookkeeping is in Pa: confirm the unit at every conversion site.
- Constants (`R_gas`, `M_earth`, `molar_mass`, `ocean_moles`) come from `calliope.constants`; a bare literal in a function body is a defect.
- Bounds on every result: T > 0; total pressure positive; partial pressures non-negative and summing to the total within the solver tolerance; mole fractions summing to 1; solubilities non-negative and finite; `log10(fO2)` finite; equilibrium constants positive (the modified-Keq form returns `log10(Keq)`, so check the sign of the exponentiation); the mantle mass `M_planet - M_core` positive.
- Mass closure: each `<e>_kg_total` is gas plus melt (the metal is in `<e>_kg_core`), and each species total is the sum of its atmosphere and melt masses. A new species updates every total it contributes to, and O carried by H2O, CO2 and SO2 stays in those totals. PROTEUS writes the element totals to `hf_row` and builds its next target from them.
- Solver intermediates: `solve.py` iterates on vectors that can turn complex, NaN or infinite. New solver code checks `np.isreal` and `np.isfinite` before the next step, keeps a numpy `RuntimeWarning` visible in the log, and does not drop an imaginary part without asserting that it is negligible.
- Defaults selected by name (the IW buffer `DEFAULT_FO2_MODEL`, the modified-Keq form, the solubility laws): a change of a default changes the pinned reference values, the docs that quote them (`docs/Explanations/cross_backend_comparison.md`, the tutorials, the `docs/Validation/` pages), the frozen fixtures of `scripts/cross_backend/`, and PROTEUS tests that pin a CALLIOPE value. `git grep` the old name and update each in the same pull request, add a guard against the sibling path to the tests, state the change in the release notes, and run the PROTEUS unit and smoke tiers against the release before PROTEUS raises its CALLIOPE pin.
- PROTEUS (`src/proteus/outgas/calliope.py`) calls `equilibrium_atmosphere`, `equilibrium_atmosphere_authoritative_O`, `get_target_from_params` and `get_target_from_pressures`, and reads `molar_mass` and `ocean_moles`; `src/proteus/outgas/lavatmos.py` uses `OxygenFugacity`. A change of a signature, a return key or a unit comes with the PROTEUS change.
- Authoritative-O mode (`planet.fO2_source = "from_O_budget"` in PROTEUS): the result carries the solved `O_kg_total` and `fO2_shift_derived`, which PROTEUS writes to `hf_row['fO2_shift_IW_derived']` next to the input `fO2_shift_IW`. PROTEUS compares the solved O budget with the user budget `O_kg_user_ic` and stops the run above 50 % divergence. Do not replace `O_kg_total` with a clipped or fallback value without a flag in the result, and do not change the input `fO2_shift_IW`.
- A config dictionary or object that carries user input is not changed after initialisation; use local variables.

## Adding a physics source

1. Create `src/calliope/<file>.py` and `tests/test_<file>.py` with the module-level `pytestmark`, and add the file to `PHYSICS_SOURCES` in `tools/check_test_quality.py`.
2. Add at least one `physics_invariant` test and one `reference_pinned` test.
3. Add `docs/Validation/<file>.md` with a subsection per case (reference, tolerance, test id, command), and list the case in `docs/Validation/index.md` and in the Validation nav of `mkdocs.yml`.

## Review

Check each change against these points and the physics and coupling contract above. `.github/copilot-instructions.md` repeats this checklist for tools that read only that file; change it together with this section.

- Bounds and mass closure hold on every returned result, including the metal reservoir and the authoritative-O path.
- Units at each bar and Pa boundary; constants from `calliope.constants`.
- A changed default comes with its updated reference values, docs and sibling-path guards.
- New solver code guards its intermediates against complex, NaN and infinite values.
- The return keys PROTEUS reads keep their names and units; `O_kg_total` and `fO2_shift_derived` are the solved values.
- Tests follow `tests/AGENTS.md`.
