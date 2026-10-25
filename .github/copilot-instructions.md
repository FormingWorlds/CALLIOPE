# CALLIOPE short instructions

The full instructions are in `AGENTS.md` (repository rules and review points) and `tests/AGENTS.md` (test rules). This file repeats what a review or a chat needs when it cannot read them.

## Commands

```bash
pip install -e ".[develop]"
pytest -m "(unit or smoke) and not skip"
ruff check --fix src/ tests/ tools/ scripts/ && ruff format src/ tests/ tools/ scripts/
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check
python tools/agents/check_agents_md.py
pre-commit run --all-files
```

## Review checklist

- T stays positive; partial pressures and solubilities stay non-negative and finite, mole fractions sum to 1, `log10(fO2)` is finite, and the mantle mass stays positive.
- Mass closure: each species total is its atmosphere plus melt mass; each `<e>_kg_total` is gas plus melt, with the metal in `<e>_kg_core`; O carried by H2O, CO2 and SO2 stays in those totals.
- Units: T in K, partial pressures in bar, masses in kg, `fO2_shift_IW` in dex, solubilities in ppmw of the melt, `P_eq` in GPa; constants from `calliope.constants`, never a retyped literal.
- A changed default (IW buffer `DEFAULT_FO2_MODEL`, modified-Keq form, solubility law) comes with its updated reference values, docs and sibling-path test guards.
- New solver code checks `np.isreal` and `np.isfinite` on its intermediates.
- PROTEUS calls `equilibrium_atmosphere`, `equilibrium_atmosphere_authoritative_O`, `get_target_from_params` and `get_target_from_pressures` (`src/proteus/outgas/calliope.py`) and reads `O_kg_total` and `fO2_shift_derived`; a change of a signature, a key or a unit comes with the PROTEUS change.
- Tests: `src/calliope/<file>.py` is tested in `tests/test_<file>.py`; each file has `pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<s>)]`; each test has a docstring, at least 2 assertions, an edge case and the error path; floats are compared with a tolerance; a test that asserts a physical invariant carries `physics_invariant`, a test against a published or analytical value also carries `reference_pinned`, and pinned values have sign, scale and wrong-formula guards.
- Commit messages, pull-request text, comments, test names and CI job names describe the change: no tool attribution, no internal plan labels, no em or en dashes; an inline comment block has at most 4 lines.

This checklist is a subset of the review points of `AGENTS.md` (its Review section and its physics and coupling contract) and `tests/AGENTS.md`; change them together.
