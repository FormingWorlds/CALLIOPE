# Validation

## What this tab checks

Validation means that we run the code on a case where another source gives the right answer, and check that the code gives that answer within a stated tolerance.
The other source is a published equation, table or measurement, an analytical result, or a second, independent code.
This tab has one page per process of CALLIOPE, in the menu on the left: chemistry, oxygen fugacity, solubility, the equilibrium solve, the metal reservoir and the mantle mass.
The table below lists every case once, so you can find a case and open its page.

**How to read the table.** Each row is one case. "Process" links to the page of the case, "Case" says what is checked, and "Test id" names the test that runs it. "Reference" is the source of the right answer, and "Status" says which test tier runs the case and whether it has a plot.

**How the case pages are laid out.** On each page, every case says why the check exists, what it tests, how we know the right answer and why it matters, with the test details in a box. Where a case has worked numbers, a second box gives them.

Each case is a test that has the `reference_pinned` marker, so the CI fails when the code leaves its tolerance. Each page names the test ids and the references of its process.

This command runs every case:

```console
pytest -m reference_pinned tests
```

## Cases by process

| Process | Case | Test id | Reference | Status |
|---|---|---|---|---|
| [Chemistry](chemistry.md) | H2O = H2 + 0.5 O2 at 2000 K | `tests/test_chemistry.py::test_modified_keq_janaf_H2_matches_closed_form_at_2000K_with_oneill` | NIST-JANAF tables, closed-form constant | unit tier; plot pair |
| [Chemistry](chemistry.md) | H2S and NH3 constants over 700-5000 K | `tests/test_chemistry.py::test_janaf_constant_matches_janaf_rows_over_fit_range` | NIST-JANAF rows | unit tier; plot pair |
| [Chemistry](chemistry.md) | H2S and NH3 fit coefficients | `tests/test_chemistry.py::test_janaf_coefficients_are_least_squares_fit_of_janaf_rows` | least-squares fit of the NIST-JANAF rows | unit tier; plot pair |
| [Chemistry](chemistry.md) | Equilibrium constants and reaction quotients | `tests/test_cross_code_consistency.py::test_calliope_equilibrium_constants_and_reaction_quotients_vs_nasa_glenn` | NASA Glenn thermodynamic data (cross-code check) | integration tier, nightly cross-code job; needs atmodeller; plot pair |
| [Oxygen fugacity](oxygen_fugacity.md#hirschmann-buffer-at-1473-k) | Hirschmann IW buffer at 1473 K and 1 bar | `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_hirschmann_value_at_1473K_matches_published_figure` | Hirschmann (2021), Table 1 and Fig. 5 | unit tier; plot pair |
| [Oxygen fugacity](oxygen_fugacity.md) | Low-pressure slope of the Hirschmann buffer | `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_hirschmann_low_pressure_slope_matches_volume_change` | volume change of Fe + 1/2 O2 = FeO (analytical) | unit tier; plot pair |
| [Oxygen fugacity](oxygen_fugacity.md) | Fischer IW buffer at 2000 K | `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_fischer_value_at_2000K_cross_check_and_regression` | Fischer et al. (2011); O'Neill and Eggins (2002) | unit tier; plot pair |
| [Solubility](solubility.md) | H2O in peridotite | `tests/test_solubility.py::TestSolubilityH2O::test_peridotite_default_matches_sossi_2023_fit` | Sossi et al. (2023) | unit tier; plot pair |
| [Solubility](solubility.md) | S2 with the Earth-mantle FeO | `tests/test_solubility.py::TestSolubilityS2_xFeO::test_default_call_matches_gaillard_2022_earth_mantle_value` | Gaillard et al. (2022) | unit tier; plot pair |
| [Solubility](solubility.md) | Noble gas Henry constants | `tests/test_solubility.py::TestSolubilityNobleGas::test_every_gas_constant_pinned_to_independent_literal`, `::test_helium_constant_matches_jambon_1986` | Jambon et al. (1986) | unit tier; plot pair |
| [Solubility](solubility.md) | Noble gas constants in atmodeller | `tests/test_solubility.py::TestSolubilityNobleGas::test_noble_constants_cross_check_against_atmodeller` | atmodeller (cross-code check) | unit tier; skips without atmodeller; plot pair |
| [Solve](solve.md) | Forward and authoritative-O round trip | `tests/test_solve.py::test_round_trip_self_consistency_at_earth_fiducial` | self-consistency of the 2 entry points | smoke tier; plot pair |
| [Solve](solve.md) | Noble gas partition | `tests/test_noble_gases.py::test_partition_ratio_matches_henry_versus_hydrostatic_column` | Henry constant against the hydrostatic column (analytical) | smoke tier; plot pair |
| [Solve](solve.md) | Gas-only grid of 72 cases | `tests/test_cross_code_consistency.py::test_calliope_vs_atmodeller_72_case_gas_only_grid` | atmodeller (cross-code check) | integration tier, nightly cross-code job; needs atmodeller; plot pair |
| [Metal reservoir](metal.md#gas-free-core-fraction) | Gas-free core fraction | `tests/test_metal.py::test_gas_free_core_fraction_curve`, `tests/test_metal.py::TestMetalReservoir::test_gas_free_core_fraction_is_dx_over_one_plus_dx` | mass balance of melt and metal (analytical) | unit tier; plot pair |
| [Metal reservoir](metal.md#equal-d-is-a-larger-melt) | Equal D is a larger melt | `tests/test_metal_solve.py::test_metal_with_equal_d_is_a_larger_melt`, `::test_equal_d_identity_plot` | identity of the element balance (analytical) | smoke tier; plot pair |
| [Metal reservoir](metal.md#o-budget-with-metal-o) | O budget with metal O | `tests/test_metal_solve.py::test_authoritative_o_budget_excludes_the_metal_o` | round trip of the buffered and the authoritative-O solve | smoke tier |
| [Metal reservoir](metal.md#equal-d-is-a-larger-melt-with-an-o-budget) | Equal D with an O budget | `tests/test_metal_solve.py::test_authoritative_o_with_equal_d_is_a_larger_melt` | identity of the element balance (analytical) | smoke tier |
| [Metal reservoir](metal.md#tagawa-et-al-2021-h-runs) | H runs of Tagawa et al. (2021) | `tests/test_metal.py::test_tagawa_runs_lie_within_one_sigma_of_eq2` | Tagawa et al. (2021), SI Table 1 | unit tier; plot pair |
| [Metal reservoir](metal.md#tsutsumi-et-al-2025-single-stage-cores) | Single-stage cores of Tsutsumi et al. (2025) | `tests/test_metal.py::TestBasis::test_tsutsumi_single_stage_core`, `tests/test_metal.py::test_tsutsumi_core_values_plot` | Tsutsumi et al. (2025), SI Table 1 | unit tier; plot pair |
| [Metal reservoir](metal.md#published-law-values) | Published law values | `tests/test_metal.py::TestLaw::test_tagawa_eq2_values`, `::test_fischer_eq1_single_stage_value`, `::test_suer_lines`, `tests/test_metal.py::TestExchange::test_ideal_exchange_mole_fractions` | Tagawa et al. (2021), Fischer et al. (2020), Suer et al. (2023) | unit tier |
| [Structure](structure.md) | Earth core mass fraction | `tests/test_structure.py::test_calculate_mantle_mass_recovers_wang_2018_earth_core_fraction` | Wang et al. (2018) | unit tier; plot pair |

The CI job of the pull requests runs the unit and smoke tiers. The nightly workflow runs all tiers, and its cross-code job installs atmodeller and runs `tests/test_cross_code_consistency.py`. That file also holds `test_required_env_guards`, which has the marker of the file but checks only the `CHEM_CROSS_CODE_REQUIRED` setting.

## Plots written by the tests

A case with a plot pair shows the tested values against their reference in a light and a dark version. The test writes the pair with the same numbers it asserts, through the `validation_plot` fixture of proteus-mpl. This command writes all pairs to `docs/assets/validation/`:

```console
pytest -m reference_pinned tests --validation-plots
```

Without `--validation-plots` the tests write nothing: each compares the data of its plot with the committed pair, and a missing or stale file fails the test.
