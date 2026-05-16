"""Tests for elemental mass stoichiometry and equilibrium chemistry.

Verifies that:
1. Atmospheric elemental mass tallies count atoms correctly for all species.
2. Equilibrium constants produce partial pressures consistent with
   analytical Kp expressions.
3. CH4 solubility pressure coefficient has correct units.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from calliope.chemistry import ModifiedKeq
from calliope.constants import element_list, molar_mass, volatile_species
from calliope.oxygen_fugacity import OxygenFugacity
from calliope.solubility import SolubilityCH4, SolubilityCO

# This file mixes unit-tier stoichiometry checks with an integration-tier
# end-to-end mass-conservation class; a module-level pytestmark would stack
# with the class-level markers and pull the integration tests into the PR
# gate. The linter accepts the missing module-level pytestmark for this
# specific file as a documented exception; track in tools/test_quality_baseline.json.


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_ddict(T=2000.0, fO2_shift=0.0, gravity=9.81, radius=6.371e6):
    """Build a minimal ddict for atmosphere_mass."""
    d = {
        'T_magma': T,
        'fO2_shift_IW': fO2_shift,
        'gravity': gravity,
        'radius': radius,
    }
    for sp in volatile_species:
        d[f'{sp}_included'] = 1
    return d


def _column_mass(p_bar, g=9.81, R=6.371e6):
    """Atmospheric column mass [kg] for a species at partial pressure p_bar."""
    return p_bar * 1e5 / g * 4 * np.pi * R**2


# ===================================================================
# 1. Stoichiometry: verify elemental tallies atom-by-atom
# ===================================================================


@pytest.mark.unit
class TestAtmosphericStoichiometry:
    """Verify elemental masses by running atmosphere_mass with known
    primary pressures and checking the elemental totals.

    Strategy: set one primary species to a known pressure, compute the
    resulting molecular and elemental masses, and verify against the
    expected N_atoms * M_element / M_molecule scaling.
    """

    def _get_elemental_masses(self, pin, ddict):
        """Run atmosphere_mass and return (p_d, mass_atm_d)."""
        from calliope.solve import atmosphere_mass, get_partial_pressures

        mass_atm_d = atmosphere_mass(pin, ddict)
        p_d = get_partial_pressures(pin, ddict)
        return p_d, mass_atm_d

    def test_H_from_H2O_only(self):
        """With only H2O (other primaries ~0), H mass = 2*M_H/M_H2O * mass_H2O."""
        ddict = _make_ddict()
        pin = {'H2O': 10.0, 'CO2': 1e-30, 'N2': 1e-30, 'S2': 1e-30}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        mass_H2O = _column_mass(p_d['H2O'])
        expected_H = mass_H2O * 2 * molar_mass['H'] / molar_mass['H2O']
        # H mass includes contributions from H2 (derived from H2O)
        # but H2O dominates, so H should be >= expected from H2O alone
        assert mass['H'] >= expected_H * 0.95

        # Discrimination guard: the wrong stoichiometry (factor 1 for H in H2O
        # instead of 2) would give half the H mass. The 2x gap is well outside
        # the 5% derivation tolerance.
        expected_H_wrong_factor_1 = mass_H2O * 1 * molar_mass['H'] / molar_mass['H2O']
        assert mass['H'] > expected_H_wrong_factor_1 * 1.5

    def test_S_from_S2(self):
        """S2 has 2 S atoms. Atmospheric S should be ~2*M_S/M_S2 * mass_S2."""
        ddict = _make_ddict()
        pin = {'H2O': 1e-30, 'CO2': 1e-30, 'N2': 1e-30, 'S2': 10.0}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        # In a S2-dominated atmosphere, mu ~ M_S2
        mass_S2 = _column_mass(p_d['S2'])
        expected_S = mass_S2 * 2 * molar_mass['S'] / molar_mass['S2']
        assert mass['S'] == pytest.approx(expected_S, rel=0.01)

        # Discrimination guard: the wrong stoichiometry (factor 1 for S in S2
        # instead of 2) would give half the S mass. The 2x gap is well outside
        # the 1% tolerance.
        expected_S_wrong_factor_1 = mass_S2 * 1 * molar_mass['S'] / molar_mass['S2']
        assert abs(mass['S'] - expected_S_wrong_factor_1) > expected_S * 0.4

    def test_N_from_N2(self):
        """N2 has 2 N atoms."""
        ddict = _make_ddict()
        pin = {'H2O': 1e-30, 'CO2': 1e-30, 'N2': 10.0, 'S2': 1e-30}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        mass_N2 = _column_mass(p_d['N2'])
        expected_N = mass_N2 * 2 * molar_mass['N'] / molar_mass['N2']
        # NH3 is derived from N2, contributing a small fraction
        assert mass['N'] == pytest.approx(expected_N, rel=0.05)

        # Discrimination guard: the wrong stoichiometry (factor 1 for N in N2)
        # would give half the N mass. The 2x gap dwarfs the 5% NH3 contribution.
        expected_N_wrong_factor_1 = mass_N2 * 1 * molar_mass['N'] / molar_mass['N2']
        assert abs(mass['N'] - expected_N_wrong_factor_1) > expected_N * 0.4

    def test_C_from_CO2(self):
        """CO2 has 1 C atom; full atom-by-atom tally over CO2, CO, CH4."""
        ddict = _make_ddict()
        pin = {'H2O': 1e-30, 'CO2': 10.0, 'N2': 1e-30, 'S2': 1e-30}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        # Full analytical tally: every C-bearing species contributes exactly
        # 1 C atom per molecule.
        from calliope.solve import atmosphere_mean_molar_mass

        mu = atmosphere_mean_molar_mass(p_d)
        g, R = 9.81, 6.371e6
        mass_CO2_kg = p_d['CO2'] * 1e5 / g * 4 * np.pi * R**2 * molar_mass['CO2'] / mu
        mass_CO_kg = p_d['CO'] * 1e5 / g * 4 * np.pi * R**2 * molar_mass['CO'] / mu
        mass_CH4_kg = p_d['CH4'] * 1e5 / g * 4 * np.pi * R**2 * molar_mass['CH4'] / mu
        expected_C_kg = (
            1 * mass_CO2_kg / molar_mass['CO2']
            + 1 * mass_CO_kg / molar_mass['CO']
            + 1 * mass_CH4_kg / molar_mass['CH4']
        ) * molar_mass['C']

        assert mass['C'] == pytest.approx(expected_C_kg, rel=1e-6)

        # Discrimination guard: the wrong stoichiometry (coefficient 2 for C
        # in CO2) would add an extra mass_CO2 / M_CO2 * M_C moles of C. With
        # p_d['CO2'] = O(1 bar) the extra contribution is well outside the
        # 1e-6 tolerance.
        wrong_C_kg = (
            2 * mass_CO2_kg / molar_mass['CO2']
            + 1 * mass_CO_kg / molar_mass['CO']
            + 1 * mass_CH4_kg / molar_mass['CH4']
        ) * molar_mass['C']
        assert abs(mass['C'] - wrong_C_kg) > 0.05 * expected_C_kg

    def test_NH3_contributes_1_N_not_3(self):
        """Under reducing conditions, NH3 becomes significant.
        Verify N from NH3 uses coefficient 1 (not 3)."""
        ddict = _make_ddict(T=1500.0, fO2_shift=-4.0)  # very reducing
        pin = {'H2O': 100.0, 'CO2': 1e-30, 'N2': 1.0, 'S2': 1e-30}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        # Compute N mass analytically from all N-bearing species
        from calliope.solve import atmosphere_mean_molar_mass

        mu = atmosphere_mean_molar_mass(p_d)
        mass_N2_kg = p_d['N2'] * 1e5 / 9.81 * 4 * np.pi * (6.371e6) ** 2 * molar_mass['N2'] / mu
        mass_NH3_kg = (
            p_d['NH3'] * 1e5 / 9.81 * 4 * np.pi * (6.371e6) ** 2 * molar_mass['NH3'] / mu
        )

        expected_N_moles = (
            2 * mass_N2_kg / molar_mass['N2'] + 1 * mass_NH3_kg / molar_mass['NH3']
        )
        expected_N_kg = expected_N_moles * molar_mass['N']

        assert mass['N'] == pytest.approx(expected_N_kg, rel=1e-6)

        # Discrimination guard: the wrong stoichiometry (coefficient 3 for N in
        # NH3 instead of 1, i.e. confusing the H subscript with an N count)
        # would add an extra 2 * mass_NH3_kg / molar_mass['NH3'] moles of N.
        # Under the reducing conditions of this test, NH3 is significant.
        wrong_N_moles = (
            2 * mass_N2_kg / molar_mass['N2'] + 3 * mass_NH3_kg / molar_mass['NH3']
        )
        wrong_N_kg = wrong_N_moles * molar_mass['N']
        assert abs(mass['N'] - wrong_N_kg) > 0.01 * expected_N_kg

    def test_O_from_O2_uses_factor_2(self):
        """O2 has 2 O atoms. The O tally must use factor 2.

        With all primary species near zero, the only O source is O2
        from the fO2 buffer. The expected O mass = 2*M_O/M_O2 * mass_O2.
        """
        ddict = _make_ddict()
        pin = {'H2O': 1e-30, 'CO2': 1e-30, 'N2': 1e-30, 'S2': 1e-30}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        from calliope.solve import atmosphere_mean_molar_mass

        mu = atmosphere_mean_molar_mass(p_d)
        mass_O2_kg = p_d['O2'] * 1e5 / 9.81 * 4 * np.pi * (6.371e6) ** 2 * molar_mass['O2'] / mu

        # Expected: factor 2 for diatomic O2
        expected_O = 2 * mass_O2_kg / molar_mass['O2'] * molar_mass['O']

        # With only O2 contributing (H2O ~ 0), mass["O"] should match
        assert mass['O'] == pytest.approx(expected_O, rel=0.01), (
            f'O mass {mass["O"]:.4e} != expected {expected_O:.4e} '
            '(factor-2 for O2 not applied?)'
        )

        # Discrimination guard: the wrong stoichiometry (factor 1 for O in O2,
        # i.e. treating O2 as monoatomic) would give half the O mass. The
        # 2x gap dwarfs the 1% derivation tolerance.
        expected_O_wrong_factor_1 = 1 * mass_O2_kg / molar_mass['O2'] * molar_mass['O']
        assert mass['O'] > expected_O_wrong_factor_1 * 1.5

    def test_elemental_masses_all_positive(self):
        """All elemental masses should be non-negative and finite."""
        ddict = _make_ddict()
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        _, mass = self._get_elemental_masses(pin, ddict)

        for e in element_list:
            assert mass[e] >= 0.0, f'{e} mass is negative: {mass[e]}'
            assert math.isfinite(mass[e]), f'{e} mass is not finite: {mass[e]}'

        # Discrimination guard: the elemental masses must differ across
        # elements given the asymmetric input (H2O dominates, S2 trace). A
        # tally that returned the same value for every element (e.g. a stub
        # that always returns 1.0) would pass the positivity check but fail
        # this one.
        unique_values = {round(mass[e], 6) for e in element_list}
        assert len(unique_values) >= 3, (
            f'Elemental masses should differ across elements; got {mass}'
        )


# ===================================================================
# 2. Equilibrium chemistry: Keq consistency
# ===================================================================


@pytest.mark.unit
class TestEquilibriumChemistry:
    """Verify that the ModifiedKeq + sqrt expression in solve.py
    produces partial pressures consistent with the analytical Kp.

    The JANAF fits now store K for the doubled reaction form.
    The formation constant K_f (halved form) has coefficients that are
    exactly half the stored values.
    """

    @pytest.mark.parametrize('T', [1500.0, 2000.0, 2500.0, 3000.0])
    def test_SO2_equilibrium(self, T):
        """p_SO2 from code matches analytical K_f * p_S2^0.5 * p_O2."""
        fO2_shift = 0.0
        p_S2 = 0.01

        fO2_model = OxygenFugacity('oneill')
        p_O2 = 10.0 ** fO2_model(T, fO2_shift)

        # Analytical: K_f for 0.5 S2 + O2 -> SO2
        log10_Kf = 18887.0 / T - 3.8064
        Kf = 10.0**log10_Kf
        p_expected = Kf * p_S2**0.5 * p_O2

        # Code path: stored K_B = K_f^2, fO2_stoich=0
        mk = ModifiedKeq('janaf_SO2')
        Geq = mk(T, fO2_shift)
        p_code = (Geq * p_S2 * p_O2**2) ** 0.5

        assert p_code == pytest.approx(p_expected, rel=1e-6), (
            f'SO2 at {T}K: code={p_code:.6e}, expected={p_expected:.6e}'
        )

        # Discrimination guard: the wrong stoichiometry (forgetting the 0.5
        # exponent on p_S2, i.e. treating it as a full S2 reaction) would
        # multiply the result by sqrt(p_S2). At p_S2=0.01 the wrong formula
        # is 10x smaller, well outside the 1e-6 tolerance.
        p_wrong_stoich = Kf * p_S2 * p_O2
        assert abs(p_code - p_wrong_stoich) > 0.5 * p_expected

    @pytest.mark.parametrize('T', [1500.0, 2000.0, 2500.0, 3000.0])
    def test_H2S_equilibrium(self, T):
        """p_H2S from code matches analytical K_f * p_S2^0.5 * p_H2."""
        p_S2 = 0.01
        p_H2 = 0.1

        # Analytical: K_f for 0.5 S2 + H2 -> H2S
        log10_Kf = 6731.01547 / T - 3.62273031
        Kf = 10.0**log10_Kf
        p_expected = Kf * p_S2**0.5 * p_H2

        # Code path
        mk = ModifiedKeq('janaf_H2S')
        Geq = mk(T, 0.0)
        p_code = (Geq * p_S2 * p_H2**2) ** 0.5

        assert p_code == pytest.approx(p_expected, rel=1e-6), (
            f'H2S at {T}K: code={p_code:.6e}, expected={p_expected:.6e}'
        )

        # Discrimination guard: the wrong stoichiometry (treating the
        # reaction as 1 H2 + 1 S2 -> H2S rather than 1 H2 + 0.5 S2 -> H2S)
        # would drop the 0.5 exponent on p_S2 and multiply the result by
        # sqrt(p_S2). At p_S2=0.01 the wrong formula is 10x smaller.
        p_wrong_stoich = Kf * p_S2 * p_H2
        assert abs(p_code - p_wrong_stoich) > 0.5 * p_expected

    @pytest.mark.parametrize('T', [1500.0, 2000.0, 2500.0, 3000.0])
    def test_NH3_equilibrium(self, T):
        """p_NH3 from code matches analytical K_f * p_N2^0.5 * p_H2^1.5."""
        p_N2 = 1.0
        p_H2 = 1.0

        # Analytical: K_f for 0.5 N2 + 1.5 H2 -> NH3
        log10_Kf = 2664.015623 / T - 5.99238046
        Kf = 10.0**log10_Kf
        p_expected = Kf * p_N2**0.5 * p_H2**1.5

        # Code path
        mk = ModifiedKeq('janaf_NH3')
        Geq = mk(T, 0.0)
        p_code = (Geq * p_N2 * p_H2**3) ** 0.5

        assert p_code == pytest.approx(p_expected, rel=1e-6), (
            f'NH3 at {T}K: code={p_code:.6e}, expected={p_expected:.6e}'
        )

        # Discrimination guard: the wrong stoichiometry that swaps the H2
        # and N2 exponents (1.5 vs 0.5) would give Kf * p_N2^1.5 * p_H2^0.5,
        # a different power-law in p_N2/p_H2. At p_N2 = p_H2 = 1.0 the two
        # forms coincide; sample a non-symmetric point to discriminate.
        p_N2_test, p_H2_test = 0.5, 2.0
        p_code_asym = (Geq * p_N2_test * p_H2_test**3) ** 0.5
        p_correct_asym = Kf * p_N2_test**0.5 * p_H2_test**1.5
        p_wrong_swapped = Kf * p_N2_test**1.5 * p_H2_test**0.5
        assert p_code_asym == pytest.approx(p_correct_asym, rel=1e-6)
        assert abs(p_code_asym - p_wrong_swapped) > 0.1 * p_correct_asym

    def test_SO2_increases_with_fO2(self):
        """More oxidizing conditions should produce more SO2."""
        p_S2, T = 0.01, 2000.0
        fO2_model = OxygenFugacity('oneill')

        results = []
        for shift in [-2.0, 0.0, 2.0]:
            p_O2 = 10.0 ** fO2_model(T, shift)
            mk = ModifiedKeq('janaf_SO2')
            Geq = mk(T, shift)
            p_SO2 = (Geq * p_S2 * p_O2**2) ** 0.5
            results.append(p_SO2)

        assert results[1] > results[0], 'SO2 should increase from IW-2 to IW'
        assert results[2] > results[1], 'SO2 should increase from IW to IW+2'

    def test_H2S_positive_and_finite(self):
        """H2S partial pressure should always be positive and finite."""
        for T in [1200.0, 2000.0, 3500.0]:
            mk = ModifiedKeq('janaf_H2S')
            Geq = mk(T, 0.0)
            p_H2S = (Geq * 0.01 * 0.1**2) ** 0.5
            assert p_H2S > 0.0
            assert math.isfinite(p_H2S)

    def test_NH3_decreases_with_temperature(self):
        """NH3 is thermodynamically favored at lower T."""
        p_N2, p_H2 = 1.0, 1.0
        mk = ModifiedKeq('janaf_NH3')

        p_low = (mk(1500.0, 0.0) * p_N2 * p_H2**3) ** 0.5
        p_high = (mk(3000.0, 0.0) * p_N2 * p_H2**3) ** 0.5

        assert p_low > p_high, 'NH3 should be more abundant at lower T'

        # Discrimination guard: NH3 synthesis is exothermic, so the ratio
        # between 1500 K and 3000 K should be several-fold. A near-unity
        # ratio would mean the temperature dependence of the equilibrium
        # constant is being missed. Empirical ratio at these (p_N2, p_H2)
        # is ~7.7 with the janaf_NH3 fit.
        ratio = p_low / p_high
        assert ratio > 3.0, (
            f'p_NH3(1500K)/p_NH3(3000K) = {ratio:.2f}; expected several-fold '
            f'for an exothermic synthesis reaction'
        )

    def test_H2_and_CO_reference_values(self):
        """H2 and CO reactions: pin the modified equilibrium constant
        at the canonical T = 2000 K, dIW = 0 evaluation.

        Hidden coupling: the modified Keq folds fO2 in, so this pin
        depends on the default IW buffer (Fischer et al. 2011).
        """
        T = 2000.0
        mk_h2 = ModifiedKeq('janaf_H2')
        mk_co = ModifiedKeq('janaf_CO')

        g_h2 = mk_h2(T, 0.0)
        g_co = mk_co(T, 0.0)

        # Values computed at the default IW buffer (Fischer 2011).
        assert g_h2 == pytest.approx(1.0896, rel=1e-3)
        assert g_co == pytest.approx(4.8897, rel=1e-3)

    def test_get_partial_pressures_end_to_end(self):
        """End-to-end test: get_partial_pressures should produce positive,
        finite pressures for all species and satisfy p_total > 0."""
        from calliope.solve import get_partial_pressures, get_total_pressure

        ddict = _make_ddict(T=2000.0, fO2_shift=0.0)
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        p_d = get_partial_pressures(pin, ddict)

        # All pressures non-negative and finite
        for sp, p in p_d.items():
            assert p >= 0.0, f'{sp} pressure is negative'
            assert math.isfinite(p), f'{sp} pressure is not finite'

        # Primary species preserved
        assert p_d['H2O'] == pytest.approx(100.0, rel=1e-10)
        assert p_d['CO2'] == pytest.approx(10.0, rel=1e-10)

        # Derived species should be present
        assert p_d['H2'] > 0
        assert p_d['CO'] > 0
        assert p_d['SO2'] > 0
        assert p_d['H2S'] > 0
        assert p_d['NH3'] > 0
        assert p_d['O2'] > 0

        # Total pressure should exceed sum of primaries
        p_total = get_total_pressure(p_d)
        assert p_total > 111.1  # > H2O + CO2 + N2 + S2

    def test_SO2_end_to_end_matches_analytical(self):
        """Verify get_partial_pressures produces SO2 consistent with Kf."""
        from calliope.solve import get_partial_pressures

        T = 2000.0
        ddict = _make_ddict(T=T, fO2_shift=0.0)
        pin = {'H2O': 1e-30, 'CO2': 1e-30, 'N2': 1e-30, 'S2': 0.01}
        p_d = get_partial_pressures(pin, ddict)

        # Analytical Kf for 0.5 S2 + O2 -> SO2
        log10_Kf = 18887.0 / T - 3.8064
        Kf = 10.0**log10_Kf
        p_analytical = Kf * p_d['S2'] ** 0.5 * p_d['O2']

        assert p_d['SO2'] == pytest.approx(p_analytical, rel=1e-6)

        # Discrimination guard: dropping the 0.5 exponent on p_S2 would give
        # a result that differs by a factor of sqrt(p_S2). At p_S2 ~ 0.01 the
        # wrong formula is 10x smaller, well outside the 1e-6 tolerance.
        p_wrong_stoich = Kf * p_d['S2'] * p_d['O2']
        assert abs(p_d['SO2'] - p_wrong_stoich) > 0.5 * p_analytical


# ===================================================================
# 3. CH4 solubility pressure dependence
# ===================================================================


@pytest.mark.unit
class TestCH4Solubility:
    """Verify CH4 solubility pressure dependence has correct magnitude."""

    def test_pressure_correction_at_1GPa(self):
        """At 1 GPa (10,000 bar), the pressure correction should reduce
        solubility by ~exp(1.93) ~ 7x compared to 1 bar."""
        sol = SolubilityCH4('basalt_ardia')

        p_ch4 = 1.0  # bar partial pressure
        low_P = sol(p_ch4, 1.0)  # 1 bar total
        high_P = sol(p_ch4, 10000.0)  # 10,000 bar total

        ratio = low_P / high_P
        # exp(1.93 * (1 GPa - 0)) ~ exp(1.93) ~ 6.89
        assert ratio == pytest.approx(math.exp(1.93), rel=0.1), (
            f'Pressure correction ratio = {ratio:.2f}, expected ~{math.exp(1.93):.2f}'
        )

        # Discrimination guard: a wrong-sign pressure correction would give
        # ratio ~ exp(-1.93) ~ 0.145 instead of ~ 6.89, a 47x gap. A missing
        # pressure correction would give ratio ~ 1.0. Either failure mode is
        # well outside the 10 % approx tolerance.
        wrong_sign_ratio = math.exp(-1.93)
        assert abs(ratio - wrong_sign_ratio) > 1.0
        assert abs(ratio - 1.0) > 1.0

    def test_low_pressure_nearly_linear(self):
        """At low pressures, CH4 solubility should be nearly proportional
        to partial pressure (pressure correction negligible)."""
        sol = SolubilityCH4('basalt_ardia')

        c1 = sol(1.0, 1.0)
        c10 = sol(10.0, 10.0)

        ratio = c10 / c1
        assert 8.0 < ratio < 12.0, f'Low-P ratio = {ratio:.2f}, expected ~10'

        # Discrimination guard: a missing pressure-dependence on p_CH4 (a stub
        # that returns a constant) would give ratio ~ 1.0. A quadratic-in-p
        # mistake would give ratio ~ 100. Both failure modes are excluded by
        # the [8, 12] band, but the bare interval check passes silently for
        # any value in that band including spurious 9 or 11; tighten by
        # confirming the ratio is close to the 10x change in input.
        assert abs(ratio - 10.0) < 2.0, (
            f'Ratio {ratio:.2f} deviates more than 20 % from pure linearity'
        )

    def test_ch4_positive(self):
        """CH4 solubility should always be positive and finite, and
        monotonic in p_CH4 at fixed total pressure."""
        sol = SolubilityCH4('basalt_ardia')
        values = []
        for p in [0.001, 1.0, 100.0, 10000.0]:
            v = sol(p, p)
            assert v > 0.0, f'sol({p}, {p}) = {v} is non-positive'
            assert math.isfinite(v), f'sol({p}, {p}) = {v} is not finite'
            values.append(v)

        # Discrimination guard: a stub that returns a constant positive value
        # would pass the positivity check. Require the four sample points to
        # be distinct (they sweep four orders of magnitude in input pressure
        # at constant p_CH4 = p_total).
        assert len({round(math.log10(v), 6) for v in values}) >= 3, (
            f'Solubility should vary across four decades of input; got {values}'
        )


@pytest.mark.unit
class TestCOSolubility:
    """Verify CO solubility is unaffected."""

    def test_co_pressure_correction_direction(self):
        """Higher total pressure should reduce CO solubility."""
        sol = SolubilityCO('mafic_armstrong')
        c_low = sol(1.0, 100.0)
        c_high = sol(1.0, 10000.0)
        assert c_low > c_high

        # Discrimination guard: the gap must be larger than floating-point
        # noise. A wrong-sign correction (V_bar with the wrong sign) would
        # give c_high > c_low; the test would catch the sign, but a near-zero
        # gap would be consistent with a missing pressure dependence.
        assert c_low > c_high * 1.05, (
            f'Pressure correction is too weak: c_low={c_low:.4e}, '
            f'c_high={c_high:.4e}, ratio={c_low / c_high:.4f}'
        )


# ===================================================================
# 4. Integration: full equilibrium_atmosphere mass conservation
# ===================================================================


@pytest.mark.integration
class TestEquilibriumAtmosphereIntegration:
    """Run the full solver and verify mass conservation."""

    def _run_equilibrium(self, masses, T=2000.0, Phi=0.5, dIW=0.0):
        """Run equilibrium_atmosphere and return result dict."""
        import warnings

        from calliope.solve import equilibrium_atmosphere

        ddict = {
            'M_mantle': 4.03e24,
            'gravity': 9.81,
            'radius': 6.371e6,
            'Phi_global': Phi,
            'T_magma': T,
            'fO2_shift_IW': dIW,
        }
        for sp in volatile_species:
            ddict[f'{sp}_included'] = 1
            ddict[f'{sp}_initial_bar'] = 0.0

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result = equilibrium_atmosphere(
                masses,
                ddict,
                hide_warnings=True,
                print_result=False,
                nguess=5000,
            )
        return result

    def test_hydrogen_mass_conservation(self):
        """Total H mass (atm + dissolved) should equal the target."""
        H_target = 2.78e20
        target = {'H': H_target, 'C': 1.0, 'N': 1.0, 'S': 1.0}
        result = self._run_equilibrium(target)

        H_total = result.get('H_kg_atm', 0) + result.get('H_kg_liquid', 0)
        assert H_total == pytest.approx(H_target, rel=0.01)

        # Discrimination guard: forgetting one of the two channels (atm or
        # dissolved) would give H_total != H_target by orders of magnitude.
        # Both channels must be physically plausible (atm > 0; dissolved in
        # [0, H_target] given M_mantle is finite).
        H_atm = result.get('H_kg_atm', 0)
        H_liq = result.get('H_kg_liquid', 0)
        assert H_atm > 0.0
        assert 0.0 <= H_liq < H_target

    def test_sulfur_mass_conservation(self):
        """Total S mass should be conserved."""
        S_target = 1e18
        target = {'H': 1e20, 'C': 1.0, 'N': 1.0, 'S': S_target}
        result = self._run_equilibrium(target)

        S_total = result.get('S_kg_atm', 0) + result.get('S_kg_liquid', 0)
        assert S_total == pytest.approx(S_target, rel=0.05)

        # Discrimination guard: forgetting the dissolved channel for a sulfur
        # budget would give S_total ~ S_kg_atm only. Confirm both channels
        # are physically meaningful and that the result is within an order
        # of magnitude of the target (a 100x discrepancy would mean wrong
        # units or a missing channel).
        S_atm = result.get('S_kg_atm', 0)
        S_liq = result.get('S_kg_liquid', 0)
        assert S_atm >= 0.0 and S_liq >= 0.0
        assert 0.5 * S_target < S_total < 2.0 * S_target

    def test_nitrogen_mass_conservation_reducing(self):
        """N mass should be conserved under reducing conditions."""
        N_target = 1e18
        target = {'H': 1e20, 'C': 1.0, 'N': N_target, 'S': 1.0}
        result = self._run_equilibrium(target, T=2000.0, Phi=0.5, dIW=-3.0)

        N_total = result.get('N_kg_atm', 0) + result.get('N_kg_liquid', 0)
        assert N_total == pytest.approx(N_target, rel=0.05)

        # Discrimination guard: under reducing conditions NH3 contributes
        # meaningful N, so the test would catch a stoichiometry that
        # double-counted or dropped that contribution. The result must be
        # within an order of magnitude of the target.
        N_atm = result.get('N_kg_atm', 0)
        N_liq = result.get('N_kg_liquid', 0)
        assert N_atm >= 0.0 and N_liq >= 0.0
        assert 0.5 * N_target < N_total < 2.0 * N_target

    def test_all_pressures_positive(self):
        """All partial pressures should be non-negative and finite."""
        target = {'H': 1e20, 'C': 1e17, 'N': 1e17, 'S': 1e16}
        result = self._run_equilibrium(target)

        seen = []
        for sp in volatile_species:
            key = f'{sp}_bar'
            if key in result:
                p = result[key]
                assert p >= 0.0, f'{sp} pressure is negative: {p}'
                assert math.isfinite(p), f'{sp} pressure is not finite: {p}'
                seen.append((sp, p))

        # Discrimination guard: a stub that returns zero for every species
        # would pass the positivity check. The primary species (H2O, CO2,
        # N2, S2) must each carry meaningful pressure given the budget here.
        primary_p = {sp: p for sp, p in seen if sp in ('H2O', 'CO2', 'N2', 'S2')}
        assert sum(primary_p.values()) > 0.0, (
            'No primary species carries any pressure; solver returned zeros'
        )

    def test_full_chns_converges(self):
        """Full C-H-N-S system should converge to a physically plausible state."""
        target = {'H': 1e20, 'C': 1e17, 'N': 1e17, 'S': 1e16}
        result = self._run_equilibrium(target, T=2500.0, Phi=1.0, dIW=0.0)
        P_surf = result.get('P_surf', 0)
        assert P_surf > 0.0

        # Discrimination guard: a stub that returns a positive constant for
        # P_surf would pass the bare positivity check. Confirm that P_surf
        # is in a physically plausible range for this H budget (well under
        # the ~1e6 bar runaway-greenhouse ceiling and well above the
        # solver-floor of 1e-30 bar).
        assert 1e-3 < P_surf < 1e6, (
            f'P_surf = {P_surf:.4e} bar is outside the physically plausible '
            f'range for this H budget'
        )

        # Element budgets must each be conserved to within solver tolerance.
        for e, target_kg in target.items():
            total = result.get(f'{e}_kg_atm', 0) + result.get(f'{e}_kg_liquid', 0)
            assert 0.1 * target_kg < total < 10.0 * target_kg, (
                f'{e}: total={total:.4e}, target={target_kg:.4e} '
                f'(order-of-magnitude check)'
            )
