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

    def test_S_from_S2(self):
        """S2 has 2 S atoms. Atmospheric S should be ~2*M_S/M_S2 * mass_S2."""
        ddict = _make_ddict()
        pin = {'H2O': 1e-30, 'CO2': 1e-30, 'N2': 1e-30, 'S2': 10.0}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        # In a S2-dominated atmosphere, mu ~ M_S2
        mass_S2 = _column_mass(p_d['S2'])
        expected_S = mass_S2 * 2 * molar_mass['S'] / molar_mass['S2']
        assert mass['S'] == pytest.approx(expected_S, rel=0.01)

    def test_N_from_N2(self):
        """N2 has 2 N atoms."""
        ddict = _make_ddict()
        pin = {'H2O': 1e-30, 'CO2': 1e-30, 'N2': 10.0, 'S2': 1e-30}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        mass_N2 = _column_mass(p_d['N2'])
        expected_N = mass_N2 * 2 * molar_mass['N'] / molar_mass['N2']
        # NH3 is derived from N2, contributing a small fraction
        assert mass['N'] == pytest.approx(expected_N, rel=0.05)

    def test_C_from_CO2(self):
        """CO2 has 1 C atom."""
        ddict = _make_ddict()
        pin = {'H2O': 1e-30, 'CO2': 10.0, 'N2': 1e-30, 'S2': 1e-30}
        p_d, mass = self._get_elemental_masses(pin, ddict)

        mass_CO2 = _column_mass(p_d['CO2'])
        expected_C = mass_CO2 * 1 * molar_mass['C'] / molar_mass['CO2']
        # CO and CH4 are derived from CO2
        assert mass['C'] >= expected_C * 0.90

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

    def test_elemental_masses_all_positive(self):
        """All elemental masses should be non-negative."""
        ddict = _make_ddict()
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        _, mass = self._get_elemental_masses(pin, ddict)

        for e in element_list:
            assert mass[e] >= 0.0, f'{e} mass is negative: {mass[e]}'


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

    def test_H2_and_CO_unchanged(self):
        """H2 and CO reactions should be unaffected by the chemistry changes."""
        T = 2000.0
        mk_h2 = ModifiedKeq('janaf_H2')
        mk_co = ModifiedKeq('janaf_CO')

        g_h2 = mk_h2(T, 0.0)
        g_co = mk_co(T, 0.0)

        # Precomputed values from before the fix (must not change)
        assert g_h2 == pytest.approx(1.469, rel=1e-2)
        assert g_co == pytest.approx(6.581, rel=1e-2)

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

    def test_low_pressure_nearly_linear(self):
        """At low pressures, CH4 solubility should be nearly proportional
        to partial pressure (pressure correction negligible)."""
        sol = SolubilityCH4('basalt_ardia')

        c1 = sol(1.0, 1.0)
        c10 = sol(10.0, 10.0)

        ratio = c10 / c1
        assert 8.0 < ratio < 12.0, f'Low-P ratio = {ratio:.2f}, expected ~10'

    def test_ch4_positive(self):
        """CH4 solubility should always be positive."""
        sol = SolubilityCH4('basalt_ardia')
        for p in [0.001, 1.0, 100.0, 10000.0]:
            assert sol(p, p) > 0.0


@pytest.mark.unit
class TestCOSolubility:
    """Verify CO solubility is unaffected."""

    def test_co_pressure_correction_direction(self):
        """Higher total pressure should reduce CO solubility."""
        sol = SolubilityCO('mafic_armstrong')
        c_low = sol(1.0, 100.0)
        c_high = sol(1.0, 10000.0)
        assert c_low > c_high


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

    def test_sulfur_mass_conservation(self):
        """Total S mass should be conserved."""
        S_target = 1e18
        target = {'H': 1e20, 'C': 1.0, 'N': 1.0, 'S': S_target}
        result = self._run_equilibrium(target)

        S_total = result.get('S_kg_atm', 0) + result.get('S_kg_liquid', 0)
        assert S_total == pytest.approx(S_target, rel=0.05)

    def test_nitrogen_mass_conservation_reducing(self):
        """N mass should be conserved under reducing conditions."""
        N_target = 1e18
        target = {'H': 1e20, 'C': 1.0, 'N': N_target, 'S': 1.0}
        result = self._run_equilibrium(target, T=2000.0, Phi=0.5, dIW=-3.0)

        N_total = result.get('N_kg_atm', 0) + result.get('N_kg_liquid', 0)
        assert N_total == pytest.approx(N_target, rel=0.05)

    def test_all_pressures_positive(self):
        """All partial pressures should be non-negative."""
        target = {'H': 1e20, 'C': 1e17, 'N': 1e17, 'S': 1e16}
        result = self._run_equilibrium(target)

        for sp in volatile_species:
            key = f'{sp}_bar'
            if key in result:
                assert result[key] >= 0.0, f'{sp} pressure is negative'

    def test_full_chns_converges(self):
        """Full C-H-N-S system should converge."""
        target = {'H': 1e20, 'C': 1e17, 'N': 1e17, 'S': 1e16}
        result = self._run_equilibrium(target, T=2500.0, Phi=1.0, dIW=0.0)
        assert result.get('P_surf', 0) > 0.0
