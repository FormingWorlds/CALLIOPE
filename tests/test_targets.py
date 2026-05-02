"""Tests for the elemental-budget and pressure-budget target builders.

`get_target_from_params` and `get_target_from_pressures` are the two
public entry points PROTEUS uses to build the H/C/N/S target inventory
that `equilibrium_atmosphere` then solves against. Neither is exercised
by the rest of the suite. This module covers both, including the four
short-circuit branches in `get_target_from_pressures` (low p_total
raises, no-S/no-C/no-N zero-pin shortcuts).
"""

from __future__ import annotations

import pytest

from calliope.constants import molar_mass, ocean_moles, volatile_species
from calliope.solve import (
    equilibrium_atmosphere,
    get_target_from_params,
    get_target_from_pressures,
)

pytestmark = pytest.mark.unit


def _base_ddict(
    M_mantle: float = 4.03e24,
    gravity: float = 9.81,
    radius: float = 6.371e6,
    Phi_global: float = 1.0,
    T_magma: float = 2500.0,
    fO2_shift_IW: float = 0.0,
) -> dict:
    """Minimal ddict for both target builders, with all species included
    and all initial pressures zero.
    """
    d = {
        'M_mantle': M_mantle,
        'gravity': gravity,
        'radius': radius,
        'Phi_global': Phi_global,
        'T_magma': T_magma,
        'fO2_shift_IW': fO2_shift_IW,
    }
    for sp in volatile_species:
        d[f'{sp}_included'] = 1
        d[f'{sp}_initial_bar'] = 0.0
    return d


# ---------------------------------------------------------------------------
# get_target_from_params
# ---------------------------------------------------------------------------


class TestGetTargetFromParams:
    """get_target_from_params reads four "user dial" entries and turns
    them into kg-of-element targets. Each test pins one dial against
    its analytical formula; the four formulas use distinct multipliers,
    so any cross-wiring (e.g. C using the N formula) would fail the
    discriminating ratio assertions.
    """

    def test_earth_like_inventory_matches_analytical(self):
        """1 ocean H, CH=0.1, N=2 ppmw, S=200 ppmw on Earth mantle.
        Discriminating: pick H/C/N/S where the implied kg quantities
        differ by 5-7 orders of magnitude, so swapping any two
        formulas would be obvious from absolute values, not just
        relative ratios.
        """
        ddict = _base_ddict()
        ddict.update(
            hydrogen_earth_oceans=1.0,
            CH_ratio=0.1,
            nitrogen_ppmw=2.0,
            sulfur_ppmw=200.0,
        )
        target = get_target_from_params(ddict)

        H_expected = 1.0 * ocean_moles * molar_mass['H2']
        C_expected = 0.1 * H_expected
        N_expected = 2.0e-6 * ddict['M_mantle']
        S_expected = 200.0e-6 * ddict['M_mantle']

        assert target['H'] == pytest.approx(H_expected, rel=1e-12)
        assert target['C'] == pytest.approx(C_expected, rel=1e-12)
        assert target['N'] == pytest.approx(N_expected, rel=1e-12)
        assert target['S'] == pytest.approx(S_expected, rel=1e-12)

        # Discriminating: orders of magnitude must be distinct so any
        # accidental row swap is visible. Earth-like values give
        # H ~ 1.5e20, C ~ 1.5e19, N ~ 8e18, S ~ 8e20 kg.
        assert target['H'] > target['C']  # H from 1 ocean > C at CH=0.1
        assert target['S'] > target['H']  # 200 ppmw S on a 4e24 mantle
        assert target['C'] > target['N']  # 1.5e19 vs 8e18

    def test_zero_oceans_zero_h_propagates_to_c(self):
        """Edge: zero hydrogen oceans => zero H => zero C (because C is
        defined as CH * H, even with CH != 0)."""
        ddict = _base_ddict()
        ddict.update(
            hydrogen_earth_oceans=0.0,
            CH_ratio=0.5,
            nitrogen_ppmw=2.0,
            sulfur_ppmw=200.0,
        )
        target = get_target_from_params(ddict)

        assert target['H'] == pytest.approx(0.0, abs=1e-30)
        # C = CH * H = 0.5 * 0 = 0, even though CH_ratio is nonzero.
        # A pure read of `C = CH_ratio * 1` would give 0.5; the test
        # discriminates between "C from CH only" and "C from CH * H".
        assert target['C'] == pytest.approx(0.0, abs=1e-30)
        # N and S are independent of the H pathway and stay nonzero
        assert target['N'] > 0.0
        assert target['S'] > 0.0

    def test_negative_sulfur_propagates_unchecked(self):
        """Physically-unreasonable input contract: the function does
        not validate sign, so negative sulfur_ppmw flows straight
        through to a negative S target. Pinning this so a future
        contributor who adds validation breaks the test loudly.
        """
        ddict = _base_ddict()
        ddict.update(
            hydrogen_earth_oceans=1.0,
            CH_ratio=0.1,
            nitrogen_ppmw=2.0,
            sulfur_ppmw=-1.0,
        )
        target = get_target_from_params(ddict)

        assert target['S'] < 0.0
        assert target['S'] == pytest.approx(-1.0e-6 * ddict['M_mantle'], rel=1e-12)

    def test_doubling_mantle_mass_doubles_n_and_s_only(self):
        """Discriminating: M_mantle scales N and S (per-mass-fraction)
        but not H (set in oceans) or C (CH * H, so independent of
        M_mantle). A bug that multiplied H or C by M_mantle would
        change those entries when we double the mantle.
        """
        d1 = _base_ddict(M_mantle=4.03e24)
        d1.update(hydrogen_earth_oceans=1.0, CH_ratio=0.1, nitrogen_ppmw=2.0, sulfur_ppmw=200.0)
        t1 = get_target_from_params(d1)

        d2 = _base_ddict(M_mantle=8.06e24)
        d2.update(hydrogen_earth_oceans=1.0, CH_ratio=0.1, nitrogen_ppmw=2.0, sulfur_ppmw=200.0)
        t2 = get_target_from_params(d2)

        assert t2['H'] == pytest.approx(t1['H'], rel=1e-12)
        assert t2['C'] == pytest.approx(t1['C'], rel=1e-12)
        assert t2['N'] == pytest.approx(2.0 * t1['N'], rel=1e-12)
        assert t2['S'] == pytest.approx(2.0 * t1['S'], rel=1e-12)


# ---------------------------------------------------------------------------
# get_target_from_pressures
# ---------------------------------------------------------------------------


class TestGetTargetFromPressures:
    """get_target_from_pressures inverts equilibrium_atmosphere: given
    initial primary partial pressures, it computes the implied total
    elemental masses by summing atmospheric column mass + dissolved
    Henry's-law mass.

    Three short-circuits exist (S, C, N pin to zero when no
    corresponding species is present in the initial atmosphere).
    Test each.
    """

    def test_low_total_pressure_raises(self):
        """Sum of initial pressures below 1e-3 bar must raise.
        Hits the `if p_tot < 1.0e-3: raise` short-circuit (line ~347).
        """
        ddict = _base_ddict()
        ddict['H2O_initial_bar'] = 1e-5
        ddict['CO2_initial_bar'] = 1e-5
        ddict['N2_initial_bar'] = 1e-5
        ddict['S2_initial_bar'] = 1e-5

        with pytest.raises(Exception, match='Initial surface pressure too low'):
            get_target_from_pressures(ddict)

    def test_no_sulfur_zeros_S(self):
        """S2/SO2/H2S all at zero initial pressure => target['S'] = 0
        without consulting solubility/atmosphere mass.
        Hits line ~356.
        """
        ddict = _base_ddict()
        ddict['H2O_initial_bar'] = 220.0
        ddict['CO2_initial_bar'] = 100.0
        ddict['N2_initial_bar'] = 1.0
        ddict['S2_initial_bar'] = 0.0
        ddict['SO2_initial_bar'] = 0.0
        ddict['H2S_initial_bar'] = 0.0
        target = get_target_from_pressures(ddict)

        assert target['S'] == pytest.approx(0.0, abs=1e-30)
        # Discriminating: H, C, N populated normally
        assert target['H'] > 0.0
        assert target['C'] > 0.0
        assert target['N'] > 0.0

    def test_no_carbon_zeros_C(self):
        """CO2/CO at zero => target['C'] = 0. Hits line ~363."""
        ddict = _base_ddict()
        ddict['H2O_initial_bar'] = 220.0
        ddict['CO2_initial_bar'] = 0.0
        ddict['CO_initial_bar'] = 0.0
        ddict['N2_initial_bar'] = 1.0
        ddict['S2_initial_bar'] = 0.01
        # CO2 = 0 => CO2 dissolution path evaluates log10(0); the helper
        # masks the resulting -inf/0 contribution but numpy still emits
        # a divide-by-zero RuntimeWarning. Expected.
        with pytest.warns(RuntimeWarning, match='divide by zero'):
            target = get_target_from_pressures(ddict)

        assert target['C'] == pytest.approx(0.0, abs=1e-30)
        assert target['H'] > 0.0
        assert target['S'] > 0.0

    def test_no_nitrogen_zeros_N(self):
        """N2/NH3 at zero => target['N'] = 0. Hits line ~370."""
        ddict = _base_ddict()
        ddict['H2O_initial_bar'] = 220.0
        ddict['CO2_initial_bar'] = 100.0
        ddict['N2_initial_bar'] = 0.0
        ddict['NH3_initial_bar'] = 0.0
        ddict['S2_initial_bar'] = 0.01
        target = get_target_from_pressures(ddict)

        assert target['N'] == pytest.approx(0.0, abs=1e-30)
        assert target['H'] > 0.0
        assert target['C'] > 0.0

    def test_round_trip_recovers_pressures(self):
        """Property: given initial pressures P0, build target = T(P0),
        then solve equilibrium_atmosphere with the target back. The
        recovered pressures must approximate P0 within solver tolerance,
        because P0 already satisfies the mass-balance equations
        (it was used to build the target in the first place).

        Warm-start with p_guess=P0 so the inner fsolve converges in
        a single iteration and the test runs in milliseconds.
        """
        ddict = _base_ddict()
        P0 = {'H2O': 220.0, 'CO2': 100.0, 'N2': 1.0, 'S2': 0.01}
        for sp, p in P0.items():
            ddict[f'{sp}_initial_bar'] = p

        target = get_target_from_pressures(ddict)

        # Sanity: the four element targets are nontrivially nonzero
        # (each pressure couples to its element via column mass +
        # solubility, both > 0).
        for e in ('H', 'C', 'N', 'S'):
            assert target[e] > 0.0, f'target[{e}] is zero, round-trip ill-posed'

        result = equilibrium_atmosphere(
            target,
            ddict,
            p_guess=P0,
            print_result=False,
            nguess=10,  # warm start should converge fast
        )

        # Recover within 1% — the target was built from these very
        # pressures, so the solver lands here exactly modulo float
        # round-off and Powell-hybrid tolerance.
        assert result['H2O_bar'] == pytest.approx(P0['H2O'], rel=0.01)
        assert result['CO2_bar'] == pytest.approx(P0['CO2'], rel=0.01)
        assert result['N2_bar'] == pytest.approx(P0['N2'], rel=0.01)
        assert result['S2_bar'] == pytest.approx(P0['S2'], rel=0.01)

    def test_unphysical_negative_initial_pressure_raises_or_propagates(self):
        """Negative initial pressure in one slot. Function does not
        validate; it adds it to p_tot and continues. If the negative is
        large enough to push p_tot below 1e-3, the low-p check raises.
        """
        ddict = _base_ddict()
        # Net total -100 bar => well below 1e-3, must raise
        ddict['H2O_initial_bar'] = -100.0
        ddict['CO2_initial_bar'] = 0.001
        ddict['N2_initial_bar'] = 0.001
        ddict['S2_initial_bar'] = 0.001

        with pytest.raises(Exception, match='Initial surface pressure too low'):
            get_target_from_pressures(ddict)
