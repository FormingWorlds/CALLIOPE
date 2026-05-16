"""Tests for partial-species inclusion paths.

Covers the `is_included(...) == False` branches of `get_partial_pressures`,
`atmosphere_mass`, and `dissolved_mass`, plus the `SolubilityN2.libourel`
solubility law that the rest of the suite never instantiates.

The existing `test_stoichiometry.py::_make_ddict` helper hardcodes every
species to `included=1`. This file uses its own builder that accepts an
`included` override dict so the False branches actually fire.
"""

from __future__ import annotations

import math

import pytest

from calliope.constants import molar_mass, volatile_species
from calliope.solubility import SolubilityN2
from calliope.solve import (
    atmosphere_mass,
    dissolved_mass,
    get_partial_pressures,
)

pytestmark = pytest.mark.unit


def _make_ddict(
    included: dict | None = None,
    T: float = 2000.0,
    fO2_shift: float = 0.0,
    gravity: float = 9.81,
    radius: float = 6.371e6,
    M_mantle: float = 4.03e24,
    Phi_global: float = 0.5,
) -> dict:
    """Build a ddict with selective species exclusion.

    Parameters
    ----------
    included
        Mapping of species name to 0/1. Anything not listed defaults to 1.
    """
    overrides = included or {}
    d = {
        'T_magma': T,
        'fO2_shift_IW': fO2_shift,
        'gravity': gravity,
        'radius': radius,
        'M_mantle': M_mantle,
        'Phi_global': Phi_global,
    }
    for sp in volatile_species:
        d[f'{sp}_included'] = overrides.get(sp, 1)
    return d


# ---------------------------------------------------------------------------
# get_partial_pressures: False branches of the is_included gating
# ---------------------------------------------------------------------------


class TestGetPartialPressuresExclusions:
    """Verify that each species-inclusion flag flips the right secondary
    species off, and only that one (or its dependents)."""

    def test_h2_excluded_zeros_h2_and_h_dependents(self):
        """With H2 excluded, every H-bearing reduced species must be zero,
        regardless of its own inclusion flag, because they all multiply
        through `p_d['H2']` in the equilibrium expressions.
        """
        ddict = _make_ddict(included={'H2': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        p_d = get_partial_pressures(pin, ddict)

        # H2 directly suppressed
        assert p_d['H2'] == pytest.approx(0.0, abs=1e-30)

        # CH4 gated on `is_included('H2') and is_included('CH4')` so it
        # vanishes even with CH4_included=1 set above.
        assert p_d['CH4'] == pytest.approx(0.0, abs=1e-30)

        # H2S has the same combined gate; NH3 has no flag check on H2 but
        # multiplies through p_d['H2']**3 so still vanishes.
        assert p_d['H2S'] == pytest.approx(0.0, abs=1e-30)
        assert p_d['NH3'] == pytest.approx(0.0, abs=1e-30)

        # Discriminating: O2 unaffected (independent of H), CO unaffected
        # (function of CO2 only). If our flag flip wrongly disabled CO,
        # this assert would fail.
        assert p_d['O2'] > 0.0
        assert p_d['CO'] > 0.0

        # Primaries preserved
        assert p_d['H2O'] == pytest.approx(100.0, rel=1e-12)
        assert p_d['CO2'] == pytest.approx(10.0, rel=1e-12)

    def test_so2_h2s_excluded_zeros_only_those(self):
        """SO2 and H2S off; S2 primary still passes through."""
        ddict = _make_ddict(included={'SO2': 0, 'H2S': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        p_d = get_partial_pressures(pin, ddict)

        assert p_d['SO2'] == pytest.approx(0.0, abs=1e-30)
        assert p_d['H2S'] == pytest.approx(0.0, abs=1e-30)

        # Discriminating: S2 primary not zeroed, downstream sulfur
        # secondaries have been the only thing knocked out.
        assert p_d['S2'] == pytest.approx(0.1, rel=1e-12)

    def test_nh3_excluded_zeros_nh3_only(self):
        """NH3 off; N2 primary unchanged."""
        ddict = _make_ddict(included={'NH3': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        p_d = get_partial_pressures(pin, ddict)

        assert p_d['NH3'] == pytest.approx(0.0, abs=1e-30)
        assert p_d['N2'] == pytest.approx(1.0, rel=1e-12)

        # Discriminating: H2S unaffected even though both depend on H2
        assert p_d['H2S'] > 0.0

    def test_unphysical_negative_partial_pressure_raises_or_clips(self):
        """The function clips negative outputs to 0 via the explicit
        non-negative-real clip at the bottom of `_get_partial_pressures`.
        Feed it a primary pressure that drives a derived species negative
        (impossible at any real fO2 but the clip path must still hold).
        """
        import warnings as _warnings

        ddict = _make_ddict()
        pin = {'H2O': -5.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        # Whether intermediate steps emit a RuntimeWarning, drop into a
        # NaN path, or produce a complex sqrt-of-negative depends on the
        # buffer; the clip must absorb all three. Suppress any warnings
        # and pin the contract that matters: every output a non-negative
        # real.
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore', RuntimeWarning)
            p_d = get_partial_pressures(pin, ddict)

        for sp, p in p_d.items():
            assert isinstance(p, float), f'{sp} = {p!r} is not a real float'
            assert p >= 0.0, f'{sp} = {p} broke the non-negative clip'


# ---------------------------------------------------------------------------
# atmosphere_mass: False branches in the elemental-tally if-blocks
# ---------------------------------------------------------------------------


class TestAtmosphereMassExclusions:
    """When secondary species are excluded, the elemental tallies must
    not include their contributions.
    """

    def test_h_tally_excludes_h2_ch4_h2s_nh3_when_off(self):
        """H is fully attributable to H2O when every H-bearing reduced
        species is excluded. The expected H mass is exactly
        2 * M_H * mass_H2O / M_H2O.
        """
        ddict_off = _make_ddict(included={'H2': 0, 'CH4': 0, 'H2S': 0, 'NH3': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        mass_off = atmosphere_mass(pin, ddict_off)

        # mass_atm_d['H2O'] is computed inside atmosphere_mass; recompute
        # the analytical column mass (modulo the mu-correction, which we
        # invert by reading mass_atm_d['H2O'] back).
        expected_H = 2 * mass_off['H2O'] / molar_mass['H2O'] * molar_mass['H']
        assert mass_off['H'] == pytest.approx(expected_H, rel=1e-9)

        # Discrimination guard: with all reduced H-bearing species ON, the
        # H tally is strictly higher because H2, CH4, H2S, and NH3 each
        # carry additional H atoms. The exact-match-to-H2O equality above
        # would not hold.
        ddict_on = _make_ddict(included={'H2': 1, 'CH4': 1, 'H2S': 1, 'NH3': 1})
        mass_on = atmosphere_mass(pin, ddict_on)
        assert mass_on['H'] > mass_off['H'], (
            f'H tally with reduced species ON ({mass_on["H"]:.4e}) is not '
            f'greater than with them OFF ({mass_off["H"]:.4e})'
        )

    def test_c_tally_excludes_co_ch4_when_off(self):
        """C tally collapses to mass_CO2 / M_CO2 when CO and CH4 off."""
        ddict_off = _make_ddict(included={'CO': 0, 'CH4': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        mass_off = atmosphere_mass(pin, ddict_off)

        expected_C = mass_off['CO2'] / molar_mass['CO2'] * molar_mass['C']
        assert mass_off['C'] == pytest.approx(expected_C, rel=1e-9)

        # Discrimination guard: with CO and CH4 ON, the C tally is strictly
        # higher (CO and CH4 each contribute additional C atoms).
        ddict_on = _make_ddict(included={'CO': 1, 'CH4': 1})
        mass_on = atmosphere_mass(pin, ddict_on)
        assert mass_on['C'] > mass_off['C']

    def test_s_tally_excludes_so2_h2s_when_off(self):
        """S tally collapses to 2 * mass_S2 / M_S2 when SO2 and H2S off."""
        ddict_off = _make_ddict(included={'SO2': 0, 'H2S': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        mass_off = atmosphere_mass(pin, ddict_off)

        expected_S = 2 * mass_off['S2'] / molar_mass['S2'] * molar_mass['S']
        assert mass_off['S'] == pytest.approx(expected_S, rel=1e-9)

        # Discrimination guard: with SO2 and H2S ON, the S tally is strictly
        # higher (each contributes one additional S atom per molecule).
        ddict_on = _make_ddict(included={'SO2': 1, 'H2S': 1})
        mass_on = atmosphere_mass(pin, ddict_on)
        assert mass_on['S'] > mass_off['S']

    def test_zero_pressure_inputs_yield_nonneg_masses(self):
        """Edge case: every primary at numerical zero. Output masses must
        be non-negative and finite, never NaN. The `max(0.0, ...)` clip
        at the end of atmosphere_mass is the contract.
        """
        ddict = _make_ddict()
        pin = {'H2O': 0.0, 'CO2': 0.0, 'N2': 0.0, 'S2': 0.0}
        mass = atmosphere_mass(pin, ddict)

        for k, v in mass.items():
            assert math.isfinite(v), f'{k} non-finite at zero input'
            assert v >= 0.0, f'{k} = {v} broke clip at zero input'


# ---------------------------------------------------------------------------
# dissolved_mass: CO and CH4 explicit zero branches
# ---------------------------------------------------------------------------


class TestDissolvedMassExclusions:
    """When CO or CH4 is excluded, dissolved_mass must store an explicit
    0.0 (lines 233 and 241), not a missing key. The downstream element
    tallies must skip the species cleanly.
    """

    def test_co_excluded_writes_zero(self):
        ddict = _make_ddict(included={'CO': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        m = dissolved_mass(pin, ddict)

        # Explicit zero, not absent
        assert 'CO' in m
        assert m['CO'] == pytest.approx(0.0, abs=1e-30)

        # Discriminating: CO2 still dissolves (gating is on CO only)
        assert m['CO2'] > 0.0

    def test_ch4_excluded_writes_zero(self):
        ddict = _make_ddict(included={'CH4': 0})
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        m = dissolved_mass(pin, ddict)

        assert 'CH4' in m
        assert m['CH4'] == pytest.approx(0.0, abs=1e-30)

        # CO2 and CO untouched
        assert m['CO2'] > 0.0
        assert m['CO'] > 0.0

    def test_zero_phi_zeros_all_dissolved(self):
        """Edge: with zero melt fraction, every dissolved mass is zero.

        `dissolved_mass` only writes the species that can dissolve in
        a silicate melt (H2O, CO2, CO, CH4, N2, S2), never O2, SO2,
        H2S, NH3, or H2. Iterate only over those keys.
        """
        ddict_zero = _make_ddict(Phi_global=0.0)
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        m_zero = dissolved_mass(pin, ddict_zero)
        dissolved_species = {'H2O', 'CO2', 'CO', 'CH4', 'N2', 'S2'}
        for sp in dissolved_species:
            assert m_zero[sp] == pytest.approx(0.0, abs=1e-30), (
                f'{sp} nonzero at Phi=0'
            )

        # Discrimination guard: at Phi=1.0 (fully molten) the same species
        # carry nonzero dissolved mass. A stub that hard-coded 0.0 for
        # every dissolved entry would also pass the Phi=0 loop above.
        ddict_full = _make_ddict(Phi_global=1.0)
        m_full = dissolved_mass(pin, ddict_full)
        nonzero_at_full = sum(
            1 for sp in dissolved_species if m_full[sp] > 0.0
        )
        assert nonzero_at_full >= 4, (
            f'Only {nonzero_at_full} of {len(dissolved_species)} species '
            f'dissolved at Phi=1.0; the zero-Phi check would be vacuous'
        )

    def test_unphysical_negative_mantle_mass_propagates(self):
        """Pin the contract: negative M_mantle is not validated here, it
        flows straight through into a negative dissolved mass. Document
        that surprise so a future caller sees it.
        """
        ddict = _make_ddict(M_mantle=-1.0e24)
        pin = {'H2O': 100.0, 'CO2': 10.0, 'N2': 1.0, 'S2': 0.1}
        m = dissolved_mass(pin, ddict)
        # Not all dissolved fields are clipped — only the per-element
        # tallies at the bottom are. The per-species mass for H2O can go
        # negative when M_mantle is negative.
        assert m['H2O'] < 0.0  # documents non-validation
        # But element-level entries are clipped to 0
        assert m['H'] == pytest.approx(0.0, abs=1e-20)


# ---------------------------------------------------------------------------
# SolubilityN2.libourel
# ---------------------------------------------------------------------------


class TestSolubilityN2Libourel:
    """The libourel mode of SolubilityN2 is a Henry's-law linear law:
    ppmw = 0.0611 * p. Discriminating values: pick pressures where p,
    p^0.5, and p^2 give visibly different results to rule out the
    plausible wrong powers.
    """

    def test_libourel_linear_in_p(self):
        sol = SolubilityN2('libourel')
        # Discriminating: at p=1, all of p^0.5/p/p^2 give 1.0. At p=100
        # they differ by a factor of 10 (sqrt) or 10000 (square), so the
        # ratio test below distinguishes the correct power.
        c1 = sol(1.0)
        c100 = sol(100.0)

        assert c1 == pytest.approx(0.0611, rel=1e-12)
        assert c100 == pytest.approx(0.0611 * 100.0, rel=1e-12)
        # Distinguishes linear from sqrt (10) or quadratic (10000)
        assert c100 / c1 == pytest.approx(100.0, rel=1e-12)

    def test_libourel_zero_p(self):
        """Edge: zero pressure should give zero dissolved concentration."""
        sol = SolubilityN2('libourel')
        assert sol(0.0) == pytest.approx(0.0, abs=1e-30)

        # Discrimination guard: a stub that always returned 0 would pass
        # the bare zero-input check. Confirm the call path is genuinely
        # linear by checking that a small positive p gives a small positive
        # output and a larger p gives a proportionally larger one.
        small = sol(1e-6)
        large = sol(1.0)
        assert small > 0.0
        assert large > small * 100.0  # linear should give exactly 1e6x

    def test_libourel_negative_p_propagates_unchecked(self):
        """The Henry's-law power-law is not guarded against p < 0;
        document that it returns a negative concentration for negative
        input. A caller that passes a sub-zero pressure has a bigger
        problem than this function silently signing the output.
        """
        sol = SolubilityN2('libourel')
        out = sol(-10.0)
        # 0.0611 * (-10)^1 = -0.611
        assert out == pytest.approx(-0.611, rel=1e-12)
        assert out < 0.0  # explicit unphysical-input contract

    def test_libourel_distinguishable_from_dasgupta(self):
        """At a representative MO state, libourel and dasgupta give
        different ppmw (otherwise the test_dissolved_mass coverage of
        the dasgupta branch would silently double as a libourel test).
        """
        sol_lib = SolubilityN2('libourel')
        sol_das = SolubilityN2('dasgupta')

        p_N2 = 1.0  # bar
        p_tot = 100.0  # bar
        T = 2000.0
        fO2_shift = 0.0

        c_lib = sol_lib(p_N2)
        c_das = sol_das(p_N2, p_tot, T, fO2_shift)

        # Both positive and finite
        assert c_lib > 0.0 and math.isfinite(c_lib)
        assert c_das > 0.0 and math.isfinite(c_das)

        # And actually distinct (relative difference > 10%) so the two
        # code paths are not accidentally measuring the same number.
        rel_diff = abs(c_lib - c_das) / max(c_lib, c_das)
        assert rel_diff > 0.1, (
            f'libourel ({c_lib:.3e}) and dasgupta ({c_das:.3e}) '
            'agree too closely; one of the two code paths may be a no-op'
        )
