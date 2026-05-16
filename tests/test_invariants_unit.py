"""Unit-tier invariants for CALLIOPE's solubility laws and atmospheric
mass tallies.

This file holds the unit-tier subset of the invariants previously
co-located in `test_invariants.py`. They test individual physics laws
(`SolubilityS2.gaillard`, `SolubilityN2.dasgupta`, `SolubilityN2.libourel`)
and the atmospheric-mass stoichiometry (`_atmosphere_mass`) directly,
without invoking the multi-species solver. The smoke-tier invariants
(which do invoke the solver) stay in `test_invariants.py`.

The four classes pulled here:

- `TestCO2AtomCounting`: closed-form C tally from the CO2 column.
- `TestS2SolubilityMonotonicity`: Gaillard monotonicity with redox.
- `TestN2SolubilityMonotonicity`: Dasgupta monotonicity and Libourel
  redox-independence.
- `TestDasguptaReducingEdgeLaw`: Dasgupta law finite at strongly-reducing
  conditions.
"""

from __future__ import annotations

import math

import pytest

from calliope.constants import molar_mass, volatile_species
from calliope.solubility import SolubilityN2, SolubilityS2
from calliope.solve import _atmosphere_mass

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _ddict(T: float = 1800.0, Phi: float = 1.0, dIW: float = 4.0) -> dict:
    """Realistic ddict with every volatile species included."""
    d = {
        'M_mantle': 4.03e24,
        'gravity': 9.81,
        'radius': 6.371e6,
        'Phi_global': Phi,
        'T_magma': T,
        'fO2_shift_IW': dIW,
    }
    for sp in volatile_species:
        d[f'{sp}_included'] = 1
        d[f'{sp}_initial_bar'] = 0.0
    return d


# ===========================================================================
# CO2 stoichiometric atom counting
# ===========================================================================


class TestCO2AtomCounting:
    """At known p_CO2, the C contribution from CO2 in the atmosphere
    is exactly (12/44) * CO2_kg_atm via stoichiometric atom counting."""

    @pytest.mark.parametrize('p_CO2_bar', [0.1, 1.0, 10.0, 100.0])
    def test_C_atom_count_from_CO2_only(self, p_CO2_bar):
        """Single-species CO2 atmosphere at p_CO2 in {0.1, 1, 10, 100} bar:
        the C tally equals (12/44) * CO2 column mass to within rel=1e-12,
        with a wrong-factor-2 discrimination guard."""
        # Single-species atmosphere with only CO2 to isolate C contribution
        ddict = _ddict(T=1800.0, Phi=0.0, dIW=0.0)
        for sp in volatile_species:
            ddict[f'{sp}_included'] = 0
        ddict['CO2_included'] = 1
        # Set all pressures to zero, then CO2
        p_d = {s: 0.0 for s in volatile_species}
        p_d['CO2'] = p_CO2_bar
        mass_atm = _atmosphere_mass(p_d, 0.0, ddict)
        # Atomic-C contribution expected from stoichiometry: 12.011 / 44.01
        # times the CO2 column mass
        expected_C_kg = mass_atm['CO2'] * molar_mass['C'] / molar_mass['CO2']
        assert mass_atm['C'] == pytest.approx(expected_C_kg, rel=1e-12)

        # Discrimination guard: the wrong stoichiometry (factor 2 for C in
        # CO2, i.e. treating CO2 as having 2 C atoms) would double the C
        # tally, well outside the rel=1e-12 tolerance.
        wrong_C_kg = 2.0 * mass_atm['CO2'] * molar_mass['C'] / molar_mass['CO2']
        assert abs(mass_atm['C'] - wrong_C_kg) > 0.5 * expected_C_kg

    def test_zero_CO2_pressure_gives_zero_C(self):
        """Sad-path: at p_CO2 = 0 the C tally from CO2 channel is zero."""
        ddict = _ddict(T=1800.0, Phi=0.0, dIW=0.0)
        for sp in volatile_species:
            ddict[f'{sp}_included'] = 0
        ddict['CO2_included'] = 1
        p_d = {s: 0.0 for s in volatile_species}
        mass_atm = _atmosphere_mass(p_d, 0.0, ddict)
        assert mass_atm.get('C', 0.0) == pytest.approx(0.0, abs=1e-30)

        # Discrimination guard: confirm the CO2 column mass is also zero.
        # A stub that returned zero only for the C key (but nonzero CO2
        # column mass) would pass the bare C == 0 check, hiding a real
        # stoichiometry bug elsewhere.
        assert mass_atm.get('CO2', 0.0) == pytest.approx(0.0, abs=1e-30)


# ===========================================================================
# S2 Gaillard monotonicity with redox
# ===========================================================================


class TestS2SolubilityMonotonicity:
    """Gaillard sulfide-saturated solubility increases as fO2 decreases
    at fixed p_S2 and T (the ``+0.5 ln(p_S2/fO2)`` term carries the
    redox dependence directly). Tested at the solubility-function level
    to isolate the law from the multi-species solver feedback loops."""

    @pytest.mark.parametrize('p_S2_bar', [0.01, 0.1, 1.0])
    @pytest.mark.parametrize('T', [1500.0, 1800.0, 2200.0])
    def test_gaillard_strictly_decreasing_with_oxidation(self, T, p_S2_bar):
        """Gaillard S2 solubility strictly decreases at every adjacent dIW
        step across [-4, +4], with a span check requiring the endpoint
        ratio to exceed 10x."""
        S2 = SolubilityS2('gaillard')
        dIWs = [-4.0, -2.0, 0.0, +2.0, +4.0]
        values = [S2.gaillard(p_S2_bar, T, dIW) for dIW in dIWs]
        for i in range(len(values) - 1):
            assert values[i] > values[i + 1], (
                f'Gaillard ppmw at dIW={dIWs[i]} ({values[i]:.4e}) '
                f'not greater than at dIW={dIWs[i+1]} ({values[i+1]:.4e}) '
                f'(T={T}, p_S2={p_S2_bar})'
            )

        # Discrimination guard: the monotonicity span across the 8-dex
        # dIW range should be meaningful (the +0.5 ln(p_S2/fO2) term gives
        # at least a 10x change in ppmw across dIWs in [-4, +4]). A near-
        # constant function would pass the > check trivially at each step.
        assert values[0] / values[-1] > 10.0, (
            f'Gaillard span too small: values[-4]/values[+4] = '
            f'{values[0] / values[-1]:.2f}, expected > 10x'
        )

    def test_negative_pressure_returns_zero(self):
        """Sad-path: at p_S2 < 1e-20 bar the implementation returns 0
        to avoid log(0). Verify the floor behaves correctly."""
        S2 = SolubilityS2('gaillard')
        assert S2.gaillard(1e-30, 1800.0, 0.0) == pytest.approx(0.0, abs=1e-30)

        # Discrimination guard: a stub that always returned 0 would pass.
        # Confirm a normal p_S2 of 0.1 bar gives a nonzero, finite value
        # so the floor branch is genuinely distinct from the main path.
        normal = S2.gaillard(0.1, 1800.0, 0.0)
        assert normal > 0.0
        assert math.isfinite(normal)


# ===========================================================================
# N2 Dasgupta monotonicity and Libourel redox-independence
# ===========================================================================


class TestN2SolubilityMonotonicity:
    """Dasgupta N2 solubility increases as fO2 decreases at fixed
    p_N2, p_tot, T (the ``-1.6 dIW`` term in the reduced-N branch
    drives this). Tested at the solubility-function level."""

    @pytest.mark.parametrize('p_N2_bar', [0.1, 1.0, 10.0])
    @pytest.mark.parametrize('T', [1500.0, 1800.0, 2200.0])
    def test_dasgupta_monotonic_with_oxidation(self, T, p_N2_bar):
        """Dasgupta N2 solubility is monotonically non-increasing as dIW
        rises from -6 to +4, with a span check requiring the endpoint
        ratio to exceed 10x."""
        N2 = SolubilityN2('dasgupta')
        dIWs = [-6.0, -4.0, -2.0, 0.0, +2.0, +4.0]
        values = [N2.dasgupta(p_N2_bar, p_N2_bar, T, dIW) for dIW in dIWs]
        # Dasgupta has two terms (reduced-N and molecular N2); the
        # reduced-N term dominates at low fO2 and falls steeply with
        # increasing fO2, so the total is monotonically decreasing.
        for i in range(len(values) - 1):
            assert values[i] >= values[i + 1], (
                f'Dasgupta ppmw at dIW={dIWs[i]} ({values[i]:.4e}) '
                f'less than at dIW={dIWs[i+1]} ({values[i+1]:.4e}) '
                f'(T={T}, p_N2={p_N2_bar})'
            )

        # Discrimination guard: the reduced-N branch carries an exp(-1.6 dIW)
        # factor, so the span across dIW in [-6, +4] should be at least
        # exp(1.6 * 10) ~ 9e6 in the reduced-N contribution. The molecular
        # N2 term puts a floor under the high-dIW end, but the span should
        # still be > 10x. A near-constant function would pass the >= check.
        assert values[0] / values[-1] > 10.0, (
            f'Dasgupta span too small: values[-6]/values[+4] = '
            f'{values[0] / values[-1]:.2f}, expected > 10x'
        )

    def test_libourel_alternative_is_redox_independent(self):
        """Sad-path / contrast: the Libourel linear Henry's law has no
        fO2 term, so it does NOT show the Dasgupta monotonicity. This
        documents that ``libourel`` and ``dasgupta`` are distinct laws
        with distinct calibration footprints."""
        N2 = SolubilityN2('libourel')
        # Libourel takes only p_N2; no fO2 dependence
        assert N2.libourel(1.0) == N2.libourel(1.0)  # trivially deterministic
        # Two different p_N2 give different values (linear)
        assert N2.libourel(2.0) == pytest.approx(2.0 * N2.libourel(1.0))


# ===========================================================================
# Dasgupta plateau behaviour at strongly-reducing fO2 (law-only)
# ===========================================================================


class TestDasguptaReducingEdgeLaw:
    """Direct check that the Dasgupta law itself stays numerically
    well-behaved at strongly-reducing conditions, independent of the
    multi-species solver."""

    @pytest.mark.parametrize('dIW', [-3.0, -4.0, -6.0, -8.0])
    def test_dasgupta_law_finite_at_reducing_edge(self, dIW):
        """The reduced-N branch grows as exp(-1.6 dIW); at dIW=-8 the
        value reaches 10^5 ppmw at p_N2=1 bar but is still
        numerically representable."""
        N2 = SolubilityN2('dasgupta')
        val = N2.dasgupta(1.0, 1.0, 1800.0, dIW)
        assert math.isfinite(val)
        assert val > 0.0
