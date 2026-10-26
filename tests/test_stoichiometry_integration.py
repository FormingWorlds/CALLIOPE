"""Integration-tier mass-conservation tests for `equilibrium_atmosphere`.

Run the full equilibrium-chemistry solver end-to-end and verify that the
per-element mass closure holds across the H/C/N/S inventories. These
tests carry a real solver call per test and live in the nightly tier;
the unit-tier stoichiometry checks (atom-by-atom tallies, equilibrium
constant identities, CH4 solubility pressure dependence) live in
`test_stoichiometry.py`.

The split is needed because pytest stacks module-level and class-level
markers; a single module-level pytestmark on a mixed-tier file would
pull the integration tests into the PR gate's unit selection.
"""

from __future__ import annotations

import inspect
import logging
import math
import re

import pytest

from calliope.constants import volatile_species
from calliope.solve import FALLBACK_ATTEMPTS, TRUNC_MASS, equilibrium_atmosphere

# 4 x the slowest CI duration of the module (rule in docs/Explanations/testing.md)
pytestmark = [pytest.mark.integration, pytest.mark.timeout(780)]


_DEFAULTS = inspect.signature(equilibrium_atmosphere).parameters
RTOL, ATOL = _DEFAULTS['rtol'].default, _DEFAULTS['atol'].default
NGUESS = 5000


class TestEquilibriumAtmosphereIntegration:
    """Run the full solver and verify mass conservation."""

    def _run_equilibrium(self, masses, T=2000.0, Phi=0.5, dIW=0.0):
        """Run equilibrium_atmosphere and return result dict."""
        import warnings

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
                nguess=NGUESS,
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
        """Full C-H-N-S system at S/H 3e-3 converges to a root that closes every element.

        S/H 3e-3 lies well above the ratio below which the restarts may not reach a root
        that keeps S (roughly 3e-4 for this case).
        """
        target = {'H': 1e20, 'C': 1e17, 'N': 1e17, 'S': 3e17}
        result = self._run_equilibrium(target, T=2500.0, Phi=1.0, dIW=0.0)
        assert not result['fallback']
        # Hirschmann buffer; the same solve on Fischer gives 0.35939 bar, on O'Neill 0.5503
        assert result['P_surf'] == pytest.approx(0.361194, rel=1e-4)
        for e, target_kg in target.items():
            assert abs(result[f'{e}_res']) <= max(target_kg * RTOL, TRUNC_MASS), e

    def test_trace_sulfur_far_below_hydrogen_returns_the_warned_fallback(self, caplog):
        """S 1e15 kg under H 1e20 kg (S/H 1e-5): the restarts reach no root that closes S.

        The solve returns the fallback after FALLBACK_ATTEMPTS more attempts; its warning
        names S alone, and the returned state passes the scalar CHNOS mass gate. The size of
        the S residual depends on round-off and is not pinned.
        """
        target = {'H': 1e20, 'C': 1e17, 'N': 1e17, 'S': 1e15}
        with caplog.at_level('WARNING', logger='fwl.calliope.solve'):
            result = self._run_equilibrium(target, T=2500.0, Phi=1.0, dIW=0.0)
        assert result['fallback'] is True
        warning = [m for m in caplog.messages if 'No root closes' in m]
        assert len(warning) == 1
        assert re.findall(r'(?:for|,) ([A-Z][a-z]?) \(', warning[0]) == ['S']
        attempts = int(re.search(r'after (\d+) attempts', warning[0]).group(1))
        # the cap ends the loop FALLBACK_ATTEMPTS after the first fallback, before nguess
        assert FALLBACK_ATTEMPTS < attempts < NGUESS
        # scalar CHNOS gate of the solve: largest budget x rtol + atol + TRUNC_MASS
        gate = max(target.values()) * RTOL + ATOL + TRUNC_MASS
        assert max(abs(result[f'{e}_res']) for e in 'HCNS') <= gate
        assert result['P_surf'] == pytest.approx(0.361192, rel=1e-4)

    def test_installation_page_example_prints_the_stated_pressure(self, capsys, caplog):
        """The solver example of docs/How-to/installation.md converges and prints 0.36 bar."""
        from pathlib import Path

        page = Path(__file__).resolve().parents[1] / 'docs' / 'How-to' / 'installation.md'
        blocks = re.findall(r'```python\n(.*?)```', page.read_text(), re.S)
        example = [b for b in blocks if 'equilibrium_atmosphere' in b]
        assert len(example) == 1
        namespace = {}
        with caplog.at_level('WARNING', logger='fwl.calliope.solve'):
            exec(example[0], namespace)
        assert capsys.readouterr().out.strip() == 'Surface pressure: 0.36 bar'
        assert not namespace['result']['fallback']
        assert namespace['result']['fO2_model'] == 'hirschmann'  # Fischer also prints 0.36
        assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []
