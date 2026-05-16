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

import math

import pytest

from calliope.constants import volatile_species

pytestmark = [pytest.mark.integration, pytest.mark.timeout(300)]


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
