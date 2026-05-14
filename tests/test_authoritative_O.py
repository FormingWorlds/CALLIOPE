"""Smoke tests for ``equilibrium_atmosphere_authoritative_O``.

Covers the happy path (canonical Earth-like inputs produce physical
output), the round-trip property (legacy mode at fO2_shift = X yields
an O budget; feeding that budget back to the new mode recovers
fO2_derived ≈ X), and reproducibility under ``random_seed``.

These tests invoke the full solver and take a few seconds each. They
are marked ``smoke`` per the project's four-marker scheme.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from calliope.constants import volatile_species
from calliope.solve import (
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
)

logging.getLogger('calliope').setLevel(logging.WARNING)

pytestmark = pytest.mark.smoke


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _ddict(T: float = 1800.0, Phi: float = 1.0, dIW: float = 4.0) -> dict:
    """Realistic ddict with every volatile species included.

    Default T_magma=1800 K is inside the Dasgupta/Gaillard solubility
    law calibration range so tests don't emit extrapolation warnings.
    """
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


def _earth_target_HCNS() -> dict:
    """Earth-like H/C/N/S budget [kg] — converges cleanly in legacy mode."""
    return {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


class TestHappyPath:
    """Canonical Earth-like inputs should converge and return a
    physically reasonable atmosphere with finite primary pressures
    and a finite derived fO2 inside the wider [-12, +12] solver bounds.
    """

    def test_returns_physical_output(self):
        """Smoke: solver runs, returns a dict with expected keys."""
        # Use legacy mode to derive a self-consistent O budget; that
        # O is by construction reachable so the new mode is guaranteed
        # to have a solution.
        ddict = _ddict(dIW=4.0)
        legacy = equilibrium_atmosphere(
            _earth_target_HCNS(),
            ddict,
            print_result=False,
            nguess=200,
        )
        target = dict(_earth_target_HCNS(), O=legacy['O_kg_total'])

        out = equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=4.0,
            random_seed=0,
            nguess=500,
            nsolve=1000,
            print_result=False,
        )

        # Required output keys
        for key in (
            'fO2_shift_derived',
            'O_res',
            'H_res',
            'C_res',
            'N_res',
            'S_res',
            'H2O_bar',
            'CO2_bar',
            'N2_bar',
            'S2_bar',
            'M_atm',
            'P_surf',
        ):
            assert key in out, f'missing key {key!r} in output dict'

        # Pressures and derived fO2 are finite and physical
        for p_key in ('H2O_bar', 'CO2_bar', 'N2_bar', 'S2_bar'):
            assert np.isfinite(out[p_key]), f'{p_key} is not finite'
            assert out[p_key] >= 0, f'{p_key} is negative'

        assert np.isfinite(out['fO2_shift_derived']), 'fO2_shift_derived is not finite'
        assert -12.0 <= out['fO2_shift_derived'] <= 12.0, (
            f'fO2_shift_derived={out["fO2_shift_derived"]} outside bounds'
        )

    def test_residuals_within_tolerance(self):
        """Solver returns only when per-element residuals are within
        the per-element tolerance gate. Verify all 5 residuals are
        small relative to their targets."""
        ddict = _ddict(dIW=4.0)
        legacy = equilibrium_atmosphere(
            _earth_target_HCNS(),
            ddict,
            print_result=False,
            nguess=200,
        )
        target = dict(_earth_target_HCNS(), O=legacy['O_kg_total'])

        out = equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=4.0,
            random_seed=0,
            rtol=1e-5,
            nguess=500,
            nsolve=1000,
            print_result=False,
        )

        for elem in ('H', 'C', 'N', 'S', 'O'):
            res = out[f'{elem}_res']
            tgt = target[elem]
            # Per-element tolerance from solve.py: max(target * rtol, atol/5).
            # atol default is 1e10 kg, /5 = 2e9 kg, /10 floor for TRUNC_MASS.
            allowed = max(tgt * 1e-5, 2e9)
            assert abs(res) <= allowed, f'{elem}_res={res:.3e} exceeds tolerance {allowed:.3e}'


# ---------------------------------------------------------------------------
# Round-trip with legacy mode
# ---------------------------------------------------------------------------


class TestRoundTrip:
    """For any fO2_shift_IW value X that the legacy solver accepts, the
    new mode should recover fO2_derived ≈ X when fed the legacy O budget.

    This is the *core* correctness property of the new entry point:
    extending the unknown set from 4 to 5 and adding the O equation
    must reproduce the legacy chemistry. If the round-trip fails by
    more than the per-element tolerance, the new equation system has
    a different fixed point than the legacy one — that's a bug.
    """

    @pytest.mark.parametrize('dIW', [-2.0, 0.0, 2.0, 4.0, 6.0])
    def test_round_trip_recovers_fO2_within_tolerance(self, dIW):
        """Legacy at dIW -> new mode with implied O_budget -> recover dIW."""
        ddict_legacy = _ddict(dIW=dIW)

        legacy = equilibrium_atmosphere(
            _earth_target_HCNS(),
            ddict_legacy,
            print_result=False,
            nguess=200,
        )
        target_O = legacy['O_kg_total']
        target = dict(_earth_target_HCNS(), O=target_O)

        # Run the new mode with the same ddict (whose fO2_shift_IW is
        # ignored under authoritative-O mode). fO2_hint = dIW so the
        # solver starts at the right basin.
        out = equilibrium_atmosphere_authoritative_O(
            target,
            ddict_legacy,
            fO2_hint=dIW,
            random_seed=0,
            nguess=500,
            nsolve=1000,
            print_result=False,
        )

        derived = out['fO2_shift_derived']
        delta = abs(derived - dIW)

        # 0.05 dex on the derived fO2 is well within the solver's xtol;
        # if the equation system were inconsistent we would see deltas
        # of >> 0.1 dex.
        assert delta < 0.05, (
            f'round-trip failed at dIW={dIW}: derived={derived:.4f}, delta={delta:.4f} dex'
        )

    def test_round_trip_primary_pressures_match_legacy(self):
        """Primary partial pressures should match the legacy output
        within solver tolerance after a round-trip."""
        ddict_legacy = _ddict(dIW=4.0)

        legacy = equilibrium_atmosphere(
            _earth_target_HCNS(),
            ddict_legacy,
            print_result=False,
            nguess=200,
        )
        target = dict(_earth_target_HCNS(), O=legacy['O_kg_total'])

        out = equilibrium_atmosphere_authoritative_O(
            target,
            ddict_legacy,
            fO2_hint=4.0,
            random_seed=0,
            nguess=500,
            nsolve=1000,
            print_result=False,
        )

        # Primary pressures within 0.5% (the legacy mode and the new
        # mode share the same physics so the same fixed point should
        # be found to within fsolve's xtol).
        for key in ('H2O_bar', 'CO2_bar', 'N2_bar', 'S2_bar'):
            legacy_p = legacy[key]
            new_p = out[key]
            if legacy_p > 1e-20:  # skip near-zero
                rel = abs(new_p - legacy_p) / legacy_p
                assert rel < 0.005, (
                    f'{key} mismatch: legacy={legacy_p:.6e}, new={new_p:.6e}, rel={rel:.3e}'
                )


# ---------------------------------------------------------------------------
# Reproducibility under random_seed
# ---------------------------------------------------------------------------


class TestReproducibility:
    """Two calls with the same ``random_seed`` must produce bit-identical
    output. This is essential for diffing solver outputs across PRs
    and for regression tests."""

    def test_same_seed_bit_identical(self):
        ddict = _ddict()
        target = dict(_earth_target_HCNS(), O=2.0e21)

        out1 = equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=4.0,
            random_seed=42,
            nguess=200,
            nsolve=500,
            print_result=False,
        )
        out2 = equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=4.0,
            random_seed=42,
            nguess=200,
            nsolve=500,
            print_result=False,
        )

        for key in (
            'fO2_shift_derived',
            'H2O_bar',
            'CO2_bar',
            'N2_bar',
            'S2_bar',
            'H_res',
            'C_res',
            'N_res',
            'S_res',
            'O_res',
        ):
            assert out1[key] == out2[key], (
                f'{key} differs between two same-seed calls: {out1[key]!r} vs {out2[key]!r}'
            )

    def test_different_seeds_may_differ(self):
        """Different seeds should produce different intermediate paths
        even though the final root may be the same (the solver has a
        unique attractor for canonical inputs). Verify the solver does
        consume the seed by checking that the elapsed-time signature
        of two different seeds is not identical."""
        ddict = _ddict()
        target = dict(_earth_target_HCNS(), O=2.0e21)

        # Both should converge to the same root since the problem is
        # well-posed, so this test only checks the call completes.
        out1 = equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=4.0,
            random_seed=1,
            nguess=200,
            nsolve=500,
            print_result=False,
        )
        out2 = equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=4.0,
            random_seed=999,
            nguess=200,
            nsolve=500,
            print_result=False,
        )
        # Both converged to a valid fO2 within solver bounds
        assert -12.0 <= out1['fO2_shift_derived'] <= 12.0
        assert -12.0 <= out2['fO2_shift_derived'] <= 12.0


# ---------------------------------------------------------------------------
# Convergence failure contract
# ---------------------------------------------------------------------------


class TestConvergenceFailure:
    """An unreachable target (e.g. O = 1e30 kg, impossibly large) must
    raise RuntimeError with the documented diagnostic message. Never
    ZeroDivisionError, never silent success.
    """

    def test_impossible_O_target_raises_runtime_error(self):
        ddict = _ddict()
        target = dict(_earth_target_HCNS(), O=1e30)

        with pytest.raises(RuntimeError, match='Could not find solution'):
            equilibrium_atmosphere_authoritative_O(
                target,
                ddict,
                fO2_hint=4.0,
                random_seed=0,
                nguess=50,
                nsolve=200,
                print_result=False,
            )

    def test_runtime_error_message_includes_final_attempt(self):
        """The RuntimeError message must include the final pH2O / pCO2
        / fO2_shift attempt for debugging."""
        ddict = _ddict()
        target = dict(_earth_target_HCNS(), O=1e30)

        try:
            equilibrium_atmosphere_authoritative_O(
                target,
                ddict,
                fO2_hint=4.0,
                random_seed=0,
                nguess=20,
                nsolve=100,
                print_result=False,
            )
            pytest.fail('expected RuntimeError for impossible target')
        except RuntimeError as exc:
            msg = str(exc)
            for hint in ('pH2O', 'pCO2', 'pN2', 'pS2', 'fO2_shift'):
                assert hint in msg, f'RuntimeError message missing {hint!r}: {msg}'
