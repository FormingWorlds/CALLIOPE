"""Input validation tests for `equilibrium_atmosphere_authoritative_O`.

Verifies that the entry point rejects bad inputs at the boundary
(raises a typed exception with a useful message) rather than letting
them propagate through fsolve and emerge as opaque ZeroDivisionError
or KeyError from deep inside the chemistry chain.

These are pure-validation tests: every case raises before any solver
call, so they run in < 100 ms each and are marked ``unit``.
"""

from __future__ import annotations

import logging

import pytest

from calliope.constants import volatile_species
from calliope.solve import equilibrium_atmosphere_authoritative_O

logging.getLogger('calliope').setLevel(logging.WARNING)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _ddict(T: float = 1800.0, Phi: float = 1.0) -> dict:
    """Realistic ddict with every volatile species included."""
    d = {
        'M_mantle': 4.03e24,
        'gravity': 9.81,
        'radius': 6.371e6,
        'Phi_global': Phi,
        'T_magma': T,
        'fO2_shift_IW': 0.0,  # ignored under authoritative-O mode
    }
    for sp in volatile_species:
        d[f'{sp}_included'] = 1
        d[f'{sp}_initial_bar'] = 0.0
    return d


def _target() -> dict:
    """Earth-like element budget [kg] with all five elements."""
    return {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20, 'O': 2.0e21}


def _base_kwargs() -> dict:
    """Minimal kwargs that skip the print log and use a small nguess."""
    return {'fO2_hint': 4.0, 'print_result': False, 'nguess': 10, 'nsolve': 100}


# ---------------------------------------------------------------------------
# target_d: required keys
# ---------------------------------------------------------------------------


class TestTargetDRequiredKeys:
    """target_d must include all five element keys (H, C, N, S, O).

    Missing keys raise KeyError; the error message names the missing
    keys so the caller can fix the input without grepping fsolve.
    """

    @pytest.mark.parametrize('missing', ['H', 'C', 'N', 'S', 'O'])
    def test_missing_element_raises_keyerror(self, missing):
        t = _target()
        del t[missing]
        with pytest.raises(KeyError, match=missing):
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **_base_kwargs())

    def test_extra_key_is_ignored(self):
        t = _target()
        t['Cl'] = 1e18  # not a tracked element
        # No exception expected at validation time. The solver itself
        # may or may not converge depending on physics, but the entry
        # point should accept extra keys silently.
        # We don't run the full solver; just check validation passes by
        # asserting fsolve actually starts.
        kwargs = _base_kwargs()
        kwargs['nguess'] = 1
        try:
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **kwargs)
        except RuntimeError:
            pass  # convergence failure is acceptable; validation passed.


class TestTargetDValueValidation:
    """Each target value must be a finite, non-negative real number."""

    @pytest.mark.parametrize('elem', ['H', 'C', 'N', 'S', 'O'])
    def test_negative_value_raises_valueerror(self, elem):
        t = _target()
        t[elem] = -1.0
        with pytest.raises(ValueError, match='non-negative'):
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **_base_kwargs())

    @pytest.mark.parametrize('elem', ['H', 'C', 'N', 'S', 'O'])
    def test_nan_value_raises_valueerror(self, elem):
        t = _target()
        t[elem] = float('nan')
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **_base_kwargs())

    @pytest.mark.parametrize('elem', ['H', 'C', 'N', 'S', 'O'])
    def test_inf_value_raises_valueerror(self, elem):
        t = _target()
        t[elem] = float('inf')
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **_base_kwargs())


# ---------------------------------------------------------------------------
# fO2_hint
# ---------------------------------------------------------------------------


class TestFO2HintValidation:
    """fO2_hint must be finite and within the solver bounds [-12, +12]."""

    def test_nan_raises_valueerror(self):
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = float('nan')
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_inf_raises_valueerror(self):
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = float('inf')
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    @pytest.mark.parametrize('hint', [-100.0, -13.0, 12.5, 100.0])
    def test_outside_bounds_raises_valueerror(self, hint):
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = hint
        with pytest.raises(ValueError, match='\\[-12, \\+12\\]'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    @pytest.mark.parametrize('hint', [-12.0, -6.0, 0.0, 4.0, 8.0, 12.0])
    def test_inside_bounds_accepted(self, hint):
        """Valid hints reach the solver loop (may not converge in 1
        attempt, but the entry-point check must not raise)."""
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = hint
        kwargs['nguess'] = 1
        try:
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)
        except RuntimeError:
            pass  # convergence failure is acceptable; validation passed.
        # No ValueError raised: validation passed.


# ---------------------------------------------------------------------------
# ddict: planet/state parameters
# ---------------------------------------------------------------------------


class TestDdictValidation:
    """M_mantle, Phi_global, T_magma must be present and physical."""

    @pytest.mark.parametrize('key', ['M_mantle', 'Phi_global', 'T_magma', 'gravity', 'radius'])
    def test_missing_required_key_raises_keyerror(self, key):
        dd = _ddict()
        del dd[key]
        with pytest.raises(KeyError, match=key):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())

    @pytest.mark.parametrize('M_mantle', [0.0, -1e24, float('nan'), float('inf')])
    def test_bad_M_mantle_raises_valueerror(self, M_mantle):
        dd = _ddict()
        dd['M_mantle'] = M_mantle
        with pytest.raises(ValueError, match='M_mantle'):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())

    @pytest.mark.parametrize('Phi', [-0.1, 1.5, 2.0, float('nan'), float('inf')])
    def test_bad_Phi_global_raises_valueerror(self, Phi):
        dd = _ddict()
        dd['Phi_global'] = Phi
        with pytest.raises(ValueError, match='Phi_global'):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())

    @pytest.mark.parametrize('T', [0.0, -300.0, float('nan'), float('inf')])
    def test_bad_T_magma_raises_valueerror(self, T):
        dd = _ddict()
        dd['T_magma'] = T
        with pytest.raises(ValueError, match='T_magma'):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())

    def test_T_magma_outside_calibration_warns(self, recwarn):
        """T_magma outside Dasgupta/Gaillard calibration emits a single
        UserWarning so the user notices the solver is extrapolating."""
        dd = _ddict(T=3000.0)
        kwargs = _base_kwargs()
        kwargs['nguess'] = 1
        try:
            equilibrium_atmosphere_authoritative_O(_target(), dd, **kwargs)
        except RuntimeError:
            pass
        msgs = [str(w.message) for w in recwarn]
        assert any('calibrated range' in m and 'extrapolate' in m for m in msgs), (
            f'expected calibration-range warning, got: {msgs}'
        )


# ---------------------------------------------------------------------------
# Solver-parameter bounds
# ---------------------------------------------------------------------------


class TestSolverParamValidation:
    """nguess and nsolve must be >= 1."""

    @pytest.mark.parametrize('nguess', [0, -1, -100])
    def test_bad_nguess_raises_valueerror(self, nguess):
        kwargs = _base_kwargs()
        kwargs['nguess'] = nguess
        with pytest.raises(ValueError, match='nguess'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    @pytest.mark.parametrize('nsolve', [0, -1])
    def test_bad_nsolve_raises_valueerror(self, nsolve):
        kwargs = _base_kwargs()
        kwargs['nsolve'] = nsolve
        with pytest.raises(ValueError, match='nsolve'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_nguess_zero_not_unbound_local_error(self):
        """nguess=0 must raise ValueError, never UnboundLocalError.

        The legacy bug was that `success` was bound inside the restart
        loop; a 0-iteration loop never bound it, and the post-loop check
        raised UnboundLocalError instead of the documented RuntimeError.
        """
        kwargs = _base_kwargs()
        kwargs['nguess'] = 0
        with pytest.raises(ValueError):
            try:
                equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)
            except UnboundLocalError:
                pytest.fail('nguess=0 raised UnboundLocalError instead of ValueError')


# ---------------------------------------------------------------------------
# p_guess
# ---------------------------------------------------------------------------


class TestPGuessValidation:
    """p_guess is optional but if supplied must be a dict with the four
    primary keys and finite values."""

    def test_p_guess_is_list_raises_typeerror(self):
        kwargs = _base_kwargs()
        kwargs['p_guess'] = [1.0, 1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match='dict'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_p_guess_empty_dict_raises_valueerror(self):
        kwargs = _base_kwargs()
        kwargs['p_guess'] = {}
        with pytest.raises(ValueError, match='missing required keys'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_p_guess_with_nan_raises_valueerror(self):
        kwargs = _base_kwargs()
        kwargs['p_guess'] = {'H2O': float('nan'), 'CO2': 1.0, 'N2': 1.0, 'S2': 1.0}
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)
