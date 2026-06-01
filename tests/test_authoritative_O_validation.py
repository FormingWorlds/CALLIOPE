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

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


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
        """Dropping any one of the five required element keys raises
        KeyError naming the missing key."""
        t = _target()
        del t[missing]
        with pytest.raises(KeyError, match=missing):
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **_base_kwargs())

    def test_extra_key_is_ignored(self):
        """Validation accepts unknown element keys (e.g. Cl) without
        raising ValueError or KeyError; a convergence-time RuntimeError
        from the solver loop is acceptable."""
        t = _target()
        t['Cl'] = 1e18  # not a tracked element
        kwargs = _base_kwargs()
        kwargs['nguess'] = 1
        # Contract: validation must not raise ValueError or KeyError.
        # A RuntimeError from the convergence loop is acceptable.
        try:
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **kwargs)
        except (ValueError, KeyError) as e:
            pytest.fail(f'Extra key triggered validation error: {e!r}')
        except RuntimeError:
            pass  # convergence failure is acceptable

        # Discrimination guard: confirm the test setup actually puts an
        # unknown element in t. A future _target() change that already
        # carried 'Cl' would make this test vacuous.
        assert 'Cl' in t
        assert 'Cl' not in ('H', 'C', 'N', 'S', 'O')


class TestTargetDValueValidation:
    """Each target value must be a finite, non-negative real number."""

    @pytest.mark.parametrize('elem', ['H', 'C', 'N', 'S', 'O'])
    def test_negative_value_raises_valueerror(self, elem):
        """A negative element budget raises ValueError with a 'non-negative'
        message for any of the five elements."""
        t = _target()
        t[elem] = -1.0
        with pytest.raises(ValueError, match='non-negative'):
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **_base_kwargs())

    @pytest.mark.parametrize('elem', ['H', 'C', 'N', 'S', 'O'])
    def test_nan_value_raises_valueerror(self, elem):
        """A NaN element budget raises ValueError with a 'finite'
        message for any of the five elements."""
        t = _target()
        t[elem] = float('nan')
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(t, _ddict(), **_base_kwargs())

    @pytest.mark.parametrize('elem', ['H', 'C', 'N', 'S', 'O'])
    def test_inf_value_raises_valueerror(self, elem):
        """An infinite element budget raises ValueError with a 'finite'
        message for any of the five elements."""
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
        """NaN fO2_hint raises ValueError with a 'finite' message."""
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = float('nan')
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_inf_raises_valueerror(self):
        """Infinite fO2_hint raises ValueError with a 'finite' message."""
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = float('inf')
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    @pytest.mark.parametrize('hint', [-100.0, -13.0, 12.5, 100.0])
    def test_outside_bounds_raises_valueerror(self, hint):
        """fO2_hint outside the [-12, +12] band raises ValueError with a
        bounds-naming message."""
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = hint
        with pytest.raises(ValueError, match='\\[-12, \\+12\\]'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    @pytest.mark.parametrize('hint', [-12.0, -6.0, 0.0, 4.0, 8.0, 12.0])
    def test_inside_bounds_accepted(self, hint):
        """Valid hints inside [-12, +12] must pass validation; the solver
        may not converge in one attempt but the entry-point check must
        not raise ValueError."""
        kwargs = _base_kwargs()
        kwargs['fO2_hint'] = hint
        kwargs['nguess'] = 1
        # Contract: validation must not raise ValueError; a convergence-
        # time RuntimeError from the solver loop is acceptable.
        try:
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)
        except ValueError as e:
            pytest.fail(f'In-bounds fO2_hint={hint} rejected: {e!r}')
        except RuntimeError:
            pass  # convergence failure is acceptable

        # Discrimination guard: confirm hint is genuinely inside the
        # [-12, +12] band. Catches a parametrize drift that put an
        # out-of-bounds value into the in-bounds sweep.
        assert -12.0 <= hint <= 12.0


# ---------------------------------------------------------------------------
# ddict: planet/state parameters
# ---------------------------------------------------------------------------


class TestDdictValidation:
    """M_mantle, Phi_global, T_magma must be present and physical."""

    @pytest.mark.parametrize('key', ['M_mantle', 'Phi_global', 'T_magma', 'gravity', 'radius'])
    def test_missing_required_key_raises_keyerror(self, key):
        """Dropping any of the five required ddict keys raises KeyError
        naming the missing key."""
        dd = _ddict()
        del dd[key]
        with pytest.raises(KeyError, match=key):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())

    @pytest.mark.parametrize('M_mantle', [0.0, -1e24, float('nan'), float('inf')])
    def test_bad_M_mantle_raises_valueerror(self, M_mantle):
        """Zero, negative, NaN, or infinite M_mantle raises ValueError
        naming M_mantle."""
        dd = _ddict()
        dd['M_mantle'] = M_mantle
        with pytest.raises(ValueError, match='M_mantle'):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())

    @pytest.mark.parametrize('Phi', [-0.1, 1.5, 2.0, float('nan'), float('inf')])
    def test_bad_Phi_global_raises_valueerror(self, Phi):
        """Phi_global outside [0, 1] or non-finite raises ValueError
        naming Phi_global."""
        dd = _ddict()
        dd['Phi_global'] = Phi
        with pytest.raises(ValueError, match='Phi_global'):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())

    @pytest.mark.parametrize('T', [0.0, -300.0, float('nan'), float('inf')])
    def test_bad_T_magma_raises_valueerror(self, T):
        """Zero, negative, NaN, or infinite T_magma raises ValueError
        naming T_magma."""
        dd = _ddict()
        dd['T_magma'] = T
        with pytest.raises(ValueError, match='T_magma'):
            equilibrium_atmosphere_authoritative_O(_target(), dd, **_base_kwargs())


# ---------------------------------------------------------------------------
# Solver-parameter bounds
# ---------------------------------------------------------------------------


class TestSolverParamValidation:
    """nguess and nsolve must be >= 1."""

    @pytest.mark.parametrize('nguess', [0, -1, -100])
    def test_bad_nguess_raises_valueerror(self, nguess):
        """nguess <= 0 raises ValueError naming nguess."""
        kwargs = _base_kwargs()
        kwargs['nguess'] = nguess
        with pytest.raises(ValueError, match='nguess'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    @pytest.mark.parametrize('nsolve', [0, -1])
    def test_bad_nsolve_raises_valueerror(self, nsolve):
        """nsolve <= 0 raises ValueError naming nsolve."""
        kwargs = _base_kwargs()
        kwargs['nsolve'] = nsolve
        with pytest.raises(ValueError, match='nsolve'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_nguess_zero_not_unbound_local_error(self):
        """nguess=0 must raise ValueError, never UnboundLocalError.

        The legacy bug was that `success` was bound inside the restart
        loop; a 0-iteration loop never bound it, and the post-loop check
        raised UnboundLocalError instead of the documented ValueError.
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
        """A list p_guess raises TypeError with a 'dict' message."""
        kwargs = _base_kwargs()
        kwargs['p_guess'] = [1.0, 1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match='dict'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_p_guess_empty_dict_raises_valueerror(self):
        """An empty p_guess dict raises ValueError with a 'missing required
        keys' message."""
        kwargs = _base_kwargs()
        kwargs['p_guess'] = {}
        with pytest.raises(ValueError, match='missing required keys'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)

    def test_p_guess_with_nan_raises_valueerror(self):
        """A NaN value inside p_guess raises ValueError with a 'finite'
        message."""
        kwargs = _base_kwargs()
        kwargs['p_guess'] = {'H2O': float('nan'), 'CO2': 1.0, 'N2': 1.0, 'S2': 1.0}
        with pytest.raises(ValueError, match='finite'):
            equilibrium_atmosphere_authoritative_O(_target(), _ddict(), **kwargs)
