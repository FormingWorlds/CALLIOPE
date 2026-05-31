"""Solver-loop branch tests for ``equilibrium_atmosphere_authoritative_O``.

The authoritative-O entry point wraps a Monte-Carlo restart loop around
fsolve / trust-constr. Between the raw solver call and accepting a root
it applies three guards that the happy-path smoke tests never trip:

- a per-element residual gate (a converged-looking root whose mass
  residual exceeds tolerance is rejected),
- a physical-box gate (a root with a derived fO2_shift outside
  [-12, +12], a negative partial pressure, or a partial pressure above
  the 1e7 bar ceiling is rejected, because the production unbounded
  fsolve does not enforce ``bounds``),
- an exception firewall (a solver call or residual evaluation that
  raises is treated as a failed attempt, not a crash).

These tests drive each guard by replacing ``opt.fsolve`` and the residual
function ``func_authoritative_O`` so the controlled root reaches the
guard deterministically, in milliseconds, without depending on whether
the real solver happens to find that corner. The accept path then runs
the real chemistry (``_get_partial_pressures`` etc.) on the controlled
root, so the output assertions exercise real physics, not the stub.

Marked ``unit``: every test stubs the iterative solver, so each runs in
well under 100 ms.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from calliope import solve as calsolve
from calliope.constants import volatile_species
from calliope.solve import (
    equilibrium_atmosphere_authoritative_O,
    get_initial_pressures_with_fO2,
)

# The module logger is 'fwl.calliope.solve'; raise it to DEBUG so the
# rejection-branch debug lines reach caplog.
_SOLVE_LOGGER = 'fwl.calliope.solve'
logging.getLogger(_SOLVE_LOGGER).setLevel(logging.DEBUG)

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


def _stub_solver(monkeypatch, sol, residual):
    """Force the iterative solver to return ``sol`` with ``residual``.

    Patches the module objects directly (not via dotted-string targets):
    ``opt`` is an aliased ``scipy.optimize`` module, so a string target
    like ``'calliope.solve.opt.fsolve'`` fails to resolve and silently
    leaves the real solver running. ``func_authoritative_O`` is replaced
    so the residual gate sees a controlled value.

    Parameters
    ----------
    sol : sequence of 5 floats
        The (H2O, CO2, N2, S2, fO2_shift) root fsolve will report.
    residual : sequence of 5 floats
        The residual the gate will evaluate at ``sol``.
    """

    def _fake_fsolve(*args, **kwargs):
        # Real signature returns (x, infodict, ier, mesg); ier == 1 is
        # nominal convergence.
        return np.asarray(sol, dtype=float), {}, 1, 'stub converged'

    monkeypatch.setattr(calsolve.opt, 'fsolve', _fake_fsolve)
    monkeypatch.setattr(
        calsolve, 'func_authoritative_O', lambda *a, **k: np.asarray(residual, dtype=float)
    )


# ---------------------------------------------------------------------------
# Cold-start helper: default-RNG branch
# ---------------------------------------------------------------------------


def test_get_initial_pressures_with_fO2_default_rng():
    """``rng=None`` falls back to the global ``np.random`` state and
    returns a five-vector whose fO2 slot is the unredrawn hint.

    Exercises the ``rng is None`` default-RNG branch. With
    ``restart=False`` the fO2 seed must be the hint verbatim (the
    redraw only happens on restart), which discriminates this path from
    the restart path that draws fO2 from Uniform(-6, +8).
    """
    np.random.seed(42)  # determinism for the global-state branch
    hint = 3.5
    x0 = get_initial_pressures_with_fO2(_target(), hint, restart=False, rng=None)

    assert len(x0) == 5, f'expected 5-vector (4 pressures + fO2), got {len(x0)}'
    pressures = x0[:4]
    assert all(np.isfinite(p) for p in pressures), 'cold-start pressures must be finite'
    # 10 ** Uniform(-12, 5) lands strictly inside (1e-12, 1e5) bar.
    assert all(1e-13 < p < 1e6 for p in pressures), (
        f'pressures out of cold-start range: {pressures}'
    )

    # Discrimination guard: restart=False must pass the hint through
    # untouched. A restart draw would generically miss 3.5 exactly.
    assert x0[4] == pytest.approx(hint), (
        'restart=False must seed fO2 from the hint, not a redraw'
    )


# ---------------------------------------------------------------------------
# Accept path
# ---------------------------------------------------------------------------


def test_accept_path_returns_derived_fo2_and_real_chemistry(monkeypatch, caplog):
    """A converged, in-box, zero-residual root is accepted and the
    derived fO2 plus real outgassed pressures flow to the output.

    Forces fsolve to return a physical root and the residual function to
    report zero residual (so the residual + box gates pass), then lets
    the real chemistry build the output dict. The derived fO2 must equal
    the root's fifth component, NOT the hint, which is the whole point of
    authoritative-O mode.
    """
    sol = [10.0, 5.0, 1.0, 0.5, 1.5]  # H2O, CO2, N2, S2 [bar]; fO2_shift = 1.5
    _stub_solver(monkeypatch, sol, residual=np.zeros(5))

    caplog.set_level(logging.INFO, logger=_SOLVE_LOGGER)
    out = equilibrium_atmosphere_authoritative_O(
        _target(),
        _ddict(),
        fO2_hint=4.0,
        opt_solver=False,
        nguess=1,
        print_result=True,
    )

    # Derived fO2 is the solver output, not the hint.
    assert out['fO2_shift_derived'] == pytest.approx(1.5)
    # Discrimination guard: a buggy implementation that echoed the hint
    # would return 4.0; the gap is 2.5, far outside any rounding. (This
    # is also the symptom of the fsolve stub failing to apply.)
    assert abs(out['fO2_shift_derived'] - 4.0) > 1.0

    # The accepted primary pressures flow through the real chemistry into
    # the output unchanged (H2O is a primary, so p_d['H2O'] == sol[0]).
    assert out['H2O_bar'] == pytest.approx(10.0, rel=1e-9)
    # P_surf sums every species' partial pressure, so it is at least the
    # single H2O contribution and strictly positive.
    assert out['P_surf'] >= out['H2O_bar'] > 0.0

    # The fifth residual (O mass balance) is reported alongside the four
    # legacy ones.
    assert 'O_res' in out and np.isfinite(out['O_res'])
    # print_result=True logs the solve header and the derived-fO2 line.
    assert 'authoritative-O mode' in caplog.text
    assert 'Derived fO2_shift' in caplog.text


def test_p_guess_tiny_pressure_collapses_upper_bound(monkeypatch):
    """A sub-1e-10 bar entry in ``p_guess`` collapses that slot's upper
    bound to 1.0 before the solve, matching the legacy degenerate-slot
    guard.

    Exercises the ``p_guess`` ub-collapse loop. The accept path then
    confirms the run still completes and returns the controlled root, so
    the collapse does not corrupt an otherwise-valid solve.
    """
    sol = [1e-12, 5.0, 1.0, 0.5, 2.0]
    _stub_solver(monkeypatch, sol, residual=np.zeros(5))

    p_guess = {'H2O': 1e-12, 'CO2': 5.0, 'N2': 1.0, 'S2': 0.5}
    out = equilibrium_atmosphere_authoritative_O(
        _target(),
        _ddict(),
        fO2_hint=2.0,
        p_guess=p_guess,
        opt_solver=False,
        nguess=1,
        print_result=False,
    )

    # The solve completes through the collapsed-bound path and returns
    # the controlled root's derived fO2.
    assert out['fO2_shift_derived'] == pytest.approx(2.0)
    # Discrimination guard: the collapse must not have rewritten the
    # derived value toward the 1.0 bound it imposes on the pressure slot.
    assert out['fO2_shift_derived'] != pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Rejection gates: each must exhaust the restart budget and raise.
# ---------------------------------------------------------------------------


def test_root_with_out_of_box_fo2_is_rejected(monkeypatch, caplog):
    """A converged root whose derived fO2 lies outside [-12, +12] is
    rejected, not returned, because unbounded fsolve does not enforce
    the physical box.

    The single restart is consumed by the rejection, so the loop raises
    RuntimeError. The debug log names the out-of-box value so an operator
    can see why an apparently-converged solve was discarded.
    """
    sol = [10.0, 5.0, 1.0, 0.5, 20.0]  # fO2_shift = 20 > +12 bound
    _stub_solver(monkeypatch, sol, residual=np.zeros(5))

    caplog.set_level(logging.DEBUG, logger=_SOLVE_LOGGER)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere_authoritative_O(
            _target(), _ddict(), fO2_hint=4.0, opt_solver=False, nguess=1, print_result=False
        )

    assert 'outside' in caplog.text, 'out-of-box rejection must be logged'
    # Discrimination guard: confirm the root really was out of the box, so
    # the rejection path (not some other failure) fired.
    assert not (-12.0 <= sol[4] <= 12.0)


def test_root_with_negative_pressure_is_rejected(monkeypatch, caplog):
    """A converged root with a negative partial pressure is rejected.

    fsolve can return a slightly negative pressure on a degenerate slot;
    the box gate rejects it rather than passing a non-physical state into
    the chemistry. With one restart the loop then raises RuntimeError.
    """
    sol = [-1.0, 5.0, 1.0, 0.5, 2.0]  # H2O = -1 bar, in-box fO2
    _stub_solver(monkeypatch, sol, residual=np.zeros(5))

    caplog.set_level(logging.DEBUG, logger=_SOLVE_LOGGER)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere_authoritative_O(
            _target(), _ddict(), fO2_hint=2.0, opt_solver=False, nguess=1, print_result=False
        )

    assert 'negative partial pressure' in caplog.text
    # Discrimination guard: the fO2 is in-box, so ONLY the negative-pressure
    # branch (not the out-of-box branch) can have rejected this root.
    assert -12.0 <= sol[4] <= 12.0 and min(sol[:4]) < 0.0


def test_root_above_pressure_ceiling_is_rejected(monkeypatch, caplog):
    """A converged root with a partial pressure above the 1e7 bar ceiling
    is rejected, because unbounded fsolve does not enforce the upper box.

    A mass budget can drive fsolve to a root that closes the element
    balance yet places one partial pressure above the documented 1e7 bar
    bound that trust-constr's ``ub`` array would have enforced. The box
    gate rejects it; with one restart the loop raises RuntimeError instead
    of returning the non-physical state.
    """
    sol = [2.0e7, 5.0, 1.0, 0.5, 2.0]  # H2O = 2e7 bar > 1e7 ceiling, in-box fO2
    _stub_solver(monkeypatch, sol, residual=np.zeros(5))

    caplog.set_level(logging.DEBUG, logger=_SOLVE_LOGGER)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere_authoritative_O(
            _target(), _ddict(), fO2_hint=2.0, opt_solver=False, nguess=1, print_result=False
        )

    assert 'ceiling' in caplog.text, 'above-ceiling rejection must be logged'
    # Discrimination guard: fO2 is in-box and every pressure is non-negative,
    # so ONLY the upper-ceiling branch (not the out-of-box or negative-pressure
    # branches) can have rejected this root.
    assert -12.0 <= sol[4] <= 12.0 and min(sol[:4]) >= 0.0 and max(sol[:4]) > 1e7


def test_root_with_excess_residual_is_rejected(monkeypatch, caplog):
    """A root that fsolve flags converged but whose mass residual blows
    past tolerance is rejected by the per-element residual gate.

    The controlled residual exceeds the kg-scale tolerance for every
    element, so the worst-element check trips and the attempt fails. The
    single restart is consumed and the loop raises RuntimeError.
    """
    sol = [10.0, 5.0, 1.0, 0.5, 2.0]
    # Residual far above max(target*rtol, TRUNC_MASS) for every element.
    _stub_solver(monkeypatch, sol, residual=np.full(5, 1e30))

    caplog.set_level(logging.DEBUG, logger=_SOLVE_LOGGER)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere_authoritative_O(
            _target(), _ddict(), fO2_hint=2.0, opt_solver=False, nguess=1, print_result=False
        )

    assert 'residual' in caplog.text
    # Discrimination guard: 1e30 kg dwarfs the largest element tolerance
    # (Earth O budget 2e21 kg * rtol 1e-5 = 2e16 kg), so the gate must trip.
    assert 1e30 > max(1e10, 2.0e21 * 1e-5)


def test_residual_evaluation_exception_is_caught(monkeypatch, caplog):
    """A residual evaluation that raises is treated as a failed attempt,
    not a crash propagated to the caller.

    fsolve reports convergence, but the post-solve residual recompute
    raises ValueError. The firewall converts that into a rejected attempt;
    with one restart the loop ends in RuntimeError, never a bare
    ValueError escaping the function.
    """
    sol = [10.0, 5.0, 1.0, 0.5, 2.0]

    def _fake_fsolve(*a, **k):
        return np.asarray(sol, dtype=float), {}, 1, 'stub converged'

    def _raising_residual(*a, **k):
        raise ValueError('residual blew up at trial point')

    monkeypatch.setattr(calsolve.opt, 'fsolve', _fake_fsolve)
    monkeypatch.setattr(calsolve, 'func_authoritative_O', _raising_residual)

    caplog.set_level(logging.DEBUG, logger=_SOLVE_LOGGER)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere_authoritative_O(
            _target(), _ddict(), fO2_hint=2.0, opt_solver=False, nguess=1, print_result=False
        )

    # The firewall logged the caught exception; the bare ValueError did
    # not escape (pytest.raises above already enforces RuntimeError).
    assert 'residual evaluation raised' in caplog.text


def test_solver_call_exception_is_caught(monkeypatch, caplog):
    """A solver call that raises ZeroDivisionError is caught and retried,
    not propagated.

    fsolve itself raises mid-iteration; the firewall marks the attempt
    failed and the loop redraws. With one restart it ends in RuntimeError.
    """

    def _raising_fsolve(*a, **k):
        raise ZeroDivisionError('residual divided by zero at trial point')

    monkeypatch.setattr(calsolve.opt, 'fsolve', _raising_fsolve)

    caplog.set_level(logging.DEBUG, logger=_SOLVE_LOGGER)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere_authoritative_O(
            _target(), _ddict(), fO2_hint=2.0, opt_solver=False, nguess=1, print_result=False
        )

    assert 'restarting' in caplog.text
