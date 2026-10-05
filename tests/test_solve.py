"""Tests for `src/calliope/solve.py`.

Exercises the public API of the equilibrium-atmosphere solver:
`equilibrium_atmosphere` (the conventional forward call with
user-supplied `fO2_shift_IW`) and `equilibrium_atmosphere_authoritative_O`
(the inverse call with user-supplied `O_kg_total`).

`solve.py` is the largest CALLIOPE source file (>1200 LOC) and its
test surface is split across this file and several **topical
cross-cutting** files for readability:

- `tests/test_authoritative_O.py` and the two siblings
  (`test_authoritative_O_monotonicity.py`,
  `test_authoritative_O_validation.py`) cover the authoritative-O
  entry point's contract, monotonicity properties, and input validation.
- `tests/test_equilibrium_paths.py` covers the forward solver's
  behaviour on multi-species compositions.
- `tests/test_partial_species.py` covers the partial-species (only-
  some-elements-included) branches.
- `tests/test_stoichiometry.py` covers stoichiometric ratios across
  the published reactions.
- `tests/test_targets.py` covers the target-element-budget
  computation that feeds both entry points.
- `tests/test_invariants.py` covers per-element / per-species
  closure invariants.
- `tests/test_invariants_hypothesis.py` covers the property-based
  fuzz tests at the slow tier.

This file is the **primary per-source test file** required by the
1:1 mirroring rule. It contains the reference-pinned anchor (round-
trip self-consistency at the Earth fiducial) plus a small set of
physics_invariant smoke tests that exercise both entry points.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

import calliope.solve as solve_mod
from calliope.constants import volatile_species
from calliope.solve import (
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
)

logging.getLogger('calliope').setLevel(logging.WARNING)

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]


def _earth_ddict(T: float = 1800.0, Phi: float = 1.0, dIW: float = 2.0) -> dict:
    """Earth-like input dict with every volatile species included.

    `T = 1800 K` is inside the Dasgupta / Gaillard solubility calibration
    range so the solver does not emit extrapolation warnings.
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
    """Earth-like H / C / N / S budget [kg] used by the legacy and the
    authoritative-O entry points; converges cleanly at fO2 = IW + 2."""
    return {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_round_trip_self_consistency_at_earth_fiducial():
    """Forward solve at `fO2 = IW + 2` and inverse solve from the
    resulting O budget recover the original `fO2_shift_IW_derived`
    within 0.05 dex.

    Anchor type: cross-implementation cross-check. CALLIOPE implements
    two distinct entry points into the equilibrium-chemistry solver:

    - `equilibrium_atmosphere` (legacy): user supplies `fO2_shift_IW`,
      solver returns species kg.
    - `equilibrium_atmosphere_authoritative_O` (authoritative-O entry point): user supplies
      the target `O_kg_total`, solver inverts to find the `fO2_shift_IW`
      that matches it.

    The round-trip is the contract: feeding the forward-mode `O_kg_total`
    output into the authoritative-O entry point must recover the original
    `fO2_shift_IW` within the documented solver tolerance.

    Discrimination guard: a regression that broke either the forward
    O mass-balance or the inverse bisection would lose the round-trip
    within 0.1 dex. The 0.05 dex envelope is half of that, so a
    coefficient-only bug would fail loudly.

    Hidden coupling: the pin uses the Earth-fiducial input
    (`T_magma = 1800 K`, `Phi = 1.0`, all volatile species included)
    with the default Fischer 2011 IW buffer.
    """
    fO2_shift = 2.0
    ddict = _earth_ddict(dIW=fO2_shift)
    target = _earth_target_HCNS()

    # Forward solve at fO2 = IW + 2 to get the resulting O_kg_total.
    legacy_result = equilibrium_atmosphere(
        target,
        ddict,
        print_result=False,
        nguess=200,
    )
    O_kg_total_forward = legacy_result['O_kg_total']
    assert O_kg_total_forward > 0
    # Scale guard: O mass at Earth-fiducial inputs is order 1e20 kg.
    assert 1e18 < O_kg_total_forward < 1e22

    # Inverse solve: target the forward-mode O budget; recover fO2_shift.
    target_with_O = dict(target, O=O_kg_total_forward)
    inverse_result = equilibrium_atmosphere_authoritative_O(
        target_with_O,
        ddict,
        fO2_hint=fO2_shift,
        random_seed=0,
        nguess=500,
        nsolve=1000,
        print_result=False,
    )
    fO2_derived = inverse_result['fO2_shift_derived']
    # Round-trip: |fO2_derived - 2.0| < 0.05 dex (well within solver xtol).
    assert fO2_derived == pytest.approx(fO2_shift, abs=0.05)


@pytest.mark.physics_invariant
def test_equilibrium_atmosphere_mass_closure_at_earth_fiducial():
    """Per-element mass closure: `sum(species_kg_total)` recovers the
    input H, C, N, S budgets within solver tolerance.

    The forward solver must conserve every input element; a regression
    that lost a species in the post-solve aggregation would violate
    this. Tolerance `rel=1e-3` is chosen for the SciPy nonlinear-solver
    convergence floor of `1e-8` on the residuals (the species sums are
    `1e20` kg-scale, so `1e-8 * 1e20 = 1e12` absolute, well under the
    1e-3 relative).
    """
    ddict = _earth_ddict(dIW=2.0)
    target = _earth_target_HCNS()
    result = equilibrium_atmosphere(target, ddict, print_result=False, nguess=200)

    # Per-element closure: sum of species masses that contain element E
    # must equal the input budget for E.
    # Hydrogen-bearing: H2O, H2, CH4, H2S, NH3.
    h_recovered = (
        result['H2O_kg_total'] * 2 / 18.015
        + result['H2_kg_total'] * 2 / 2.016
        + result['CH4_kg_total'] * 4 / 16.04
        + result['H2S_kg_total'] * 2 / 34.08
        + result['NH3_kg_total'] * 3 / 17.03
    ) * 1.008  # convert moles-of-H back to kg
    assert h_recovered == pytest.approx(target['H'], rel=1e-3)
    # Carbon-bearing: CO2, CO, CH4. Independent invariant on the C
    # budget; a regression that lost only one of the H species would
    # not necessarily corrupt the C closure, and vice versa.
    c_recovered = (
        result['CO2_kg_total'] * 1 / 44.01
        + result['CO_kg_total'] * 1 / 28.01
        + result['CH4_kg_total'] * 1 / 16.04
    ) * 12.011  # convert moles-of-C back to kg
    assert c_recovered == pytest.approx(target['C'], rel=1e-3)


@pytest.mark.parametrize('dIW', [-2.0, 0.0, 4.0])
def test_equilibrium_atmosphere_returns_positive_O_kg_total(dIW):
    """`O_kg_total` is positive at every dIW in the realistic
    {-2, 0, +4} range.

    Boundedness check: the solver must not return zero or negative
    oxygen mass for any physically valid input. A regression that
    introduced a clip-to-zero on a negative intermediate would fail
    this; a sign flip on the O mass-balance equation would also fail.
    """
    ddict = _earth_ddict(dIW=dIW)
    target = _earth_target_HCNS()
    result = equilibrium_atmosphere(target, ddict, print_result=False, nguess=200)
    O_kg = result['O_kg_total']
    assert O_kg > 0
    # Scale guard: O mass stays within [1e16, 1e23] kg over the dIW
    # range; the upper bound catches a unit-conversion bug, the lower
    # bound catches a clip-to-near-zero regression.
    assert 1e16 < O_kg < 1e23


# Low-H, C-rich inventory [kg]: besides the physical root the CHNOS residual
# has a root with pH2O near -1.7e4 bar that closes H through CH4 and H2S.
_LOW_H = {'H': 1.0e19, 'C': 5.0e19, 'N': 1.0e17, 'S': 1.0e18}
_NO_N = {'H': 1.0e20, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}
_COLD = dict(
    xtol=1e-6,
    rtol=1e-4,
    atol=1e16,
    nguess=1000,
    nsolve=3000,
    print_result=False,
    opt_solver=False,
)


def _spy_roots(monkeypatch):
    """Record every root fsolve returns, in call order, and whether it converged."""
    roots, converged, real = [], [], solve_mod.opt.fsolve

    def spy(*args, **kwargs):
        out = real(*args, **kwargs)
        roots.append(np.array(out[0]))
        converged.append(out[2] == 1)
        return out

    monkeypatch.setattr(solve_mod.opt, 'fsolve', spy)
    return roots, converged


@pytest.mark.physics_invariant
@pytest.mark.parametrize('seed', [0, 3, 42])
def test_low_h_cold_start_rejects_the_negative_water_root(monkeypatch, seed):
    """A low-H, C-rich cold start at IW+2 and 1500 K returns the physical root.

    Seeds 0, 3 and 42 reach the root with negative pH2O, which closes the mass
    balance at 18.58 bar with no H2O or H2 in the atmosphere. It is rejected,
    and the solver goes on to the physical root at 24.82 bar.
    """
    roots, converged = _spy_roots(monkeypatch)
    np.random.seed(seed)
    r = equilibrium_atmosphere(
        dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **_COLD
    )
    assert any(ok and root[0] < -1.0e3 for root, ok in zip(roots[:-1], converged[:-1]))
    assert roots[-1][0] > 0.0
    assert r['H2O_bar'] > 0.0
    assert r['P_surf'] == pytest.approx(24.8219, rel=1e-4)
    assert r['H_kg_total'] == pytest.approx(_LOW_H['H'], rel=1e-3)


@pytest.mark.physics_invariant
def test_absent_element_root_keeps_its_inert_negative_primary(monkeypatch):
    """With no N, fsolve leaves pN2 negative along a flat residual direction.

    Every N species is clipped to zero, so the root still closes the mass
    balance and is accepted.
    """
    roots, converged = _spy_roots(monkeypatch)
    np.random.seed(0)
    r = equilibrium_atmosphere(
        dict(_NO_N), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **_COLD
    )
    assert roots[-1][2] < 0.0
    np.testing.assert_array_equal([r['N2_bar'], r['NH3_bar']], 0.0)
    assert r['H2O_bar'] == pytest.approx(0.178874, rel=1e-4)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('target', 'seed'),
    [
        ({'H': 0.0, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}, 0),
        ({'H': 0.0, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}, 2),
        ({'H': 1.0e15, 'C': 5.0e19, 'N': 1.0e17, 'S': 1.0e18}, 0),
        ({'H': 1.0e15, 'C': 5.0e19, 'N': 1.0e17, 'S': 1.0e18}, 1),
    ],
    ids=['no_H_no_N-0', 'no_H_no_N-2', 'trace_H-0', 'trace_H-1'],
)
def test_sub_gate_hydrogen_reports_no_phantom_species(monkeypatch, target, seed):
    """An H budget below the scalar mass gate (1.5e16 kg here) yields a state
    free of H-bearing species from a negative pH2O, C and S close on their own
    tolerance rather than on the scalar gate, and H stays at its budget.

    These seeds reach a converged root with negative pH2O, whose pH2**2 forms
    CH4 and H2S.
    """
    roots, converged = _spy_roots(monkeypatch)
    np.random.seed(seed)
    r = equilibrium_atmosphere(
        dict(target), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **_COLD
    )
    assert any(ok and root[0] < 0.0 for root, ok in zip(roots, converged))
    assert r['CH4_bar'] < 1.0e-20
    assert r['H2S_bar'] < 1.0e-12
    assert r['NH3_bar'] < 1.0e-12
    for e in 'CS':
        assert abs(r[e + '_res']) <= max(target[e] * _COLD['rtol'], solve_mod.TRUNC_MASS), e
    assert r['H_kg_total'] == pytest.approx(target['H'], rel=1e-3, abs=solve_mod.TRUNC_MASS)


def test_unjudgeable_clipped_root_is_rejected(monkeypatch):
    """A root with a negative primary whose clipped residual is NaN is a
    failed attempt: with one attempt the solver raises."""
    raw = np.array([-1.0, 5.0, 1.0, 0.5])
    seen = []

    def _residual(x, *args):
        seen.append(np.array(x))
        return [0.0] * 4 if x[0] < 0.0 else [float('nan')] * 4

    monkeypatch.setattr(solve_mod.opt, 'fsolve', lambda *a, **k: (raw, {}, 1, 'stub'))
    monkeypatch.setattr(solve_mod, 'func', _residual)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere(
            dict(_LOW_H),
            _earth_ddict(T=1500.0, dIW=2.0),
            p_guess=None,
            **{**_COLD, 'nguess': 1},
        )
    np.testing.assert_array_equal(seen[-1][0], 0.0)
    assert seen[-1][1] == pytest.approx(5.0, rel=1e-12)


def test_clipped_root_must_close_each_noble_gas(monkeypatch):
    """The per-element check of a root with a negative primary covers the active
    noble gases: a clipped state that closes CHNOS but not He is rejected."""
    raw = np.array([1.0, 5.0, -3.0, 0.5, 2.0])  # N2 = -3 bar, He = 2 bar
    seen = []

    def _residual(x, *args):
        seen.append(np.array(x))
        return [0.0] * 5 if x[2] < 0.0 else [0.0, 0.0, 0.0, 0.0, 1.0e18]

    monkeypatch.setattr(solve_mod.opt, 'fsolve', lambda *a, **k: (raw, {}, 1, 'stub'))
    monkeypatch.setattr(solve_mod, 'func', _residual)
    ddict = dict(_earth_ddict(T=1500.0, dIW=2.0), He_included=1)
    target = dict(_LOW_H, He=1.0e17)
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere(
            target,
            ddict,
            p_guess={'H2O': 1.0, 'CO2': 5.0, 'N2': 1.0, 'S2': 0.5, 'He': 2.0},
            **{**_COLD, 'nguess': 1},
        )
    np.testing.assert_array_equal(seen[-1][2], 0.0)
    assert seen[-1][4] == pytest.approx(2.0, rel=1e-12)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('seed', 'p_surf'),
    [(0, 24.910292949351646), (1, 24.91031882820358), (2, 24.910318448637597)],
)
def test_trace_sulfur_with_absent_nitrogen_solves_as_before(monkeypatch, seed, p_surf):
    """An S budget below the scalar mass gate with no N leaves pN2 and pS2
    negative at the accepted root, both inert. The clip moves no mass there,
    so the solve returns the root it returns without the sign check (pinned
    P_surf) rather than holding trace S to a closure it never reaches."""
    roots, _ = _spy_roots(monkeypatch)
    np.random.seed(seed)
    target = {'H': 1.0e20, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e14}
    r = equilibrium_atmosphere(target, _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **_COLD)
    assert roots[-1][3] < 0.0
    assert r['P_surf'] == pytest.approx(p_surf, rel=1e-12)
    np.testing.assert_array_equal([r['N2_bar'], r['S2_bar']], 0.0)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('seed', 'p_surf', 'h2o'),
    [(0, 9.284216977332548, 0.40219572788976693), (7, 9.284231983250477, 0.40219618681983893)],
)
def test_all_positive_root_is_unchanged(monkeypatch, seed, p_surf, h2o):
    """At the Earth fiducial every accepted root has positive primaries, so
    the sign check does not run and the result is the root fsolve returned
    (pinned values)."""
    roots, _ = _spy_roots(monkeypatch)
    np.random.seed(seed)
    r = equilibrium_atmosphere(_earth_target_HCNS(), _earth_ddict(), p_guess=None, **_COLD)
    assert np.all(roots[-1] > 0.0)
    assert r['P_surf'] == pytest.approx(p_surf, rel=1e-12)
    assert r['H2O_bar'] == pytest.approx(h2o, rel=1e-12)


def _stub_buffered(monkeypatch, raw, residual):
    """Make fsolve return ``raw`` and replace the residual function by ``residual``."""
    seen = []

    def _residual(x, *args):
        seen.append(np.array(x))
        return residual(x)

    monkeypatch.setattr(solve_mod.opt, 'fsolve', lambda *a, **k: (np.array(raw), {}, 1, 'stub'))
    monkeypatch.setattr(solve_mod, 'func', _residual)
    return seen


@pytest.mark.parametrize('n2', [1.0e-30, -1.0e-30])
def test_verdict_does_not_depend_on_an_inert_sign(monkeypatch, n2):
    """A root whose N residual (5e15 kg) passes the scalar gate (1.5e16 kg) but
    not the N tolerance (1e13 kg) is accepted whatever the sign of an inert
    pN2, because clipping that pN2 changes no residual."""
    seen = _stub_buffered(monkeypatch, [1.0, 5.0, n2, 0.5], lambda x: [0.0, 0.0, 5.0e15, 0.0])
    r = equilibrium_atmosphere(
        dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 1}
    )
    assert r['H2O_bar'] == pytest.approx(1.0, rel=1e-12)
    assert r['N2_bar'] == pytest.approx(max(n2, 0.0), abs=1e-40)
    # The clipped state is evaluated only for the negative sign.
    assert len(seen) == (3 if n2 < 0.0 else 2)


@pytest.mark.parametrize(('factor', 'accepted'), [(0.5, True), (2.0, False)])
def test_clip_may_worsen_an_element_by_its_tolerance_only(monkeypatch, factor, accepted):
    """A clip that worsens the N residual by half the N tolerance
    (rtol * 1e17 kg) is accepted; by twice the tolerance it is rejected."""
    worse = factor * _LOW_H['N'] * _COLD['rtol']
    _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0, 0.0, 0.0 if x[0] < 0.0 else worse, 0.0],
    )
    args = (dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0))
    kw = dict(p_guess=None, **{**_COLD, 'nguess': 1})
    if accepted:
        r = equilibrium_atmosphere(*args, **kw)
        np.testing.assert_array_equal([r['H2O_bar'], r['H2_bar']], 0.0)
        assert r['CO2_bar'] == pytest.approx(5.0, rel=1e-12)
    else:
        with pytest.raises(RuntimeError, match='Could not find solution'):
            equilibrium_atmosphere(*args, **kw)
        assert worse > _LOW_H['N'] * _COLD['rtol']


def test_clip_that_closes_an_element_is_accepted(monkeypatch):
    """A clip that brings an element closer to its budget is accepted, as for
    a dry inventory where the negative pH2O's phantom CH4 and H2S are dropped:
    the raw N residual (2x its tolerance, inside the scalar gate) falls to 0."""
    tol_n = _LOW_H['N'] * _COLD['rtol']
    seen = _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0, 0.0, 2.0 * tol_n if x[0] < 0.0 else 0.0, 0.0],
    )
    r = equilibrium_atmosphere(
        dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 1}
    )
    np.testing.assert_array_equal([r['H2O_bar'], seen[-1][0]], 0.0)
    assert r['CO2_bar'] == pytest.approx(5.0, rel=1e-12)


def _authoritative_o(target, seed, dIW=2.0):
    """Authoritative-O solve at 1500 K with the O budget of a seeded buffered solve."""
    ddict = _earth_ddict(T=1500.0, dIW=dIW)
    np.random.seed(1)
    ref = equilibrium_atmosphere(dict(target), ddict, p_guess=None, **_COLD)
    tgt = dict(target, O=ref['O_kg_total'])
    return equilibrium_atmosphere_authoritative_O(
        tgt, ddict, fO2_hint=dIW, random_seed=seed, **_COLD
    ), tgt


@pytest.mark.physics_invariant
@pytest.mark.parametrize(('dIW', 'seed'), [(2.0, 42), (-2.0, 4)])
def test_authoritative_o_accepts_an_inert_negative_primary(monkeypatch, dIW, seed):
    """The authoritative-O path accepts a root whose pN2 is negative when N is
    absent, instead of restarting until fsolve lands on pN2 near zero, and
    recovers the buffered fO2 and the O budget. The fO2 offset is not a
    pressure and is never clipped, so the reducing case keeps it at -2.
    """
    roots, converged = _spy_roots(monkeypatch)
    r, tgt = _authoritative_o(_NO_N, seed=seed, dIW=dIW)
    assert roots[-1][2] < -1.0e-6
    assert r['fO2_shift_derived'] == pytest.approx(dIW, abs=1e-4)
    assert r['O_kg_total'] == pytest.approx(tgt['O'], rel=1e-6)
    np.testing.assert_array_equal(r['N2_bar'], 0.0)


@pytest.mark.physics_invariant
def test_authoritative_o_low_h_returns_the_physical_root():
    """The low-H inventory on the authoritative-O path returns positive pH2O
    and the surface pressure of the physical buffered root. This guards the
    sign check against rejecting the physical root; the authoritative-O solve
    does not reach the negative-water root for this inventory.
    """
    r, tgt = _authoritative_o(_LOW_H, seed=42)
    assert r['H2O_bar'] > 0.0
    assert r['P_surf'] == pytest.approx(24.8219, rel=1e-4)
    assert r['fO2_shift_derived'] == pytest.approx(2.0, abs=1e-4)
