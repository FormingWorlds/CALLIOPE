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

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(360)]


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
    with the default Hirschmann 2021 IW buffer.
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
def test_ddict_buffer_reaches_every_fo2_channel():
    """With ddict['fO2_model'] = 'oneill', the O2 partial pressure, the H2,
    CO, CH4 and SO2 modified equilibrium constants and the S2 solubility all
    reference the O'Neill & Eggins (2002) IW buffer. At 2500 K O'Neill sits
    0.713 dex below the default Hirschmann (2021) buffer, so a channel left
    on the default fails by far more than the tolerances. The Dasgupta N2
    law takes the shift itself as its Delta-IW, so it is not converted to
    another buffer. An unknown buffer name, or None, in ddict raises
    ValueError instead of falling back to the default.
    """
    from calliope.chemistry import ModifiedKeq
    from calliope.oxygen_fugacity import OxygenFugacity
    from calliope.solubility import SolubilityN2, SolubilityS2

    T, dIW = 2500.0, 2.0
    ddict = dict(_earth_ddict(T=T, dIW=dIW), fO2_model='oneill')
    pin = {'H2O': 10.0, 'CO2': 50.0, 'N2': 1.0, 'S2': 0.1}
    p = solve_mod.get_partial_pressures(pin, ddict)

    assert np.log10(p['O2']) == pytest.approx(OxygenFugacity('oneill')(T, dIW), abs=1e-12)
    assert p['H2'] / p['H2O'] == pytest.approx(
        ModifiedKeq('janaf_H2', 'oneill')(T, dIW), rel=1e-12
    )
    assert p['CO'] / p['CO2'] == pytest.approx(
        ModifiedKeq('janaf_CO', 'oneill')(T, dIW), rel=1e-12
    )
    # CH4 carries one O2 in its modified constant, so the default buffer
    # would move this ratio by the full 0.713 dex.
    assert p['CH4'] / (p['CO2'] * p['H2'] ** 2) == pytest.approx(
        ModifiedKeq('schaefer_CH4', 'oneill')(T, dIW), rel=1e-12
    )
    assert p['SO2'] ** 2 / (p['S2'] * p['O2'] ** 2) == pytest.approx(
        ModifiedKeq('janaf_SO2', 'oneill')(T, dIW), rel=1e-10
    )
    prefactor = 1e-6 * ddict['M_mantle'] * ddict['Phi_global']
    assert solve_mod.dissolved_mass(pin, ddict)['S2'] == pytest.approx(
        prefactor * SolubilityS2(fO2_model='oneill')(pin['S2'], T, dIW), rel=1e-12
    )
    p_tot = solve_mod.get_total_pressure(p)
    n2_raw = SolubilityN2('dasgupta')(pin['N2'], p_tot, T, dIW)
    n2_mass = solve_mod.dissolved_mass(pin, ddict)['N2']
    assert n2_mass == pytest.approx(prefactor * n2_raw, rel=1e-12)
    # Wrong-law guard: converting the shift to Hirschmann's Delta-IW (dIW - 0.713)
    # would grow the reduced-N term by exp(1.6 * 0.713) = 3.13.
    offset = OxygenFugacity('hirschmann')(T) - OxygenFugacity('oneill')(T)
    n2_converted = SolubilityN2('dasgupta')(pin['N2'], p_tot, T, dIW - offset)
    assert n2_converted / n2_raw > 1.01

    # Wrong-buffer guard: without the key the solve uses Hirschmann, whose
    # O2 pressure at the same shift is 0.713 dex higher.
    p_default = solve_mod.get_partial_pressures(pin, _earth_ddict(T=T, dIW=dIW))
    assert np.log10(p_default['O2'] / p['O2']) == pytest.approx(0.713, abs=2e-3)

    for bad in ('typo', None):
        with pytest.raises(ValueError, match='Unknown IW buffer'):
            solve_mod.get_partial_pressures(pin, dict(ddict, fO2_model=bad))
        with pytest.raises(ValueError, match='Unknown IW buffer'):
            solve_mod.dissolved_mass(pin, dict(ddict, fO2_model=bad))


@pytest.mark.physics_invariant
def test_fo2_model_argument_sets_the_buffer_of_the_whole_solve():
    """equilibrium_atmosphere(fO2_model='oneill') solves on the O'Neill
    buffer and reaches the same state as ddict['fO2_model'] = 'oneill',
    without writing the key into the caller's ddict. Its O2 pressure sits
    IW_Hirschmann - IW_ONeill (0.160 dex at 1800 K) below the default solve
    at the same shift. A buffer that contradicts ddict, or an unknown name,
    is refused before any solving.
    """
    from calliope.oxygen_fugacity import OxygenFugacity

    T, dIW = 1800.0, 2.0
    ddict = _earth_ddict(T=T, dIW=dIW)
    target = _earth_target_HCNS()
    kw = dict(print_result=False, nguess=200)

    np.random.seed(0)
    by_arg = equilibrium_atmosphere(target, ddict, fO2_model='oneill', **kw)
    assert 'fO2_model' not in ddict
    np.random.seed(0)
    by_key = equilibrium_atmosphere(target, dict(ddict, fO2_model='oneill'), **kw)
    np.testing.assert_array_equal(
        [by_arg[s + '_bar'] for s in volatile_species],
        [by_key[s + '_bar'] for s in volatile_species],
    )
    assert np.log10(by_arg['O2_bar']) == pytest.approx(
        OxygenFugacity('oneill')(T, dIW), abs=1e-10
    )
    for e in 'HCNS':
        assert abs(by_arg[e + '_res']) <= 1e-5 * target[e] + 1e10, e

    # Wrong-buffer guard against the default solve at the same shift.
    np.random.seed(0)
    default = equilibrium_atmosphere(target, ddict, **kw)
    offset = OxygenFugacity()(T) - OxygenFugacity('oneill')(T)
    assert offset == pytest.approx(0.160, abs=2e-3)
    assert np.log10(default['O2_bar'] / by_arg['O2_bar']) == pytest.approx(offset, abs=1e-10)

    with pytest.raises(ValueError, match='contradicts'):
        equilibrium_atmosphere(
            target, dict(ddict, fO2_model='fischer'), fO2_model='oneill', **kw
        )
    with pytest.raises(ValueError, match='Unknown IW buffer'):
        equilibrium_atmosphere(target, ddict, fO2_model='fisher', **kw)
    with pytest.raises(ValueError, match='Unknown IW buffer'):
        equilibrium_atmosphere_authoritative_O(
            dict(target, O=1.0e20), ddict, fO2_model='fisher', **kw
        )


@pytest.mark.physics_invariant
def test_authoritative_o_reports_the_shift_against_the_chosen_buffer():
    """The authoritative-O solve recovers an absolute fO2; fO2_model sets only
    the buffer the shift is reported against. An O budget made on O'Neill at
    IW+2 round-trips to +2 on O'Neill, and the same budget reported against
    the default Hirschmann buffer gives 2 - (IW_Hirschmann - IW_ONeill),
    1.840 at 1800 K, with the same physical state.
    """
    from calliope.oxygen_fugacity import OxygenFugacity

    T, dIW = 1800.0, 2.0
    ddict = _earth_ddict(T=T, dIW=dIW)
    target = _earth_target_HCNS()
    np.random.seed(0)
    fwd = equilibrium_atmosphere(
        target, ddict, fO2_model='oneill', print_result=False, nguess=200
    )
    tgt = dict(target, O=fwd['O_kg_total'])
    kw = dict(
        fO2_hint=dIW,
        p_guess={s: fwd[s + '_bar'] for s in ('H2O', 'CO2', 'N2', 'S2')},
        random_seed=0,
        nguess=500,
        nsolve=1000,
        print_result=False,
    )
    on_oneill = equilibrium_atmosphere_authoritative_O(tgt, ddict, fO2_model='oneill', **kw)
    on_default = equilibrium_atmosphere_authoritative_O(tgt, ddict, **kw)

    offset = OxygenFugacity()(T) - OxygenFugacity('oneill')(T)
    assert on_oneill['fO2_shift_derived'] == pytest.approx(dIW, abs=1e-4)
    assert on_default['fO2_shift_derived'] == pytest.approx(dIW - offset, abs=1e-4)
    # Guard: the offset is far above the 1e-4 tolerance, so reporting
    # against the wrong buffer cannot pass either pin.
    assert offset == pytest.approx(0.160, abs=2e-3)
    assert on_default['O2_bar'] == pytest.approx(on_oneill['O2_bar'], rel=1e-3)
    assert on_default['P_surf'] == pytest.approx(on_oneill['P_surf'], rel=1e-4)


@pytest.mark.physics_invariant
def test_authoritative_o_dissolved_n_follows_the_shift_on_each_buffer():
    """At 3000 K and IW-2 the Hirschmann (2021) and Fischer (2011) buffers sit
    0.954 and 1.063 dex above O'Neill & Eggins (2002). The same O budget then
    gives nearly the same absolute fO2 (within 0.03 dex) on all three buffers,
    but a shift lower by about the offset. The Dasgupta N2 law takes that shift
    as its Delta-IW, so its reduced-N term grows and the melt holds 3.6x
    (Hirschmann) and 4.1x (Fischer) the N it holds on O'Neill: 2.50e18 and
    2.85e18 kg against 7.02e17 kg of the 8e18 kg budget, with P_surf 2.9 and
    3.5 % lower. Every element still closes on every buffer.
    """
    from calliope.oxygen_fugacity import OxygenFugacity

    T, dIW = 3000.0, -2.0
    ddict = _earth_ddict(T=T, dIW=dIW)
    target = _earth_target_HCNS()
    np.random.seed(0)
    fwd = equilibrium_atmosphere(
        target, ddict, fO2_model='oneill', print_result=False, nguess=200
    )
    tgt = dict(target, O=fwd['O_kg_total'])
    kw = dict(
        fO2_hint=dIW,
        p_guess={s: fwd[s + '_bar'] for s in ('H2O', 'CO2', 'N2', 'S2')},
        random_seed=0,
        nguess=500,
        nsolve=1000,
        print_result=False,
    )
    runs = {
        m: equilibrium_atmosphere_authoritative_O(tgt, ddict, fO2_model=m, **kw)
        for m in ('oneill', 'hirschmann', 'fischer')
    }

    ref = runs['oneill']
    assert ref['fO2_shift_derived'] == pytest.approx(dIW, abs=1e-4)
    assert ref['N2_kg_liquid'] == pytest.approx(7.02e17, rel=1e-2)
    expected = {'hirschmann': (2.50e18, 0.971), 'fischer': (2.85e18, 0.965)}
    for m, r in runs.items():
        for e in 'HCNS':
            assert abs(r[e + '_res']) <= 1e-5 * target[e] + 1e10, (m, e)
        assert np.log10(r['O2_bar'] / ref['O2_bar']) == pytest.approx(0.0, abs=0.03), m
        if m == 'oneill':
            continue
        offset = OxygenFugacity(m)(T) - OxygenFugacity('oneill')(T)
        # The shift falls by nearly the offset; N2 dissolution takes up the rest.
        assert r['fO2_shift_derived'] == pytest.approx(dIW - offset, abs=0.05), m
        n_liq, p_ratio = expected[m]
        assert r['N2_kg_liquid'] == pytest.approx(n_liq, rel=1e-2), m
        assert r['P_surf'] / ref['P_surf'] == pytest.approx(p_ratio, rel=1e-3), m
        # Guard: a buffer-converted N2 law would keep the O'Neill N within 1e-3.
        assert r['N2_kg_liquid'] / ref['N2_kg_liquid'] > 3.0, m
    assert OxygenFugacity()(T) - OxygenFugacity('oneill')(T) == pytest.approx(0.954, abs=2e-3)


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


@pytest.fixture(autouse=True)
def _fresh_empty_budget_warnings(monkeypatch):
    """Each test sees the fallback warning as if in a new process."""
    monkeypatch.setattr(solve_mod, '_FALLBACK_WARNED', set())


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


def _stub_buffered(monkeypatch, raw, residual):
    """Make fsolve return ``raw`` and replace the residual function by ``residual``."""
    seen = []

    def _residual(x, *args):
        seen.append(np.array(x))
        return residual(x)

    monkeypatch.setattr(solve_mod.opt, 'fsolve', lambda *a, **k: (np.array(raw), {}, 1, 'stub'))
    monkeypatch.setattr(solve_mod, 'func', _residual)
    return seen


def _one_attempt():
    """Single-attempt solve of the low-H inventory at 1500 K and IW+2."""
    return equilibrium_atmosphere(
        dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 1}
    )


@pytest.mark.physics_invariant
def test_low_h_cold_start_rejects_the_negative_water_root(monkeypatch, caplog):
    """A low-H, C-rich cold start at IW+2 and 1500 K returns the physical root.

    Some cold starts reach the root with negative pH2O, which closes the mass
    balance with no H2O or H2 in the atmosphere. It is rejected, and the solver
    goes on to the physical root at 25.05 bar with the Hirschmann (2021) IW
    buffer (24.82 bar with the Fischer 2011 buffer). Which seeds reach
    the negative root depends on the platform's floating point (seed 6 does on
    macOS arm64, not on Linux x86_64), so every seed must return the physical
    root and at least one of them must have rejected the negative one.
    """
    roots, converged = _spy_roots(monkeypatch)
    rejected = 0
    for seed in (4, 6, 7):
        roots.clear()
        converged.clear()
        np.random.seed(seed)
        r = equilibrium_atmosphere(
            dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **_COLD
        )
        rejected += any(ok and root[0] < -1.0e3 for root, ok in zip(roots[:-1], converged[:-1]))
        assert roots[-1][0] > 0.0, seed
        assert r['H2O_bar'] > 0.0, seed
        assert r['P_surf'] == pytest.approx(25.0485, rel=1e-4), seed
        assert r['H_kg_total'] == pytest.approx(_LOW_H['H'], rel=1e-3), seed
    assert 'returning the last root' not in caplog.text
    assert rejected >= 1


@pytest.mark.physics_invariant
def test_absent_element_root_is_accepted_and_reported_clipped(monkeypatch):
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
        ({'H': 1.0e15, 'C': 5.0e19, 'N': 1.0e17, 'S': 1.0e18}, 1),
    ],
    ids=['no_H_no_N-0', 'no_H_no_N-2', 'trace_H-1'],
)
def test_sub_gate_hydrogen_reports_no_phantom_species(monkeypatch, target, seed):
    """An H budget below the scalar mass gate (1.5e16 kg here) yields a state
    free of H-bearing species from a negative pH2O, C and S close on their own
    tolerance rather than on the scalar gate, and H is kept. The scalar gate does
    not bound a 1e15 kg H budget, so H is pinned to 10 %: dropping the phantom
    species of a negative pH2O would remove it whole.

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
    assert r['H_kg_total'] == pytest.approx(target['H'], rel=0.1, abs=solve_mod.TRUNC_MASS)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(('fO2_model', 'seed'), [('hirschmann', 0), ('oneill', 3)])
def test_root_that_empties_hydrogen_is_refused_and_h_is_kept(
    monkeypatch, caplog, fO2_model, seed
):
    """A 1e15 kg H budget lies below the scalar gate (1.5e16 kg), so the gate does
    not see a root that loses it. These seeds converge on a root with pH2O near
    -74 bar (Hirschmann) or -83 bar (O'Neill), which loses the whole H budget once
    clipped. It is refused, and the solve goes on to a root that keeps H, with no
    fallback warning and no phantom CH4, H2S or NH3. C and S close on their own
    tolerance. Which other roots a seed reaches depends on the platform's floating
    point (on macOS arm64 both seeds also reach a root that has lost its H before
    any clip, on Linux x86_64 only the O'Neill one does), so the test asserts the
    outcome; the refusal of an unclipped emptied root is pinned by
    test_emptied_root_is_refused_for_a_later_root_that_keeps_the_budget.
    """
    target = {'H': 1.0e15, 'C': 5.0e19, 'N': 1.0e17, 'S': 1.0e18}
    roots, converged = _spy_roots(monkeypatch)
    np.random.seed(seed)
    r = equilibrium_atmosphere(
        dict(target),
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess=None,
        fO2_model=fO2_model,
        **_COLD,
    )
    refused = [root for root, ok in zip(roots[:-1], converged[:-1]) if ok]
    assert any(root[0] < -1.0 for root in refused)
    assert roots[-1][0] > 0.0
    assert 'returning the last root' not in caplog.text
    # The scalar gate does not bound H here, so H is pinned to 10 %, as in
    # test_sub_gate_hydrogen_reports_no_phantom_species; an emptied state gives 0 kg.
    assert r['H_kg_total'] == pytest.approx(target['H'], rel=0.1)
    assert r['CH4_bar'] < 1.0e-20
    assert r['H2S_bar'] < 1.0e-12
    assert r['NH3_bar'] < 1.0e-12
    for e in 'CS':
        assert abs(r[e + '_res']) <= max(target[e] * _COLD['rtol'], solve_mod.TRUNC_MASS), e


def test_unjudgeable_clipped_root_is_rejected(monkeypatch):
    """A root with a negative primary whose clipped residual is NaN is a
    failed attempt: with one attempt the solver raises."""
    seen = _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0] * 4 if x[0] < 0.0 else [float('nan')] * 4,
    )
    with pytest.raises(RuntimeError, match='Could not find solution'):
        _one_attempt()
    np.testing.assert_array_equal(seen[-1][0], 0.0)
    assert seen[-1][1] == pytest.approx(5.0, rel=1e-12)


def test_clipped_root_must_close_each_noble_gas(monkeypatch, caplog):
    """The check of a root with a negative primary covers the active noble
    gases: a clip that worsens He beyond its tolerance is rejected, and since
    the clipped state also fails the He gate it is no fallback either."""
    seen = _stub_buffered(
        monkeypatch,
        [1.0, 5.0, -3.0, 0.5, 2.0],  # N2 = -3 bar, He = 2 bar
        lambda x: [0.0] * 5 if x[2] < 0.0 else [0.0, 0.0, 0.0, 0.0, 1.0e18],
    )
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere(
            dict(_LOW_H, He=1.0e17),
            dict(_earth_ddict(T=1500.0, dIW=2.0), He_included=1),
            p_guess={'H2O': 1.0, 'CO2': 5.0, 'N2': 1.0, 'S2': 0.5, 'He': 2.0},
            **{**_COLD, 'nguess': 1},
        )
    assert 'returning the last root' not in caplog.text
    np.testing.assert_array_equal(seen[-1][2], 0.0)
    assert seen[-1][4] == pytest.approx(2.0, rel=1e-12)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('seed', 'p_surf'),
    [(0, 25.13491096061388), (1, 25.134909250760263), (2, 25.134908903155235)],
)
def test_trace_sulfur_with_absent_nitrogen_returns_through_the_fallback(caplog, seed, p_surf):
    """An S budget below the scalar mass gate with no N: no attempt in 1000
    reaches a root that keeps S, so every root that passes the gate is refused
    and the last one returns through the fallback, whose warning names S. pN2
    and pS2 are reported clipped to zero (P_surf to rel 1e-5, the xtol noise)
    and H, C close."""
    np.random.seed(seed)
    target = {'H': 1.0e20, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e14}
    r = equilibrium_atmosphere(target, _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **_COLD)
    assert 'after 1000 attempts' in caplog.text
    assert 'tolerance for S (1e+14 kg, tolerance 1e+10 kg, budget 1e+14 kg)' in caplog.text
    np.testing.assert_array_equal([r['N2_bar'], r['S2_bar'], r['S_kg_total']], 0.0)
    assert r['P_surf'] == pytest.approx(p_surf, rel=1e-5)
    for e in 'HC':
        assert abs(r[e + '_res']) <= target[e] * _COLD['rtol'], e


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('seed', 'p_surf', 'h2o'),
    [(0, 9.397713354566907, 0.4022191760949055), (7, 9.397713513933748, 0.4022193071737796)],
)
def test_all_positive_root_is_unchanged(monkeypatch, caplog, seed, p_surf, h2o):
    """At the Earth fiducial every accepted root has positive primaries, so
    the clip changes nothing and the result is the root fsolve returned. The
    pins hold to rel 1e-5, the solver's xtol noise across platforms."""
    roots, _ = _spy_roots(monkeypatch)
    np.random.seed(seed)
    r = equilibrium_atmosphere(_earth_target_HCNS(), _earth_ddict(), p_guess=None, **_COLD)
    assert np.all(roots[-1] > 0.0)
    np.testing.assert_array_equal(
        [r[s + '_bar'] for s in ('H2O', 'CO2', 'N2', 'S2')], roots[-1]
    )
    assert r['P_surf'] == pytest.approx(p_surf, rel=1e-5)
    assert r['H2O_bar'] == pytest.approx(h2o, rel=1e-5)
    assert 'No root closes' not in caplog.text


@pytest.mark.parametrize('n2', [1.0e-30, -1.0e-30])
def test_verdict_does_not_depend_on_an_inert_sign(monkeypatch, caplog, n2):
    """A root whose N residual (5e15 kg) passes the scalar gate (1.5e16 kg) but
    not the N tolerance (1e13 kg) is accepted whatever the sign of an inert
    pN2, because clipping that pN2 changes no residual."""
    seen = _stub_buffered(monkeypatch, [1.0, 5.0, n2, 0.5], lambda x: [0.0, 0.0, 5.0e15, 0.0])
    r = _one_attempt()
    assert r['H2O_bar'] == pytest.approx(1.0, rel=1e-12)
    np.testing.assert_array_equal(r['N2_bar'], max(n2, 0.0))
    assert 'returning the last root' not in caplog.text
    # The clipped state is evaluated only for the negative sign.
    assert len(seen) == (3 if n2 < 0.0 else 2)


@pytest.mark.parametrize(('factor', 'accepted'), [(0.5, True), (2.0, False)])
def test_clip_may_worsen_an_element_by_its_tolerance_only(
    monkeypatch, caplog, factor, accepted
):
    """A clip that worsens the N residual by half the N tolerance
    (rtol * 1e17 kg) is accepted; by twice the tolerance it is rejected, and
    with no attempt left it is returned through the warned fallback."""
    worse = factor * _LOW_H['N'] * _COLD['rtol']
    _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0, 0.0, 0.0 if x[0] < 0.0 else worse, 0.0],
    )
    r = _one_attempt()
    np.testing.assert_array_equal([r['H2O_bar'], r['H2_bar']], 0.0)
    assert r['CO2_bar'] == pytest.approx(5.0, rel=1e-12)
    assert ('tolerance for N (' in caplog.text) == (not accepted)


def test_clip_that_closes_an_element_is_accepted(monkeypatch, caplog):
    """A clip that brings an element closer to its budget is accepted, as for
    a dry inventory where the negative pH2O's phantom CH4 and H2S are dropped:
    the raw N residual (2x its tolerance, inside the scalar gate) falls to 0."""
    tol_n = _LOW_H['N'] * _COLD['rtol']
    seen = _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0, 0.0, 2.0 * tol_n if x[0] < 0.0 else 0.0, 0.0],
    )
    r = _one_attempt()
    np.testing.assert_array_equal([r['H2O_bar'], seen[-1][0]], 0.0)
    assert 'returning the last root' not in caplog.text
    assert r['CO2_bar'] == pytest.approx(5.0, rel=1e-12)


@pytest.mark.parametrize('nguess', [2, 5])
def test_fallback_never_returns_a_state_outside_the_mass_gate(caplog, nguess):
    """With few attempts the low-H cold start can end on the negative-water
    root only. Clipped, that root loses its H and most of its C and S, so it
    fails the scalar gate and is no fallback: the solve raises or finds the
    physical root, never the 6.9 bar state."""
    outcomes = []
    for seed in (3, 6, 7):
        np.random.seed(seed)
        try:
            r = equilibrium_atmosphere(
                dict(_LOW_H),
                _earth_ddict(T=1500.0, dIW=2.0),
                p_guess=None,
                **{**_COLD, 'nguess': nguess},
            )
            outcomes.append(r['P_surf'])
        except RuntimeError:
            outcomes.append(None)
    assert all(p is None or p == pytest.approx(25.0485, rel=1e-4) for p in outcomes)
    assert 'returning the last root' not in caplog.text


@pytest.mark.parametrize(('second_in_gate', 'co2'), [(True, 7.0), (False, 5.0)])
def test_fallback_is_the_last_gate_passing_root(monkeypatch, caplog, second_in_gate, co2):
    """Of two attempts rejected by the sign check, the solver falls back to the
    later one whose clipped state passes the mass gate: the second if it does,
    else the first."""
    roots = iter([np.array([-1.0, 5.0, 1.0, 0.5]), np.array([-1.0, 7.0, 1.0, 0.5])])
    monkeypatch.setattr(solve_mod.opt, 'fsolve', lambda *a, **k: (next(roots), {}, 1, 'stub'))
    worse = 2.0 * _LOW_H['N'] * _COLD['rtol']
    beyond = 1.0e17 if not second_in_gate else worse

    def _residual(x, *a):
        if x[0] < 0.0:
            return [0.0, 0.0, 0.0, 0.0]
        return [0.0, 0.0, worse if x[1] < 6.0 else beyond, 0.0]

    monkeypatch.setattr(solve_mod, 'func', _residual)
    r = equilibrium_atmosphere(
        dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 2}
    )
    assert 'after 2 attempts' in caplog.text
    assert r['CO2_bar'] == pytest.approx(co2, rel=1e-12)


def test_clipped_state_outside_the_mass_gate_is_no_fallback(monkeypatch, caplog):
    """A raw root inside the scalar gate whose clipped state falls outside it
    (2e16 kg against 1.5e16 kg) is rejected and kept as no fallback: with no
    root that passed the gate, the solve raises without the fallback warning."""
    _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0 if x[0] < 0.0 else 2.0e16, 0.0, 0.0, 0.0],
    )
    with pytest.raises(RuntimeError, match='Could not find solution'):
        _one_attempt()
    assert 'returning the last root' not in caplog.text
    assert 2.0e16 > max(_LOW_H.values()) * _COLD['rtol'] + _COLD['atol']


def test_accepted_clipped_state_must_pass_the_mass_gate(monkeypatch):
    """A clip that worsens N by less than its tolerance (1e12 kg against 1e13 kg)
    but moves the clipped state just outside the scalar gate (1.5e16 kg) is not
    accepted: the clipped state passes the same gate as any accepted root."""
    gate = max(_LOW_H.values()) * _COLD['rtol'] + _COLD['atol'] + solve_mod.TRUNC_MASS
    raw_n, clip_n = gate - 5.0e11, gate + 5.0e11
    _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0, 0.0, raw_n if x[0] < 0.0 else clip_n, 0.0],
    )
    with pytest.raises(RuntimeError, match='Could not find solution'):
        _one_attempt()
    assert clip_n - raw_n < _LOW_H['N'] * _COLD['rtol']


def test_accepted_clipped_state_must_pass_the_noble_gate(monkeypatch, caplog):
    """A clip that worsens He by less than its tolerance (7e12 kg against 1e13 kg)
    but leaves the clipped He residual above that tolerance is not accepted,
    and the clipped state is no fallback either."""
    _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5, 2.0],
        lambda x: [0.0, 0.0, 0.0, 0.0, 6.0e12 if x[0] < 0.0 else 1.3e13],
    )
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere(
            dict(_LOW_H, He=1.0e17),
            dict(_earth_ddict(T=1500.0, dIW=2.0), He_included=1),
            p_guess={'H2O': 1.0, 'CO2': 5.0, 'N2': 1.0, 'S2': 0.5, 'He': 2.0},
            **{**_COLD, 'nguess': 1},
        )
    assert 'returning the last root' not in caplog.text
    assert 1.3e13 - 6.0e12 < 1.0e17 * _COLD['rtol']


def test_sub_gate_budget_worsened_by_the_clip_is_rejected(monkeypatch, caplog):
    """A clip that worsens an element beyond its tolerance is rejected even when
    that element's budget lies below the scalar gate; with no attempt left the
    clipped root, which passes the gate, returns through the warned fallback."""
    target = dict(_LOW_H, H=1.0e15)
    _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0 if x[0] < 0.0 else 1.0e15, 0.0, 0.0, 0.0],
    )
    r = equilibrium_atmosphere(
        target, _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 1}
    )
    np.testing.assert_array_equal([r['H2O_bar'], r['H2_bar'], r['CH4_bar']], 0.0)
    assert 'tolerance for H (1e+15 kg' in caplog.text
    assert 1.0e15 < max(target.values()) * _COLD['rtol'] + _COLD['atol']


def test_trace_hydrogen_without_opt_solver_returns_through_the_fallback(caplog):
    """With the library tolerances and nguess and fsolve only, no attempt keeps a
    trace H budget (1e12 kg): every root that passes the mass gate loses it. The
    last of them returns through the fallback, whose warning names H with its
    residual, tolerance and budget, and S, which that root misses by 1.59e13 kg
    against a 1e13 kg tolerance."""
    np.random.seed(1)
    r = equilibrium_atmosphere(
        {'H': 1.0e12, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18},
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess=None,
        xtol=1e-8,
        rtol=1e-5,
        atol=1e10,
        nguess=7500,
        nsolve=1500,
        print_result=False,
        opt_solver=False,
    )
    assert 'after 7500 attempts' in caplog.text
    assert 'tolerance for H (1e+12 kg, tolerance 1e+07 kg, budget 1e+12 kg)' in caplog.text
    assert 'S (1.59e+13 kg, tolerance 1e+13 kg' in caplog.text
    assert r['H_kg_total'] == pytest.approx(0.0, abs=1.0e12 * 1e-5)


_SUB_GATE_N = dict(_LOW_H, N=1.0e15)


_PRIMARIES, _STUB_ROOT = ('H2O', 'CO2', 'N2', 'S2'), [1.0, 5.0, 1.0, 0.5]


def _empty_budget_solve(monkeypatch, residual, target=_SUB_GATE_N):
    """One-attempt solve of an all-positive stub root with the given residual."""
    _stub_buffered(monkeypatch, _STUB_ROOT, lambda x: residual)
    return equilibrium_atmosphere(
        dict(target), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 1}
    )


@pytest.mark.parametrize(
    ('n_res', 'refused'),
    [(-1.0e15, True), (-0.5e15, False), (1.0e15, False)],
    ids=['emptied', 'half_lost', 'overfilled'],
)
def test_refusal_needs_the_whole_budget_lost(monkeypatch, caplog, n_res, refused):
    """A state whose N residual is minus the N budget (1e15 kg, inside the scalar
    gate) is refused, and with no attempt left it returns through the fallback,
    whose warning names N; a state that loses half of it or overfills it is
    accepted."""
    r = _empty_budget_solve(monkeypatch, [0.0, 0.0, n_res, 0.0])
    assert ('tolerance for N (1e+15 kg' in caplog.text) == refused
    np.testing.assert_array_equal([r[g + '_bar'] for g in _PRIMARIES], _STUB_ROOT)


def test_emptied_root_is_refused_for_a_later_root_that_keeps_the_budget(monkeypatch, caplog):
    """An all-positive root that empties a 1e15 kg H budget passes the scalar
    gate (1.5e16 kg) and has nothing to clip, yet is refused; the solve goes on
    to the next attempt, whose root closes every element, and returns it with
    no fallback warning."""
    roots = iter([np.array([1.0, 5.0, 1.0, 0.5]), np.array([2.0, 5.0, 1.0, 0.5])])
    monkeypatch.setattr(solve_mod.opt, 'fsolve', lambda *a, **k: (next(roots), {}, 1, 'stub'))
    target = dict(_LOW_H, H=1.0e15)
    monkeypatch.setattr(
        solve_mod, 'func', lambda x, *a: [-target['H'] if x[0] < 1.5 else 0.0, 0.0, 0.0, 0.0]
    )
    r = equilibrium_atmosphere(
        target, _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 2}
    )
    assert r['H2O_bar'] == pytest.approx(2.0, rel=1e-12)
    assert 'No root closes' not in caplog.text
    assert target['H'] < max(target.values()) * _COLD['rtol'] + _COLD['atol']


def test_budget_within_the_mass_floor_may_be_emptied(monkeypatch, caplog):
    """A 5 kg N budget lost whole lies within the 10 kg floor of the N tolerance,
    so the state is accepted rather than refused."""
    r = _empty_budget_solve(monkeypatch, [0.0, 0.0, -5.0, 0.0], target=dict(_LOW_H, N=5.0))
    assert 'No root closes' not in caplog.text
    np.testing.assert_array_equal([r[g + '_bar'] for g in _PRIMARIES], _STUB_ROOT)
    assert 5.0 < solve_mod.TRUNC_MASS


def test_zero_budget_is_never_empty(monkeypatch, caplog):
    """An element with no budget and no residual is not refused as emptied."""
    r = _empty_budget_solve(monkeypatch, [0.0, 0.0, 0.0, 0.0], target=dict(_LOW_H, N=0.0))
    assert 'No root closes' not in caplog.text
    np.testing.assert_array_equal([r[g + '_bar'] for g in _PRIMARIES], _STUB_ROOT)


def test_fallback_warns_once_per_element(monkeypatch, caplog):
    """A second solve that falls back on the same emptied budget logs at debug,
    not warning."""
    caplog.set_level(logging.DEBUG, logger='fwl.calliope.solve')
    for _ in range(2):
        _empty_budget_solve(monkeypatch, [0.0, 0.0, -1.0e15, 0.0])
    records = [r for r in caplog.records if 'No root closes' in r.getMessage()]
    assert [r.levelno for r in records] == [logging.WARNING, logging.DEBUG]
    assert all(
        'N (1e+15 kg, tolerance 1e+11 kg, budget 1e+15 kg)' in r.getMessage() for r in records
    )


def test_fallback_that_empties_a_budget_warns_once(monkeypatch, caplog):
    """A fallback state that empties a budget gives the fallback warning only."""
    _stub_buffered(
        monkeypatch,
        [-1.0, 5.0, 1.0, 0.5],
        lambda x: [0.0, 0.0, -1.0e15 if x[0] >= 0.0 else 0.0, 0.0],
    )
    equilibrium_atmosphere(
        dict(_SUB_GATE_N),
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess=None,
        **{**_COLD, 'nguess': 1},
    )
    warned = [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert len(warned) == 1
    assert 'tolerance for N (1e+15 kg' in warned[0].getMessage()


_LIBRARY = dict(xtol=1e-8, rtol=1e-5, atol=1e10, nguess=50, nsolve=3000, print_result=False)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('target', 'seed'),
    [
        ({'H': 1.0e12, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}, 0),
        ({'H': 1.0e12, 'C': 5.0e19, 'N': 1.0e14, 'S': 1.0e18}, 0),
        ({'H': 1.0e13, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}, 1),
    ],
    ids=['H1e12_N0', 'H1e12_N1e14', 'H1e13_N0'],
)
def test_trace_hydrogen_falls_back_to_the_clipped_root(caplog, target, seed):
    """At library tolerances with few attempts, a trace H budget is closed only
    by roots whose negative pH2O carries the H in CH4 and H2S. Clipping such a
    root empties the H budget, far beyond the H tolerance, so it is rejected
    (the phantom H2S also takes S beyond its tolerance); the solve returns the
    last one clipped, with a warning listing H and S, and no H-bearing species."""
    np.random.seed(seed)
    r = equilibrium_atmosphere(
        dict(target),
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess=None,
        opt_solver=False,
        **_LIBRARY,
    )
    assert 'tolerance for H (' in caplog.text
    assert ', S (' in caplog.text
    np.testing.assert_array_equal(
        [r[s + '_bar'] for s in ('H2O', 'H2', 'CH4', 'H2S', 'NH3')], 0.0
    )
    assert r['P_surf'] == pytest.approx(25.0274, rel=1e-4)
    # The returned clipped state still passes the scalar mass gate.
    gate = max(target.values()) * _LIBRARY['rtol'] + _LIBRARY['atol'] + solve_mod.TRUNC_MASS
    assert max(abs(r[e + '_res']) for e in 'HCNS') <= gate


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
    absent, instead of restarting until fsolve lands on a non-negative pN2, and
    recovers the buffered fO2 and the O budget. The fO2 offset is not a
    pressure and is never clipped, so the reducing case keeps it at -2.
    """
    roots, converged = _spy_roots(monkeypatch)
    r, tgt = _authoritative_o(_NO_N, seed=seed, dIW=dIW)
    assert roots[-1][2] < 0.0
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
    assert r['P_surf'] == pytest.approx(25.0485, rel=1e-4)
    assert r['fO2_shift_derived'] == pytest.approx(2.0, abs=1e-4)


@pytest.mark.parametrize('T', [1800, np.int64(1800)], ids=['python_int', 'numpy_int64'])
def test_integer_T_magma_solves_like_the_float_value(T):
    """An integer T_magma reaches every equilibrium constant and gives the same
    seeded Earth atmosphere, H2S and NH3 included, as T_magma = 1800.0."""
    np.random.seed(0)
    ref = equilibrium_atmosphere(_earth_target_HCNS(), _earth_ddict(), p_guess=None, **_COLD)
    np.random.seed(0)
    r = equilibrium_atmosphere(_earth_target_HCNS(), _earth_ddict(T=T), p_guess=None, **_COLD)
    assert r['H2S_bar'] > 0.0
    assert r['NH3_bar'] > 0.0
    for key in ('P_surf', 'H2S_bar', 'NH3_bar'):
        assert r[key] == pytest.approx(ref[key], rel=1e-12), key
