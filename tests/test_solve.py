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
import re

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
    CO and CH4 modified equilibrium constants, the SO2 pressure (through pO2)
    and the S2 solubility all reference the O'Neill & Eggins (2002) IW buffer.
    At 2500 K O'Neill sits 0.713 dex below the default Hirschmann (2021)
    buffer, so a channel left on the default fails by far more than the
    tolerances. The Dasgupta N2
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
    # janaf_SO2 carries no O2 (fO2_stoich = 0), so SO2 follows the buffer only
    # through pO2: pin it to the O'Neill fO2 itself; on Hirschmann it is 0.713 dex higher.
    so2_k = ModifiedKeq('janaf_SO2', 'oneill')(T, dIW)
    assert so2_k == pytest.approx(ModifiedKeq('janaf_SO2', 'hirschmann')(T, dIW), rel=1e-15)
    so2_oneill = np.sqrt(so2_k * pin['S2']) * 10.0 ** OxygenFugacity('oneill')(T, dIW)
    assert p['SO2'] == pytest.approx(so2_oneill, rel=1e-10)
    so2_default = np.sqrt(so2_k * pin['S2']) * 10.0 ** OxygenFugacity('hirschmann')(T, dIW)
    assert np.log10(so2_default / so2_oneill) == pytest.approx(0.713, abs=2e-3)
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
    """The authoritative-O solve recovers an absolute fO2, and fO2_model sets
    the buffer the shift is reported against. An O budget made on O'Neill at
    IW+2 round-trips to +2 on O'Neill, and the same budget reported against
    the default Hirschmann buffer gives 2 - (IW_Hirschmann - IW_ONeill),
    1.840 at 1800 K. The buffer is not only a label: the Dasgupta N2 law takes
    the shift itself, so the state can change with it. At this T and shift the
    O2 pressure agrees to 1e-3 and P_surf to 1e-4; at 3000 K and IW-2 it does
    not (test_authoritative_o_dissolved_n_follows_the_shift_on_each_buffer).
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
# First guesses [bar] at 1500 K and IW+2 whose first attempt converges on a root with
# negative pH2O, also from 0.5 to 2 times their values; seeded restarts follow it.
_NEG_H2O_GUESS = dict(zip(('H2O', 'CO2', 'N2', 'S2'), (4.0e4, 0.019, 3.7, 4.6e-12)))
_NEG_H2O_GUESS_NO_H = dict(zip(('H2O', 'CO2', 'N2', 'S2'), (1.6, 520.0, 53.0, 3.8e-9)))


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


def _fallback_attempts(caplog):
    """Attempt count in the single fallback warning of a solve."""
    found = re.findall(r'No root closes every element after (\d+) attempts', caplog.text)
    assert len(found) == 1, found
    return int(found[0])


def _fallback_elements(caplog):
    """Elements the single fallback warning of a solve names as over tolerance."""
    _fallback_attempts(caplog)
    return set(re.findall(r'(\w+) \([^()]* kg, tolerance ', caplog.text))


def _over_tolerance(r, target, rtol):
    """Elements whose returned residual exceeds max(rtol * budget, TRUNC_MASS)."""
    return {
        e for e in target if abs(r[e + '_res']) > max(target[e] * rtol, solve_mod.TRUNC_MASS)
    }


def _one_attempt():
    """Single-attempt solve of the low-H inventory at 1500 K and IW+2."""
    return equilibrium_atmosphere(
        dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 1}
    )


@pytest.mark.physics_invariant
def test_low_h_cold_start_rejects_the_negative_water_root(monkeypatch, caplog):
    """A low-H, C-rich cold start at IW+2 and 1500 K returns the physical root.

    A cold start from _NEG_H2O_GUESS converges first on the root with negative
    pH2O (near -2.8e4 bar), which closes the mass balance with no H2O or H2 in
    the atmosphere. It is rejected, and whichever seeded restarts follow, the
    solver goes on to the physical root at 25.05 bar with the Hirschmann (2021)
    IW buffer (25.23 bar with the O'Neill 2002 buffer).
    """
    roots, converged = _spy_roots(monkeypatch)
    for seed in (4, 6, 7):
        roots.clear()
        converged.clear()
        np.random.seed(seed)
        r = equilibrium_atmosphere(
            dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=_NEG_H2O_GUESS, **_COLD
        )
        assert converged[0] and roots[0][0] < -1.0e3, seed
        assert roots[-1][0] > 0.0, seed
        assert r['H2O_bar'] > 0.0, seed
        assert r['P_surf'] == pytest.approx(25.0485, rel=1e-4), seed
        assert abs(r['P_surf'] - 25.233) > 0.05, seed  # O'Neill 2002 buffer
        assert r['H_kg_total'] == pytest.approx(_LOW_H['H'], rel=1e-3), seed
    assert 'No root closes' not in caplog.text


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
    ('target', 'guess', 'seed'),
    [
        ({'H': 0.0, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}, _NEG_H2O_GUESS_NO_H, 0),
        ({'H': 0.0, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}, _NEG_H2O_GUESS_NO_H, 2),
        ({'H': 1.0e15, 'C': 5.0e19, 'N': 1.0e17, 'S': 1.0e18}, _NEG_H2O_GUESS, 9),
    ],
    ids=['no_H_no_N-0', 'no_H_no_N-2', 'trace_H-9'],
)
def test_sub_gate_hydrogen_reports_no_phantom_species(monkeypatch, target, guess, seed):
    """An H budget below the scalar mass gate (1.5e16 kg here) yields a state
    free of H-bearing species from a negative pH2O, and H, C and S close on
    their own tolerance rather than on the scalar gate: dropping the phantom
    species of a negative pH2O would remove a 1e15 kg H budget whole.

    The first guess converges on a root with negative pH2O, whose pH2**2 forms
    CH4 and H2S; seeded restarts follow.
    """
    roots, converged = _spy_roots(monkeypatch)
    np.random.seed(seed)
    r = equilibrium_atmosphere(
        dict(target), _earth_ddict(T=1500.0, dIW=2.0), p_guess=dict(guess), **_COLD
    )
    assert converged[0] and roots[0][0] < 0.0
    assert r['fallback'] is False
    assert r['CH4_bar'] < 1.0e-20
    assert r['H2S_bar'] < 1.0e-12
    assert r['NH3_bar'] < 1.0e-12
    for e in 'HCS':
        assert abs(r[e + '_res']) <= max(target[e] * _COLD['rtol'], solve_mod.TRUNC_MASS), e


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('fO2_model', 'seed', 'p_h2o'), [('hirschmann', 8, -111.8), ('oneill', 8, -139.8)]
)
def test_root_that_empties_hydrogen_is_refused_and_h_is_kept(
    monkeypatch, caplog, fO2_model, seed, p_h2o
):
    """A 1e15 kg H budget lies below the scalar gate (1.5e16 kg), so the gate does
    not see a root that loses it. From _NEG_H2O_GUESS the first attempt converges
    on a root with pH2O near -112 bar (Hirschmann) or -140 bar (O'Neill), which
    loses the whole H budget once clipped. It is refused, and the seeded
    restarts go on to a root that keeps H, within FALLBACK_ATTEMPTS, with no
    fallback and no phantom CH4, H2S or NH3. H, C and S close on their own
    tolerance. The refusal of an unclipped emptied root is pinned by
    test_emptied_root_is_refused_for_a_later_root_that_keeps_the_budget.
    """
    target = {'H': 1.0e15, 'C': 5.0e19, 'N': 1.0e17, 'S': 1.0e18}
    roots, converged = _spy_roots(monkeypatch)
    np.random.seed(seed)
    r = equilibrium_atmosphere(
        dict(target),
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess=dict(_NEG_H2O_GUESS),
        fO2_model=fO2_model,
        **_COLD,
    )
    assert converged[0] and roots[0][0] == pytest.approx(p_h2o, abs=0.1)
    # The other buffer's root lies 28 bar away.
    assert abs(roots[0][0] - (-139.8 if p_h2o > -120.0 else -111.8)) > 20.0
    assert roots[-1][0] > 0.0
    assert 'No root closes' not in caplog.text
    assert r['fallback'] is False
    assert r['CH4_bar'] < 1.0e-20
    assert r['H2S_bar'] < 1.0e-12
    assert r['NH3_bar'] < 1.0e-12
    for e in 'HCS':
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
    assert 'No root closes' not in caplog.text
    np.testing.assert_array_equal(seen[-1][2], 0.0)
    assert seen[-1][4] == pytest.approx(2.0, rel=1e-12)


@pytest.mark.physics_invariant
@pytest.mark.parametrize('seed', [0, 1, 2])
def test_trace_sulfur_with_absent_nitrogen_returns_through_the_fallback(caplog, seed):
    """An S budget below the scalar mass gate with no N: no root keeps S, so every
    root that passes the gate is refused. FALLBACK_ATTEMPTS after the first of
    them the best returns through the fallback, well before the 1000 attempts
    run out, and its warning names S. pN2 and pS2 are reported clipped to zero
    and H, C close. Which refused state ranks best depends on the seed, and the
    states differ within rtol (1e-4): over seeds 0 to 7, P_surf spans 25.1335
    to 25.1352 bar, so it is pinned to rel 1e-4."""
    np.random.seed(seed)
    target = {'H': 1.0e20, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e14}
    r = equilibrium_atmosphere(target, _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **_COLD)
    attempts = _fallback_attempts(caplog)
    assert solve_mod.FALLBACK_ATTEMPTS < attempts < _COLD['nguess']
    assert 'tolerance for S (1e+14 kg, tolerance 1e+10 kg, budget 1e+14 kg)' in caplog.text
    np.testing.assert_array_equal([r['N2_bar'], r['S2_bar'], r['S_kg_total']], 0.0)
    assert r['P_surf'] == pytest.approx(25.1349, rel=1e-4)
    assert abs(r['P_surf'] - 25.318) > 0.05  # O'Neill 2002 buffer
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
    assert abs(r['P_surf'] - 9.272) > 0.05  # O'Neill 2002 buffer
    assert r['H2O_bar'] == pytest.approx(h2o, rel=1e-5)
    assert 'No root closes' not in caplog.text
    assert r['fallback'] is False


@pytest.mark.parametrize('n2', [1.0e-30, -1.0e-30])
def test_verdict_does_not_depend_on_an_inert_sign(monkeypatch, caplog, n2):
    """A root whose N residual (5e15 kg) passes the scalar gate (1.5e16 kg) but
    not the N tolerance (1e13 kg) is refused whatever the sign of an inert pN2,
    because clipping that pN2 changes no residual: with no attempt left it
    returns through the fallback, whose warning names N."""
    seen = _stub_buffered(monkeypatch, [1.0, 5.0, n2, 0.5], lambda x: [0.0, 0.0, 5.0e15, 0.0])
    r = _one_attempt()
    assert r['H2O_bar'] == pytest.approx(1.0, rel=1e-12)
    np.testing.assert_array_equal(r['N2_bar'], max(n2, 0.0))
    assert 'tolerance for N (5e+15 kg, tolerance 1e+13 kg' in caplog.text
    # The raw root is judged first; the clipped state follows it only for the
    # negative sign (later evaluations are the polish and the reported state).
    assert seen[0][2] == n2
    assert (seen[1][2] == pytest.approx(0.0, abs=1.0e-40)) == (n2 < 0.0)


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
    assert 'No root closes' not in caplog.text
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
    assert 'No root closes' not in caplog.text
    # From _NEG_H2O_GUESS the single attempt is that root, so the solve raises.
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere(
            dict(_LOW_H),
            _earth_ddict(T=1500.0, dIW=2.0),
            p_guess=dict(_NEG_H2O_GUESS),
            **{**_COLD, 'nguess': 1},
        )


_TOL_N, _TOL_S = _LOW_H['N'] * _COLD['rtol'], _LOW_H['S'] * _COLD['rtol']


_N_TWICE = [0.0, 0.0, 2.0 * _TOL_N, 0.0]


@pytest.mark.parametrize(
    ('first', 'second', 'co2'),
    [
        (_N_TWICE, [0.0, 0.0, 1.5 * _TOL_N, 0.0], 7.0),
        (_N_TWICE, [0.0, 0.0, 3.0 * _TOL_N, 0.0], 5.0),
        (_N_TWICE, [0.0, 0.0, 2.0 * _TOL_N, 0.0], 5.0),
        (_N_TWICE, [0.0, 0.0, 1.1 * _TOL_N, 1.1 * _TOL_S], 5.0),
        (_N_TWICE, [0.0, 0.0, 1.0e17, 0.0], 5.0),
        (
            [0.0, 0.0, 2.0 * _TOL_N, 1.05 * _TOL_S],
            [0.0, 0.0, 1.9 * _TOL_N, 1.9 * _TOL_S],
            7.0,
        ),
    ],
    ids=['closer', 'further', 'tied', 'more_elements', 'outside_gate', 'worst_not_sum'],
)
def test_fallback_is_the_best_gate_passing_root(monkeypatch, caplog, first, second, co2):
    """Of two attempts whose clipped states miss closure, the first misses N by
    twice its tolerance. The solver falls back to the second only if it misses
    fewer elements, or as many by a smaller largest ratio of residual to
    tolerance, and passes the mass gate; a tie keeps the first. With N and S
    missed by both, the second wins on its largest ratio (1.9 against 2.0)
    although its ratios sum higher (3.8 against 3.05)."""
    roots = iter([np.array([-1.0, 5.0, 1.0, 0.5]), np.array([-1.0, 7.0, 1.0, 0.5])])
    monkeypatch.setattr(solve_mod.opt, 'fsolve', lambda *a, **k: (next(roots), {}, 1, 'stub'))

    def _residual(x, *a):
        if x[0] < 0.0:
            return [0.0, 0.0, 0.0, 0.0]
        return first if x[1] < 6.0 else second

    monkeypatch.setattr(solve_mod, 'func', _residual)
    r = equilibrium_atmosphere(
        dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), p_guess=None, **{**_COLD, 'nguess': 2}
    )
    assert 'after 2 attempts' in caplog.text
    assert r['CO2_bar'] == pytest.approx(co2, rel=1e-12)
    kept = second if co2 > 6.0 else first
    assert r['N_res'] == pytest.approx(kept[2], rel=1e-12)
    assert r['S_res'] == pytest.approx(kept[3], rel=1e-12)


def test_fallback_search_stops_fallback_attempts_after_the_first(monkeypatch, caplog):
    """Once a state missing closure is kept, the solve makes FALLBACK_ATTEMPTS
    more attempts, counted from the first such state even when a better one
    follows, and returns the best; with no such state it uses every attempt
    and raises. Attempts 3 and 50 of the first solve are negative-pH2O roots
    whose clipped states miss N by twice and by 1.5 times its tolerance; every
    other root misses H by 1e17 kg, outside the 1.5e16 kg mass gate."""
    calls, candidate_at = [], {3: 5.0, 50: 7.0}

    def _fsolve(*a, **k):
        calls.append(1)
        h2o = -1.0 if len(calls) in candidate_at else 1.0e3
        return np.array([h2o, candidate_at.get(len(calls), 5.0), 1.0, 0.5]), {}, 1, 'stub'

    def _residual(x, *a):
        if x[0] > 0.0:
            return [1.0e17, 0.0, 0.0, 0.0]
        n_res = (1.5 if x[1] > 6.0 else 2.0) * _TOL_N if x[0] == 0.0 else 0.0
        return [0.0, 0.0, n_res, 0.0]

    monkeypatch.setattr(solve_mod.opt, 'fsolve', _fsolve)
    monkeypatch.setattr(solve_mod, 'func', _residual)
    nguess = solve_mod.FALLBACK_ATTEMPTS + 50
    solve = dict(p_guess=None, **{**_COLD, 'nguess': nguess})
    r = equilibrium_atmosphere(dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), **solve)
    # docs/Explanations/mass_balance.md quotes 100 and the seeds of a solvable
    # H 1e15 kg budget it loses; a change of the limit must update that page.
    assert solve_mod.FALLBACK_ATTEMPTS == 100
    assert 50 < 3 + solve_mod.FALLBACK_ATTEMPTS < 50 + solve_mod.FALLBACK_ATTEMPTS
    assert len(calls) == 3 + solve_mod.FALLBACK_ATTEMPTS
    assert _fallback_attempts(caplog) == len(calls)
    assert r['N_res'] == pytest.approx(1.5 * _TOL_N, rel=1e-12)
    assert r['CO2_bar'] == pytest.approx(7.0, rel=1e-12)
    calls.clear()
    candidate_at.clear()
    with pytest.raises(RuntimeError, match='Could not find solution'):
        equilibrium_atmosphere(dict(_LOW_H), _earth_ddict(T=1500.0, dIW=2.0), **solve)
    assert len(calls) == nguess


# Library tolerances for H 1e12 kg, C 5e19 kg, no N and S 1e18 kg: the per-element
# tolerances and the scalar mass gate the polish is judged against.
_POLISH_TOL = np.array([1.0e7, 5.0e14, solve_mod.TRUNC_MASS, 1.0e13])
_POLISH_GATE = 5.0e19 * 1e-5 + 1e10 + solve_mod.TRUNC_MASS


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('s_root', 'at_root', 'nsolve', 'polished', 'converged'),
    [
        (0.7, {}, 1500, True, True),
        (0.7, {}, 2, True, False),
        (-0.1, {}, 1500, False, True),
        (0.7, {0: -1.0e17}, 1500, False, True),
        (0.7, {2: 100.0}, 1500, False, True),
        (None, {}, 1500, False, False),
    ],
    ids=['closes_s', 'stops_at_maxfev', 'negative_root', 'leaves_gate', 'no_better', 'flat'],
)
def test_polished_fallback_is_kept_only_when_valid_and_better(
    monkeypatch, s_root, at_root, nsolve, polished, converged
):
    """A clipped fallback with pH2O and pN2 at zero has lost its H (1e12 kg) and
    misses S by 1.59e13 kg against a 1e13 kg tolerance, as a negative-pH2O root
    does once its phantom H2S is dropped. Re-solving pCO2 and pS2 with the zero
    primaries held at zero gives a state that closes S, and it replaces the
    fallback. The polished state is judged by its own residuals, not by the
    solver's verdict: a re-solve stopped by maxfev (nsolve = 2) after it has
    reached the root reports failure and is kept too. The fallback is kept as it is
    when the polished state has a negative pS2, lies outside the scalar mass
    gate (H at 1e17 kg), misses as many elements by as much (N at 100 kg for S),
    or does not move (an S residual flat in pS2)."""
    state = np.array([0.0, 6.0, 0.0, 0.7 - 1.59e-5])
    # S residual [kg] is -1.59e13 at the fallback's pS2 and crosses zero at s_root
    # (flat when s_root is None); C closes at pCO2 = 6 bar, H stays lost.
    slope = 0.0 if s_root is None else 1.59e13 / (s_root - state[3])

    def _residual(x, *a):
        res = [-1.0e12, (x[1] - 6.0) * 1.0e19, 0.0, (x[3] - state[3]) * slope - 1.59e13]
        if abs(x[3] - 0.7) < 1.0e-6:
            for i, v in at_root.items():
                res[i] = v
        return res

    verdicts, real_root = [], solve_mod.opt.root

    def _root(*a, **k):
        result = real_root(*a, **k)
        verdicts.append(bool(result.success))
        return result

    monkeypatch.setattr(solve_mod, 'func', _residual)
    monkeypatch.setattr(solve_mod.opt, 'root', _root)
    resid = np.abs(np.asarray(_residual(state)))
    assert solve_mod._fallback_rank(resid, _POLISH_TOL)[0] == 2
    out, out_resid = solve_mod._polish_fallback(
        state, resid, {}, {}, _POLISH_TOL, _POLISH_GATE, 1e-8, nsolve
    )
    assert verdicts == [converged]
    assert np.all(out >= 0.0)
    np.testing.assert_array_equal(out[[0, 2]], 0.0)
    np.testing.assert_allclose(out_resid, np.abs(_residual(out)), rtol=1e-12)
    if polished:
        assert out[3] == pytest.approx(0.7, abs=1.0e-6)
        assert out_resid[3] <= _POLISH_TOL[3]
        assert solve_mod._fallback_rank(out_resid, _POLISH_TOL) == (1, 1.0e5)
    else:
        np.testing.assert_array_equal(out, state)
        np.testing.assert_array_equal(out_resid, resid)


@pytest.mark.parametrize(
    ('he_at_root', 'polished'), [(0.0, True), (100.0, False)], ids=['he_closes', 'he_misses']
)
def test_polished_fallback_must_close_each_noble_gas(monkeypatch, he_at_root, polished):
    """With He active and held at zero, the re-solve closes S and halves the H
    loss, so it ranks better (H at 5e11 kg over its 1e7 kg tolerance is the
    worst miss) and keeps the CHNOS elements inside the scalar gate. It is kept
    when He still closes; a He residual of 100 kg against its 10 kg tolerance
    fails the per-gas gate and the fallback is returned as it is."""
    tol = np.append(_POLISH_TOL, solve_mod.TRUNC_MASS)
    state = np.array([0.0, 6.0, 0.0, 0.7 - 1.59e-5, 0.0])

    def _residual(x, *a):
        at_root = abs(x[3] - 0.7) < 1.0e-6
        h, he = (-5.0e11, he_at_root) if at_root else (-1.0e12, 0.0)
        return [h, (x[1] - 6.0) * 1.0e19, 0.0, (x[3] - 0.7) * 1.0e18, he]

    monkeypatch.setattr(solve_mod, 'func', _residual)
    resid = np.abs(np.asarray(_residual(state)))
    out, out_resid = solve_mod._polish_fallback(
        state, resid, {}, {}, tol, _POLISH_GATE, 1e-8, 1500
    )
    if polished:
        assert out[3] == pytest.approx(0.7, abs=1.0e-6)
        assert solve_mod._fallback_rank(out_resid, tol) == (1, 5.0e4)
    else:
        np.testing.assert_array_equal(out, state)
        np.testing.assert_array_equal(out_resid, resid)
    # Even with the He miss the polished state (H and He over, worst ratio 5e4)
    # outranks the fallback, so only the per-gas gate can reject it.
    assert (2, 5.0e4) < solve_mod._fallback_rank(resid, tol) == (2, 1.0e5)


def test_polish_leaves_a_fallback_with_no_zero_primary_untouched(monkeypatch):
    """With every primary positive there is nothing to hold at zero, so the
    fallback is returned as it is and the residual is not evaluated again; with
    every primary at zero there is nothing to re-solve."""
    calls = []
    monkeypatch.setattr(solve_mod, 'func', lambda x, *a: calls.append(x) or [0.0] * 4)
    resid = np.array([1.0e12, 0.0, 0.0, 1.59e13])
    for state in (np.array([1.0e-3, 6.0, 1.0e-9, 0.7]), np.zeros(4)):
        out, out_resid = solve_mod._polish_fallback(
            state, resid, {}, {}, _POLISH_TOL, _POLISH_GATE, 1e-8, 1500
        )
        assert out is state
        assert out_resid is resid
    assert calls == []


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
    assert 'No root closes' not in caplog.text
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
    assert 'No root closes' not in caplog.text
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
    """With the library tolerances and fsolve only, no attempt keeps a trace H
    budget (1e12 kg): every root that passes the mass gate loses it, and most
    are clipped negative-pH2O roots whose dropped H2S misses S by 1.59e13 kg
    against a 1e13 kg tolerance. The solve returns the closest of them
    FALLBACK_ATTEMPTS after the first, well before the 500 attempts run out,
    re-solved with pH2O held at zero so that C and S close again. It is flagged,
    and its warning names H alone, with its residual, tolerance and budget."""
    target = {'H': 1.0e12, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}
    np.random.seed(0)
    r = equilibrium_atmosphere(
        dict(target),
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess=None,
        xtol=1e-8,
        rtol=1e-5,
        atol=1e10,
        nguess=500,
        nsolve=1500,
        print_result=False,
        opt_solver=False,
    )
    assert solve_mod.FALLBACK_ATTEMPTS < _fallback_attempts(caplog) < 500
    assert r['fallback'] is True
    assert 'tolerance for H (1e+12 kg, tolerance 1e+07 kg, budget 1e+12 kg)' in caplog.text
    assert _fallback_elements(caplog) == _over_tolerance(r, target, 1e-5) == {'H'}
    assert abs(r['S_res']) <= 1.0e18 * 1e-5
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
    [(-1.0e15, True), (-0.5e15, True), (1.0e15, True), (1.01e11, True), (0.99e11, False)],
    ids=['emptied', 'half_lost', 'overfilled', 'just_over', 'just_under'],
)
def test_any_element_beyond_its_tolerance_is_refused(monkeypatch, caplog, n_res, refused):
    """A 1e15 kg N budget lies inside the scalar gate (1.5e16 kg), which would
    pass any of these states. A state whose N residual exceeds the N tolerance
    (1e11 kg), by losing, halving or overfilling the budget or by 1 %, is refused,
    and with no attempt left it returns through the fallback, whose warning names
    N; a state 1 % inside the tolerance is accepted without a warning."""
    r = _empty_budget_solve(monkeypatch, [0.0, 0.0, n_res, 0.0])
    assert ('No root closes' in caplog.text) == refused
    assert (f'tolerance for N ({abs(n_res):.3g} kg' in caplog.text) == refused
    np.testing.assert_array_equal([r[g + '_bar'] for g in _PRIMARIES], _STUB_ROOT)


def test_emptied_root_is_refused_for_a_later_root_that_keeps_the_budget(monkeypatch, caplog):
    """An all-positive root that empties a 1e15 kg H budget passes the scalar
    gate (1.5e16 kg) and has nothing to clip, yet is refused; the solve goes on
    to the next attempt, whose root closes every element, and returns it with
    no fallback warning and the fallback flag unset, although a fallback was
    stored on the way."""
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
    assert r['fallback'] is False
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


def test_fallback_warns_on_every_call(monkeypatch, caplog):
    """A second solve that falls back on the same emptied budget warns again,
    so each step of a coupled run that falls back is visible, and is flagged
    again. A solve in between that closes every element, with N 1 % inside its
    tolerance, neither warns nor is flagged."""
    caplog.set_level(logging.DEBUG, logger='fwl.calliope.solve')
    flags = [
        _empty_budget_solve(monkeypatch, n_res)['fallback']
        for n_res in (
            [0.0, 0.0, -1.0e15, 0.0],
            [0.0, 0.0, 0.99e11, 0.0],
            [0.0, 0.0, -1.0e15, 0.0],
        )
    ]
    assert flags == [True, False, True]
    records = [r for r in caplog.records if 'No root closes' in r.getMessage()]
    assert [r.levelno for r in records] == [logging.WARNING, logging.WARNING]
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
    ],
    ids=['H1e12_N0', 'H1e12_N1e14'],
)
def test_trace_hydrogen_falls_back_to_the_clipped_root(caplog, target, seed):
    """At library tolerances, a 1e12 kg H budget is closed only by roots whose
    negative pH2O carries the H in CH4 and H2S (H 1e13 kg is not used: a root
    that closes every element exists there). Clipping such a
    root empties the H budget, far beyond the H tolerance, so it is rejected, and
    the dropped phantom H2S takes S beyond its tolerance (1.59e13 kg against
    1e13 kg for 1e12 kg of H). The solve returns the best of them clipped and
    re-solved with pH2O held at zero, which closes C, N and S again: the warning
    lists H alone, and no H-bearing species is reported."""
    np.random.seed(seed)
    r = equilibrium_atmosphere(
        dict(target),
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess=None,
        opt_solver=False,
        **_LIBRARY,
    )
    assert r['fallback'] is True
    assert _fallback_elements(caplog) == _over_tolerance(r, target, _LIBRARY['rtol']) == {'H'}
    assert r['H_kg_total'] == pytest.approx(0.0, abs=target['H'] * _LIBRARY['rtol'])
    np.testing.assert_array_equal(
        [r[s + '_bar'] for s in ('H2O', 'H2', 'CH4', 'H2S', 'NH3')], 0.0
    )
    assert r['P_surf'] == pytest.approx(25.0274, rel=1e-4)
    assert abs(r['P_surf'] - 25.212) > 0.05  # O'Neill 2002 buffer
    # The returned clipped state still passes the scalar mass gate.
    gate = max(target.values()) * _LIBRARY['rtol'] + _LIBRARY['atol'] + solve_mod.TRUNC_MASS
    assert max(abs(r[e + '_res']) for e in 'HCNS') <= gate


@pytest.mark.physics_invariant
def test_library_defaults_warn_for_every_element_a_trace_state_misses(caplog):
    """With the library defaults (opt_solver on) and a 1e12 kg H budget at
    1500 K and IW+2, no root keeps H, and a returned state can miss a trace
    element by orders of magnitude while passing the scalar mass gate
    (5e14 kg). Such a state is never returned silently: it is flagged, and the
    warning names every element whose residual exceeds its own tolerance and no
    other. The first guess (pH2O 800 bar, pCO2 1600 bar) reaches a gate-passing
    root that loses H on the first attempt (from 0.5 to 1.3 times it as well),
    so 5 attempts keep this fast; N2 and S2 start above 1e-10 bar so the
    solver bounds stay the library ones."""
    target = {'H': 1.0e12, 'C': 5.0e19, 'N': 0.0, 'S': 1.0e18}
    r = equilibrium_atmosphere(
        dict(target),
        _earth_ddict(T=1500.0, dIW=2.0),
        p_guess={'H2O': 800.0, 'CO2': 1600.0, 'N2': 1.0e-9, 'S2': 1.0e-9},
        nguess=5,
        print_result=False,
    )
    assert _fallback_attempts(caplog) == 5
    assert r['fallback'] is True
    over = _over_tolerance(r, target, 1e-5)
    assert _fallback_elements(caplog) == over
    # The re-solve with pH2O held at zero closes the S the clip took away.
    assert over == {'H'}
    # The scalar gate alone would pass this state: H misses by its whole budget,
    # 1e5 times its tolerance, yet sits 500 times inside the gate.
    gate = 5.0e19 * 1e-5 + 1e10 + solve_mod.TRUNC_MASS
    assert 1.0e12 * 1e-5 < abs(r['H_res']) < gate / 100.0
    assert r['H_kg_total'] == pytest.approx(0.0, abs=1.0e12 * 1e-5)


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
    assert abs(r['P_surf'] - 25.233) > 0.05  # O'Neill 2002 buffer
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
