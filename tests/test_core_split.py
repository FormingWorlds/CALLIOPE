"""Tests for `src/calliope/core_split.py`: the core and mantle split of a bulk inventory.

Contract clauses: at fixed mass-basis D the core holds D x / (1 + D x) of each element,
x = M_metal / M_silicate; the element ledger bulk = core + mantle holds to 1e-12 on every
split with core and mantle >= 0; a law outside its range is evaluated
there (no clamp) and flagged with its distance per axis; O and Si are refused as bulk
elements; every input outside its domain raises ValueError naming it. Hypothesis tests
use ``derandomize=True``. Commands: docs/How-to/build_tests.md.
"""

from __future__ import annotations

import json
import logging
import math
import platform
from pathlib import Path

import pytest

from calliope import core_split, metal
from calliope.core_split import DEFAULT_LAWS, split_core_mantle
from calliope.metal import get_law, log10_D, mass_D, register_law

hypothesis = pytest.importorskip('hypothesis')
from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Earth-like masses [kg] and bulk inventories [kg]
M_METAL, M_SIL = 1.94e24, 4.03e24
BULK = {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}
DECL = {
    'basis': 'mass',
    'P_range_GPa': (0.0, 100.0),
    'T_range_K': (1000.0, 5000.0),
    'dIW_range': (-5.0, 0.0),
    'source': 'constant law of the tests',
    'doi': '10.0000/none',
}


@pytest.fixture(autouse=True)
def _own_registry(monkeypatch):
    monkeypatch.setattr(metal, '_LAWS', dict(metal._LAWS))
    monkeypatch.setattr(core_split, '_WARNED', set())


def _const(name, log10_d):
    register_law('N', name, lambda *a: log10_d, **DECL)
    return name


def _ledger(bulk, core, mantle):
    for el, m in bulk.items():
        assert core[el] >= 0.0 and mantle[el] >= 0.0
        assert abs(core[el] + mantle[el] - m) <= 1e-12 * m


def _x_fe(w):
    """x'_Fe of a metal of light-element mass fractions w (H and C left out of the basis)."""
    n = {el: v / metal._mass(el) for el, v in w.items()}
    n_fe = (1.0 - sum(w.values())) / metal._mass('Fe')
    return n_fe / (n_fe + sum(v for el, v in n.items() if el not in ('H', 'C')))


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_constant_d_gives_the_closed_form_core_fraction(validation_plot, legend_clear):
    """Mass balance: a constant D puts D x / (1 + D x) of the bulk into the core."""
    x = M_METAL / M_SIL
    ds = [0.1, 1.0, 10.0, 100.0, 1000.0]
    f = []
    for k, d in enumerate(ds):
        core, mantle, _ = split_core_mantle(
            {'N': 1e19}, M_METAL, M_SIL, dIW=-2.0, laws={'N': _const(f'c{k}', math.log10(d))}
        )
        _ledger({'N': 1e19}, core, mantle)
        f.append(core['N'] / 1e19)
        # the mantle mass fraction is the bulk over M_silicate + D M_metal
        assert mantle['N'] / M_SIL == pytest.approx(1e19 / (M_SIL + d * M_METAL), rel=1e-12)
    assert f == pytest.approx([d * x / (1 + d * x) for d in ds], rel=1e-12)
    # D x without the melt back-reaction overshoots by more than 10 % from D = 1 on
    assert all(d * x - fi > 0.1 * fi for fi, d in zip(f[1:], ds[1:]))

    def figure():
        import matplotlib.pyplot as plt

        grid = [10 ** (k / 20) for k in range(-25, 66)]
        fig, ax = plt.subplots()
        ax.plot(grid, [d * x / (1 + d * x) for d in grid], label=r'$Dx/(1+Dx)$')
        ax.plot(ds, f, 'o', label='split_core_mantle')
        ax.set(
            xscale='log',
            xlabel='D (metal / silicate, mass)',
            ylabel='core fraction of the bulk',
        )
        ax.minorticks_off()
        ax.legend(loc='upper left')
        fig.tight_layout()
        return fig

    legend_clear(figure())
    validation_plot(figure, 'core_split_fraction')


@pytest.mark.physics_invariant
def test_default_laws_split_an_earth_like_bulk():
    """Default laws at 60 GPa, 4000 K, dIW -1.5: ledger closes, metal composition is a fixed point."""
    core, mantle, flags = split_core_mantle(
        BULK, M_METAL, M_SIL, dIW=-1.5, w_fixed={'O': 0.02, 'Si': 0.03}
    )
    _ledger(BULK, core, mantle)
    assert set(core) == set(BULK) and DEFAULT_LAWS == {
        'H': 'tagawa2021_eq2',
        'C': 'fischer2020_eq1',
        'N': 'suer2023_line',
        'S': 'suer2017_eq7',
    }
    # the D of the returned metal and melt reproduce the returned split
    w = {'O': 0.02, 'Si': 0.03, **{el: core[el] / M_METAL for el in BULK}}
    w_melt = {el: mantle[el] / M_SIL for el in ('H', 'C')}
    for el in BULK:
        d = mass_D(get_law(el, DEFAULT_LAWS[el]), 60.0, 4000.0, -1.5, w, w_melt)
        assert core[el] / M_METAL == pytest.approx(d * mantle[el] / M_SIL, rel=1e-10)
    # S follows Suer 2017, not the Suer 2023 line of the placeholder file (D 10 to 173)
    alt3, _, _ = split_core_mantle(BULK, M_METAL, M_SIL, dIW=-3.0, laws={'S': 'suer2023_line'})
    ref3, _, _ = split_core_mantle(BULK, M_METAL, M_SIL, dIW=-3.0)
    assert abs(alt3['S'] - ref3['S']) > 0.01 * BULK['S']
    # at 60 GPa and 4000 K the S law is inside its range, the H law (dIW -1.33 to -0.88) is not
    assert flags['S']['out_of_range'] == 0 and flags['H']['out_of_range'] == 1


@pytest.mark.physics_invariant
def test_out_of_range_is_evaluated_unclamped_and_flagged(caplog):
    """300 GPa, 5500 K (a CMB point): D from the law at that point, flag 1 and the distances."""
    with caplog.at_level(logging.WARNING, logger='fwl.calliope.core_split'):
        core, mantle, flags = split_core_mantle(
            {'S': 1e20}, M_METAL, M_SIL, dIW=-1.5, P_eq=300.0, T_eq=5500.0
        )
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-1.5, P_eq=300.0, T_eq=5500.0)
    assert {k: v for k, v in flags['S'].items() if k != 'D_residual'} == {
        'out_of_range': 1,
        'P_GPa': 209.0,
        'T_K': 1400.0,
        'dIW': 0.0,
    }
    assert 0.0 <= flags['S']['D_residual'] <= 1e-12
    w = {'S': core['S'] / M_METAL}
    d = mass_D(get_law('S', 'suer2017_eq7'), 300.0, 5500.0, -1.5, w, {})
    assert core['S'] / mantle['S'] == pytest.approx(d * M_METAL / M_SIL, rel=1e-10)
    # a clamp to (91 GPa, 4100 K) gives a different D
    clamped = mass_D(get_law('S', 'suer2017_eq7'), 91.0, 4100.0, -1.5, w, {})
    assert abs(math.log10(clamped / d)) > 0.1
    assert sum('outside its range' in r.message for r in caplog.records) == 1


def test_zero_metal_and_zero_bulk():
    """M_metal 0 leaves every element in the mantle; a bulk of 0 gives 0 in both."""
    core, mantle, flags = split_core_mantle(BULK, 0.0, M_SIL, dIW=-2.0)
    assert core == {el: 0.0 for el in BULK} and mantle == BULK
    assert all(flags[el]['D_residual'] == 0.0 for el in BULK)
    core, mantle, _ = split_core_mantle({'S': 0.0}, M_METAL, M_SIL, dIW=-2.0)
    assert core == {'S': 0.0} and mantle == {'S': 0.0}


@pytest.mark.parametrize('el', ['O', 'Si'])
def test_exchange_elements_are_refused_as_bulk(el):
    """O and Si enter only through w_fixed: their bulk is rock, outside the volatile inventory."""
    with pytest.raises(ValueError, match=f'{el}.*w_fixed'):
        split_core_mantle({el: 1e20}, M_METAL, M_SIL, dIW=-2.0)


BAD = [True, 'x', None, math.nan, math.inf, -math.inf, -1.0]
OK = {'bulk_kg': {'S': 1e20}, 'M_metal': M_METAL, 'M_silicate': M_SIL, 'dIW': -2.0}


@pytest.mark.parametrize('value', BAD)
@pytest.mark.parametrize('field', ['bulk_kg', 'M_metal', 'M_silicate', 'P_eq', 'T_eq'])
def test_bad_input_is_refused_by_name(field, value):
    """A boolean, string, None, NaN, inf or negative value raises ValueError naming the field."""
    with pytest.raises(ValueError, match=field):
        split_core_mantle(**{**OK, field: {'S': value} if field == 'bulk_kg' else value})
    assert core_split._WARNED == set()


@pytest.mark.parametrize(
    'field, value',
    [('M_silicate', 0.0), ('P_eq', 0.0), ('T_eq', 0.0), ('dIW', math.nan), ('dIW', True)],
)
def test_zero_silicate_pressure_temperature_and_bad_diw_are_refused(field, value):
    """M_silicate, P_eq and T_eq must be above 0; dIW must be a finite number."""
    with pytest.raises(ValueError, match=field):
        split_core_mantle(**{**OK, field: value})
    assert core_split._WARNED == set()


@pytest.mark.parametrize(
    'w_fixed',
    [{'O': 1.0}, {'O': 0.6, 'Si': 0.5}, {'O': -0.1}, {'Si': True}, {'S': 0.1}, {'Fe': 0.1}],
)
def test_bad_fixed_metal_is_refused(w_fixed):
    """w_fixed: elements outside the bulk in [0, 1) with a sum below 1, and not an element split."""
    with pytest.raises(ValueError, match='w_fixed'):
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-2.0, w_fixed=w_fixed)


def test_unknown_element_or_law_is_refused():
    """A bulk element without a law, or a law name that is not registered, raises by name."""
    with pytest.raises(ValueError, match='Ar'):
        split_core_mantle({'Ar': 1e15}, M_METAL, M_SIL, dIW=-2.0)
    with pytest.raises(metal.PartitionTableError, match='no law'):
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-2.0, laws={'S': 'none'})


def test_fixed_point_that_does_not_settle_raises(monkeypatch):
    """Plain passes that swing, then relaxed passes from the start metal; past the limits RuntimeError."""
    seen = []
    real = core_split.mass_D

    def spy(entry, p, t, d_iw, w, *rest):
        seen.append(w['S'])
        return real(entry, p, t, d_iw, w, *rest)

    monkeypatch.setattr(core_split, 'mass_D', spy)
    monkeypatch.setattr(core_split, 'MAX_PASSES', 3)
    monkeypatch.setattr(core_split, 'MAX_RELAXED_PASSES', 4)
    with pytest.raises(RuntimeError, match='4 relaxed passes'):
        split_core_mantle({'S': 2e22}, 1e22, 4e24, dIW=-1.5)
    assert len(seen) == 3 + 2 * 4
    assert seen[3] == seen[7] == 0.0  # each relaxed sweep starts again from the start metal


def test_broken_ledger_raises(monkeypatch):
    """A split whose core and mantle do not sum to the bulk raises RuntimeError naming the element."""
    real = core_split._masses

    def leak(bulk, d, m_metal, m_sil):
        core, mantle = real(bulk, d, m_metal, m_sil)
        return core, {el: 1.001 * v for el, v in mantle.items()}

    monkeypatch.setattr(core_split, '_masses', leak)
    with pytest.raises(RuntimeError, match='S ledger'):
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-2.0)


@pytest.mark.physics_invariant
@settings(derandomize=True, max_examples=200, deadline=None)
@given(
    frac=st.floats(0.0, 0.5),
    log10_d=st.floats(-6.0, 6.0),
    m_metal=st.floats(1e15, 1e24),
    m_sil=st.floats(1e15, 1e24),
)
def test_ledger_closes_over_the_domain(frac, log10_d, m_metal, m_sil):
    """Ledger: bulk = core + mantle to 1e-12, core and mantle >= 0, for D 1e-6 to 1e6 and masses 1e15 to 1e24 kg.

    The bulk is ``frac`` M_metal with frac up to 0.5, so the metal keeps its iron; 2 ulp of 0
    cover the round-off of a subnormal bulk.
    """
    bulk, name = frac * m_metal, f'd{log10_d!r}'
    if ('N', name) not in metal._LAWS:
        _const(name, log10_d)
    core, mantle, _ = split_core_mantle({'N': bulk}, m_metal, m_sil, dIW=-2.0, laws={'N': name})
    assert core['N'] >= 0.0
    assert mantle['N'] >= 0.0
    assert abs(core['N'] + mantle['N'] - bulk) <= 1e-12 * bulk + 2 * math.ulp(0.0)


def test_registered_callable_is_used_by_the_split():
    """A registered law for S replaces the default when named in laws, with the D it returns."""
    register_law('S', 'test_d100', lambda *a: 2.0, **DECL)
    core, mantle, flags = split_core_mantle(
        {'S': 1e20}, M_METAL, M_SIL, dIW=-2.0, laws={'S': 'test_d100'}
    )
    assert core['S'] / M_METAL == pytest.approx(100.0 * mantle['S'] / M_SIL, rel=1e-12)
    assert flags['S']['out_of_range'] == 0
    assert log10_D(get_law('S', 'test_d100'), 1.0, 2000.0, 0.0) == pytest.approx(2.0, rel=1e-12)


S_ENTRY = get_law('S', 'suer2017_eq7')
PTE = metal.PartitionTableError


@pytest.mark.parametrize(
    'laws, err, word',
    [
        ({'S': 'suer2023_line', 's': 'none'}, ValueError, 'outside bulk_kg'),
        ([], ValueError, 'must be a dict'),
        (['S'], ValueError, 'must be a dict'),
        ({'S': None}, PTE, 'needs a law entry'),
        ({'S': lambda *a: 1.0}, PTE, 'needs a law entry'),
        ({'S': 'N-entry'}, PTE, 'needs a law entry'),
        ({'S': {'element': 'S', 'form': 'log10_linear'}}, PTE, "'quantity'"),
        ({'S': {'element': 'S', 'form': 'callable'}}, PTE, 'func'),
        ({'S': {**S_ENTRY, 'a': '1'}}, PTE, "'a' must be a finite number"),
        ({'S': {**S_ENTRY, 'T_range_K': 5}}, PTE, "'T_range_K' must be"),
        ({'S': {**S_ENTRY, 'T_range_K': {1873.0: 0, 4100.0: 0}}}, PTE, "'T_range_K' must be"),
    ],
    ids=[
        'typo-key',
        'empty-list',
        'list',
        'none',
        'bare-callable',
        'other-element-entry',
        'entry-without-coefficients',
        'callable-without-func',
        'string-coefficient',
        'scalar-range',
        'dict-range',
    ],
)
def test_laws_argument_is_checked(laws, err, word):
    """A bad laws argument or a bad law raises its own error type naming the fault."""
    if laws == {'S': 'N-entry'}:
        laws = {'S': get_law('N', 'suer2023_line')}
    with pytest.raises(err, match=word) as info:
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-1.5, laws=laws)
    assert type(info.value) is err
    assert core_split._WARNED == set()


def test_overflow_of_d_times_metal_raises_a_named_error():
    """log10 D 290 with 1e24 kg of metal overflows: ValueError naming the element, no NaN."""
    with pytest.raises(ValueError, match="'N'.*not finite"):
        split_core_mantle({'N': 1e19}, 1e24, M_SIL, dIW=-2.0, laws={'N': _const('huge', 290.0)})
    assert core_split._WARNED == set()


def test_bulk_beyond_the_float_range_raises_and_a_bulk_above_the_metal_mass_settles():
    """A bulk of 10**400 kg raises ValueError naming it; S above the metal mass finds its root past FeS."""
    with pytest.raises(ValueError, match='bulk_kg'):
        split_core_mantle({'S': 10**400}, M_METAL, M_SIL, dIW=-1.5)
    core, mantle, flags = split_core_mantle({'S': 5e24}, M_METAL, M_SIL, dIW=-1.5)
    assert core['S'] / M_METAL == pytest.approx(
        _bisect_s_root(5e24, M_METAL, M_SIL, -1.5), rel=1e-9
    )
    assert flags['metal']['beyond_FeS'] == 1


def test_each_law_warns_once_with_rounded_distances(caplog):
    """A law warns once; a law of other ranges under the same name warns again, distances rounded."""
    register_law('N', 'mine', lambda *a: 1.0, **DECL)
    with caplog.at_level(logging.WARNING, logger='fwl.calliope.core_split'):
        for _ in range(2):
            split_core_mantle({'H': 1e19}, M_METAL, M_SIL, dIW=-1.5)  # H dIW -1.33..-0.88
            split_core_mantle({'N': 1e19}, M_METAL, M_SIL, dIW=-6.0, laws={'N': 'mine'})
        register_law(
            'N', 'mine', lambda *a: 1.0, replace=True, **{**DECL, 'dIW_range': (-5.5, 0)}
        )
        split_core_mantle({'N': 1e19}, M_METAL, M_SIL, dIW=-6.0, laws={'N': 'mine'})
    msgs = [r.message for r in caplog.records if 'outside its range' in r.message]
    assert len(msgs) == 3
    assert "'dIW': 0.17}" in msgs[0] and '0.1699' not in msgs[0]
    assert "'dIW': 0.5}" in msgs[2]


@pytest.mark.physics_invariant
def test_s_law_takes_its_feo_from_diw_and_the_metal():
    """X_FeO of the S law is X'_Fe 10^(dIW/2): D_S rises with dIW, the melt oxides do not set it."""
    bulk = {'S': 1e15}  # trace S, so the metal stays Fe with the fixed O and Si
    fixed = {'O': 0.02, 'Si': 0.03}
    share = []
    for diw in (-3.0, -1.5, -0.5):
        core, mantle, _ = split_core_mantle(bulk, M_METAL, M_SIL, dIW=diw, w_fixed=fixed)
        w = {**fixed, 'S': core['S'] / M_METAL}
        x_feo = _x_fe(w) * 10 ** (diw / 2)
        n = {el: v / metal._mass(el) for el, v in w.items()}
        n_tot = sum(n.values()) + (1.0 - sum(w.values())) / metal._mass('Fe')
        x_o = n['O'] / n_tot  # x' and x agree without H and C
        log_d = -3.50 + 3000 / 4000 + 33 * 60 / 4000 + math.log10(x_feo) + 5.372
        log_d += 14 * math.log10(1 - x_o)
        d_mass = 10**log_d * n_tot / metal.silicate_units_per_kg()
        assert core['S'] / M_METAL == pytest.approx(d_mass * mantle['S'] / M_SIL, rel=1e-9)
        share.append(core['S'] / 1e15)
    assert share[0] < share[1] < share[2]
    # 1.5 dex of dIW is 0.75 dex of X_FeO: a D_S ratio of 10^0.75 between -3 and -1.5
    assert (share[1] / (1 - share[1])) / (share[0] / (1 - share[0])) == pytest.approx(
        10**0.75, rel=1e-9
    )


@pytest.mark.physics_invariant
def test_melt_oxides_and_nbo_t_reach_the_laws():
    """NBO/T moves D_C by 10^(0.118 dNBO/T); twice the melt FeO moves D_S only through the units."""
    base = split_core_mantle({'C': 1e15, 'S': 1e15}, M_METAL, M_SIL, dIW=-1.5)
    low = split_core_mantle({'C': 1e15, 'S': 1e15}, M_METAL, M_SIL, dIW=-1.5, nbo_t=0.0)
    ratio = (low[0]['C'] / low[1]['C']) / (base[0]['C'] / base[1]['C'])
    assert ratio == pytest.approx(10 ** (0.118 * 2.6), rel=1e-9)
    ox = dict(metal.KLB1, FeO=16.2, MgO=39.22 - 8.1)
    rich = split_core_mantle({'S': 1e15}, M_METAL, M_SIL, dIW=-1.5, oxides=ox)
    s_ratio = (rich[0]['S'] / rich[1]['S']) / (base[0]['S'] / base[1]['S'])
    units = metal.silicate_units_per_kg() / metal.silicate_units_per_kg(ox)
    assert s_ratio == pytest.approx(units, rel=1e-9)
    assert abs(s_ratio - 1.0) > 1e-3  # the units differ, so the check can fail


@pytest.mark.physics_invariant
def test_metal_beyond_fes_is_flagged_not_clamped():
    """A bulk S that makes the metal richer in S than FeS sets beyond_FeS; an Earth-like one does not."""
    _, _, earth = split_core_mantle(BULK, M_METAL, M_SIL, dIW=-1.5, w_fixed={'O': 0.02})
    assert earth['metal']['beyond_FeS'] == 0
    assert 0.0 < earth['metal']['x_light'] < 0.2
    core, mantle, flags = split_core_mantle({'S': 1.5e22}, 1e22, M_SIL, dIW=-1.5)
    assert flags['metal']['beyond_FeS'] == 1
    assert flags['metal']['x_S'] > 0.5
    assert abs(core['S'] + mantle['S'] - 1.5e22) <= 1e-12 * 1.5e22


@pytest.mark.physics_invariant
@settings(derandomize=True, max_examples=60, deadline=None)
@given(
    frac=st.floats(0.0, 0.05),
    diw=st.floats(-3.0, -0.5),
    ratio=st.floats(0.05, 2.0),
)
def test_ledger_closes_with_the_default_laws(frac, diw, ratio):
    """Ledger with metal feedback: the default H, C, N, S laws, dIW -3 to -0.5, M_metal/M_silicate 0.05 to 2.

    Each element is frac M_metal (frac up to 0.05), so the eps terms and mole bases feed
    back on D through the metal composition.
    """
    m_metal = ratio * M_SIL
    bulk = {el: frac * m_metal for el in ('H', 'C', 'N', 'S')}
    core, mantle, flags = split_core_mantle(bulk, m_metal, M_SIL, dIW=diw, w_fixed={'O': 0.02})
    w = {'O': 0.02, **{el: core[el] / m_metal for el in bulk}}
    w_melt = {el: mantle[el] / M_SIL for el in ('H', 'C')}
    for el, m in bulk.items():
        assert abs(core[el] + mantle[el] - m) <= 1e-12 * m + 2 * math.ulp(0.0)
        if m > 0.0:  # the D of the returned metal reproduces the split
            d = mass_D(get_law(el, DEFAULT_LAWS[el]), 60.0, 4000.0, diw, w, w_melt)
            assert core[el] / m_metal == pytest.approx(d * mantle[el] / M_SIL, rel=1e-8)


def test_diw_above_the_ideal_mixing_limit_raises_by_name():
    """X_FeO = x'_Fe 10^(dIW/2) cannot pass 1: dIW above -2 log10 x'_Fe raises a ValueError on dIW."""
    for diw in (0.3, 700.0):
        with pytest.raises(ValueError, match='dIW.*limit') as info:
            split_core_mantle({'S': 1e15}, M_METAL, M_SIL, dIW=diw, w_fixed={'O': 0.02})
        assert type(info.value) is ValueError
    core, _, _ = split_core_mantle({'S': 1e15}, M_METAL, M_SIL, dIW=0.0, w_fixed={'O': 0.02})
    assert core['S'] > 0.99e15


def test_callable_in_the_split_gets_the_melt_as_given():
    """A callable law in the split receives the melt cation fractions of oxides, summing to 1."""
    seen = []
    register_law('S', 'spy', lambda *a: seen.append(a) or 1.0, **DECL)
    split_core_mantle({'S': 1e15}, M_METAL, M_SIL, dIW=-1.5, laws={'S': 'spy'})
    x_sil = seen[-1][4]
    assert sum(x_sil.values()) == pytest.approx(1.0, rel=1e-12)
    assert x_sil['FeO'] == pytest.approx(metal.cation_fraction(metal.KLB1, 'FeO'), rel=1e-12)


@pytest.mark.physics_invariant
def test_beyond_fes_compares_s_with_fe_also_with_o_and_si():
    """With O and Si in the metal, S/Fe > 1 sets beyond_FeS although x_S < 0.5; S/Fe < 1 does not."""
    fixed = {'O': 0.05, 'Si': 0.05}
    out = []
    for bulk_s in (3e22, 5e22):  # S/Fe 0.65 and 1.26 at dIW -1
        core, _, flags = split_core_mantle({'S': bulk_s}, 1e22, M_SIL, dIW=-1.0, w_fixed=fixed)
        w_s = core['S'] / 1e22
        n = {'S': w_s / 32.06, 'O': 0.05 / 15.999, 'Si': 0.05 / 28.085}
        n_fe = (1 - w_s - 0.1) / 55.845
        assert flags['metal']['x_S'] == pytest.approx(
            n['S'] / (n_fe + sum(n.values())), rel=1e-3
        )
        assert flags['metal']['S_over_Fe'] == pytest.approx(n['S'] / n_fe, rel=1e-3)
        out.append(flags['metal'])
    assert out[0]['beyond_FeS'] == 0 and out[0]['S_over_Fe'] < 1.0
    assert out[1]['beyond_FeS'] == 1 and out[1]['x_S'] < 0.5 < out[1]['S_over_Fe']


def _limit(w):
    """-2 log10 x'_Fe of a metal of light-element mass fractions w (no H, no C)."""
    n = sum(v / metal._mass(el) for el, v in w.items())
    n_fe = (1.0 - sum(w.values())) / metal._mass('Fe')
    return -2.0 * math.log10(n_fe / (n_fe + n))


def test_diw_limit_is_that_of_the_start_metal():
    """The first pass has no S in the metal, so dIW may reach -2 log10 x'_Fe of w_fixed, not more."""
    lim = _limit({'O': 0.02})
    assert 0.0 < lim < 0.2
    core, _, _ = split_core_mantle(
        {'S': 1e15}, M_METAL, M_SIL, dIW=lim - 1e-9, w_fixed={'O': 0.02}
    )
    assert core['S'] > 0.99e15
    with pytest.raises(ValueError, match='ideal-mixing limit'):
        split_core_mantle({'S': 1e15}, M_METAL, M_SIL, dIW=lim + 1e-9, w_fixed={'O': 0.02})


def test_very_reduced_diw_and_a_metal_with_iron_at_round_off_raise_by_name():
    """dIW -700 underflows X_FeO: ValueError naming dIW; Fe at round-off gives finite flags."""
    with pytest.raises(ValueError, match='dIW -700.*underflow'):
        split_core_mantle({'S': 1e15}, M_METAL, M_SIL, dIW=-700.0)
    fixed = {'O': 0.5, 'Si': 0.5 - 1e-15}
    core, mantle, flags = split_core_mantle({'N': 0.0}, M_METAL, M_SIL, dIW=-2.0, w_fixed=fixed)
    assert math.isfinite(flags['metal']['S_over_Fe'])
    assert core['N'] == 0.0 and mantle['N'] == 0.0


@pytest.mark.physics_invariant
def test_s_share_of_the_tutorial_planet_against_diw():
    """The S core shares quoted on the Explanation page: 0.871, 0.974, 0.992 at dIW -3, -1.5, -0.5."""
    bulk = {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}
    share = []
    for diw in (-3.0, -1.5, -0.5):
        core, _, _ = split_core_mantle(
            bulk, M_METAL, M_SIL, dIW=diw, w_fixed={'O': 0.02, 'Si': 0.03}
        )
        share.append(core['S'] / bulk['S'])
    assert share == pytest.approx([0.8714, 0.9744, 0.9918], abs=5e-5)
    assert share[0] < share[1] < share[2]


@pytest.mark.parametrize(
    'w_fixed, word',
    [({'O': 1.0}, 'must be below 1'), ({'O': 0.5, 'Si': 0.5}, 'leaves no iron')],
)
def test_fixed_metal_messages_name_the_rule(w_fixed, word):
    """An element at 1 names the per-element rule; a sum at 1 names the missing iron."""
    with pytest.raises(ValueError, match=word):
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-2.0, w_fixed=w_fixed)
    assert core_split._WARNED == set()


def test_unknown_law_name_is_named_with_its_element():
    """laws={'H': 'nope'} raises PartitionTableError that names laws['H'] and the name."""
    with pytest.raises(PTE, match=r"laws\['H'\] = 'nope'") as info:
        split_core_mantle({'H': 1e15}, M_METAL, M_SIL, dIW=-1.5, laws={'H': 'nope'})
    assert type(info.value) is PTE


def test_subnormal_bulk_closes_within_two_ulp():
    """D 1 and M_metal = M_silicate halve a bulk of 5e-324 kg: both parts round to 0, within 2 ulp."""
    core, mantle, _ = split_core_mantle(
        {'N': 5e-324}, 1e20, 1e20, dIW=-2.0, laws={'N': _const('one', 0.0)}
    )
    assert core['N'] + mantle['N'] < 5e-324
    assert abs(core['N'] + mantle['N'] - 5e-324) <= 2 * math.ulp(0.0)


def test_iron_moles_come_from_the_iron_mass():
    """With 1e-13 of Fe by mass the S/Fe flag uses the Fe moles of the Fe mass, not a difference."""
    register_law('S', 'd1', lambda *a: 0.0, **DECL)
    fixed = {'O': 1 - 1e-12}
    core, _, flags = split_core_mantle(
        {'S': 1e9}, 1e22, 1e22, dIW=-2.0, laws={'S': 'd1'}, w_fixed=fixed
    )
    w_s = core['S'] / 1e22
    n_fe = (1.0 - (fixed['O'] + w_s)) / metal._mass('Fe')  # the order of the code's sum
    assert flags['metal']['S_over_Fe'] == pytest.approx(
        w_s / metal._mass('S') / n_fe, rel=1e-12
    )
    assert flags['metal']['beyond_FeS'] == 0  # S/Fe about 0.09


def test_unknown_range_is_flagged_out_of_range():
    """A callable law with an empty T range is always out of range: flag 1 and distance NaN."""
    register_law('N', 'no_t', lambda *a: 1.0, **{**DECL, 'T_range_K': ()})
    _, _, flags = split_core_mantle({'N': 1e19}, M_METAL, M_SIL, dIW=-2.0, laws={'N': 'no_t'})
    assert flags['N']['out_of_range'] == 1
    assert math.isnan(flags['N']['T_K'])


def test_ledger_undershoot_raises(monkeypatch):
    """A mantle 1e-9 short of the bulk breaks the ledger as an excess does."""
    real = core_split._masses

    def short(bulk, d, m_metal, m_sil):
        core, mantle = real(bulk, d, m_metal, m_sil)
        return core, {el: (1 - 1e-9) * v for el, v in mantle.items()}

    monkeypatch.setattr(core_split, '_masses', short)
    with pytest.raises(RuntimeError, match='S ledger'):
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-1.5)
    assert core_split._WARNED == set()


SPLIT_REFERENCE = Path(__file__).parent / 'data' / 'split_reference.json'


def _bisect_s_root(bulk, m_metal, m_sil, diw):
    """Root of w = D(w) w_sil(w) M/M_metal for an S-only metal, by bisection on (0, 1)."""
    law = get_law('S', 'suer2017_eq7')

    def excess(w):
        d = mass_D(law, 60.0, 4000.0, diw, {'S': w}, {})
        return d * m_metal / (m_sil + d * m_metal) * bulk / m_metal - w

    lo, hi = 0.0, 1.0 - 1e-12
    assert excess(lo) > 0.0 > excess(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if excess(mid) > 0.0 else (lo, mid)
    return 0.5 * (lo + hi)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'bulk, m_metal, m_sil',
    [(2e22, 1e22, 4e24), (3e22, 1e22, 4e24), (1e22, 1e22, 1e22)],
    ids=['s2e22', 's3e22', 'small-silicate'],
)
def test_s_rich_metal_poor_split_finds_the_root(bulk, m_metal, m_sil):
    """Where plain substitution oscillates, the split returns the bisection root past FeS."""
    root = _bisect_s_root(bulk, m_metal, m_sil, -1.5)
    core, mantle, flags = split_core_mantle({'S': bulk}, m_metal, m_sil, dIW=-1.5)
    assert core['S'] / m_metal == pytest.approx(root, rel=1e-9)
    assert abs(core['S'] + mantle['S'] - bulk) <= 1e-12 * bulk
    assert flags['metal']['beyond_FeS'] == 1


def test_no_metal_returns_the_bulk_at_any_diw():
    """M_metal 0 uses no D: an oxidised dIW past the FeO limit still gives core 0, mantle = bulk."""
    core, mantle, _ = split_core_mantle({'S': 1e20}, 0.0, M_SIL, dIW=0.5)
    assert core == {'S': 0.0}
    assert mantle == {'S': 1e20}


def test_converging_cases_equal_the_stored_split():
    """32 splits that converge by plain substitution (tutorial, trace S, a grid) give the stored outputs.

    Equal to the bit on the platform that stored them; elsewhere rel 1e-12 (libm round-off).
    """
    stored = json.loads(SPLIT_REFERENCE.read_text())
    meta = stored.pop('meta')
    same = (platform.system(), platform.machine()) == (meta['platform'], meta['machine'])
    assert len(stored) == 32  # every stored case is compared
    for name, ref in stored.items():
        c = ref['case']
        core, mantle, flags = split_core_mantle(
            c['bulk'],
            c['m_metal'],
            c['m_sil'],
            dIW=c['diw'],
            P_eq=c.get('p', 60.0),
            T_eq=c.get('t', 4000.0),
            w_fixed=c['w_fixed'],
        )
        for part, got in (('core', core), ('mantle', mantle), ('metal', flags['metal'])):
            want = ref[part]
            assert got == (want if same else pytest.approx(want, rel=1e-12, abs=0.0)), name


@pytest.mark.physics_invariant
def test_s_law_ranges_and_their_flags():
    """Suer 2017: 1-91 GPa, 1873-4100 K, dIW -1.7 to -1.0; dIW -3 and -0.5 lie 1.3 and 0.5 outside."""
    s = get_law('S', 'suer2017_eq7')
    assert (s['P_range_GPa'], s['T_range_K'], s['dIW_range']) == (
        (1.0, 91.0),
        (1873.0, 4100.0),
        (-1.7, -1.0),
    )
    for diw, dist in ((-3.0, 1.3), (-0.5, 0.5)):
        _, _, flags = split_core_mantle(
            {'S': 1e15}, M_METAL, M_SIL, dIW=diw, w_fixed={'O': 0.02}
        )
        assert flags['S']['out_of_range'] == 1
        assert flags['S']['dIW'] == pytest.approx(dist, abs=1e-12)


def _scripted(g_of_w):
    """An evaluate of _relaxed_fixed_point whose metal S is g_of_w(w_S); records each w_S."""
    seen = []

    def evaluate(w, w_melt):
        seen.append(w['S'])
        g = {'S': g_of_w(w['S'])}
        return {'S': g['S']}, {'S': 1.0 - g['S']}, g, {}

    return evaluate, seen


def test_relaxed_passes_damp_a_rising_map():
    """g = 0.5 w + 0.1 rises with s = 0.5: the second move is 1/(1 + 0.5) of its step, 0.1 -> 0.1333."""
    evaluate, seen = _scripted(lambda w: 0.5 * w + 0.1)
    core, _, w = core_split._relaxed_fixed_point(evaluate, {'S': 0.0}, {})
    assert seen[:3] == pytest.approx([0.0, 0.1, 0.1 + 0.05 / 1.5], rel=1e-12)  # full step: 0.15
    assert w['S'] == pytest.approx(0.2, rel=1e-11)


def test_relaxed_step_takes_at_most_half_of_the_iron_left():
    """g = 10 asks for 10 of S: the first move halves until it is at most 0.5, so w_S = 0.3125."""
    evaluate, seen = _scripted(lambda w: 0.4 if w > 0.0 else 10.0)
    core_split._relaxed_fixed_point(evaluate, {'S': 0.0}, {})
    assert seen[1] == pytest.approx(
        0.3125, rel=1e-15
    )  # 10 halved 5 times; 0.625 for a cap of 1
    assert seen[:1] == [0.0]


@pytest.mark.parametrize(
    'bulk, m_metal, m_sil, diw, w_fixed',
    [
        ({'C': 1.5e22}, 1e22, 1e18, -1.5, None),
        ({'N': 0.959e22}, 1e22, 1.26e16, -5.0, {'O': 0.02, 'Si': 0.03}),
    ],
    ids=['carbon-fills-the-metal', 'nitrogen-beyond-the-metal'],
)
def test_bulk_whose_iron_runs_out_raises_by_name(bulk, m_metal, m_sil, diw, w_fixed):
    """A bulk whose passes leave the metal without iron raises the named ValueError, not another error."""
    with pytest.raises(ValueError, match='no iron left') as info:
        split_core_mantle(bulk, m_metal, m_sil, dIW=diw, w_fixed=w_fixed)
    assert type(info.value) is ValueError


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'bulk, m_metal, m_sil, diw',
    [(1.5e22, 1e22, 1e18, -1.5), (1e23, 1e22, 1e21, -0.5)],
    ids=['iron-1e-6', 'iron-2e-5'],
)
def test_near_pure_s_root_on_the_float_grid(bulk, m_metal, m_sil, diw):
    """Roots with 1e-6 to 2e-5 of iron, where g falls by about 1e6 per unit w, match bisection."""
    root = _bisect_s_root(bulk, m_metal, m_sil, diw)
    core, mantle, flags = split_core_mantle({'S': bulk}, m_metal, m_sil, dIW=diw)
    assert core['S'] / m_metal == pytest.approx(root, rel=1e-9)
    assert abs(core['S'] + mantle['S'] - bulk) <= 1e-12 * bulk
    assert flags['metal']['beyond_FeS'] == 1


@pytest.mark.physics_invariant
def test_coupled_s_rich_split_reaches_the_damped_reference():
    """S, C, H and N at 2.5 % iron (dIW -3, O 0.1): the root of a damped plain iteration (lam 0.05).

    Reference from ~1000 passes of w <- w + 0.05 (g(w) - w) to a relative step below 1e-13.
    """
    bulk = {'S': 8e21, 'C': 2e21, 'H': 2e21, 'N': 1e21}
    ref = {
        'S': 0.413987228997558,
        'C': 0.16177301812214018,
        'H': 0.1997038436197129,
        'N': 0.09970598652514251,
    }
    core, mantle, _ = split_core_mantle(bulk, 1e22, 1e20, dIW=-3.0, w_fixed={'O': 0.1})
    assert {el: core[el] / 1e22 for el in bulk} == pytest.approx(ref, rel=1e-9)
    assert all(abs(core[el] + mantle[el] - m) <= 1e-12 * m for el, m in bulk.items())


def test_steep_step_above_the_grid_tolerance_is_not_accepted():
    """A move within 2 ulp does not end the passes while the step is 2e-6 of g: the next pass settles on w."""
    w1 = 0.5 + math.ulp(0.5)
    script = iter([w1, w1 - 1e-6])

    def evaluate(w, w_melt):
        g_s = next(script, w['S'])
        return {}, {}, {'A': 0.2, 'S': g_s}, {}

    _, _, w, _ = core_split._relaxed_passes(evaluate, {'A': 0.1, 'S': 0.5}, {}, False)
    assert w['A'] == pytest.approx(0.2, rel=1e-15)
    assert abs(w['S'] - w1) <= 2 * math.ulp(w1)


def test_relaxed_passes_wait_for_the_melt_basis_to_settle():
    """g settled on the start, the melt H new on pass 1 and moved on pass 2: pass 3 must confirm it."""
    calls = []
    melts = iter([{'H': 0.1}, {'H': 0.2}])

    def evaluate(w, w_melt):
        calls.append(dict(w_melt))
        return {}, {}, dict(w), next(melts, {'H': 0.2})

    _, _, w, _ = core_split._relaxed_passes(evaluate, {'S': 0.3}, {}, False)
    assert calls == [{}, {'H': 0.1}, {'H': 0.2}]
    assert w == {'S': 0.3}


def test_common_step_settles_a_split_that_the_per_element_passes_do_not():
    """Sweep case 1870: the per-element passes run first and do not settle; the smallest common lam does."""
    calls = []
    real = core_split.mass_D

    def count(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    m = 1.4793984719876185e24
    bulk = {
        'H': 0.0013584026541234068 * m,
        'C': 0.004212317293360554 * m,
        'S': 1.4239059823940543 * m,
    }
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(core_split, 'mass_D', count)
        core, mantle, _ = split_core_mantle(
            bulk, m, 7.608546877446469e18, dIW=-4.682711923437374, w_fixed={'O': 0.1}
        )
    assert len(calls) > len(bulk) * (core_split.MAX_PASSES + core_split.MAX_RELAXED_PASSES)
    assert core['S'] == pytest.approx(1.264557613090235e24, rel=1e-9)
    assert all(core[el] + mantle[el] == pytest.approx(bulk[el], rel=1e-12) for el in bulk)


def test_metal_with_no_iron_outside_h_by_round_off_raises_by_name():
    """w_H = 1 - 2**-53 leaves exactly 0 mol outside H: the named ValueError, not a division by 0."""
    with pytest.raises(ValueError, match='no iron left outside H and C'):
        mass_D(get_law('H', 'tagawa2021_eq2'), 60.0, 4000.0, -2.0, {'H': 1 - 2**-53}, {})


@pytest.mark.parametrize(
    'bulk, m_metal, m_sil, diw, w_fixed',
    [
        (
            {'H': 8.234751575583354e23, 'N': 8.141738774031332e22},
            3.0594828320755087e24,
            2.436509988289884e16,
            -0.8658248542625239,
            {'O': 0.1},
        ),
        (
            {'H': 9.153087750914148e23, 'N': 4.6589689694365475e19, 'S': 1.2152201420659063e20},
            2.3336975758573758e24,
            2.129139491150588e16,
            -2.3214007087220727,
            {'O': 0.02},
        ),
    ],
    ids=['sweep-1532-H', 'sweep-1721-S'],
)
def test_plain_split_with_a_small_silicate_settles_the_melt_basis(
    bulk, m_metal, m_sil, diw, w_fixed
):
    """The melt H sets the silicate units of D: at the returned split D reproduces the core of every element."""
    core, mantle, _ = split_core_mantle(bulk, m_metal, m_sil, dIW=diw, w_fixed=w_fixed)
    w = {**w_fixed, **{el: core[el] / m_metal for el in bulk}}
    w_melt = {el: mantle[el] / m_sil for el in ('H', 'C') if el in bulk}
    for el in bulk:
        d = mass_D(get_law(el, DEFAULT_LAWS[el]), 60.0, 4000.0, diw, w, w_melt)
        assert core[el] / m_metal == pytest.approx(d * mantle[el] / m_sil, rel=1e-9), el
    assert all(core[el] + mantle[el] == pytest.approx(bulk[el], rel=1e-12) for el in bulk)


@pytest.mark.parametrize(
    'bulk, m_metal, m_sil, diw, w_fixed',
    [
        (
            {'S': 1.0152954463428923e23},
            5.53739704857299e22,
            1.5858640020848288e17,
            -1.4382133911099835,
            {'O': 0.0},
        ),
        (
            {'C': 9.68644358404661e22, 'S': 4.157184340766592e17},
            4.485866476890163e22,
            4.377813186132594e18,
            -0.661571639511006,
            {'O': 0.0},
        ),
        (
            {'N': 3.9377888855367536e21, 'S': 6.698036413445549e21},
            4.933646958854052e21,
            1.5066317385381114e17,
            -0.6260226487946895,
            {'O': 0.1},
        ),
    ],
    ids=[
        'sweep-359-root-at-fe-1.5e-8',
        'sweep-2268-root-at-fe-1.8e-5',
        'sweep-849-no-root-on-a-grid',
    ],
)
def test_iron_poor_split_that_does_not_settle_raises_runtime_error(
    bulk, m_metal, m_sil, diw, w_fixed
):
    """S- or C-rich, iron-poor splits whose relaxed passes do not settle raise RuntimeError by name.

    359 and 2268 have a root with iron (Fe 1.5e-8 and 1.8e-5 by mass, by bisection) that the relaxed
    passes do not settle on; 849 has no root on a grid scan, which does not prove that none exists.
    """
    calls = []
    real = core_split.mass_D

    def count(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(core_split, 'mass_D', count)
        with pytest.raises(RuntimeError, match='did not settle in 2 x 1000 relaxed passes'):
            split_core_mantle(bulk, m_metal, m_sil, dIW=diw, w_fixed=w_fixed)
    assert len(calls) >= len(bulk) * 2 * core_split.MAX_RELAXED_PASSES


def test_plain_passes_check_the_melt_basis_of_the_first_pass():
    """g settled on the start but the first melt H is new: the plain passes take a second pass to check it."""
    calls = []

    def evaluate(w, w_melt):
        calls.append(dict(w_melt))
        return {}, {}, dict(w), {'H': 0.1}

    _, _, w = core_split._plain_fixed_point(evaluate, {'S': 0.3}, {'O': 0.02})
    assert calls == [{}, {'H': 0.1}]
    assert w == {'O': 0.02, 'S': 0.3}


def test_grid_rule_accepts_a_move_of_up_to_2_ulp():
    """A steep falling g after a 2-ulp step, root 1 ulp below w: pass 2 brackets it and returns the metal at w."""
    w0 = 0.5
    w1 = math.nextafter(math.nextafter(w0, 1.0), 1.0)
    root = math.nextafter(w1, 0.0)
    k = 1e-12 / (w1 - root)
    calls = []

    def evaluate(w, w_melt):
        calls.append(w['S'])
        g_s = w1 if len(calls) == 1 else w['S'] - k * (w['S'] - root)
        return {'A': 0.2, 'S': g_s}, {'A': 0.8, 'S': 1.0 - g_s}, {'A': 0.2, 'S': g_s}, {}

    core, mantle, w, _ = core_split._relaxed_passes(evaluate, {'A': 0.1, 'S': w0}, {}, False)
    step = 4 * math.ulp(w1)
    assert calls == [w0, w1, w1 - step, w1 + step]
    assert w['S'] == w1
    assert core['S'] == pytest.approx(w1, rel=1e-15, abs=0.0)  # M_metal 1 in this map
    assert core['S'] + mantle['S'] == pytest.approx(1.0, rel=1e-15, abs=0.0)


def test_grid_rule_of_the_common_step_judges_each_element_at_its_own_lam():
    """A steep A sets a tiny common lam; B (slope -0.5) is 4.5e-10 off its root and must keep moving."""
    a0 = 0.5
    a1 = math.nextafter(math.nextafter(a0, 1.0), 1.0)
    script_a = iter([a1, a1 - 5e-10])

    calls = []

    def evaluate(w, w_melt):
        calls.append(1)
        return {}, {}, {'A': next(script_a, w['A']), 'B': 0.2 - 0.5 * (w['B'] - 0.2)}, {}

    _, _, w, _ = core_split._relaxed_passes(
        evaluate, {'A': a0, 'B': 0.2 * (1 + 6e-10)}, {}, True
    )
    assert w['B'] == pytest.approx(0.2, rel=1e-13, abs=0.0)
    assert w['A'] == pytest.approx(a1, rel=1e-9)
    assert len(calls) == 6  # judged at the common lam, B takes one more pass


def _d_residuals(bulk, m_metal, m_sil, diw, w_fixed, core, mantle):
    """|w - D(w, w_melt) w_sil| / w of each element at a returned split."""
    w = {**w_fixed, **{el: core[el] / m_metal for el in bulk}}
    w_melt = {el: mantle[el] / m_sil for el in ('H', 'C') if el in bulk}
    return {
        el: abs(
            w[el]
            - mass_D(get_law(el, DEFAULT_LAWS[el]), 60.0, 4000.0, diw, w, w_melt)
            * mantle[el]
            / m_sil
        )
        / w[el]
        for el in bulk
    }


@pytest.mark.parametrize(
    'bulk, m_metal, m_sil, diw, ref',
    [
        (
            {'H': 1.9363249777556582e23, 'C': 8.163503547150277e22},
            1.84443083841283e23,
            1.7931639051646671e18,
            -3.6689582794162297,
            {'H': 0.4316511365896101, 'C': 0.44260285373610686},
        ),
        (
            {
                'S': 1.6351025151692482e24,
                'H': 4.026924495995128e22,
                'C': 2.0258800008695458e24,
                'N': 3.113694670982211e23,
            },
            2.6777700956688316e24,
            4.869342993603356e19,
            -0.9403058010523511,
            {
                'S': 0.5580179768265563,
                'H': 0.014249695242000643,
                'C': 0.2417249413734452,
                'N': 0.11627933728438677,
            },
        ),
    ],
    ids=['G1-iron-floor-before-the-common-step', 'G3-root-with-5-percent-iron'],
)
def test_relaxed_split_reaches_the_damped_root_with_iron(bulk, m_metal, m_sil, diw, ref):
    """Iron-poor H-C and S-C-H-N splits: the root that a damped iteration (lam 0.005) reaches, with D settled.

    G3 has a second root with iron 6e-7; the split prefers the one with more iron.
    """
    core, mantle, flags = split_core_mantle(bulk, m_metal, m_sil, dIW=diw, w_fixed={'O': 0.02})
    assert {el: core[el] / m_metal for el in bulk} == pytest.approx(ref, rel=1e-9)
    res = _d_residuals(bulk, m_metal, m_sil, diw, {'O': 0.02}, core, mantle)
    assert max(res.values()) <= 1e-9
    assert {el: flags[el]['D_residual'] for el in bulk} == pytest.approx(
        res, rel=1e-6, abs=1e-15
    )


@pytest.mark.parametrize(
    'iron, raises', [(5e-9, True), (2e-8, False)], ids=['iron-5e-9-raises', 'iron-2e-8-returns']
)
def test_iron_floor_is_1e_8_of_the_metal(iron, raises):
    """A map that keeps the metal at iron 5e-9 runs out of iron; one at iron 2e-8 settles there."""

    def evaluate(w, w_melt):
        return {}, {}, {'S': 1.0 - iron}, {}

    if raises:
        with pytest.raises(ValueError, match='below 1e-08 by mass'):
            core_split._relaxed_passes(evaluate, {'S': 0.0}, {}, False)
    else:
        _, _, w, _ = core_split._relaxed_passes(evaluate, {'S': 0.0}, {}, False)
        assert w['S'] == 1.0 - iron
    assert core_split.IRON_FLOOR == pytest.approx(1e-8, rel=0.0)


def test_plain_passes_check_a_carbon_only_melt_basis():
    """The melt C alone moves on pass 2: the plain passes take pass 3 to confirm it."""
    calls = []
    melts = iter([{'C': 0.1}, {'C': 0.2}])

    def evaluate(w, w_melt):
        calls.append(dict(w_melt))
        return {}, {}, dict(w), next(melts, {'C': 0.2})

    _, _, w = core_split._plain_fixed_point(evaluate, {'S': 0.3}, {})
    assert calls == [{}, {'C': 0.1}, {'C': 0.2}]
    assert w == {'S': 0.3}


@pytest.mark.parametrize('ulps', [3, 4])
def test_grid_rule_refuses_a_move_of_more_than_2_ulp(ulps):
    """After a 3- or 4-ulp secant step the move is above 2 ulp: no accept at w, the passes end within 4 ulp of the root."""
    w0 = 0.5
    w1 = w0
    for _ in range(ulps):
        w1 = math.nextafter(w1, 1.0)
    root = math.nextafter(w1, 0.0)
    k = 1e-12 / (w1 - root)
    calls = []

    def evaluate(w, w_melt):
        calls.append(w['S'])
        g_s = w1 if len(calls) == 1 else w['S'] - k * (w['S'] - root)
        return {'A': 0.2, 'S': g_s}, {'A': 0.8, 'S': 1.0 - g_s}, {'A': 0.2, 'S': g_s}, {}

    _, _, w, _ = core_split._relaxed_passes(evaluate, {'A': 0.1, 'S': w0}, {}, False)
    assert w['S'] != w1 and abs(w['S'] - root) <= 4 * math.ulp(root)
    assert calls[:2] == [w0, w1] and w1 + 4 * math.ulp(w1) not in calls


def test_iron_floor_can_stop_a_split_that_has_a_root_with_iron():
    """Small-silicate sweep 1078 (H-N): the passes drain the iron to the floor and raise ValueError.

    A damped iteration (lam 0.005) reaches a root with Fe 2.87e-2 by mass; the relaxed passes do not.
    """
    with pytest.raises(ValueError, match='a pass leaves iron below 1e-08') as info:
        split_core_mantle(
            {'H': 8.691548542853415e21, 'N': 1.2021745831093811e19},
            8.768672306066827e21,
            2.125733426050056e18,
            dIW=-3.5303822081080702,
            w_fixed={'O': 0.1},
        )
    assert type(info.value) is ValueError


@pytest.mark.parametrize(
    'bulk, m_metal, m_sil, diw, w_fixed',
    [
        (
            {'H': 1.026913291017055e16, 'N': 2.0745447651584106e17, 'S': 4.908478118313359e21},
            4.1883763504334535e21,
            7.758815031595311e17,
            -0.7860190639804525,
            {'O': 0.02},
        ),
        (
            {'N': 1.385645895332543e20, 'S': 6.831823027319011e24},
            2.9383194256916654e24,
            1.9975254165041137e18,
            -1.2177369244970537,
            {'O': 0.1},
        ),
    ],
    ids=['broad-1695', 'small-silicate-786'],
)
def test_common_passes_that_do_not_settle_after_an_iron_out_raise_the_iron_error(
    bulk, m_metal, m_sil, diw, w_fixed
):
    """The per-element passes run out of iron and the common passes do not settle: the iron ValueError stands."""
    with pytest.raises(ValueError, match='a pass leaves iron below 1e-08') as info:
        split_core_mantle(bulk, m_metal, m_sil, dIW=diw, w_fixed=w_fixed)
    assert type(info.value) is ValueError


def test_return_at_the_evaluated_metal_keeps_every_mantle_non_negative():
    """A grid accept with a mantle below core (w - g) / g is not returned at w: the passes go on to the root."""
    w0 = 0.5
    w1 = math.nextafter(math.nextafter(w0, 1.0), 1.0)
    root = math.nextafter(w1, 0.0)
    k = 1e-12 / (w1 - root)
    calls = []

    def evaluate(w, w_melt):
        calls.append(w['S'])
        g_s = w1 if len(calls) == 1 else w['S'] - k * (w['S'] - root)
        return {'A': 0.2, 'S': g_s}, {'A': 0.8, 'S': 1e-13}, {'A': 0.2, 'S': g_s}, {}

    core, mantle, w, _ = core_split._relaxed_passes(evaluate, {'A': 0.1, 'S': w0}, {}, False)
    assert min(mantle.values()) >= 0.0
    assert w['S'] == pytest.approx(root, rel=1e-15, abs=0.0)


def test_bracket_uses_the_melt_basis_of_the_evaluated_metal():
    """The root moves with the melt H: the bracket at w takes the melt that evaluated w, not the newer one."""
    w0 = 0.5
    w1 = math.nextafter(math.nextafter(w0, 1.0), 1.0)
    root0 = math.nextafter(w1, 0.0)
    k = 1e-12 / (w1 - root0)
    calls = []

    def evaluate(w, w_melt):
        calls.append(w['S'])
        root = root0 + 0.01 * (w_melt.get('H', 0.5) - 0.5)
        g_s = w1 if len(calls) == 1 else w['S'] - k * (w['S'] - root)
        melt = {'H': 0.5 if len(calls) == 1 else 0.5 + 4e-13}
        return {'A': 0.2, 'S': g_s}, {'A': 0.8, 'S': 1.0 - g_s}, {'A': 0.2, 'S': g_s}, melt

    _, _, w, _ = core_split._relaxed_passes(evaluate, {'A': 0.1, 'S': w0}, {}, False)
    step = 4 * math.ulp(w1)
    assert calls == [w0, w1, w1 - step, w1 + step]
    assert w['S'] == w1


def test_d_residual_is_taken_at_the_metal_of_the_returned_core():
    """S 1.5 x M_metal returns at an evaluated metal 1 ulp from core / M_metal: the flag uses core / M_metal."""
    got = {}
    real = core_split._relaxed_fixed_point

    def capture(*args):
        out = real(*args)
        got['w'] = out[2]
        return out

    m_metal, m_sil = 1.3e22, 1.3e18
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(core_split, '_relaxed_fixed_point', capture)
        core, mantle, flags = split_core_mantle({'S': 1.5 * m_metal}, m_metal, m_sil, dIW=-1.5)
    law = get_law('S', DEFAULT_LAWS['S'])

    def residual(w):
        x = core['S'] / m_metal
        return abs(x - mass_D(law, 60.0, 4000.0, -1.5, w, {}) * mantle['S'] / m_sil) / x

    at_core = residual({'S': core['S'] / m_metal})
    assert flags['S']['D_residual'] == pytest.approx(at_core, rel=1e-6)
    assert abs(residual(got['w']) - at_core) > 0.1 * at_core  # the 2 metals give different D


def _steep_stub(fail_at=(), zero_element=False):
    """The 2-ulp grid stub: root 1 ulp below w1; ValueError at the metals in ``fail_at``."""
    w0 = 0.5
    w1 = math.nextafter(math.nextafter(w0, 1.0), 1.0)
    root = math.nextafter(w1, 0.0)
    k = 1e-12 / (w1 - root)
    calls = []

    def evaluate(w, w_melt):
        calls.append(w['S'])
        if w['S'] in fail_at:
            raise ValueError('law outside its domain')
        g_s = w1 if len(calls) == 1 else w['S'] - k * (w['S'] - root)
        g = {'A': 0.2, 'S': g_s, **({'B': 0.0} if zero_element else {})}
        return (
            dict(g),
            {'A': 0.8, 'S': 1.0 - g_s, **({'B': 0.0} if zero_element else {})},
            g,
            {},
        )

    start = {'A': 0.1, 'S': w0, **({'B': 0.0} if zero_element else {})}
    return evaluate, start, calls, w1, root


@pytest.mark.parametrize('sides', ['both', 'upper'])
def test_bracket_end_that_raises_is_no_bracket(sides):
    """A ValueError at a bracket end reads as no sign change: no return at w, the passes go on to the root."""
    w1 = math.nextafter(math.nextafter(0.5, 1.0), 1.0)
    step = 4 * math.ulp(w1)
    fail_at = (w1 - step, w1 + step) if sides == 'both' else (w1 + step,)
    evaluate, start, calls, w1, root = _steep_stub(fail_at)
    _, _, w, _ = core_split._relaxed_passes(evaluate, start, {}, False)
    assert w['S'] != w1
    assert w['S'] == pytest.approx(root, rel=1e-15, abs=0.0)


def test_return_at_w_takes_a_zero_mantle():
    """An element with no mass has rest 0 at w: the grid accept returns at w all the same."""
    evaluate, start, calls, w1, root = _steep_stub(zero_element=True)
    core, mantle, w, _ = core_split._relaxed_passes(evaluate, start, {}, False)
    step = 4 * math.ulp(w1)
    assert calls == [0.5, w1, w1 - step, w1 + step]
    assert (core['B'], mantle['B'], w['S']) == (0.0, 0.0, w1)


_G3 = (
    {
        'S': 1.6351025151692482e24,
        'H': 4.026924495995128e22,
        'C': 2.0258800008695458e24,
        'N': 3.113694670982211e23,
    },
    2.6777700956688316e24,
    4.869342993603356e19,
    -0.9403058010523511,
)


@pytest.mark.parametrize('start', range(1, 6))
def test_iron_poor_split_with_2_roots_returns_the_root_with_more_iron(start):
    """The S-C-H-N split of G3, perturbed at round-off: the root with 5 % iron, not the one with 6e-7."""
    np = pytest.importorskip('numpy')
    eps = np.random.default_rng(3).uniform(-1, 1, (start + 1, 7))[start] * 4e-16
    bulk = {el: m * (1 + e) for (el, m), e in zip(_G3[0].items(), eps)}
    m_metal, m_sil, diw = (v * (1 + e) for v, e in zip(_G3[1:], eps[4:]))
    core, mantle, flags = split_core_mantle(bulk, m_metal, m_sil, dIW=diw, w_fixed={'O': 0.02})
    iron = 0.98 - sum(core.values()) / m_metal
    assert iron == pytest.approx(0.04972804927351793, rel=1e-9)
    assert (
        max(_d_residuals(bulk, m_metal, m_sil, diw, {'O': 0.02}, core, mantle).values()) <= 1e-9
    )
    assert flags['metal']['beyond_FeS'] == 1


def _relaxed_out(iron, at_w=True):
    """A relaxed return with metal S 1 - iron: beyond FeS (S/Fe > 1 in moles) for iron < 0.635."""
    return {'S': 0.3}, {'S': 0.7}, {'S': 1.0 - iron}, at_w


@pytest.mark.parametrize(
    'iron, at_w, damped, want',
    [
        (0.1, False, 0.2, 'out'),
        (0.637, True, 0.7, 'out'),
        (0.634, True, 0.7, 'damped'),
        (0.5, True, 0.6, 'damped'),
        (0.1, True, 0.2, 'damped'),
        (0.1, True, 0.1, 'out'),
        (0.1, True, 0.05, 'out'),
        (0.1, True, None, 'out'),
        (0.1, True, 'iron-out', 'out'),
        (0.1, True, 'law-error', 'out'),
    ],
    ids=[
        'settled-return',
        'below-FeS',
        'just-beyond-FeS',
        'S-over-Fe-1.7',
        'more-iron',
        'same-iron',
        'less-iron',
        'damped-unsettled',
        'damped-iron-out',
        'damped-law-error',
    ],
)
def test_more_iron_takes_the_damped_root_only_with_more_iron(iron, at_w, damped, want):
    """Damped passes run after a return at w beyond FeS; their root replaces it only with more iron."""
    calls = []

    def passes(*args):
        calls.append(args)
        if damped == 'iron-out':
            raise core_split._IronOut('no iron')
        if damped == 'law-error':
            raise ValueError('law outside its domain')
        return None if damped is None else ({'S': 0.1}, {'S': 0.9}, {'S': 1.0 - damped}, False)

    out = _relaxed_out(iron, at_w)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(core_split, '_relaxed_passes', passes)
        got = core_split._more_iron('ev', {'S': 0.0}, {'O': 0.02}, out)
    ran = at_w and iron < 0.635
    assert calls == ([('ev', {'S': 0.0}, {'O': 0.02}, True, 0.02, 20000)] if ran else [])
    assert got == (out[:3] if want == 'out' else ({'S': 0.1}, {'S': 0.9}, {'S': 1.0 - damped}))
    assert (core_split.DAMPED_LAM, core_split.MAX_DAMPED_PASSES) == (0.02, 20000)


@pytest.mark.parametrize('first', ['per-element', 'common'])
def test_relaxed_fixed_point_sends_each_return_to_more_iron(first):
    """A return at w beyond FeS from the per-element or the common passes goes to the damped passes."""

    def passes(evaluate, start, fixed, common, lam_all=None, n=None):
        if lam_all is not None:
            return {'S': 0.1}, {'S': 0.9}, {'S': 0.8}, False
        if first == 'per-element' or common:
            return _relaxed_out(0.01)
        return None

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(core_split, '_relaxed_passes', passes)
        got = core_split._relaxed_fixed_point('ev', {'S': 0.0}, {})
    assert got == ({'S': 0.1}, {'S': 0.9}, {'S': 0.8})


def test_relaxed_passes_take_a_common_lam_and_a_pass_limit():
    """lam_all moves every element by lam_all of its step; ``passes`` caps the passes."""
    calls = []

    def evaluate(w, w_melt):
        calls.append(dict(w))
        return {}, {}, {'S': 0.5, 'C': 0.25}, {}

    assert core_split._relaxed_passes(evaluate, {'S': 0.0, 'C': 0.0}, {}, True, 0.02, 2) is None
    assert calls == [{'S': 0.0, 'C': 0.0}, {'S': 0.01, 'C': 0.005}]


def test_relaxed_passes_mark_a_return_at_w():
    """The 4th item of a relaxed return is True at the evaluated metal, False at a settled g."""
    evaluate, start, calls, w1, root = _steep_stub()
    assert core_split._relaxed_passes(evaluate, start, {}, False)[3] is True

    def settled(w, w_melt):
        return {'S': 0.3}, {'S': 0.7}, {'S': 0.3}, {}

    assert core_split._relaxed_passes(settled, {'S': 0.0}, {}, False)[3] is False
