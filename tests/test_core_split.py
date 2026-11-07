"""Tests for `src/calliope/core_split.py`: the core and mantle split of a bulk inventory.

Contract clauses: at fixed mass-basis D the core holds D x / (1 + D x) of each element,
x = M_metal / M_silicate; the element ledger bulk = core + mantle holds to 1e-12 on every
split (invariant I3) with core and mantle >= 0; a law outside its range is evaluated
there (no clamp) and flagged with its distance per axis; O and Si are refused as bulk
elements; every input outside its domain raises ValueError naming it. Hypothesis tests
use ``derandomize=True``. Commands: docs/How-to/build_tests.md.
"""

from __future__ import annotations

import logging
import math

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
    alt, _, _ = split_core_mantle(
        BULK,
        M_METAL,
        M_SIL,
        dIW=-1.5,
        laws={'S': 'suer2023_line'},
        w_fixed={'O': 0.02, 'Si': 0.03},
    )
    assert abs(alt['S'] - core['S']) > 0.01 * BULK['S']
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
    assert flags['S'] == {'out_of_range': 1, 'P_GPa': 209.0, 'T_K': 1400.0, 'dIW': 0.0}
    d = mass_D(
        get_law('S', 'suer2017_eq7'), 300.0, 5500.0, -1.5, {'S': core['S'] / M_METAL}, {}
    )
    assert core['S'] / mantle['S'] == pytest.approx(d * M_METAL / M_SIL, rel=1e-10)
    # a clamp to (91 GPa, 4100 K) gives a different D
    clamped = mass_D(
        get_law('S', 'suer2017_eq7'), 91.0, 4100.0, -1.5, {'S': core['S'] / M_METAL}, {}
    )
    assert abs(math.log10(clamped / d)) > 0.1
    assert sum('outside its range' in r.message for r in caplog.records) == 1


def test_zero_metal_and_zero_bulk():
    """M_metal 0 leaves every element in the mantle; a bulk of 0 gives 0 in both."""
    core, mantle, _ = split_core_mantle(BULK, 0.0, M_SIL, dIW=-2.0)
    assert core == {el: 0.0 for el in BULK} and mantle == BULK
    core, mantle, _ = split_core_mantle({'S': 0.0}, M_METAL, M_SIL, dIW=-2.0)
    assert core == {'S': 0.0} and mantle == {'S': 0.0}


@pytest.mark.parametrize('el', ['O', 'Si'])
def test_exchange_elements_are_refused_as_bulk(el):
    """O and Si enter only through w_fixed: their bulk is rock, outside the volatile inventory."""
    with pytest.raises(ValueError, match=f'{el}.*w_fixed'):
        split_core_mantle({el: 1e20}, M_METAL, M_SIL, dIW=-2.0)


BAD = [True, 'x', None, math.nan, math.inf, -math.inf]


@pytest.mark.parametrize('value', BAD + [-1.0])
@pytest.mark.parametrize('field', ['bulk', 'M_metal', 'M_silicate', 'P_eq', 'T_eq'])
def test_bad_input_is_refused_by_name(field, value):
    """A boolean, string, None, NaN, inf or negative value raises ValueError naming the field."""
    kw = {'P_eq': 60.0, 'T_eq': 4000.0}
    args = [{'S': 1e20}, M_METAL, M_SIL]
    if field == 'bulk':
        args[0] = {'S': value}
    elif field in ('M_metal', 'M_silicate'):
        args[1 if field == 'M_metal' else 2] = value
    else:
        kw[field] = value
    with pytest.raises(ValueError, match=field if field != 'bulk' else 'S'):
        split_core_mantle(*args, dIW=-2.0, **kw)


@pytest.mark.parametrize(
    'field, value',
    [('M_silicate', 0.0), ('P_eq', 0.0), ('T_eq', 0.0), ('dIW', math.nan), ('dIW', True)],
)
def test_zero_silicate_pressure_temperature_and_bad_diw_are_refused(field, value):
    """M_silicate, P_eq and T_eq must be above 0; dIW must be a finite number."""
    kw = {'dIW': -2.0, 'M_silicate': M_SIL, field: value}
    with pytest.raises(ValueError, match=field):
        split_core_mantle({'S': 1e20}, M_METAL, kw.pop('M_silicate'), **kw)


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
    """A D that alternates with the metal composition never settles: RuntimeError after 50 passes."""
    calls = []

    def flip(*a, **k):
        calls.append(1)
        return 10.0 if len(calls) % 2 else 20.0

    monkeypatch.setattr(core_split, 'mass_D', flip)
    with pytest.raises(RuntimeError, match='50 passes'):
        split_core_mantle({'S': 1e20}, M_METAL, M_SIL, dIW=-2.0)
    assert len(calls) == 50


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
    """I3: bulk = core + mantle to 1e-12, core and mantle >= 0, for D 1e-6 to 1e6 and masses 1e15 to 1e24 kg.

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
