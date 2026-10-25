"""Tests for the metal reservoir in `src/calliope/solve.py`.

Contract: with metal_mass > 0 the H, C, N and S balances carry a core term
D_i M_metal w_i(melt). With the same mass-basis D for all four elements, no O or Si
in the metal and dIW a parameter, the solve equals the solve without metal and a
melt mass larger by D M_metal (an identity of the element balance); the gas
reservoir matters, so the gas-free core fraction Dx/(1 + Dx) misses the result.
Every element closes over gas, melt and core, and metal_mass 0 is the solve
without metal. The O of the metal comes from the silicate FeO and never enters the
volatile O budget: in the authoritative-O solve gas and melt O close on the budget
with metal O present, in 'fixed' and 'partition' mode alike.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from calliope.constants import volatile_species
from calliope.metal import load_partition_table
from calliope.solve import equilibrium_atmosphere

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

ELEMS = ('H', 'C', 'N', 'S')
P_GUESS = {'H2O': 50.0, 'CO2': 20.0, 'N2': 1.0, 'S2': 1.0}


def _ddict(m_mantle=4.03e24):
    d = {
        'M_mantle': m_mantle,
        'gravity': 9.81,
        'radius': 6.371e6,
        'Phi_global': 1.0,
        'T_magma': 1800.0,
        'fO2_shift_IW': 2.0,
    }
    for sp in volatile_species:
        d[f'{sp}_included'] = 1
        d[f'{sp}_initial_bar'] = 0.0
    return d


def _equal_d_table(d):
    """Placeholder table with a constant mass-basis D for H, C, N and S."""
    table = load_partition_table()
    for el in ELEMS:
        table[el] = dict(
            table[el], basis='mass', a=math.log10(d), b=0.0, c=0.0, d=0.0, e=0.0, eps={}
        )
    return table


def _metal(d, m_metal):
    return {
        'metal_mass': m_metal,
        'partition_table': _equal_d_table(d),
        'P_eq': 50.0,
        'T_eq': 3500.0,
        'dIW_eq': -2.0,
        'dIW_eq_mode': 'parameter',
        'O_mode_core': 'fixed',
        'w_O_core': 0.0,
        'Si_mode_core': 'fixed',
        'w_Si_core': 0.0,
    }


TARGET = {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}


def _solve(ddict, target=TARGET, p_guess=P_GUESS):
    """One seeded solve, so a restart draws the same pressures for both runs compared."""
    np.random.seed(0)
    return equilibrium_atmosphere(
        target,
        ddict,
        p_guess=p_guess,
        print_result=False,
        nguess=200,
        xtol=1e-12,
        rtol=1e-10,
        atol=1.0,
        opt_solver=False,
    )


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
@pytest.mark.parametrize('d', [1.0, 30.0])
def test_metal_with_equal_d_is_a_larger_melt(d):
    """Equal D for H, C, N, S: metal solve = no-metal solve with M_melt + D M_metal."""
    m_metal = 1.9e24
    with_metal = _solve({**_ddict(), **_metal(d, m_metal)})
    larger = _solve(_ddict(4.03e24 + d * m_metal))
    for sp in ('H2O', 'CO2', 'N2', 'S2', 'H2', 'CO', 'CH4'):
        assert with_metal[sp + '_bar'] == pytest.approx(larger[sp + '_bar'], rel=1e-6)
    for el in ELEMS:
        melt, core = with_metal[el + '_kg_liquid'], with_metal[el + '_kg_core']
        assert core == pytest.approx(d * m_metal * melt / 4.03e24, rel=1e-12)
        assert melt + core == pytest.approx(larger[el + '_kg_liquid'], rel=1e-6)
        closure = with_metal[el + '_kg_atm'] + melt + core
        assert closure == pytest.approx(TARGET[el], rel=1e-4)
    # the gas holds the rest: core = (budget - gas) Dx/(1 + Dx), short of the gas-free value
    x = m_metal / 4.03e24
    gas_free = TARGET['H'] * d * x / (1 + d * x)
    gas_share = with_metal['H_kg_atm'] / TARGET['H']
    assert with_metal['H_kg_core'] / gas_free == pytest.approx(1 - gas_share, abs=2e-5)
    assert gas_share > 0


def test_zero_metal_mass_is_the_solve_without_metal():
    """metal_mass 0 gives the solve without metal, bit for bit, and zero core masses."""
    plain = _solve(_ddict())
    zero = _solve({**_ddict(), **_metal(30.0, 0.0)})
    for sp in ('H2O', 'CO2', 'N2', 'S2'):
        assert zero[sp + '_bar'] == plain[sp + '_bar']
    assert all(zero[el + '_kg_core'] == 0 for el in ELEMS)


@pytest.mark.physics_invariant
def test_partition_mode_puts_exchange_o_and_si_in_the_metal():
    """O and Si 'partition': the converged metal holds the exchange x_O, x_Si of its dIW."""
    from calliope.metal import exchange_metal, molar_mass_rock

    ddict = {**_ddict(), **_metal(30.0, 1.9e24)}
    ddict.update(O_mode_core='partition', Si_mode_core='partition', T_eq=3500.0)
    out = _solve(ddict)
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    x = exchange_metal(ddict['partition_table'], 50.0, 3500.0, -2.0)
    masses = {'O': 0.015999, 'Si': molar_mass_rock['Si']}
    n_total = (1 - sum(w.values())) / molar_mass_rock['Fe'] + sum(
        w[el]
        / (masses.get(el) or {'H': 0.001008, 'C': 0.012011, 'N': 0.014007, 'S': 0.03206}[el])
        for el in w
    )
    for el in ('O', 'Si'):
        assert w[el] / masses[el] / n_total == pytest.approx(x[el], rel=1e-4)
    assert 0 < w['O'] < 0.2 and 0 < w['Si'] < 0.2
    assert out['O_kg_core'] == pytest.approx(w['O'] * 1.9e24, rel=1e-12)


@pytest.mark.physics_invariant
def test_diw_from_feo_matches_the_melt_and_metal():
    """dIW_eq_mode 'from_FeO': dIW = 2 log10(x_FeO / x_Fe) on the KLB-1 melt and the final metal."""
    from calliope.metal import KLB1, molar_mass_rock, oxide_mole_fraction

    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'dIW_eq_mode': 'from_FeO', 'w_Si_core': 0.05}
    out = _solve(ddict)
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    m = {'H': 0.001008, 'C': 0.012011, 'N': 0.014007, 'S': 0.03206, 'O': 0.015999}
    m['Si'] = molar_mass_rock['Si']
    n_fe = (1 - sum(w.values())) / molar_mass_rock['Fe']
    x_fe = n_fe / (n_fe + sum(w[el] / m[el] for el in w))
    expected = 2 * math.log10(oxide_mole_fraction(KLB1, 'FeO') / x_fe)
    assert out['dIW_eq'] == pytest.approx(expected, abs=1e-4)
    # KLB-1 with 8.1 wt% FeO sits near IW-2.3 against nearly pure Fe metal
    assert -2.6 < out['dIW_eq'] < -2.0


def test_unknown_metal_modes_raise():
    """An O or Si mode or a dIW mode outside the known ones is refused before solving."""
    for key, value in (('O_mode_core', 'free'), ('Si_mode_core', 1), ('dIW_eq_mode', 'buffer')):
        with pytest.raises(ValueError) as exc:
            _solve({**_ddict(), **_metal(30.0, 1.9e24), key: value})
        assert key in str(exc.value) or repr(value) in str(exc.value)


def _solve_auth(ddict, target, p_guess=None):
    """One seeded authoritative-O solve at the tight acceptance of `_solve`."""
    from calliope.solve import equilibrium_atmosphere_authoritative_O

    try:
        return equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=2.0,
            p_guess=p_guess or dict(P_GUESS, fO2_shift_IW=2.0),
            print_result=False,
            nguess=200,
            xtol=1e-12,
            rtol=1e-10,
            atol=1.0,
            opt_solver=False,
            random_seed=0,
        )
    except Exception as exc:
        pytest.fail(f'authoritative-O solve with metal raised {type(exc).__name__}: {exc}')


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_authoritative_o_budget_excludes_the_metal_o():
    """Authoritative-O with 0.1 wt% O in the metal returns the fO2 of the buffered solve.

    The O budget is the gas and melt O of the buffered metal solve at IW+2; the metal O is
    rock O from FeO, so the derived fO2 must be IW+2 and gas + melt O must equal the budget.
    """
    m_metal, w_o = 1.9e24, 1e-3
    buffered = _solve({**_ddict(), **_metal(30.0, m_metal)})
    target = dict(TARGET, O=buffered['O_kg_total'])
    out = _solve_auth({**_ddict(), **_metal(30.0, m_metal), 'w_O_core': w_o}, target)
    assert out['fO2_shift_derived'] == pytest.approx(2.0, abs=1e-6)
    for sp in ('H2O', 'CO2', 'N2', 'S2', 'H2', 'CO'):
        assert out[sp + '_bar'] == pytest.approx(buffered[sp + '_bar'], rel=1e-6)
    assert out['O_kg_core'] > 0
    assert out['O_kg_core'] == pytest.approx(w_o * m_metal, rel=1e-12)
    gas_melt_o = out['O_kg_atm'] + out['O_kg_liquid']
    assert abs(gas_melt_o - target['O']) <= target['O'] * 1e-10 + 10.0
    for el in ELEMS:
        closure = out[el + '_kg_atm'] + out[el + '_kg_liquid'] + out[el + '_kg_core']
        assert abs(closure - target[el]) <= target[el] * 1e-10 + 10.0, el
    # metal O (1.9e21 kg) is 22 times gas and melt O (8.75e19 kg): charged to the budget it
    # would move fO2 far from IW+2
    assert out['O_kg_core'] > 10 * buffered['O_kg_total']
    assert out['H_kg_core'] > 0.5 * TARGET['H']


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_authoritative_o_with_equal_d_is_a_larger_melt():
    """Authoritative-O, equal D = 30, no metal O: the solve of a melt larger by D M_metal.

    The larger melt also holds the O of its dissolved species, so its O budget is larger
    by D M_metal times the O mass fraction of the melt of the metal solve.
    """
    d, m_metal, m_melt = 30.0, 1.9e24, 4.03e24
    metal = {**_ddict(), **_metal(d, m_metal)}
    target = dict(TARGET, O=_solve(metal)['O_kg_total'])
    with_metal = _solve_auth(metal, target)
    o_extra = d * m_metal * with_metal['O_kg_liquid'] / m_melt
    larger = _solve_auth(_ddict(m_melt + d * m_metal), dict(target, O=target['O'] + o_extra))
    for sp in ('H2O', 'CO2', 'N2', 'S2', 'H2', 'CO', 'CH4'):
        assert with_metal[sp + '_bar'] == pytest.approx(larger[sp + '_bar'], rel=1e-6)
    assert with_metal['fO2_shift_derived'] == pytest.approx(
        larger['fO2_shift_derived'], abs=1e-6
    )
    for el in ELEMS:
        melt, core = with_metal[el + '_kg_liquid'], with_metal[el + '_kg_core']
        assert core == pytest.approx(d * m_metal * melt / m_melt, rel=1e-12)
        assert melt + core == pytest.approx(larger[el + '_kg_liquid'], rel=1e-6)
    assert with_metal['O_kg_core'] == 0
    # the core holds a large share of H, so the identity is not the trivial no-metal solve
    assert with_metal['H_kg_core'] > 0.5 * TARGET['H']


@pytest.mark.physics_invariant
def test_partition_o_far_above_the_o_budget_solves():
    """Partition-mode metal O (3.6 wt%) is 800 times the volatile O budget and still solves."""
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': load_partition_table()}
    ddict.update(O_mode_core='partition', Si_mode_core='partition')
    buffered = _solve(ddict)
    target = dict(TARGET, O=buffered['O_kg_total'])
    out = _solve_auth(ddict, target)
    assert out['fO2_shift_derived'] == pytest.approx(2.0, abs=1e-6)
    assert out['O_kg_core'] > 100 * target['O']
    assert out['O_kg_core'] == pytest.approx(buffered['O_kg_core'], rel=1e-6)
    assert abs(out['O_kg_atm'] + out['O_kg_liquid'] - target['O']) <= target['O'] * 1e-10 + 10


@pytest.mark.parametrize('authoritative', [False, True])
def test_metal_without_melt_is_refused_by_name(authoritative):
    """metal_mass > 0 with Phi_global 0 raises the named melt error before any solve."""
    from calliope.solve import equilibrium_atmosphere_authoritative_O

    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'Phi_global': 0.0}
    try:
        if authoritative:
            equilibrium_atmosphere_authoritative_O(
                dict(TARGET, O=1e21), ddict, print_result=False, nguess=3
            )
        else:
            equilibrium_atmosphere(TARGET, ddict, print_result=False, nguess=3)
        raised = None
    except Exception as err:
        raised = err
    assert type(raised) is ValueError, raised
    assert 'needs melt' in str(raised) and 'Phi_global' in str(raised)


@pytest.mark.physics_invariant
def test_trace_element_converges_to_the_relative_stop():
    """C at w about 1e-6 in the metal: the relative stop holds it to the tol 1e-12 result.

    The absolute stop of 1e-6 in mass fraction left the C core mass 2.0e-5 off.
    """
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': load_partition_table()}
    ddict.update(O_mode_core='partition', Si_mode_core='partition')
    kw = dict(p_guess=P_GUESS, print_result=False, nguess=200, opt_solver=False, rtol=1e-5)
    np.random.seed(0)
    loose = equilibrium_atmosphere(TARGET, ddict, **kw)
    np.random.seed(0)
    tight = equilibrium_atmosphere(TARGET, {**ddict, 'tol_w_core': 1e-12}, **kw)
    assert 1e-7 < loose['w_C_core'] < 1e-5
    for el in ELEMS:
        assert loose[el + '_kg_core'] == pytest.approx(tight[el + '_kg_core'], rel=3e-6), el


@pytest.mark.physics_invariant
@pytest.mark.parametrize('authoritative', [False, True])
def test_metal_solve_with_a_noble_gas(authoritative):
    """He in the solve: the metal passes keep the He pressure in the warm start; He closes."""
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'He_included': 1}
    target = dict(TARGET, He=1e16)
    try:
        out = _solve(ddict, target, dict(P_GUESS, He=1e-3))
    except Exception as exc:
        pytest.fail(f'metal solve with He raised {type(exc).__name__}: {exc}')
    if authoritative:
        target['O'] = out['O_kg_total']
        out = _solve_auth(ddict, target, dict(P_GUESS, fO2_shift_IW=2.0, He=1e-3))
        assert out['fO2_shift_derived'] == pytest.approx(2.0, abs=1e-6)
    he = out['He_kg_atm'] + out['He_kg_liquid']
    assert abs(he - 1e16) <= 1e16 * 1e-10 + 10.0
    assert out['He_bar'] > 0
    for el in ELEMS:
        closure = out[el + '_kg_atm'] + out[el + '_kg_liquid'] + out[el + '_kg_core']
        assert closure == pytest.approx(target[el], rel=1e-4), el
    assert 'He_kg_core' not in out


def _stub_solve(core_of_call, close=True):
    """A solve stand-in whose H core mass follows core_of_call(n) on its n-th call."""
    calls = []

    def solve(target_d, ddict, **kwargs):
        calls.append(1)
        out = {s + '_bar': 1.0 for s in ('H2O', 'CO2', 'N2', 'S2')}
        for el in ELEMS:
            core = core_of_call(len(calls)) if el == 'H' else 0.0
            out[el + '_kg_core'] = core
            out[el + '_kg_liquid'] = 1e18
            out[el + '_kg_atm'] = target_d[el] - 1e18 - core if close else 0.0
        return out

    return solve, calls


def test_fixed_point_that_does_not_settle_raises():
    """A metal H that alternates between 2 values never settles: RuntimeError after 30 passes."""
    from calliope.solve import _with_metal

    solve, calls = _stub_solve(lambda n: 1e19 * (n % 2))
    kwargs = {'rtol': 1e-5, 'atol': 1e10, 'p_guess': None}
    with pytest.raises(RuntimeError, match='did not converge'):
        _with_metal(solve, TARGET, {**_ddict(), **_metal(30.0, 1.9e24)}, kwargs)
    assert len(calls) == 30


def test_state_that_breaks_closure_raises():
    """A returned state whose gas, melt and core miss a budget fails the closure check."""
    from calliope.solve import _with_metal

    solve, calls = _stub_solve(lambda n: 1e19, close=False)
    kwargs = {'rtol': 1e-5, 'atol': 1e10, 'p_guess': None}
    with pytest.raises(RuntimeError, match='closure'):
        _with_metal(solve, TARGET, {**_ddict(), **_metal(30.0, 1.9e24)}, kwargs)
    assert len(calls) == 2


def _offset_solve(miss):
    """A solve stand-in whose gas, melt and core miss each budget by miss[el] kg."""

    def solve(target_d, ddict, **kwargs):
        out = {s + '_bar': 1.0 for s in ('H2O', 'CO2', 'N2', 'S2')}
        for el in ELEMS + ('O',):
            out[el + '_kg_core'] = 0.0
            out[el + '_kg_liquid'] = 1e18
            out[el + '_kg_atm'] = target_d.get(el, 0.0) - 1e18 + miss.get(el, 0.0)
        out['fO2_shift_derived'] = 2.0
        return out

    return solve


# budgets at rtol 1e-5, atol 1e10: the buffered gate is 8e20 rtol + atol = 8.0e15 kg for
# every element; with an O budget each element has its own gate, 8e13 kg for N, 1e17 kg for O
@pytest.mark.parametrize(
    ('o_budget', 'miss', 'closes'),
    [
        (False, {'N': 1e15}, True),
        (False, {'N': 1e16}, False),
        (True, {'O': 5e16}, True),
        (True, {'O': 2e17}, False),
        (True, {'N': 1e15}, False),
        (True, {'O': -1.9e21}, False),
    ],
    ids=[
        'buffered-N-in',
        'buffered-N-out',
        'O-in',
        'O-out',
        'O-budget-N-out',
        'metal-O-charged-to-budget',
    ],
)
def test_closure_check_uses_the_gate_of_the_solve(o_budget, miss, closes):
    """The fixed-point closure check accepts what the solve accepts and refuses the rest."""
    from calliope.solve import _with_metal

    kwargs = {'rtol': 1e-5, 'atol': 1e10, 'p_guess': None}
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'w_O_core': 1e-3}
    target = dict(TARGET, O=1e22)
    try:
        out = _with_metal(_offset_solve(miss), target, ddict, kwargs, o_budget=o_budget)
        raised = None
    except Exception as exc:
        raised = exc
    if closes:
        assert raised is None, raised
        assert out['O_kg_core'] == pytest.approx(1.9e21, rel=1e-12)
    else:
        assert isinstance(raised, RuntimeError), raised
        assert f'{next(iter(miss))} closure' in str(raised)


def test_out_of_range_and_ideal_activities_are_flagged_and_warned_once(caplog, monkeypatch):
    """Flags for every partitioned element; one warning per element and process.

    At P_eq 50 GPa, T_eq 3500 K, dIW -2 the H law (dIW -1.33 to -0.88) is outside its
    range and O and Si have no T or dIW range, so all 3 are flagged; C (dIW -3.9 to -0.53)
    is inside.
    """
    import calliope.solve as solve_mod

    monkeypatch.setattr(solve_mod, '_METAL_WARNED', set())
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': load_partition_table()}
    ddict.update(O_mode_core='partition', Si_mode_core='partition')
    with caplog.at_level('WARNING', logger='fwl.calliope.solve'):
        out = _solve(ddict)
        first = [r.getMessage() for r in caplog.records]
        caplog.clear()
        _solve(ddict)
        second = [r.getMessage() for r in caplog.records]
    flags = {el: out[el + '_D_out_of_range'] for el in ('H', 'C', 'O', 'Si')}
    assert flags == {'H': 1, 'C': 0, 'O': 1, 'Si': 1}
    for el in ('H', 'O', 'Si'):
        assert sum(f'{el} partition law' in m and 'outside' in m for m in first) == 1, el
    assert not any('C partition law' in m for m in first)
    assert sum('ideal activities' in m for m in first) == 2
    assert second == []
    # fixed O and Si use no law, so they carry no flag
    fixed = _solve({**_ddict(), **_metal(30.0, 1.9e24)})
    assert 'O_D_out_of_range' not in fixed and 'Si_D_out_of_range' not in fixed


@pytest.mark.parametrize('value', [math.nan, -1.9e24, 'x', None], ids=repr)
@pytest.mark.parametrize('authoritative', [False, True])
def test_bad_metal_mass_is_refused_by_both_solves(value, authoritative):
    """A NaN, negative or non-number metal_mass raises ValueError; it never solves without metal."""
    from calliope.solve import equilibrium_atmosphere_authoritative_O

    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'metal_mass': value}
    try:
        if authoritative:
            equilibrium_atmosphere_authoritative_O(
                dict(TARGET, O=1e21), ddict, print_result=False, nguess=3
            )
        else:
            equilibrium_atmosphere(TARGET, ddict, print_result=False, nguess=3)
        raised = None
    except Exception as err:
        raised = err
    assert isinstance(raised, ValueError), raised
    assert 'metal_mass' in str(raised)


def test_table_without_a_needed_element_is_refused():
    """The solve needs laws for H, C, N, S and for each partitioned element of O and Si."""
    from calliope.metal import PartitionTableError

    table = {k: v for k, v in load_partition_table().items() if k != 'C'}
    for extra in ({}, {'Si_mode_core': 'partition'}):
        t = table if not extra else {k: v for k, v in table.items() if k != 'Si'} | {'C': 1}
        try:
            _solve({**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': t, **extra})
            raised = None
        except Exception as err:
            raised = err
        assert isinstance(raised, PartitionTableError), raised
        assert ("'C'" if not extra else "'Si'") in str(raised)


@pytest.mark.physics_invariant
def test_placeholder_laws_are_evaluated_on_the_converged_metal_and_melt():
    """Real placeholder laws, half-molten mantle, own melt and NBO/T: D matches its inputs.

    At convergence each reported D equals the law on the reported metal and melt
    compositions, and each core mass equals D M_metal w_melt with the melt of
    Phi_global M_mantle; P and T swapped, NBO/T 2.6, the KLB-1 melt or a melt
    without its volatiles each give a D outside the tolerance.
    """
    from calliope.metal import KLB1, mass_D

    oxides = {**KLB1, 'FeO': 10.0, 'MgO': KLB1['MgO'] - 1.9}
    m_metal, m_melt, nbo_t = 1.9e24, 0.5 * 4.03e24, 1.5
    ddict = {**_ddict(), **_metal(30.0, m_metal), 'Phi_global': 0.5}
    ddict.update(
        partition_table=load_partition_table(),
        silicate_oxides=oxides,
        silicate_nbo_t=nbo_t,
        w_O_core=0.02,
        w_Si_core=0.03,
        tol_w_core=1e-12,
    )
    out = _solve(ddict)
    table = ddict['partition_table']
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    w_melt = {el: out[el + '_kg_liquid'] / m_melt for el in ('H', 'C')}
    for el in ELEMS:
        d = mass_D(table[el], 50.0, 3500.0, -2.0, w, w_melt, oxides, nbo_t)
        assert out[el + '_D_metal'] == pytest.approx(d, rel=1e-6), el
        core = d * m_metal * out[el + '_kg_liquid'] / m_melt
        assert out[el + '_kg_core'] == pytest.approx(core, rel=1e-6), el
    wrong = {
        'P_T_swapped': mass_D(table['C'], 3500.0, 50.0, -2.0, w, w_melt, oxides, nbo_t),
        'nbo_t_default': mass_D(table['C'], 50.0, 3500.0, -2.0, w, w_melt, oxides, 2.6),
        'klb1_melt': mass_D(table['H'], 50.0, 3500.0, -2.0, w, w_melt, KLB1, nbo_t),
        'no_melt_volatiles': mass_D(table['H'], 50.0, 3500.0, -2.0, w, {}, oxides, nbo_t),
    }
    for name, d in wrong.items():
        el = 'C' if name in ('P_T_swapped', 'nbo_t_default') else 'H'
        assert abs(d / out[el + '_D_metal'] - 1) > 1e-5, name


def test_negative_core_mass_names_the_core_in_the_error():
    """A returned negative core mass fails the check with the core mass in the message."""
    from calliope.solve import _with_metal

    solve, _ = _stub_solve(lambda n: -1e6)
    kwargs = {'rtol': 1e-5, 'atol': 1e10, 'p_guess': None}
    try:
        _with_metal(solve, TARGET, {**_ddict(), **_metal(30.0, 1.9e24)}, kwargs)
        raised = None
    except Exception as err:
        raised = err
    assert isinstance(raised, RuntimeError), raised
    assert 'core' in str(raised) and '-1e+06' in str(raised)
