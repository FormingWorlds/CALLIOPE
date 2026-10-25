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


@pytest.mark.parametrize('authoritative', [False, True])
def test_zero_metal_mass_is_the_solve_without_metal(authoritative):
    """metal_mass 0 gives the solve without metal bit for bit, plus 4 zero core masses."""
    if authoritative:
        target = dict(TARGET, O=_solve(_ddict())['O_kg_total'])
        plain = _solve_auth(_ddict(), target)
        zero = _solve_auth({**_ddict(), **_metal(30.0, 0.0)}, target)
    else:
        plain = _solve(_ddict())
        zero = _solve({**_ddict(), **_metal(30.0, 0.0)})
    assert set(zero) - set(plain) == {el + '_kg_core' for el in ELEMS}
    assert all(zero[el + '_kg_core'] == 0 for el in ELEMS)
    nan = [k for k in plain if plain[k] != plain[k]]
    assert all(zero[k] != zero[k] for k in nan)
    assert [k for k in plain if k not in nan and zero[k] != plain[k]] == []


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
    # x' of Tagawa et al. 2021 Eq. 4: mole fractions over the metal without H and C
    n_total = (1 - sum(w.values())) / molar_mass_rock['Fe'] + sum(
        w[el] / (masses.get(el) or {'N': 0.014007, 'S': 0.03206}[el])
        for el in w
        if el not in ('H', 'C')
    )
    for el in ('O', 'Si'):
        assert w[el] / masses[el] / n_total == pytest.approx(x[el], rel=1e-4)
    assert 0 < w['O'] < 0.2 and 0 < w['Si'] < 0.2
    assert out['O_kg_core'] == pytest.approx(w['O'] * 1.9e24, rel=1e-12)


@pytest.mark.physics_invariant
def test_diw_from_feo_matches_the_melt_and_metal():
    """dIW_eq_mode 'from_FeO': dIW = 2 log10(x_FeO / x'_Fe) on the KLB-1 melt and the final metal.

    x_FeO is on the single-cation basis of the melt and x'_Fe over the metal without H and C
    (Tagawa et al. 2021 Eq. 4); the oxide-molecule x_FeO would give -2.3763.
    """
    from calliope.metal import KLB1, molar_mass_rock, silicate_units_per_kg

    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'dIW_eq_mode': 'from_FeO', 'w_Si_core': 0.05}
    out = _solve(ddict)
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    m = {'N': 0.014007, 'S': 0.03206, 'O': 0.015999, 'Si': molar_mass_rock['Si']}
    n_fe = (1 - sum(w.values())) / molar_mass_rock['Fe']
    x_fe = n_fe / (n_fe + sum(w[el] / m[el] for el in m))
    x_feo = 0.081 / (molar_mass_rock['Fe'] + 0.015999) / silicate_units_per_kg(KLB1)
    assert out['dIW_eq'] == pytest.approx(2 * math.log10(x_feo / x_fe), abs=1e-4)
    assert out['dIW_eq'] == pytest.approx(-2.3975, abs=1e-4)
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


@pytest.mark.physics_invariant
def test_user_silicate_oxides_set_the_exchange_feo_and_diw():
    """The melt oxides of the caller, not KLB-1, set x_SiO2, x_FeO of dIW and the melt FeO."""
    from calliope.metal import KLB1, exchange_metal, molar_mass_rock, silicate_units_per_kg

    oxides = {**KLB1, 'FeO': 10.0, 'SiO2': KLB1['SiO2'] + 2.0, 'MgO': KLB1['MgO'] - 3.9}
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': load_partition_table()}
    ddict.update(O_mode_core='partition', Si_mode_core='partition', dIW_eq_mode='from_FeO')
    ddict.update(silicate_oxides=oxides)
    out = _solve(ddict)
    m_feo = molar_mass_rock['Fe'] + 0.015999
    assert out['n_FeO_silicate'] == pytest.approx(4.03e24 * 0.10 / m_feo, rel=1e-4)
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    m = {'N': 0.014007, 'S': 0.03206, 'O': 0.015999, 'Si': molar_mass_rock['Si']}
    n_fe = (1 - sum(w.values())) / molar_mass_rock['Fe']
    n_free = n_fe + sum(w[el] / m[el] for el in m)
    x_feo = 0.10 / m_feo / silicate_units_per_kg(oxides)
    assert out['dIW_eq'] == pytest.approx(2 * math.log10(x_feo / (n_fe / n_free)), abs=1e-4)
    x = exchange_metal(ddict['partition_table'], 50.0, 3500.0, out['dIW_eq'], oxides)
    for el in ('O', 'Si'):
        assert w[el] / m[el] / n_free == pytest.approx(x[el], rel=1e-4), el


def test_placeholder_table_warns_once(caplog, monkeypatch):
    """A table whose meta says placeholder warns once per process; a real table does not."""
    import calliope.solve as solve_mod

    monkeypatch.setattr(solve_mod, '_METAL_WARNED', set())
    ddict = {**_ddict(), **_metal(30.0, 1.9e24)}
    real = dict(
        ddict,
        partition_table={
            **ddict['partition_table'],
            'meta': {'name': 'x', 'placeholder': False},
        },
    )
    with caplog.at_level('WARNING', logger='fwl.calliope.solve'):
        _solve(real)
        assert not [m for m in caplog.messages if 'placeholder' in m]
        _solve(ddict)
        _solve(ddict)
    assert len([m for m in caplog.messages if 'placeholder' in m]) == 1


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


# M_FeO, M_Si from the element masses (kg/mol); KLB-1 FeO 8.10 wt% in 4.03e24 kg of melt
_M_O, _M_SI, _M_FE = 0.015999, 0.0280855, 0.055845
_N_FEO_SIL = 4.03e24 * 0.081 / (_M_FE + _M_O)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('w_o', 'w_si', 'flag'),
    [(1e-3, 0.0, 0), (5e-3, 0.0, 1), (0.02, 0.03, 1), (1e-6, 0.05, 1)],
    ids=['ratio-0.026', 'ratio-0.13', 'si-returns-feo', 'trace-o-in-si-rich-metal'],
)
def test_feo_exchange_balances_o_and_si_atoms(w_o, w_si, flag, caplog, monkeypatch):
    """n_FeO_exchange = n_O - 2 n_Si of the metal; above 10 % of the melt FeO it is flagged.

    FeO = Fe + O takes 1 FeO per metal O, SiO2 + 2 Fe = Si + 2 FeO gives back 2 FeO per
    metal Si; the ratio to the KLB-1 melt FeO is 0.026, 0.13, -0.37 and -1.49 for the 4
    cases. In the last, n_O is 2e-5 of 2 n_Si.
    """
    import calliope.solve as solve_mod

    monkeypatch.setattr(solve_mod, '_METAL_WARNED', set())
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'w_O_core': w_o, 'w_Si_core': w_si}
    with caplog.at_level('WARNING', logger='fwl.calliope.solve'):
        try:
            out = _solve(ddict)
            _solve(ddict)
        except Exception as exc:
            pytest.fail(f'metal solve raised {type(exc).__name__}: {exc}')
    n_ex = w_o * 1.9e24 / _M_O - 2 * w_si * 1.9e24 / _M_SI
    assert out['n_FeO_exchange'] == pytest.approx(n_ex, rel=1e-4)
    assert out['FeO_kg_exchange'] == pytest.approx(n_ex * (_M_FE + _M_O), rel=1e-4)
    assert out['n_FeO_silicate'] == pytest.approx(_N_FEO_SIL, rel=1e-4)
    assert out['FeO_exchange_out_of_range'] == flag
    warned = [m for m in caplog.messages if '10 % of the silicate FeO' in m]
    assert len(warned) == flag
    if w_si:
        # without the Si term the exchange would be n_O alone, 2.4e24 mol of the other sign
        assert abs(out['n_FeO_exchange'] - w_o * 1.9e24 / _M_O) > abs(n_ex)


@pytest.mark.physics_invariant
def test_partition_o_budget_round_trip_on_a_partial_melt():
    """Partition O and Si, 70 % melt, IW+1: the O-budget solve returns IW+1.

    Each pass evaluates D_C on the O and Si of the same pass (the C law has eps_O = -19.5),
    and the melt FeO is that of Phi_global M_mantle, not of the whole mantle.
    """
    from calliope.metal import KLB1, mass_D

    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': load_partition_table()}
    ddict.update(O_mode_core='partition', Si_mode_core='partition', P_eq=40.0)
    ddict.update(Phi_global=0.7, fO2_shift_IW=1.0)
    buffered = _solve(ddict)
    target = dict(TARGET, O=buffered['O_kg_total'])
    out = _solve_auth(ddict, target, dict(P_GUESS, fO2_shift_IW=1.0))
    assert out['fO2_shift_derived'] == pytest.approx(1.0, abs=1e-6)
    assert abs(out['O_kg_atm'] + out['O_kg_liquid'] - target['O']) <= target['O'] * 1e-10 + 10
    assert out['O_kg_core'] == pytest.approx(buffered['O_kg_core'], rel=1e-6)
    assert out['n_FeO_silicate'] == pytest.approx(0.7 * _N_FEO_SIL, rel=1e-4)
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    w_melt = {el: out[el + '_kg_liquid'] / (0.7 * 4.03e24) for el in ('H', 'C')}
    at = (ddict['partition_table']['C'], 40.0, 3500.0, -2.0)
    d_c = mass_D(*at, w, w_melt, KLB1, 2.6)
    no_o = mass_D(*at, {**w, 'O': 0.0, 'Si': 0.0}, w_melt, KLB1, 2.6)
    assert out['C_D_metal'] == pytest.approx(d_c, rel=1e-5)
    assert d_c / no_o > 5


@pytest.mark.physics_invariant
def test_si_partition_with_diw_from_feo_settles(monkeypatch):
    """Si 'partition' with dIW from FeO, 34 % melt, 4487 K: the metal Si settles in 10 passes.

    More metal Si lowers x'_Fe, raises dIW and lowers the exchange x_Si; undamped, the Si
    of successive passes alternates about the fixed point with a ratio of about -0.6.
    """
    import calliope.solve as solve_mod
    from calliope.metal import exchange_metal

    passes, real = [0], solve_mod._with_metal

    def counting(solve, *args, **kwargs):
        def counted(*a, **k):
            passes[0] += 1
            return solve(*a, **k)

        return real(counted, *args, **kwargs)

    monkeypatch.setattr(solve_mod, '_with_metal', counting)

    ddict = {**_ddict(), **_metal(30.0, 1.6166e24), 'partition_table': load_partition_table()}
    ddict.update(P_eq=52.64, T_eq=4487.46, Phi_global=0.3386, silicate_nbo_t=3.105)
    ddict.update(Si_mode_core='partition', dIW_eq_mode='from_FeO', w_O_core=0.00534)
    try:
        out = _solve(ddict)
    except Exception as exc:
        pytest.fail(f'metal solve raised {type(exc).__name__}: {exc}')
    x_si = exchange_metal(
        ddict['partition_table'], 52.64, 4487.46, out['dIW_eq'], elements=('Si',)
    )
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    m = {'N': 0.014007, 'S': 0.03206, 'O': 0.015999, 'Si': 0.028085}
    n_free = (1 - sum(w.values())) / 0.055845 + sum(w[el] / m[el] for el in m)
    assert w['Si'] / m['Si'] / n_free == pytest.approx(x_si['Si'], rel=1e-4)
    assert -2.3 < out['dIW_eq'] < -2.0
    assert 1 <= passes[0] <= 10


@pytest.mark.parametrize(
    'extra',
    [
        {'Si_mode_core': 'partition', 'dIW_eq': -6.0},
        {'O_mode_core': 'partition', 'dIW_eq_mode': 'from_FeO', 'T_eq': 5000.0, 'feo': 17.9},
        {'O_mode_core': 'partition', 'dIW_eq': 615.9, 'T_eq': 5000.0},
    ],
    ids=['si-at-iw-minus-6', 'o-from-feo-rich-melt', 'o-exchange-overflows'],
)
def test_exchange_without_iron_raises_the_named_error(extra, request):
    """An exchange without iron that does not fall, or one that is not finite, raises ValueError.

    The message names the partitioned element and its mode field only.
    """
    from calliope.metal import KLB1

    extra = dict(extra)
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': load_partition_table()}
    if 'feo' in extra:
        feo = extra.pop('feo')
        ddict['silicate_oxides'] = {
            **KLB1,
            'FeO': feo,
            'MgO': KLB1['MgO'] - (feo - KLB1['FeO']),
        }
    ddict.update(extra)
    el = 'Si' if 'Si_mode_core' in extra else 'O'
    overflow = request.node.callspec.id == 'o-exchange-overflows'
    message = 'exchange is not finite' if overflow else 'does not fall'
    with pytest.raises(ValueError, match=message) as exc:
        _solve(ddict)
    if message == 'does not fall':
        assert f"{el}_mode_core 'partition'" in str(exc.value)
        assert 'no fixed point was reached from the start metal' in str(exc.value)
        assert ({'O', 'Si'} - {el}).pop() + '_mode_core' not in str(exc.value)


def test_stop_reads_the_whole_exchange_step(monkeypatch):
    """A damped step below tol_w_core with a whole step 500 times above it does not stop.

    The stub exchange is linear in w_Si with slope -1000 (lam 1/1001): pass 2 moves w_Si by
    0.4 tol w_Si onto the fixed point while its whole step is 400 tol w_Si, so the solve stops
    on pass 3.
    """
    import calliope.solve as solve_mod
    from calliope.metal import molar_mass_rock

    w0, tol, s = 0.05, 1e-6, -1000.0
    a = 0.4 * tol * w0
    w_star = w0 + a / (1 - s)
    targets = iter([w0 + a, w0 + a + s * a] + [w_star] * 28)
    monkeypatch.setattr(solve_mod, 'exchange_metal', lambda *args, **kw: {'Si': next(targets)})
    m_si = molar_mass_rock['Si']
    monkeypatch.setattr(solve_mod, '_metal_moles', lambda w, exclude=(): (1.0 / m_si, {}))
    solve, calls = _stub_solve(lambda n: 1e19)
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'Si_mode_core': 'partition', 'w_Si_core': w0}
    ddict['tol_w_core'] = tol
    out = solve_mod._with_metal(
        solve, TARGET, ddict, {'rtol': 1e-5, 'atol': 1e10, 'p_guess': None}
    )
    assert len(calls) == 3
    assert out['w_Si_core'] == pytest.approx(w_star, rel=1e-12)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('t_eq', 'both'), [(4800.0, False), (5000.0, True), (6500.0, False), (8000.0, False)]
)
def test_exchange_beyond_the_iron_of_pure_fe_metal_settles(t_eq, both):
    """From pure Fe metal the first exchange holds no iron (x_Si 1.10 at 4800 K); the solve settles.

    The steps of a pass stop short of a metal without iron, and the converged metal holds the
    exchange x' of its own dIW. At 6500 and 8000 K the exchange sum also exceeds 1 on later
    passes while it falls toward its fixed point.
    """
    from calliope.metal import exchange_metal

    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'partition_table': load_partition_table()}
    ddict.update(Si_mode_core='partition', dIW_eq_mode='from_FeO', T_eq=t_eq)
    if both:
        ddict['O_mode_core'] = 'partition'
    try:
        out = _solve(ddict)
    except Exception as exc:
        pytest.fail(f'metal solve raised {type(exc).__name__}: {exc}')
    elements = ('O', 'Si') if both else ('Si',)
    x = exchange_metal(ddict['partition_table'], 50.0, t_eq, out['dIW_eq'], elements=elements)
    w = {el: out['w_' + el + '_core'] for el in ('H', 'C', 'N', 'S', 'O', 'Si')}
    m = {'N': 0.014007, 'S': 0.03206, 'O': 0.015999, 'Si': 0.028085}
    n_free = (1 - sum(w.values())) / 0.055845 + sum(w[el] / m[el] for el in m)
    for el in elements:
        assert w[el] / m[el] / n_free == pytest.approx(x[el], rel=1e-4), el
    assert sum(w.values()) < 0.6


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_equal_d_identity_plot(validation_plot):
    """Equal D: partial pressures with metal against a melt larger by D M_metal, both solves.

    Buffered at D = 1 and 30 and authoritative-O at D = 30 (its O budget raised by the O of
    the extra melt): every partial pressure agrees to relative 1e-6.
    """
    sp = ('H2O', 'CO2', 'N2', 'S2', 'H2', 'CO', 'CH4')
    pairs = {}
    for d in (1.0, 30.0):
        metal = _solve({**_ddict(), **_metal(d, 1.9e24)})
        larger = _solve(_ddict(4.03e24 + d * 1.9e24))
        pairs[f'D = {d:g}'] = (larger, metal)
    metal_dict = {**_ddict(), **_metal(30.0, 1.9e24)}
    target = dict(TARGET, O=_solve(metal_dict)['O_kg_total'])
    auth = _solve_auth(metal_dict, target)
    o_extra = 30.0 * 1.9e24 * auth['O_kg_liquid'] / 4.03e24
    big = _solve_auth(_ddict(4.03e24 + 30.0 * 1.9e24), dict(target, O=target['O'] + o_extra))
    pairs['D = 30, O budget'] = (big, auth)
    for name, (ref, got) in pairs.items():
        for s in sp:
            assert got[s + '_bar'] == pytest.approx(ref[s + '_bar'], rel=1e-6), (name, s)
    # the metal changes the pressures: the solve without metal differs by far more
    plain = _solve(_ddict())
    assert abs(plain['H2O_bar'] / pairs['D = 30'][1]['H2O_bar'] - 1) > 0.5

    def figure():
        import matplotlib.pyplot as plt

        values = [r[s + '_bar'] for pair in pairs.values() for r in pair for s in sp]
        lo, hi = math.floor(math.log10(min(values))), math.ceil(math.log10(max(values)))
        edge, ticks = [10.0**lo, 10.0**hi], [10.0**k for k in range(lo, hi + 1, 4)]
        fig, ax = plt.subplots()
        for (name, (ref, got)), mk in zip(pairs.items(), ('o', 's', '+')):
            x = [ref[s + '_bar'] for s in sp]
            ax.loglog(x, [got[s + '_bar'] for s in sp], mk, mfc='none', ms=9, label=name)
        ax.plot(edge, edge, ls='--', color=plt.rcParams['text.color'], label='1:1')
        ax.set(xlabel='p, larger melt [bar]', ylabel='p, with metal [bar]')
        ax.set(xlim=edge, ylim=edge, xticks=ticks, yticks=ticks)
        ax.minorticks_off()
        ax.legend(fontsize='small')
        fig.tight_layout()
        return fig

    validation_plot(figure, 'metal_equal_d_larger_melt')
