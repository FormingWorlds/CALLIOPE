"""Tests for the metal reservoir in `src/calliope/solve.py`.

Contract: with metal_mass > 0 the H, C, N and S balances carry a core term
D_i M_metal w_i(melt). With the same mass-basis D for all four elements, no O or Si
in the metal and dIW a parameter, the solve equals the solve without metal and a
melt mass larger by D M_metal (an identity of the element balance); the gas
reservoir matters, so the gas-free core fraction Dx/(1 + Dx) misses the result.
Every element closes over gas, melt and core, and metal_mass 0 is the solve
without metal. In the authoritative-O solve the O of the metal is a reservoir of
the O budget: gas, melt and metal O close on the budget, and a metal O above the
budget is refused before the solve.
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


def _solve(ddict):
    """One seeded solve, so a restart draws the same pressures for both runs compared."""
    np.random.seed(0)
    return equilibrium_atmosphere(
        TARGET,
        ddict,
        p_guess=P_GUESS,
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


def _solve_auth(ddict, target):
    """One seeded authoritative-O solve at the tight acceptance of `_solve`."""
    from calliope.solve import equilibrium_atmosphere_authoritative_O

    try:
        return equilibrium_atmosphere_authoritative_O(
            target,
            ddict,
            fO2_hint=2.0,
            p_guess=dict(P_GUESS, fO2_shift_IW=2.0),
            print_result=False,
            nguess=200,
            xtol=1e-12,
            rtol=1e-10,
            atol=1.0,
            opt_solver=False,
            random_seed=0,
        )
    except (RuntimeError, ValueError) as exc:
        pytest.fail(f'authoritative-O solve with metal raised {type(exc).__name__}: {exc}')


# O of gas and melt in the solve without metal at fO2 IW+2 (TARGET, _ddict) is 1.22e21 kg
O_GAS_MELT = 1.2e21


@pytest.mark.physics_invariant
def test_authoritative_o_closes_the_o_budget_with_metal():
    """Authoritative-O with 0.1 wt% O in the metal: gas + melt + metal O equals the O budget."""
    m_metal, w_o = 1.9e24, 1e-3
    target = dict(TARGET, O=O_GAS_MELT + w_o * m_metal)
    out = _solve_auth({**_ddict(), **_metal(30.0, m_metal), 'w_O_core': w_o}, target)
    assert out['O_kg_core'] == pytest.approx(w_o * m_metal, rel=1e-12)
    assert out['w_O_core'] == pytest.approx(w_o, rel=1e-12)
    for el in ('H', 'C', 'N', 'S', 'O'):
        closure = out[el + '_kg_atm'] + out[el + '_kg_liquid'] + out[el + '_kg_core']
        assert abs(closure - target[el]) <= target[el] * 1e-10 + 10.0, el
    # dropping the metal O would leave a gap of 1.9e21 kg, far above the 3.1e11 kg tolerance
    assert out['O_kg_core'] > 1e6 * target['O'] * 1e-10
    # gas and melt hold the budget less the metal O
    assert out['O_kg_atm'] + out['O_kg_liquid'] == pytest.approx(O_GAS_MELT, rel=1e-6)
    assert -12 < out['fO2_shift_derived'] < 12
    assert out['H_kg_core'] > 0


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_authoritative_o_with_equal_d_is_a_larger_melt():
    """Authoritative-O, equal D = 30, no metal O: the solve of a melt larger by D M_metal.

    The larger melt also holds the O of its dissolved species, so its O budget is larger
    by D M_metal times the O mass fraction of the melt of the metal solve.
    """
    d, m_metal, m_melt = 30.0, 1.9e24, 4.03e24
    target = dict(TARGET, O=O_GAS_MELT)
    with_metal = _solve_auth({**_ddict(), **_metal(d, m_metal)}, target)
    o_extra = d * m_metal * with_metal['O_kg_liquid'] / m_melt
    larger = _solve_auth(_ddict(m_melt + d * m_metal), dict(target, O=O_GAS_MELT + o_extra))
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
    assert with_metal['H_kg_core'] > 0.1 * TARGET['H']


def test_metal_o_above_the_o_budget_is_refused():
    """A metal O larger than the O budget leaves no O for gas and melt: refused before solving."""
    from calliope.solve import equilibrium_atmosphere_authoritative_O

    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'w_O_core': 0.01}
    with pytest.raises(Exception) as exc:
        equilibrium_atmosphere_authoritative_O(
            dict(TARGET, O=O_GAS_MELT), ddict, print_result=False, nguess=3, random_seed=0
        )
    assert exc.type is ValueError
    assert 'O budget' in str(exc.value)
    # the edge: 1 kg of O left for gas and melt is not refused by the check
    from calliope.solve import _check_metal_o

    _check_metal_o(1.9e22 - 1.0, {'O': 1.9e22})
    with pytest.raises(ValueError, match='O budget'):
        _check_metal_o(1.9e22, {'O': 1.9e22})


@pytest.mark.physics_invariant
@pytest.mark.parametrize('authoritative', [False, True])
def test_metal_solve_with_a_noble_gas(authoritative):
    """He in the solve: the metal passes keep the He pressure in the warm start; He closes."""
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'He_included': 1}
    target = dict(TARGET, He=1e16)
    if authoritative:
        target['O'] = O_GAS_MELT
        out = _solve_auth(ddict, target)
    else:
        try:
            out = _solve_noble(ddict, target)
        except (RuntimeError, ValueError) as exc:
            pytest.fail(f'metal solve with He raised {type(exc).__name__}: {exc}')
    he = out['He_kg_atm'] + out['He_kg_liquid']
    assert abs(he - 1e16) <= 1e16 * 1e-10 + 10.0
    assert out['He_bar'] > 0
    for el in ELEMS:
        closure = out[el + '_kg_atm'] + out[el + '_kg_liquid'] + out[el + '_kg_core']
        assert closure == pytest.approx(target[el], rel=1e-4), el
    assert 'He_kg_core' not in out


def _solve_noble(ddict, target):
    """One seeded solve with He, at the tight acceptance of `_solve`."""
    np.random.seed(0)
    return equilibrium_atmosphere(
        target,
        ddict,
        p_guess=dict(P_GUESS, He=1e-3),
        print_result=False,
        nguess=200,
        xtol=1e-12,
        rtol=1e-10,
        atol=1.0,
        opt_solver=False,
    )


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


def test_o_budget_state_that_breaks_o_closure_raises():
    """With an O budget, a returned state whose gas, melt and metal O miss it fails closure."""
    from calliope.solve import _with_metal

    def solve(target_d, ddict, **kwargs):
        out = {s + '_bar': 1.0 for s in ('H2O', 'CO2', 'N2', 'S2')}
        for el in ELEMS:
            out[el + '_kg_core'] = 0.0
            out[el + '_kg_liquid'] = 1e18
            out[el + '_kg_atm'] = target_d[el] - 1e18
        out['O_kg_liquid'], out['O_kg_atm'] = 1e18, 1e18
        out['fO2_shift_derived'] = 2.0
        return out

    kwargs = {'rtol': 1e-5, 'atol': 1e10, 'p_guess': None}
    ddict = {**_ddict(), **_metal(30.0, 1.9e24), 'w_O_core': 1e-3}
    with pytest.raises(RuntimeError, match='O closure'):
        _with_metal(solve, dict(TARGET, O=1e22), ddict, kwargs, o_budget=True)
    # the same state with an O budget equal to gas, melt and metal O closes
    out = _with_metal(solve, dict(TARGET, O=2e18 + 1.9e21), ddict, dict(kwargs), o_budget=True)
    assert out['O_kg_core'] == pytest.approx(1.9e21, rel=1e-12)
    assert out['w_O_core'] == pytest.approx(1e-3, rel=1e-12)
