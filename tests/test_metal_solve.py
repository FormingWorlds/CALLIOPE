"""Tests for the metal reservoir in `src/calliope/solve.py` (equilibrium_atmosphere).

Contract: with metal_mass > 0 the H, C, N and S balances carry a core term
D_i M_metal w_i(melt). With the same mass-basis D for all four elements, no O or Si
in the metal and dIW a parameter, the solve equals the solve without metal and a
melt mass larger by D M_metal (an identity of the element balance); the gas
reservoir matters, so the gas-free core fraction Dx/(1 + Dx) misses the result.
Every element closes over gas, melt and core, and metal_mass 0 is the solve
without metal.
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
