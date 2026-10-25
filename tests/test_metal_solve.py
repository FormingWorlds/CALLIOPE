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


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
@pytest.mark.parametrize('d', [1.0, 30.0])
def test_metal_with_equal_d_is_a_larger_melt(d):
    """Equal D for H, C, N, S: metal solve = no-metal solve with M_melt + D M_metal."""
    m_metal = 1.9e24
    with_metal = equilibrium_atmosphere(
        TARGET,
        {**_ddict(), **_metal(d, m_metal)},
        p_guess=P_GUESS,
        print_result=False,
        nguess=200,
    )
    larger = equilibrium_atmosphere(
        TARGET, _ddict(4.03e24 + d * m_metal), p_guess=P_GUESS, print_result=False, nguess=200
    )
    for sp in ('H2O', 'CO2', 'N2', 'S2', 'H2', 'CO', 'CH4'):
        assert with_metal[sp + '_bar'] == pytest.approx(larger[sp + '_bar'], rel=1e-6)
    for el in ELEMS:
        melt, core = with_metal[el + '_kg_liquid'], with_metal[el + '_kg_core']
        assert core == pytest.approx(d * m_metal * melt / 4.03e24, rel=1e-12)
        assert melt + core == pytest.approx(larger[el + '_kg_liquid'], rel=1e-6)
        closure = with_metal[el + '_kg_atm'] + melt + core
        assert closure == pytest.approx(TARGET[el], rel=1e-4)
    # the gas holds mass: the gas-free core fraction Dx/(1 + Dx) misses the H core
    x = m_metal / 4.03e24
    gas_free = TARGET['H'] * d * x / (1 + d * x)
    assert abs(with_metal['H_kg_core'] / gas_free - 1) > 1e-3


def test_zero_metal_mass_is_the_solve_without_metal():
    """metal_mass 0 gives the solve without metal, bit for bit, and zero core masses."""
    plain = equilibrium_atmosphere(TARGET, _ddict(), p_guess=P_GUESS, print_result=False)
    zero = equilibrium_atmosphere(
        TARGET, {**_ddict(), **_metal(30.0, 0.0)}, p_guess=P_GUESS, print_result=False
    )
    for sp in ('H2O', 'CO2', 'N2', 'S2'):
        assert zero[sp + '_bar'] == plain[sp + '_bar']
    assert all(zero[el + '_kg_core'] == 0 for el in ELEMS)
