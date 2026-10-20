"""Elemental mass-balance consistency tests over the 75-case parameter grid.

Verifies stoichiometric species-to-element mass recovery against the input
elemental inventory across the 60-point buffered grid (5 fO2 x 4 T x 3
inventories), the 12-point authoritative-oxygen grid (4 T x 3 inventories),
2 volatile-inventory configurations using elemental budgets from PROTEUS
reference inputs at T = 2000 K, and a 5-noble-gas complement configuration.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pytest

import calliope.solve as solve_mod
from calliope.constants import molar_mass, noble_gases, ocean_moles, volatile_species
from calliope.solve import (
    TRUNC_MASS,
    atmosphere_mass,
    dissolved_mass,
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
    func,
    func_authoritative_O,
)

logging.getLogger('fwl.calliope').setLevel(logging.WARNING)

pytestmark = [
    pytest.mark.integration,
    pytest.mark.physics_invariant,
    pytest.mark.timeout(300),
]

SPECIES_STOICH: dict[str, dict[str, int]] = {
    'H2O': {'H': 2, 'O': 1},
    'CO2': {'C': 1, 'O': 2},
    'N2': {'N': 2},
    'S2': {'S': 2},
    'SO2': {'S': 1, 'O': 2},
    'H2S': {'H': 2, 'S': 1},
    'NH3': {'N': 1, 'H': 3},
    'H2': {'H': 2},
    'CH4': {'C': 1, 'H': 4},
    'CO': {'C': 1, 'O': 1},
    'O2': {'O': 2},
    'He': {'He': 1},
    'Ne': {'Ne': 1},
    'Ar': {'Ar': 1},
    'Kr': {'Kr': 1},
    'Xe': {'Xe': 1},
}

RTOL_DERIVED = 1.0e-5
ATOL_KG = 1.0e10


def _build_ddict(
    m_mantle: float,
    g_surf: float,
    r_int: float,
    phi_global: float,
    t_magma: float,
    fo2_shift: float,
    active_nobles: tuple[str, ...] = (),
) -> dict[str, float | int]:
    """Construct CALLIOPE parameter dictionary matching PROTEUS defaults."""
    ddict: dict[str, float | int] = {
        'M_mantle': m_mantle,
        'gravity': g_surf,
        'radius': r_int,
        'Phi_global': phi_global,
        'T_magma': t_magma,
        'fO2_shift_IW': fo2_shift,
    }
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0
    for ng in noble_gases:
        ddict[f'{ng}_included'] = 1 if ng in active_nobles else 0
        ddict[f'{ng}_initial_bar'] = 0.0
        ddict[f'{ng}_ppmw'] = 0.0
    return ddict


def _reconstruct_from_species(
    out: dict[str, float],
    r_int: float,
    g_surf: float,
    elements: tuple[str, ...],
) -> dict[str, dict[str, float]]:
    """Reconstruct per-element gas and liquid masses [kg] from species partial pressures and dissolved masses."""
    p_tot = sum(float(out.get(f'{sp}_bar', 0.0)) for sp in SPECIES_STOICH)
    mu_atm = (
        sum(float(out.get(f'{sp}_bar', 0.0)) * float(molar_mass[sp]) for sp in SPECIES_STOICH)
        / p_tot
    )
    p_to_m = 1.0e5 * 4.0 * math.pi * (r_int**2) / g_surf
    res: dict[str, dict[str, float]] = {
        el: {'gas_kg': 0.0, 'liq_kg': 0.0, 'total_kg': 0.0} for el in elements
    }
    for sp, stoich in SPECIES_STOICH.items():
        m_sp = float(molar_mass[sp])
        p_bar = float(out.get(f'{sp}_bar', 0.0))
        sp_gas_from_bar = p_bar * p_to_m * (m_sp / mu_atm)
        sp_liq = float(out.get(f'{sp}_kg_liquid', 0.0))
        for el, n_atom in stoich.items():
            if el not in res:
                continue
            w_frac_gas = 1.0 if sp in noble_gases else (n_atom * float(molar_mass[el])) / m_sp
            w_frac_liq = (
                1.0
                if sp in ('N2', 'S2', *noble_gases)
                else (n_atom * float(molar_mass[el])) / m_sp
            )
            res[el]['gas_kg'] += sp_gas_from_bar * w_frac_gas
            res[el]['liq_kg'] += sp_liq * w_frac_liq
    for el in elements:
        res[el]['total_kg'] = res[el]['gas_kg'] + res[el]['liq_kg']
    return res


def _build_grid_cases() -> list[dict]:
    """Build the 75-case parameter grid (60 buffered + 12 authoritative O + 3 reference-budget/noble)."""
    mass_ocean = float(ocean_moles * molar_mass['H2'])
    m_earth = 5.972e24
    r_earth = 6.335439e6
    cmf = 0.325
    m_mantle_default = (1.0 - cmf) * m_earth
    g_default = 9.81

    inventories: dict[str, dict[str, float]] = {}
    for label, scale in (('low', 0.1), ('default', 1.0), ('high', 10.0)):
        h_kg = scale * mass_ocean
        c_kg = 1.0 * h_kg
        n_kg = 0.5 * h_kg
        s_kg = 1.0 * h_kg
        o_stoich_kg = (
            h_kg * (molar_mass['O'] / molar_mass['H2'])
            + c_kg * (2.0 * molar_mass['O'] / molar_mass['C'])
            + s_kg * (2.0 * molar_mass['O'] / molar_mass['S'])
        )
        inventories[label] = {
            'H': h_kg,
            'C': c_kg,
            'N': n_kg,
            'S': s_kg,
            'O': o_stoich_kg,
        }

    fo2_grid = (-4.0, -2.0, 0.0, 2.0, 4.0)
    t_grid = (1500.0, 2000.0, 2500.0, 3000.0)
    inv_labels = ('low', 'default', 'high')

    cases: list[dict] = []
    for inv_label in inv_labels:
        inv = inventories[inv_label]
        for t_k in t_grid:
            for fo2 in fo2_grid:
                cases.append(
                    {
                        'case_id': f'buffered_{inv_label}_T{int(t_k)}_IW{fo2:+.0f}',
                        'mode': 'buffered',
                        'T_K': t_k,
                        'fO2_shift_IW': fo2,
                        'M_mantle': m_mantle_default,
                        'R_int': r_earth,
                        'gravity': g_default,
                        'Phi_global': 1.0,
                        'targets': {k: inv[k] for k in ('H', 'C', 'N', 'S')},
                        'nobles': (),
                    }
                )
    for inv_label in inv_labels:
        inv = inventories[inv_label]
        for t_k in t_grid:
            cases.append(
                {
                    'case_id': f'authO_{inv_label}_T{int(t_k)}',
                    'mode': 'authoritative_O',
                    'T_K': t_k,
                    'fO2_shift_IW': 2.0,
                    'M_mantle': m_mantle_default,
                    'R_int': r_earth,
                    'gravity': g_default,
                    'Phi_global': 1.0,
                    'targets': {k: inv[k] for k in ('H', 'C', 'N', 'S', 'O')},
                    'nobles': (),
                }
            )
    # Volatile budgets and fO2_shift_IW from PROTEUS reference inputs evaluated at T = 2000 K.
    cases.append(
        {
            'case_id': 'budget_tutorial_earth_T2000_IW+2',
            'mode': 'buffered',
            'T_K': 2000.0,
            'fO2_shift_IW': 2.0,
            'M_mantle': m_mantle_default,
            'R_int': r_earth,
            'gravity': g_default,
            'Phi_global': 1.0,
            'targets': {'H': 4.7e20, 'C': 2.73e20, 'N': 0.0, 'S': 0.0},
            'nobles': (),
        }
    )
    m_p_cmb = 0.5 * m_earth
    m_m_cmb = (1.0 - cmf) * m_p_cmb
    r_cmb = (0.5 ** (1.0 / 3.0)) * r_earth
    g_cmb = g_default * (0.5 ** (1.0 / 3.0))
    cases.append(
        {
            'case_id': 'budget_cmb_singleimpact_T2000_IW+1',
            'mode': 'buffered',
            'T_K': 2000.0,
            'fO2_shift_IW': 1.0,
            'M_mantle': m_m_cmb,
            'R_int': r_cmb,
            'gravity': g_cmb,
            'Phi_global': 1.0,
            'targets': {'H': 2.4e20, 'C': 1.4e20, 'N': 4.0e18, 'S': 1.0e20},
            'nobles': (),
        }
    )
    noble_mass = 1.0e-6 * m_mantle_default
    noble_targets = {k: inventories['default'][k] for k in ('H', 'C', 'N', 'S')}
    for ng in ('He', 'Ne', 'Ar', 'Kr', 'Xe'):
        noble_targets[ng] = noble_mass
    cases.append(
        {
            'case_id': 'budget_noble_gas_complement_T2000_IW+2',
            'mode': 'buffered_noble',
            'T_K': 2000.0,
            'fO2_shift_IW': 2.0,
            'M_mantle': m_mantle_default,
            'R_int': r_earth,
            'gravity': g_default,
            'Phi_global': 1.0,
            'targets': noble_targets,
            'nobles': ('He', 'Ne', 'Ar', 'Kr', 'Xe'),
        }
    )
    return cases


def _construct_zero_target_p_guess(
    target: dict[str, float],
    randomize: bool = False,
) -> dict[str, float] | None:
    """Construct an initial guess respecting zero-target elements.

    When any element among ('H', 'C', 'N', 'S') has a zero target, return a guess
    with zero pressure for species containing zero-target elements. When randomize is
    True, draw non-zero slots from the active NumPy random generator.
    """
    if not any(float(target.get(el, 1.0)) == 0.0 for el in ('H', 'C', 'N', 'S')):
        return None
    if randomize:
        return {
            'H2O': (
                float(10 ** np.random.uniform(low=-1.0, high=2.0))
                if float(target.get('H', 1.0)) > 0
                else 0.0
            ),
            'CO2': (
                float(10 ** np.random.uniform(low=0.0, high=3.0))
                if float(target.get('C', 1.0)) > 0
                else 0.0
            ),
            'N2': (
                float(10 ** np.random.uniform(low=-2.0, high=1.0))
                if float(target.get('N', 1.0)) > 0
                else 0.0
            ),
            'S2': (
                float(10 ** np.random.uniform(low=-2.0, high=1.0))
                if float(target.get('S', 1.0)) > 0
                else 0.0
            ),
        }
    return {
        'H2O': 10.0 if float(target.get('H', 1.0)) > 0 else 0.0,
        'CO2': 100.0 if float(target.get('C', 1.0)) > 0 else 0.0,
        'N2': 1.0 if float(target.get('N', 1.0)) > 0 else 0.0,
        'S2': 1.0 if float(target.get('S', 1.0)) > 0 else 0.0,
    }


def test_element_balance_across_75_case_grid():
    """Verify species-reconstructed elemental mass conservation across all 75 grid cases."""
    cases = _build_grid_cases()
    assert len(cases) == 75
    orig_fsolve = solve_mod.opt.fsolve

    for case in cases:
        case_id = case['case_id']
        ddict = _build_ddict(
            m_mantle=case['M_mantle'],
            g_surf=case['gravity'],
            r_int=case['R_int'],
            phi_global=case['Phi_global'],
            t_magma=case['T_K'],
            fo2_shift=case['fO2_shift_IW'],
            active_nobles=case['nobles'],
        )
        target = dict(case['targets'])
        np.random.seed(42)
        if case['mode'] == 'authoritative_O':
            out = equilibrium_atmosphere_authoritative_O(
                target,
                ddict,
                fO2_hint=case['fO2_shift_IW'],
                rtol=RTOL_DERIVED,
                atol=ATOL_KG,
                xtol=1.0e-8,
                random_seed=42,
                print_result=False,
                opt_solver=False,
            )
            checked_elements = ('H', 'C', 'N', 'S', 'O')
        else:
            p_guess = _construct_zero_target_p_guess(target, randomize=False)
            out = equilibrium_atmosphere(
                target,
                ddict,
                p_guess=p_guess,
                rtol=RTOL_DERIVED,
                atol=ATOL_KG,
                xtol=1.0e-8,
                print_result=False,
                opt_solver=False,
            )
            checked_elements = ('H', 'C', 'N', 'S') + case['nobles']
            assert math.isfinite(out['O_kg_total']), f'{case_id}: O_kg_total not finite'
            assert out['O_kg_total'] > 0.0, f'{case_id}: O_kg_total <= 0'

        # Verify VMR normalization and hydrostatic total mass identity on every case.
        vmr_sum = sum(float(out.get(f'{sp}_vmr', 0.0)) for sp in SPECIES_STOICH)
        assert vmr_sum == pytest.approx(1.0, rel=1.0e-12, abs=0.0), (
            f'{case_id}: vmr_sum={vmr_sum}'
        )
        expected_m_atm = (
            1.0e5
            * 4.0
            * math.pi
            * (case['R_int'] ** 2)
            / case['gravity']
            * float(out['P_surf'])
        )
        assert float(out['M_atm']) == pytest.approx(expected_m_atm, rel=1.0e-12, abs=0.0), (
            f'{case_id}: M_atm={out["M_atm"]} != {expected_m_atm}'
        )

        rec_elements = (
            checked_elements if case['mode'] == 'authoritative_O' else checked_elements + ('O',)
        )
        rec = _reconstruct_from_species(
            out, r_int=case['R_int'], g_surf=case['gravity'], elements=rec_elements
        )
        for el in rec_elements:
            assert rec[el]['gas_kg'] == pytest.approx(
                out[f'{el}_kg_atm'], rel=1.0e-11, abs=1.0
            ), (
                f'{case_id} {el}: reconstructed gas {rec[el]["gas_kg"]!r} != {out[f"{el}_kg_atm"]!r}'
            )
            assert rec[el]['liq_kg'] == pytest.approx(
                out[f'{el}_kg_liquid'], rel=1.0e-11, abs=1.0
            ), (
                f'{case_id} {el}: reconstructed liq {rec[el]["liq_kg"]!r} != {out[f"{el}_kg_liquid"]!r}'
            )
            assert out[f'{el}_kg_total'] == pytest.approx(
                out[f'{el}_kg_atm'] + out[f'{el}_kg_liquid'], rel=1.0e-12, abs=1.0
            ), f'{case_id} {el}: total != atm + liquid'

        # Solver acceptance criterion: per-element for authoritative_O and nobles; scalar CHNOS max for buffered.
        hcns_max = max(float(target[k]) for k in ('H', 'C', 'N', 'S'))
        for el in checked_elements:
            t_kg = float(target[el])
            if case['mode'] == 'authoritative_O' or el in noble_gases:
                tol_kg = max(RTOL_DERIVED * t_kg, float(TRUNC_MASS))
            else:
                tol_kg = hcns_max * RTOL_DERIVED + ATOL_KG + float(TRUNC_MASS)
            rec_tot = rec[el]['total_kg']
            assert out[f'{el}_res'] == pytest.approx(
                out[f'{el}_kg_total'] - t_kg, rel=1.0e-11, abs=1.0
            ), f'{case_id} {el}: res mismatch'
            abs_err = abs(rec_tot - t_kg)
            assert abs_err <= tol_kg, (
                f'{case_id} {el}: |recovered - target| = {abs_err!r} > tol {tol_kg!r}'
            )

        # Verify multi-seed stability across 5 seeds on 2-element, 4-element, and 5-element cases.
        if case_id in (
            'budget_tutorial_earth_T2000_IW+2',
            'budget_cmb_singleimpact_T2000_IW+1',
            'authO_default_T2000',
        ):
            fsolve_first_x0s: list[tuple[float, ...]] = []
            for seed in (0, 1, 2, 3, 4):
                first_x0 = [None]

                def _spy_fsolve(func, x0, *a, **k):
                    if first_x0[0] is None:
                        first_x0[0] = tuple(np.round(np.asarray(x0, dtype=float), 12))
                    return orig_fsolve(func, x0, *a, **k)

                solve_mod.opt.fsolve = _spy_fsolve
                try:
                    if case['mode'] == 'authoritative_O':
                        np.random.seed(0)
                        out_s = equilibrium_atmosphere_authoritative_O(
                            target,
                            ddict,
                            fO2_hint=case['fO2_shift_IW'],
                            rtol=RTOL_DERIVED,
                            atol=ATOL_KG,
                            xtol=1.0e-8,
                            random_seed=seed,
                            print_result=False,
                            opt_solver=False,
                        )
                    else:
                        np.random.seed(seed)
                        p_guess = _construct_zero_target_p_guess(target, randomize=True)
                        out_s = equilibrium_atmosphere(
                            target,
                            ddict,
                            p_guess=p_guess,
                            rtol=RTOL_DERIVED,
                            atol=ATOL_KG,
                            xtol=1.0e-8,
                            print_result=False,
                            opt_solver=False,
                        )
                finally:
                    solve_mod.opt.fsolve = orig_fsolve
                fsolve_first_x0s.append(first_x0[0])

                rec_s = _reconstruct_from_species(
                    out_s,
                    r_int=case['R_int'],
                    g_surf=case['gravity'],
                    elements=checked_elements,
                )
                for el in checked_elements:
                    tol_s = (
                        max(RTOL_DERIVED * float(target[el]), float(TRUNC_MASS))
                        if case['mode'] == 'authoritative_O'
                        else hcns_max * RTOL_DERIVED + ATOL_KG + float(TRUNC_MASS)
                    )
                    err_s = abs(rec_s[el]['total_kg'] - float(target[el]))
                    assert err_s <= tol_s, f'{case_id} seed={seed} {el}: {err_s!r} > {tol_s!r}'

                for sp in (
                    'H2O',
                    'H2',
                    'O2',
                    'CO2',
                    'CO',
                    'CH4',
                    'N2',
                    'NH3',
                    'S2',
                    'SO2',
                    'H2S',
                ):
                    p_ref = float(out[f'{sp}_bar'])
                    p_s = float(out_s[f'{sp}_bar'])
                    req_els = {
                        'H2O': ('H',),
                        'H2': ('H',),
                        'O2': ('O',) if case['mode'] == 'authoritative_O' else (),
                        'CO2': ('C',),
                        'CO': ('C',),
                        'CH4': ('C', 'H'),
                        'N2': ('N',),
                        'NH3': ('N', 'H'),
                        'S2': ('S',),
                        'SO2': ('S',),
                        'H2S': ('S', 'H'),
                    }[sp]
                    if any(float(target.get(el, 1.0)) == 0.0 for el in req_els):
                        assert p_ref == 0.0, (
                            f'{case_id} {sp} p_ref={p_ref} != 0 for zero element target'
                        )
                        assert p_s == 0.0, (
                            f'{case_id} seed={seed} {sp} p_s={p_s} != 0 for zero target'
                        )
                    else:
                        # 1e-10 bar floor: 19x-25x below atol mass resolution (Earth 1.945e-9,
                        # cmb 2.45e-9 bar, tolerance in eq_atm); authO_default_T2000 NH3
                        # at 5.13e-10 bar < N gate: rtol*target, ~1.5e-4 bar as m g/A.
                        p_floor = 1.0e-10  # bar
                        assert (p_ref > p_floor) == (p_s > p_floor), (
                            f'{case_id} seed={seed} {sp} p_bar floor mismatch: '
                            f'ref={p_ref:.4e} bar, seed={p_s:.4e} bar'
                        )
                        if p_ref > p_floor:
                            diff_dex = abs(math.log10(p_s / p_ref))
                            assert diff_dex <= 1.0e-4, (
                                f'{case_id} seed={seed} {sp} p_bar offset '
                                f'{diff_dex:.4e} dex > 1e-4 dex'
                            )
            assert len(fsolve_first_x0s) == 5, f'{case_id}: expected 5 seeded solves'
            assert len(set(fsolve_first_x0s)) == 5, (
                f'{case_id}: duplicate initial fsolve x0 across seeds {fsolve_first_x0s}'
            )

        if case_id == 'budget_tutorial_earth_T2000_IW+2':
            # Solver does not hold zero-target species below atol.
            gate_zero_kg = hcns_max * RTOL_DERIVED + ATOL_KG + float(TRUNC_MASS)
            max_n2_kg = 0.0
            fsolve_none_x0s: list[tuple[float, ...]] = []
            for seed_none in (0, 8, 32, 93, 138):
                first_x0_none = [None]

                def _spy_fsolve_none(func, x0, *a, **k):
                    if first_x0_none[0] is None:
                        first_x0_none[0] = tuple(np.round(np.asarray(x0, dtype=float), 12))
                    return orig_fsolve(func, x0, *a, **k)

                solve_mod.opt.fsolve = _spy_fsolve_none
                try:
                    np.random.seed(seed_none)
                    out_none = equilibrium_atmosphere(
                        target,
                        ddict,
                        p_guess=None,
                        rtol=RTOL_DERIVED,
                        atol=ATOL_KG,
                        xtol=1.0e-8,
                        print_result=False,
                        opt_solver=False,
                    )
                finally:
                    solve_mod.opt.fsolve = orig_fsolve
                fsolve_none_x0s.append(first_x0_none[0])

                for el_res in ('H', 'C', 'N', 'S'):
                    assert abs(out_none[f'{el_res}_res']) <= gate_zero_kg, (
                        f'{case_id} seed={seed_none} p_guess=None {el_res}_res '
                        f'{out_none[f"{el_res}_res"]!r} > gate {gate_zero_kg!r}'
                    )
                for sp in ('N2', 'NH3', 'S2', 'SO2', 'H2S'):
                    assert 0.0 <= float(out_none[f'{sp}_kg_total']), (
                        f'{case_id} seed={seed_none} p_guess=None {sp} mass negative'
                    )
                max_n2_kg = max(max_n2_kg, float(out_none['N2_kg_total']))
            assert len(fsolve_none_x0s) == 5, f'{case_id}: expected 5 seeds'
            assert len(set(fsolve_none_x0s)) == 5, (
                f'{case_id}: duplicate initial fsolve x0 across seeds {fsolve_none_x0s}'
            )
            # Pin is expected to change if the solver holds zero-target slots at 0.
            assert max_n2_kg > ATOL_KG, (
                f'{case_id} max N2 {max_n2_kg!r} not > ATOL_KG {ATOL_KG!r}'
            )


def test_species_to_element_summation_and_residual_vectors():
    """Verify atmosphere_mass, dissolved_mass, func, and func_authoritative_O against stoichiometric sums."""
    nobles = ('He', 'Ne', 'Ar', 'Kr', 'Xe')
    ddict = _build_ddict(
        m_mantle=4.0311e24,
        g_surf=9.81,
        r_int=6.335439e6,
        phi_global=1.0,
        t_magma=2200.0,
        fo2_shift=-1.5,
        active_nobles=nobles,
    )
    pin = {
        'H2O': 25.0,
        'CO2': 40.0,
        'N2': 5.0,
        'S2': 8.0,
        'He': 0.5,
        'Ne': 0.4,
        'Ar': 0.3,
        'Kr': 0.2,
        'Xe': 0.1,
    }
    m_atm = atmosphere_mass(pin, ddict)
    m_int = dissolved_mass(pin, ddict)

    all_elements = ('H', 'C', 'N', 'S', 'O') + nobles
    expected_atm = {el: 0.0 for el in all_elements}
    expected_int = {el: 0.0 for el in all_elements}
    for sp, stoich in SPECIES_STOICH.items():
        m_sp = float(molar_mass[sp])
        sp_atm = float(m_atm.get(sp, 0.0))
        sp_int = float(m_int.get(sp, 0.0))
        for el, n_atom in stoich.items():
            w_frac_gas = 1.0 if sp in noble_gases else (n_atom * float(molar_mass[el])) / m_sp
            w_frac_liq = (
                1.0
                if sp in ('N2', 'S2', *noble_gases)
                else (n_atom * float(molar_mass[el])) / m_sp
            )
            expected_atm[el] += sp_atm * w_frac_gas
            expected_int[el] += sp_int * w_frac_liq

    for el in all_elements:
        assert m_atm[el] > 0.0, f'm_atm[{el}] <= 0'
        assert m_atm[el] == pytest.approx(expected_atm[el], rel=1.0e-12, abs=0.0), (
            f'm_atm[{el}] mismatch'
        )
        assert el in m_int, f'{el} missing from m_int'
        assert m_int[el] > 0.0, f'm_int[{el}] <= 0'
        assert m_int[el] == pytest.approx(expected_int[el], rel=1.0e-12, abs=0.0), (
            f'm_int[{el}] mismatch'
        )

    # Verify secondary gas species contributions at fixed mean molar mass equal their stoichiometric sums.
    primary_sp = ('H2O', 'CO2', 'N2', 'S2', *nobles)
    secondary_sp = ('H2', 'CH4', 'H2S', 'NH3', 'CO', 'SO2', 'O2')
    for el in ('H', 'C', 'N', 'S', 'O'):
        prim_sum = sum(
            float(m_atm[sp])
            * (
                1.0
                if sp in noble_gases
                else (SPECIES_STOICH[sp].get(el, 0) * float(molar_mass[el]))
                / float(molar_mass[sp])
            )
            for sp in primary_sp
            if el in SPECIES_STOICH[sp]
        )
        sec_sum = sum(
            float(m_atm[sp])
            * (SPECIES_STOICH[sp].get(el, 0) * float(molar_mass[el]))
            / float(molar_mass[sp])
            for sp in secondary_sp
            if el in SPECIES_STOICH[sp]
        )
        assert sec_sum > 0.0, f'secondary species contribution to {el} must be positive'
        assert float(m_atm[el]) - prim_sum == pytest.approx(sec_sum, rel=1.0e-11, abs=0.0), (
            f'secondary species stoichiometric sum mismatch for {el}'
        )

    # Verify O2 contributes to m_atm['O'] at oxidizing fO2 where pO2 is large.
    ddict_ox = dict(ddict, fO2_shift_IW=10.0)
    pin_ox = dict(pin, S2=0.0)
    m_atm_ox = atmosphere_mass(pin_ox, ddict_ox)
    o_without_o2 = (
        m_atm_ox['H2O'] / molar_mass['H2O']
        + 2.0 * m_atm_ox['CO2'] / molar_mass['CO2']
        + m_atm_ox['CO'] / molar_mass['CO']
        + 2.0 * m_atm_ox['SO2'] / molar_mass['SO2']
    ) * molar_mass['O']
    assert m_atm_ox['O'] > o_without_o2 * 1.001
    assert m_atm_ox['O'] - o_without_o2 == pytest.approx(
        m_atm_ox['O2'] * (2.0 * molar_mass['O'] / molar_mass['O2']), rel=1.0e-12, abs=0.0
    )

    # Verify func and func_authoritative_O residual vectors equal (m_atm + m_int - target).
    target = {el: 0.75 * (expected_atm[el] + expected_int[el]) for el in all_elements}
    x_buf = [pin['H2O'], pin['CO2'], pin['N2'], pin['S2']] + [pin[ng] for ng in nobles]
    res_buf = func(x_buf, ddict, target)
    buf_elems = ('H', 'C', 'N', 'S') + nobles
    assert len(res_buf) == len(buf_elems)
    for idx, el in enumerate(buf_elems):
        expected_res = expected_atm[el] + expected_int[el] - target[el]
        assert res_buf[idx] == pytest.approx(expected_res, rel=1.0e-12, abs=0.0), (
            f'res_buf[{el}] mismatch'
        )
        assert abs(res_buf[idx]) > TRUNC_MASS

    # Use a distinct fO2_shift (+3.0 != ddict['fO2_shift_IW'] = -1.5) to verify x_auth[4] is used.
    fo2_auth = 3.0
    ddict_auth_ref = dict(ddict, fO2_shift_IW=fo2_auth)
    m_atm_auth = atmosphere_mass(pin, ddict_auth_ref)
    m_int_auth = dissolved_mass(pin, ddict_auth_ref)
    x_auth = [pin['H2O'], pin['CO2'], pin['N2'], pin['S2'], fo2_auth] + [
        pin[ng] for ng in nobles
    ]
    res_auth = func_authoritative_O(x_auth, ddict, target)
    auth_elems = ('H', 'C', 'N', 'S', 'O') + nobles
    assert len(res_auth) == len(auth_elems)
    assert res_auth[0] != pytest.approx(res_buf[0], rel=1.0e-3, abs=0.0)
    for idx, el in enumerate(auth_elems):
        expected_res_auth = float(m_atm_auth[el]) + float(m_int_auth[el]) - target[el]
        assert res_auth[idx] == pytest.approx(expected_res_auth, rel=1.0e-12, abs=0.0), (
            f'res_auth[{el}] mismatch at fO2_shift={fo2_auth}'
        )
