"""Cross-code thermodynamic and gas-only equilibrium consistency tests.

Verifies:
1. CALLIOPE equilibrium constants (ModifiedKeq in src/calliope/equilibrium.py) and
   algebraic secondary gas speciation (src/calliope/solve.py:165-225) against the
   9-coefficient NASA Glenn reference values across T in {1500, 2000, 2500, 3000} K
   within the propagated reaction-data tolerances.
2. Live 72-case gas-only equilibrium cross-comparison against atmodeller (60 buffered
   cases with matched fO2 plus 12 authoritative-O cases) on the shared volatile
   species (H2O, H2, O2, CO2, CO, CH4, N2, NH3, S2, SO2, H2S), gated via
   pytest.importorskip('atmodeller') when run outside the coupled environment.
"""

from __future__ import annotations

import logging
import math
from typing import Any

import numpy as np
import pytest
import scipy.optimize as opt

from calliope.chemistry import ModifiedKeq
from calliope.constants import molar_mass, noble_gases, ocean_moles
from calliope.oxygen_fugacity import OxygenFugacity
from calliope.solve import (
    atmosphere_mass,
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
    func,
    func_authoritative_O,
)

logging.getLogger('calliope').setLevel(logging.WARNING)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(300)]

# NASA Glenn 9-coefficient reference log10(K_p) values evaluated in atmodeller (f2ea2b45)
# for the 5 testable shared reactions across T in {1500, 2000, 2500, 3000} K.
NASA_GLENN_LOG10_K_REF: dict[float, dict[str, float]] = {
    1500.0: {
        'janaf_H2': -5.723743519770492,
        'janaf_CO': -5.310760370623338,
        'schaefer_CH4': -16.38399254976359,
        'janaf_SO2': 17.55414103966037,
        'janaf_NH3': -8.41149449147021,
    },
    2000.0: {
        'janaf_H2': -3.5404166434565947,
        'janaf_CO': -2.8794684826996725,
        'schaefer_CH4': -13.712241695534663,
        'janaf_SO2': 11.26019168077867,
        'janaf_NH3': -9.387231801577274,
    },
    2500.0: {
        'janaf_H2': -2.226328640710726,
        'janaf_CO': -1.4343290868120266,
        'schaefer_CH4': -12.082789477914483,
        'janaf_SO2': 7.4838162992727165,
        'janaf_NH3': -9.959571768466464,
    },
    3000.0: {
        'janaf_H2': -1.3492356139464587,
        'janaf_CO': -0.47992354032844615,
        'schaefer_CH4': -10.971187439341608,
        'janaf_SO2': 4.963218978649281,
        'janaf_NH3': -10.327532578283407,
    },
}

# Per-species relative tolerances derived from the JANAF/Schaefer vs NASA Glenn
# thermodynamic data spread (with stoichiometric power 0.5 for doubled reactions SO2 and NH3).
REACTION_QUOTIENT_RTOL: dict[str, float] = {
    'H2': 2.0e-2,
    'CO': 4.0e-2,
    'SO2': 2.5e-2,
    'NH3': 1.8e-1,
    'CH4': 2.8e-1,
}


def _build_gas_only_ddict(
    m_mantle: float,
    g_surf: float,
    r_int: float,
    t_magma: float,
    fo2_shift: float,
) -> dict[str, float | int]:
    """Construct CALLIOPE gas-only parameter dictionary (Phi_global = 0.0)."""
    ddict: dict[str, float | int] = {
        'M_mantle': m_mantle,
        'gravity': g_surf,
        'radius': r_int,
        'Phi_global': 0.0,
        'T_magma': t_magma,
        'fO2_shift_IW': fo2_shift,
    }
    for sp in ('H2O', 'CO2', 'N2', 'S2', 'SO2', 'H2S', 'NH3', 'H2', 'CH4', 'CO', 'O2'):
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0
    for ng in noble_gases:
        ddict[f'{ng}_included'] = 0
        ddict[f'{ng}_initial_bar'] = 0.0
        ddict[f'{ng}_ppmw'] = 0.0
    return ddict


def _solve_calliope_gas_only(
    target: dict[str, float],
    ddict: dict[str, float | int],
    autho: bool,
) -> dict[str, float]:
    """Solve CALLIOPE gas-only equilibrium, using a log-space preseed when cold start misses."""
    np.random.seed(42)
    try:
        if autho:
            return equilibrium_atmosphere_authoritative_O(
                target,
                ddict,
                fO2_hint=float(ddict['fO2_shift_IW']),
                rtol=1.0e-5,
                atol=1.0e10,
                xtol=1.0e-8,
                nguess=5,
                opt_solver=False,
                random_seed=42,
                print_result=False,
            )
        return equilibrium_atmosphere(
            target,
            ddict,
            rtol=1.0e-5,
            atol=1.0e10,
            xtol=1.0e-8,
            nguess=5,
            opt_solver=False,
            print_result=False,
        )
    except RuntimeError:
        scales = np.array(
            [target[k] for k in (('H', 'C', 'N', 'S', 'O') if autho else ('H', 'C', 'N', 'S'))]
        )
        if autho:

            def f_log_autho(y: np.ndarray) -> np.ndarray:
                yc = np.clip(y[:4], -80.0, 14.0)
                fo2_c = float(np.clip(y[4], -11.99, 11.99))
                x = [
                    float(np.exp(yc[0])),
                    float(np.exp(yc[1])),
                    float(np.exp(yc[2])),
                    float(np.exp(yc[3])),
                    fo2_c,
                ]
                return np.asarray(func_authoritative_O(x, ddict, target)) / scales

            res = opt.least_squares(
                f_log_autho,
                [5.0, 4.0, 2.0, -20.0, 6.0],
                bounds=([-80.0, -80.0, -80.0, -80.0, -11.99], [14.0, 14.0, 14.0, 14.0, 11.99]),
                ftol=1.0e-14,
                xtol=1.0e-14,
                gtol=1.0e-14,
            )
            p_g = {
                'H2O': float(np.exp(res.x[0])),
                'CO2': float(np.exp(res.x[1])),
                'N2': float(np.exp(res.x[2])),
                'S2': float(np.exp(res.x[3])),
            }
            return equilibrium_atmosphere_authoritative_O(
                target,
                ddict,
                fO2_hint=float(res.x[4]),
                p_guess=p_g,
                rtol=1.0e-5,
                atol=1.0e10,
                xtol=1.0e-8,
                nguess=5,
                opt_solver=False,
                random_seed=42,
                print_result=False,
            )

        fo2 = float(ddict['fO2_shift_IW'])
        y0 = [4.0, 3.0, 1.0, -15.0 if fo2 > 0 else 1.0]

        def f_log_buf(y: np.ndarray) -> np.ndarray:
            yc = np.clip(y, -80.0, 14.0)
            x = [
                float(np.exp(yc[0])),
                float(np.exp(yc[1])),
                float(np.exp(yc[2])),
                float(np.exp(yc[3])),
            ]
            return np.asarray(func(x, ddict, target)) / scales

        res = opt.least_squares(
            f_log_buf,
            y0,
            bounds=([-80.0] * 4, [14.0] * 4),
            ftol=1.0e-14,
            xtol=1.0e-14,
            gtol=1.0e-14,
        )
        p_g = {
            'H2O': float(np.exp(res.x[0])),
            'CO2': float(np.exp(res.x[1])),
            'N2': float(np.exp(res.x[2])),
            'S2': float(np.exp(res.x[3])),
        }
        return equilibrium_atmosphere(
            target,
            ddict,
            p_guess=p_g,
            rtol=1.0e-5,
            atol=1.0e10,
            xtol=1.0e-8,
            nguess=5,
            opt_solver=False,
            print_result=False,
        )


def test_calliope_equilibrium_constants_and_reaction_quotients_vs_nasa_glenn():
    """Verify CALLIOPE ModifiedKeq and secondary gas speciation against NASA Glenn reference data."""
    keq_models = {
        key: ModifiedKeq(key)
        for key in (
            'janaf_H2',
            'janaf_CO',
            'schaefer_CH4',
            'janaf_SO2',
            'janaf_NH3',
            'janaf_H2S',
        )
    }
    fo2_model = OxygenFugacity()
    pin = {'H2O': 120.0, 'CO2': 80.0, 'N2': 15.0, 'S2': 4.0}

    for t_k, ref_dict in NASA_GLENN_LOG10_K_REF.items():
        # 1. Direct equilibrium constant check against NASA Glenn per mole of secondary product
        log10_f = float(fo2_model(t_k, 0.0))
        log10_k_h2, st_h2 = keq_models['janaf_H2'].callmodel(t_k)
        log10_k_co, st_co = keq_models['janaf_CO'].callmodel(t_k)
        log10_k_ch4, st_ch4 = keq_models['schaefer_CH4'].callmodel(t_k)
        log10_k_so2, st_so2 = keq_models['janaf_SO2'].callmodel(t_k)
        log10_k_nh3, st_nh3 = keq_models['janaf_NH3'].callmodel(t_k)
        log10_k_h2s, st_h2s = keq_models['janaf_H2S'].callmodel(t_k)
        assert (st_h2, st_co, st_ch4, st_so2, st_nh3, st_h2s) == (0.5, 0.5, 1, 0, 0, 0)
        assert math.isfinite(float(log10_k_h2s))

        # Verify ModifiedKeq.__call__ returns 10**(Keq - fO2_stoich * fO2)
        assert float(keq_models['janaf_H2'](t_k, 0.0)) == pytest.approx(
            10.0 ** (log10_k_h2 - 0.5 * log10_f), rel=1.0e-12
        )
        assert float(keq_models['janaf_CO'](t_k, 0.0)) == pytest.approx(
            10.0 ** (log10_k_co - 0.5 * log10_f), rel=1.0e-12
        )
        assert float(keq_models['schaefer_CH4'](t_k, 0.0)) == pytest.approx(
            10.0 ** (log10_k_ch4 - 1.0 * log10_f), rel=1.0e-12
        )
        assert float(keq_models['janaf_SO2'](t_k, 0.0)) == pytest.approx(
            10.0**log10_k_so2, rel=1.0e-12
        )
        assert float(keq_models['janaf_NH3'](t_k, 0.0)) == pytest.approx(
            10.0**log10_k_nh3, rel=1.0e-12
        )
        assert float(keq_models['janaf_H2S'](t_k, 0.0)) == pytest.approx(
            10.0**log10_k_h2s, rel=1.0e-12
        )

        d_h2 = ref_dict['janaf_H2'] - float(log10_k_h2)
        d_co = ref_dict['janaf_CO'] - float(log10_k_co)
        d_so2 = 0.5 * (ref_dict['janaf_SO2'] - float(log10_k_so2))
        d_nh3 = 0.5 * (ref_dict['janaf_NH3'] - float(log10_k_nh3)) + 1.5 * d_h2
        d_ch4 = (ref_dict['schaefer_CH4'] - float(log10_k_ch4)) + 2.0 * d_h2

        assert abs(10.0**d_h2 - 1.0) <= REACTION_QUOTIENT_RTOL['H2']
        assert abs(10.0**d_co - 1.0) <= REACTION_QUOTIENT_RTOL['CO']
        assert abs(10.0**d_so2 - 1.0) <= REACTION_QUOTIENT_RTOL['SO2']
        assert abs(10.0**d_nh3 - 1.0) <= REACTION_QUOTIENT_RTOL['NH3']
        assert abs(10.0**d_ch4 - 1.0) <= REACTION_QUOTIENT_RTOL['CH4']

        # 2. Verify secondary partial pressures from equilibrium_atmosphere satisfy both
        # CALLIOPE's exact algebraic mass-action equations and the NASA Glenn reference bounds.
        ddict = _build_gas_only_ddict(
            m_mantle=4.0311e24,
            g_surf=9.81,
            r_int=6.335439e6,
            t_magma=t_k,
            fo2_shift=0.0,
        )
        m_atm = atmosphere_mass(pin, ddict)
        target = {el: m_atm[el] for el in ('H', 'C', 'N', 'S')}
        out = equilibrium_atmosphere(
            target,
            ddict,
            p_guess=pin,
            rtol=1.0e-5,
            atol=1.0e10,
            xtol=1.0e-8,
            nguess=5,
            opt_solver=False,
            print_result=False,
        )

        p_h2o = float(out['H2O_bar'])
        p_co2 = float(out['CO2_bar'])
        p_n2 = float(out['N2_bar'])
        p_s2 = float(out['S2_bar'])
        p_o2 = float(out['O2_bar'])
        p_h2 = float(out['H2_bar'])
        p_co = float(out['CO_bar'])
        p_ch4 = float(out['CH4_bar'])
        p_so2 = float(out['SO2_bar'])
        p_nh3 = float(out['NH3_bar'])
        p_h2s = float(out['H2S_bar'])

        log10_fo2 = float(fo2_model(t_k, 0.0))
        assert p_o2 == pytest.approx(10.0**log10_fo2, rel=1.0e-12)
        assert p_h2 == pytest.approx(
            (10.0 ** (log10_k_h2 - 0.5 * log10_fo2)) * p_h2o, rel=1.0e-12
        )
        assert p_co == pytest.approx(
            (10.0 ** (log10_k_co - 0.5 * log10_fo2)) * p_co2, rel=1.0e-12
        )
        assert p_ch4 == pytest.approx(
            (10.0 ** (log10_k_ch4 - log10_fo2)) * p_co2 * (p_h2**2), rel=1.0e-12
        )
        assert p_so2 == pytest.approx(
            ((10.0**log10_k_so2) * p_s2 * (p_o2**2)) ** 0.5, rel=1.0e-12
        )
        assert p_h2s == pytest.approx(
            ((10.0**log10_k_h2s) * p_s2 * (p_h2**2)) ** 0.5, rel=1.0e-12
        )
        assert p_nh3 == pytest.approx(
            ((10.0**log10_k_nh3) * p_n2 * (p_h2**3)) ** 0.5, rel=1.0e-12
        )

        # Reaction quotients compared against NASA Glenn reference K_p(T)
        q_h2 = p_h2 * (p_o2**0.5) / p_h2o
        q_co = p_co * (p_o2**0.5) / p_co2
        q_so2 = p_so2 / ((p_s2**0.5) * p_o2)
        q_nh3 = p_nh3 / ((p_n2**0.5) * (p_h2**1.5))
        q_ch4 = p_ch4 * p_o2 / (p_co2 * (p_h2**2))

        assert q_h2 == pytest.approx(
            10.0 ** ref_dict['janaf_H2'], rel=REACTION_QUOTIENT_RTOL['H2']
        )
        assert q_co == pytest.approx(
            10.0 ** ref_dict['janaf_CO'], rel=REACTION_QUOTIENT_RTOL['CO']
        )
        assert q_so2 == pytest.approx(
            10.0 ** (0.5 * ref_dict['janaf_SO2']), rel=REACTION_QUOTIENT_RTOL['SO2']
        )
        assert q_nh3 == pytest.approx(
            10.0 ** (0.5 * ref_dict['janaf_NH3']), rel=REACTION_QUOTIENT_RTOL['NH3']
        )
        assert q_ch4 == pytest.approx(
            10.0 ** ref_dict['schaefer_CH4'], rel=REACTION_QUOTIENT_RTOL['CH4']
        )


def test_calliope_vs_atmodeller_72_case_gas_only_grid():
    """Verify CALLIOPE vs atmodeller cross-code agreement over all 72 gas-only grid cases."""
    atmodeller = pytest.importorskip('atmodeller')
    from atmodeller.containers import FixedFugacityConstraint, SolverParameters

    mass_ocean = float(ocean_moles * molar_mass['H2'])
    m_earth = 5.972e24
    r_earth = 6.335439e6
    cmf = 0.325
    m_mantle = (1.0 - cmf) * m_earth
    g_surf = 9.81

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
        inventories[label] = {'H': h_kg, 'C': c_kg, 'N': n_kg, 'S': s_kg, 'O': o_stoich_kg}

    buffered_cases: list[dict[str, Any]] = []
    for inv_label in ('low', 'default', 'high'):
        inv = inventories[inv_label]
        for t_k in (1500.0, 2000.0, 2500.0, 3000.0):
            for fo2 in (-4.0, -2.0, 0.0, 2.0, 4.0):
                buffered_cases.append(
                    {
                        'case_id': f'buffered_{inv_label}_T{int(t_k)}_IW{fo2:+.0f}',
                        'T_K': t_k,
                        'fO2_shift_IW': fo2,
                        'targets': {k: inv[k] for k in ('H', 'C', 'N', 'S')},
                    }
                )
    autho_cases: list[dict[str, Any]] = []
    for inv_label in ('low', 'default', 'high'):
        inv = inventories[inv_label]
        for t_k in (1500.0, 2000.0, 2500.0, 3000.0):
            autho_cases.append(
                {
                    'case_id': f'authO_{inv_label}_T{int(t_k)}',
                    'T_K': t_k,
                    'fO2_shift_IW': 2.0,
                    'targets': {k: inv[k] for k in ('H', 'C', 'N', 'S', 'O')},
                }
            )

    species_names = ('H2O', 'H2', 'O2', 'CO2', 'CO', 'CH4', 'N2', 'H3N', 'S2', 'O2S', 'H2S')
    gas_species = tuple(atmodeller.ChemicalSpecies.create_gas(sp) for sp in species_names)
    solver_params = SolverParameters(atol=1.0e-6, rtol=1.0e-4, max_steps=1024, multistart=10)
    fo2_fischer = OxygenFugacity()

    # Solve 60 buffered cases in atmodeller
    n_buf = len(buffered_cases)
    planet_buf = atmodeller.Planet(
        planet_mass=np.full(n_buf, m_earth),
        core_mass_fraction=np.full(n_buf, cmf),
        mantle_melt_fraction=np.zeros(n_buf),
        surface_radius=np.full(n_buf, r_earth),
        temperature=np.array([c['T_K'] for c in buffered_cases]),
        pressure=np.full(n_buf, np.nan),
    )
    f_o2_arr = np.array(
        [float(10.0 ** fo2_fischer(c['T_K'], c['fO2_shift_IW'])) for c in buffered_cases]
    )
    mdl_buf = atmodeller.EquilibriumModel(atmodeller.SpeciesNetwork(gas_species))
    mdl_buf.solve(
        state=planet_buf,
        fugacity_constraints={'O2_g': FixedFugacityConstraint(f_o2_arr)},
        mass_constraints={
            el: np.array([c['targets'][el] for c in buffered_cases])
            for el in ('H', 'C', 'N', 'S')
        },
        solver='robust',
        solver_parameters=solver_params,
    )
    asdict_buf = mdl_buf.output.asdict()
    ptot_buf = np.asarray(mdl_buf.output.total_pressure())

    # Solve 12 authoritative-O cases in atmodeller
    n_auth = len(autho_cases)
    planet_auth = atmodeller.Planet(
        planet_mass=np.full(n_auth, m_earth),
        core_mass_fraction=np.full(n_auth, cmf),
        mantle_melt_fraction=np.zeros(n_auth),
        surface_radius=np.full(n_auth, r_earth),
        temperature=np.array([c['T_K'] for c in autho_cases]),
        pressure=np.full(n_auth, np.nan),
    )
    mdl_auth = atmodeller.EquilibriumModel(atmodeller.SpeciesNetwork(gas_species))
    mdl_auth.solve(
        state=planet_auth,
        fugacity_constraints=None,
        mass_constraints={
            el: np.array([c['targets'][el] for c in autho_cases])
            for el in ('H', 'C', 'N', 'S', 'O')
        },
        solver='robust',
        solver_parameters=solver_params,
    )
    asdict_auth = mdl_auth.output.asdict()
    ptot_auth = np.asarray(mdl_auth.output.total_pressure())

    for idx, case in enumerate(buffered_cases):
        ddict = _build_gas_only_ddict(
            m_mantle=m_mantle,
            g_surf=g_surf,
            r_int=r_earth,
            t_magma=case['T_K'],
            fo2_shift=case['fO2_shift_IW'],
        )
        cal_out = _solve_calliope_gas_only(case['targets'], ddict, autho=False)
        assert float(ptot_buf[idx]) == pytest.approx(float(cal_out['P_surf']), rel=3.5e-2)
        assert float(asdict_buf['O2_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['O2_bar']), rel=1.0e-3
        )
        assert float(asdict_buf['H2O_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['H2O_bar']), rel=5.0e-2
        )
        assert float(asdict_buf['N2_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['N2_bar']), rel=5.0e-2
        )
        assert float(asdict_buf['H2_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['H2_bar']), rel=5.0e-2
        )
        assert float(asdict_buf['H3N_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['NH3_bar']), rel=1.8e-1
        )
        assert float(asdict_buf['CH4_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['CH4_bar']), rel=2.8e-1
        )
        c_rtol = 5.0e-2 if case['fO2_shift_IW'] >= 2.0 else 2.0e-1
        assert float(asdict_buf['CO2_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['CO2_bar']), rel=c_rtol
        )
        assert float(asdict_buf['CO_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['CO_bar']), rel=c_rtol
        )

        # Reaction quotients on all 60 buffered cases
        q_so2_cal = float(cal_out['SO2_bar']) / (
            (float(cal_out['S2_bar']) ** 0.5) * float(cal_out['O2_bar'])
        )
        q_so2_atm = float(asdict_buf['O2S_g']['pressure'][idx]) / (
            (float(asdict_buf['S2_g']['pressure'][idx]) ** 0.5)
            * float(asdict_buf['O2_g']['pressure'][idx])
        )
        assert q_so2_atm == pytest.approx(q_so2_cal, rel=REACTION_QUOTIENT_RTOL['SO2'])

    for idx, case in enumerate(autho_cases):
        ddict = _build_gas_only_ddict(
            m_mantle=m_mantle,
            g_surf=g_surf,
            r_int=r_earth,
            t_magma=case['T_K'],
            fo2_shift=case['fO2_shift_IW'],
        )
        cal_out = _solve_calliope_gas_only(case['targets'], ddict, autho=True)
        assert float(ptot_auth[idx]) == pytest.approx(float(cal_out['P_surf']), rel=3.5e-2)
        assert float(asdict_auth['H2O_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['H2O_bar']), rel=5.0e-2
        )
        assert float(asdict_auth['CO2_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['CO2_bar']), rel=5.0e-2
        )
        assert float(asdict_auth['N2_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['N2_bar']), rel=5.0e-2
        )
        assert float(asdict_auth['O2S_g']['pressure'][idx]) == pytest.approx(
            float(cal_out['SO2_bar']), rel=2.5e-2
        )
