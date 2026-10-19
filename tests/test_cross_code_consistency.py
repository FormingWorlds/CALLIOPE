"""Cross-code thermodynamic and gas-only equilibrium consistency tests.

Verifies:
1. CALLIOPE equilibrium constants (ModifiedKeq) and algebraic secondary gas
   speciation against the 9-coefficient NASA Glenn reference values across
   T in {1500, 2000, 2500, 3000} K within asymmetric log10 reaction-data
   bounds, while reporting the H2S equilibrium-constant offsets.
2. Live 72-case gas-only equilibrium cross-comparison against atmodeller (60
   buffered cases with matched fO2 plus 12 authoritative-O cases) on the shared
   volatile species (H2O, H2, O2, CO2, CO, CH4, N2, NH3, S2, SO2, H2S), skipped
   via pytest.importorskip('atmodeller') when run outside the coupled
   environment. At the stoichiometric H2O + CO2 + SO2 oxygen endpoint, the
   codes agree on solved log10(fO2) within 0.05 dex at T >= 2000 K, whereas at
   1500 K the ill-conditioned titration endpoint amplifies H2S/S2 and CH4
   constant differences into a log10(fO2) difference of up to 1.15 dex (-0.30,
   -0.60, and -1.15 dex for atmodeller minus CALLIOPE across low, default, and
   high inventories).
"""

from __future__ import annotations

import logging
import math
from typing import Any

import numpy as np
import pytest
import scipy.optimize as opt

from calliope.chemistry import ModifiedKeq
from calliope.constants import molar_mass, noble_gases, ocean_moles, volatile_species
from calliope.oxygen_fugacity import OxygenFugacity
from calliope.solve import (
    TRUNC_MASS,
    atmosphere_mass,
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
    func,
    func_authoritative_O,
)

logging.getLogger('fwl.calliope').setLevel(logging.WARNING)

pytestmark = [
    pytest.mark.integration,
    pytest.mark.reference_pinned,
    pytest.mark.timeout(180),
]

# NASA Glenn 9-coefficient reference log10(K_p) values from atmodeller 1.0.2
# across T in {1500, 2000, 2500, 3000} K (rounded to 6 decimal places).
NASA_GLENN_LOG10_K_REF: dict[float, dict[str, float]] = {
    1500.0: {
        'janaf_H2': -5.723739,
        'janaf_CO': -5.310764,
        'schaefer_CH4': -16.383992,
        'janaf_SO2': 17.554142,
        'janaf_NH3': -8.411493,
        'janaf_H2S': 1.154363,
    },
    2000.0: {
        'janaf_H2': -3.540416,
        'janaf_CO': -2.879472,
        'schaefer_CH4': -13.712241,
        'janaf_SO2': 11.260191,
        'janaf_NH3': -9.387227,
        'janaf_H2S': -0.421587,
    },
    2500.0: {
        'janaf_H2': -2.226332,
        'janaf_CO': -1.434333,
        'schaefer_CH4': -12.082791,
        'janaf_SO2': 7.483816,
        'janaf_NH3': -9.959570,
        'janaf_H2S': -1.363303,
    },
    3000.0: {
        'janaf_H2': -1.349239,
        'janaf_CO': -0.479923,
        'schaefer_CH4': -10.971185,
        'janaf_SO2': 4.963220,
        'janaf_NH3': -10.327529,
        'janaf_H2S': -1.987895,
    },
}

# Asymmetric [lo_dex, hi_dex] bounds for log10(K_CALLIOPE / K_NASA_Glenn) per mole of product:
# measured signed range across 1500..3000 K plus a 0.003 to 0.010 dex margin on each side.
K_OFFSET_BOUNDS_DEX: dict[str, tuple[float, float]] = {
    'H2': (-0.011, 0.009),
    'CO': (-0.011, 0.019),
    'SO2': (0.003, 0.011),
    'NH3': (-0.021, 0.070),
    'CH4': (0.049, 0.111),
}

# Asymmetric [lo_dex, hi_dex] bounds for log10(P_CALLIOPE / P_atmodeller) on the 60 buffered cases:
# centred on the measured offset (gravity -0.0053 dex plus propagated Delta log10 K) plus ~0.008 dex margin.
BUFFERED_LOG10_BOUNDS_DEX: dict[str, tuple[float, float]] = {
    'P_surf': (-0.015, 0.005),
    'H2O': (-0.020, 0.006),
    'N2': (-0.012, 0.001),
    'H2': (-0.022, 0.010),
    'NH3': (-0.048, 0.073),
    'CH4_ox': (0.010, 0.116),
    'CH4_red': (-0.012, 0.110),
    'CO2_ox': (-0.015, 0.003),
    'CO2_red': (-0.075, 0.005),
    'CO_ox': (-0.016, 0.012),
    'CO_red': (-0.080, 0.004),
    'SO2': (-0.300, 0.390),
}

# Asymmetric [lo_dex, hi_dex] bounds for log10(P_CALLIOPE / P_atmodeller) on authoritative-O cases:
# major species reflect the -0.0053 dex gravity offset; T >= 2000 K redox species reflect Delta log10 K.
AUTHO_LOG10_BOUNDS_DEX: dict[str, tuple[float, float]] = {
    'P_surf': (-0.008, -0.002),
    'H2O': (-0.008, -0.002),
    'CO2': (-0.009, -0.002),
    'N2': (-0.008, -0.002),
    'SO2': (-0.008, -0.002),
    'O2_hiT': (-0.010, 0.050),
    'O2_1500K': (0.250, 1.200),
    'H2_hiT': (-0.030, 0.008),
    'CO_hiT': (-0.040, 0.008),
    'CH4_hiT': (0.001, 0.105),
    'NH3_hiT': (-0.012, 0.062),
}


def _assert_in_dex_bounds(
    actual_dex: float,
    bounds: tuple[float, float],
    label: str,
) -> None:
    """Assert that a signed log10 ratio [dex] lies inside the asymmetric interval [lo, hi]."""
    lo_dex, hi_dex = bounds
    assert lo_dex <= actual_dex <= hi_dex, (
        f'{label}: log10 offset {actual_dex:+.6f} dex outside [{lo_dex:+.4f}, {hi_dex:+.4f}] dex'
    )


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
    for sp in volatile_species:
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
    kw: dict[str, Any] = {
        'rtol': 1.0e-5,
        'atol': 1.0e10,
        'xtol': 1.0e-8,
        'nguess': 5,
        'opt_solver': False,
        'print_result': False,
    }
    fo2 = float(ddict['fO2_shift_IW'])
    try:
        if autho:
            return equilibrium_atmosphere_authoritative_O(
                target, ddict, fO2_hint=fo2, random_seed=42, **kw
            )
        return equilibrium_atmosphere(target, ddict, **kw)
    except RuntimeError:
        pass

    scales = np.array(
        [target[k] for k in (('H', 'C', 'N', 'S', 'O') if autho else ('H', 'C', 'N', 'S'))]
    )
    y0 = [5.0, 4.0, 2.0, -20.0, 6.0] if autho else [4.0, 3.0, 1.0, -15.0 if fo2 > 0 else 1.0]
    lo = [-80.0] * 4 + ([-11.99] if autho else [])
    hi = [14.0] * 4 + ([11.99] if autho else [])

    def f_log(y: np.ndarray) -> np.ndarray:
        yc = np.clip(y[:4], -80.0, 14.0)
        x: list[float] = [float(np.exp(v)) for v in yc]
        if autho:
            x.append(float(np.clip(y[4], -11.99, 11.99)))
            return np.asarray(func_authoritative_O(x, ddict, target)) / scales
        return np.asarray(func(x, ddict, target)) / scales

    res = opt.least_squares(
        f_log, y0, bounds=(lo, hi), ftol=1.0e-14, xtol=1.0e-14, gtol=1.0e-14
    )
    p_g = {sp: float(np.exp(res.x[i])) for i, sp in enumerate(('H2O', 'CO2', 'N2', 'S2'))}
    if autho:
        return equilibrium_atmosphere_authoritative_O(
            target, ddict, fO2_hint=float(res.x[4]), p_guess=p_g, random_seed=42, **kw
        )
    return equilibrium_atmosphere(target, ddict, p_guess=p_g, **kw)


def test_calliope_equilibrium_constants_and_reaction_quotients_vs_nasa_glenn(record_property):
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
    st_expected = {
        'janaf_H2': 0.5,
        'janaf_CO': 0.5,
        'schaefer_CH4': 1.0,
        'janaf_SO2': 0.0,
        'janaf_NH3': 0.0,
        'janaf_H2S': 0.0,
    }
    fo2_model = OxygenFugacity()
    pin = {'H2O': 120.0, 'CO2': 80.0, 'N2': 15.0, 'S2': 4.0}
    h2s_report_rows: list[str] = []

    for t_k, ref_dict in NASA_GLENN_LOG10_K_REF.items():
        t_label = f'T={int(t_k)} K'
        log10_f = float(fo2_model(t_k, 0.0))
        lk: dict[str, float] = {}
        for key, model in keq_models.items():
            val, st = model.callmodel(t_k)
            lk[key] = float(val)
            assert float(st) == pytest.approx(st_expected[key], abs=0.0), (
                f'{t_label} {key} stoich'
            )
            # Check ModifiedKeq.__call__ at fO2_shift = 0 and nonzero shifts (+2 vs -2).
            assert float(model(t_k, 0.0)) == pytest.approx(
                10.0 ** (lk[key] - st_expected[key] * log10_f), rel=1.0e-12, abs=0.0
            ), f'{t_label} {key} __call__'
            ratio_shift = float(model(t_k, 2.0)) / float(model(t_k, -2.0))
            assert ratio_shift == pytest.approx(
                10.0 ** (-4.0 * st_expected[key]), rel=1.0e-12, abs=0.0
            ), f'{t_label} {key} fO2 shift ratio'

        # Record and print the H2S equilibrium constant offset (reported, not asserted).
        d_h2s = 0.5 * (lk['janaf_H2S'] - ref_dict['janaf_H2S'])
        h2s_report_rows.append(f'{t_label}: d_log10_K_H2S={d_h2s:+.6f} dex')

        # 2. Solve gas-only equilibrium and verify reaction quotients in log10 space.
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

        q_h2 = p_h2 * (p_o2**0.5) / p_h2o
        q_co = p_co * (p_o2**0.5) / p_co2
        q_so2 = p_so2 / ((p_s2**0.5) * p_o2)
        q_nh3 = p_nh3 / ((p_n2**0.5) * (p_h2**1.5))
        q_ch4 = p_ch4 * p_o2 / (p_co2 * (p_h2**2))
        q_h2s = p_h2s / ((p_s2**0.5) * p_h2)

        assert math.log10(q_h2) == pytest.approx(lk['janaf_H2'], rel=1.0e-12, abs=0.0), (
            f'{t_label} H2 mass-action'
        )
        assert math.log10(q_co) == pytest.approx(lk['janaf_CO'], rel=1.0e-12, abs=0.0), (
            f'{t_label} CO mass-action'
        )
        assert math.log10(q_so2) == pytest.approx(
            0.5 * lk['janaf_SO2'], rel=1.0e-12, abs=0.0
        ), f'{t_label} SO2 mass-action'
        assert math.log10(q_nh3) == pytest.approx(
            0.5 * lk['janaf_NH3'], rel=1.0e-12, abs=0.0
        ), f'{t_label} NH3 mass-action'
        assert math.log10(q_ch4) == pytest.approx(lk['schaefer_CH4'], rel=1.0e-12, abs=0.0), (
            f'{t_label} CH4 mass-action'
        )
        assert math.log10(q_h2s) == pytest.approx(
            0.5 * lk['janaf_H2S'], rel=1.0e-12, abs=0.0
        ), f'{t_label} H2S mass-action'

        _assert_in_dex_bounds(
            math.log10(q_h2) - ref_dict['janaf_H2'], K_OFFSET_BOUNDS_DEX['H2'], f'{t_label} H2'
        )
        _assert_in_dex_bounds(
            math.log10(q_co) - ref_dict['janaf_CO'], K_OFFSET_BOUNDS_DEX['CO'], f'{t_label} CO'
        )
        _assert_in_dex_bounds(
            math.log10(q_so2) - 0.5 * ref_dict['janaf_SO2'],
            K_OFFSET_BOUNDS_DEX['SO2'],
            f'{t_label} SO2',
        )
        _assert_in_dex_bounds(
            math.log10(q_nh3) - 0.5 * ref_dict['janaf_NH3'],
            K_OFFSET_BOUNDS_DEX['NH3'],
            f'{t_label} NH3',
        )
        _assert_in_dex_bounds(
            math.log10(q_ch4) - ref_dict['schaefer_CH4'],
            K_OFFSET_BOUNDS_DEX['CH4'],
            f'{t_label} CH4',
        )

    h2s_summary = '; '.join(h2s_report_rows)
    record_property('h2s_k_offsets_dex', h2s_summary)
    print(f'H2S equilibrium-constant offsets vs NASA Glenn: {h2s_summary}')


def test_calliope_vs_atmodeller_72_case_gas_only_grid(record_property):
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
    sp_map = {sp.name: sp for sp in gas_species}

    # Verify NASA_GLENN_LOG10_K_REF matches atmodeller's live get_gibbs_over_RT polynomials.
    ln10 = math.log(10.0)
    for t_k, ref_dict in NASA_GLENN_LOG10_K_REF.items():
        g = {
            name: float(np.squeeze(sp.data.get_gibbs_over_RT(np.array([t_k]))))
            for name, sp in sp_map.items()
        }
        assert -(g['H2_g'] + 0.5 * g['O2_g'] - g['H2O_g']) / ln10 == pytest.approx(
            ref_dict['janaf_H2'], rel=0.0, abs=1.0e-6
        )
        assert -(g['CO_g'] + 0.5 * g['O2_g'] - g['CO2_g']) / ln10 == pytest.approx(
            ref_dict['janaf_CO'], rel=0.0, abs=1.0e-6
        )
        assert -(g['CH4_g'] + g['O2_g'] - g['CO2_g'] - 2.0 * g['H2_g']) / ln10 == pytest.approx(
            ref_dict['schaefer_CH4'], rel=0.0, abs=1.0e-6
        )
        assert -(2.0 * g['O2S_g'] - g['S2_g'] - 2.0 * g['O2_g']) / ln10 == pytest.approx(
            ref_dict['janaf_SO2'], rel=0.0, abs=1.0e-6
        )
        assert -(2.0 * g['H3N_g'] - g['N2_g'] - 3.0 * g['H2_g']) / ln10 == pytest.approx(
            ref_dict['janaf_NH3'], rel=0.0, abs=1.0e-6
        )
        assert -(2.0 * g['H2S_g'] - g['S2_g'] - 2.0 * g['H2_g']) / ln10 == pytest.approx(
            ref_dict['janaf_H2S'], rel=0.0, abs=1.0e-6
        )

    solver_params = SolverParameters(atol=1.0e-6, rtol=1.0e-4, max_steps=1024, multistart=10)
    fo2_fischer = OxygenFugacity()

    # Solve 60 buffered cases in atmodeller and assert per-case solver convergence.
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
    buf_success = np.asarray(mdl_buf.output.multi_attempt_solution.solver_success, dtype=bool)
    assert np.all(buf_success), (
        f'atmodeller buffered solve failed on {[buffered_cases[i]["case_id"] for i in np.where(~buf_success)[0]]}'
    )
    asdict_buf = mdl_buf.output.asdict()
    ptot_buf = np.asarray(mdl_buf.output.total_pressure())

    # Solve 12 authoritative-O cases in atmodeller and assert per-case solver convergence.
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
    auth_success = np.asarray(mdl_auth.output.multi_attempt_solution.solver_success, dtype=bool)
    assert np.all(auth_success), (
        f'atmodeller authO solve failed on {[autho_cases[i]["case_id"] for i in np.where(~auth_success)[0]]}'
    )
    asdict_auth = mdl_auth.output.asdict()
    ptot_auth = np.asarray(mdl_auth.output.total_pressure())

    sulfur_report_rows: list[str] = []

    for idx, case in enumerate(buffered_cases):
        cid = case['case_id']
        ddict = _build_gas_only_ddict(
            m_mantle=m_mantle,
            g_surf=g_surf,
            r_int=r_earth,
            t_magma=case['T_K'],
            fo2_shift=case['fO2_shift_IW'],
        )
        cal_out = _solve_calliope_gas_only(case['targets'], ddict, autho=False)
        hcns_max = max(float(case['targets'][k]) for k in ('H', 'C', 'N', 'S'))
        cal_gate = hcns_max * 1.0e-5 + 1.0e10 + float(TRUNC_MASS)
        for el in ('H', 'C', 'N', 'S'):
            assert abs(float(cal_out[f'{el}_res'])) <= cal_gate, f'{cid} CALLIOPE {el}_res'

        _assert_in_dex_bounds(
            math.log10(float(cal_out['P_surf']) / float(ptot_buf[idx])),
            BUFFERED_LOG10_BOUNDS_DEX['P_surf'],
            f'{cid} P_surf',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['H2O_bar']) / float(asdict_buf['H2O_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX['H2O'],
            f'{cid} H2O',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['N2_bar']) / float(asdict_buf['N2_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX['N2'],
            f'{cid} N2',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['H2_bar']) / float(asdict_buf['H2_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX['H2'],
            f'{cid} H2',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['NH3_bar']) / float(asdict_buf['H3N_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX['NH3'],
            f'{cid} NH3',
        )
        ox_suffix = 'ox' if case['fO2_shift_IW'] >= 2.0 else 'red'
        _assert_in_dex_bounds(
            math.log10(float(cal_out['CH4_bar']) / float(asdict_buf['CH4_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX[f'CH4_{ox_suffix}'],
            f'{cid} CH4',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['CO2_bar']) / float(asdict_buf['CO2_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX[f'CO2_{ox_suffix}'],
            f'{cid} CO2',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['CO_bar']) / float(asdict_buf['CO_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX[f'CO_{ox_suffix}'],
            f'{cid} CO',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['SO2_bar']) / float(asdict_buf['O2S_g']['pressure'][idx])),
            BUFFERED_LOG10_BOUNDS_DEX['SO2'],
            f'{cid} SO2',
        )

        # Reaction quotients on all 60 buffered cases (including sulfur Q_SO2).
        q_so2_cal = float(cal_out['SO2_bar']) / (
            (float(cal_out['S2_bar']) ** 0.5) * float(cal_out['O2_bar'])
        )
        q_so2_atm = float(asdict_buf['O2S_g']['pressure'][idx]) / (
            (float(asdict_buf['S2_g']['pressure'][idx]) ** 0.5)
            * float(asdict_buf['O2_g']['pressure'][idx])
        )
        _assert_in_dex_bounds(
            math.log10(q_so2_cal / q_so2_atm), K_OFFSET_BOUNDS_DEX['SO2'], f'{cid} Q_SO2'
        )

        d_h2s_case = math.log10(
            float(cal_out['H2S_bar']) / float(asdict_buf['H2S_g']['pressure'][idx])
        )
        d_s2_case = math.log10(
            float(cal_out['S2_bar']) / float(asdict_buf['S2_g']['pressure'][idx])
        )
        d_so2_case = math.log10(
            float(cal_out['SO2_bar']) / float(asdict_buf['O2S_g']['pressure'][idx])
        )
        sulfur_report_rows.append(
            f'{cid}: d_H2S={d_h2s_case:+.4f}, d_S2={d_s2_case:+.4f}, d_SO2={d_so2_case:+.4f}'
        )

    for idx, case in enumerate(autho_cases):
        cid = case['case_id']
        ddict = _build_gas_only_ddict(
            m_mantle=m_mantle,
            g_surf=g_surf,
            r_int=r_earth,
            t_magma=case['T_K'],
            fo2_shift=case['fO2_shift_IW'],
        )
        cal_out = _solve_calliope_gas_only(case['targets'], ddict, autho=True)
        for el in ('H', 'C', 'N', 'S', 'O'):
            cal_el_gate = max(float(case['targets'][el]) * 1.0e-5, float(TRUNC_MASS))
            assert abs(float(cal_out[f'{el}_res'])) <= cal_el_gate, f'{cid} CALLIOPE {el}_res'

        _assert_in_dex_bounds(
            math.log10(float(cal_out['P_surf']) / float(ptot_auth[idx])),
            AUTHO_LOG10_BOUNDS_DEX['P_surf'],
            f'{cid} P_surf',
        )
        _assert_in_dex_bounds(
            math.log10(
                float(cal_out['H2O_bar']) / float(asdict_auth['H2O_g']['pressure'][idx])
            ),
            AUTHO_LOG10_BOUNDS_DEX['H2O'],
            f'{cid} H2O',
        )
        _assert_in_dex_bounds(
            math.log10(
                float(cal_out['CO2_bar']) / float(asdict_auth['CO2_g']['pressure'][idx])
            ),
            AUTHO_LOG10_BOUNDS_DEX['CO2'],
            f'{cid} CO2',
        )
        _assert_in_dex_bounds(
            math.log10(float(cal_out['N2_bar']) / float(asdict_auth['N2_g']['pressure'][idx])),
            AUTHO_LOG10_BOUNDS_DEX['N2'],
            f'{cid} N2',
        )
        _assert_in_dex_bounds(
            math.log10(
                float(cal_out['SO2_bar']) / float(asdict_auth['O2S_g']['pressure'][idx])
            ),
            AUTHO_LOG10_BOUNDS_DEX['SO2'],
            f'{cid} SO2',
        )

        # Assert solved redox state (O2 and reduced species at T >= 2000 K; titration bound at 1500 K).
        d_o2 = math.log10(
            float(cal_out['O2_bar']) / float(asdict_auth['O2_g']['pressure'][idx])
        )
        if case['T_K'] >= 2000.0:
            _assert_in_dex_bounds(d_o2, AUTHO_LOG10_BOUNDS_DEX['O2_hiT'], f'{cid} O2')
            _assert_in_dex_bounds(
                math.log10(
                    float(cal_out['H2_bar']) / float(asdict_auth['H2_g']['pressure'][idx])
                ),
                AUTHO_LOG10_BOUNDS_DEX['H2_hiT'],
                f'{cid} H2',
            )
            _assert_in_dex_bounds(
                math.log10(
                    float(cal_out['CO_bar']) / float(asdict_auth['CO_g']['pressure'][idx])
                ),
                AUTHO_LOG10_BOUNDS_DEX['CO_hiT'],
                f'{cid} CO',
            )
            _assert_in_dex_bounds(
                math.log10(
                    float(cal_out['CH4_bar']) / float(asdict_auth['CH4_g']['pressure'][idx])
                ),
                AUTHO_LOG10_BOUNDS_DEX['CH4_hiT'],
                f'{cid} CH4',
            )
            _assert_in_dex_bounds(
                math.log10(
                    float(cal_out['NH3_bar']) / float(asdict_auth['H3N_g']['pressure'][idx])
                ),
                AUTHO_LOG10_BOUNDS_DEX['NH3_hiT'],
                f'{cid} NH3',
            )
        else:
            _assert_in_dex_bounds(d_o2, AUTHO_LOG10_BOUNDS_DEX['O2_1500K'], f'{cid} O2 1500K')

        # Gas reaction-quotient check on all 12 authoritative-O cases (including 1500 K).
        q_so2_cal = float(cal_out['SO2_bar']) / (
            (float(cal_out['S2_bar']) ** 0.5) * float(cal_out['O2_bar'])
        )
        q_so2_atm = float(asdict_auth['O2S_g']['pressure'][idx]) / (
            (float(asdict_auth['S2_g']['pressure'][idx]) ** 0.5)
            * float(asdict_auth['O2_g']['pressure'][idx])
        )
        _assert_in_dex_bounds(
            math.log10(q_so2_cal / q_so2_atm), K_OFFSET_BOUNDS_DEX['SO2'], f'{cid} Q_SO2'
        )

        d_h2s_case = math.log10(
            float(cal_out['H2S_bar']) / float(asdict_auth['H2S_g']['pressure'][idx])
        )
        d_s2_case = math.log10(
            float(cal_out['S2_bar']) / float(asdict_auth['S2_g']['pressure'][idx])
        )
        d_so2_case = math.log10(
            float(cal_out['SO2_bar']) / float(asdict_auth['O2S_g']['pressure'][idx])
        )
        sulfur_report_rows.append(
            f'{cid}: d_H2S={d_h2s_case:+.4f}, d_S2={d_s2_case:+.4f}, d_SO2={d_so2_case:+.4f}'
        )

    assert len(sulfur_report_rows) == 72
    record_property('sulfur_72_case_report_dex', '\n'.join(sulfur_report_rows))
    print('72-case sulfur species log10(CALLIOPE/atmodeller) report [dex]:')
    for row in sulfur_report_rows:
        print(f'  {row}')
