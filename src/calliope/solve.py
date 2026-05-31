from __future__ import annotations

import logging
import warnings

import numpy as np
import scipy.optimize as opt

from .chemistry import ModifiedKeq
from .constants import (
    element_list,
    molar_mass,
    ocean_moles,
    volatile_species,
)
from .oxygen_fugacity import OxygenFugacity
from .solubility import (
    SolubilityCH4,
    SolubilityCO,
    SolubilityCO2,
    SolubilityH2O,
    SolubilityN2,
    SolubilityS2,
)

log = logging.getLogger('fwl.' + __name__)

# Equilibrium-chemistry mass-balance solver. Original formulation by
# Bower et al. (2022): https://doi.org/10.3847/PSJ/ac5fb1

# Floor [kg] added to the residual tolerance so benign solutions are
# not rejected for sub-10-kg mass-balance noise.
TRUNC_MASS = 1e1

# Solver bounds, in one place so the cold-start guess range, the trust-constr
# box, and the fO2-hint validation cannot drift apart.
#
# Pressure cold-start draw is log-uniform over [P_GUESS_MIN_BAR, P_GUESS_MAX_BAR].
# The trust-constr box and the acceptance gate allow pressures up to
# P_CEILING_BAR, which is well above any realistic magma-ocean surface pressure;
# the wider box lets the solver explore from a poor cold start without escaping
# to non-physical territory. Sub-Neptune surface pressures can exceed the
# default guess maximum, so the guess helpers accept an optional ``p_max`` to
# widen the cold-start range without touching the box.
P_GUESS_MIN_BAR = 1.0e-12
P_GUESS_MAX_BAR = 1.0e5
P_CEILING_BAR = 1.0e7

# fO2-shift cold-start redraw range (log10 IW offset) and the hard solver box.
# The redraw range spans reducing-mantle to highly-oxidized; the wider hard box
# gives trust-constr room on poor cold starts.
FO2_GUESS_MIN = -6.0
FO2_GUESS_MAX = 8.0
FO2_HARD_MIN = -12.0
FO2_HARD_MAX = 12.0

# Surface pressure below which volatile mixing ratios are reported as zero,
# rather than dividing P_surf into a denormal and amplifying floating-point
# noise into spurious mixing ratios.
P_SURF_FLOOR_BAR = 1.0e-30


def is_included(gas, ddict):
    return bool(ddict[gas + '_included'] > 0)


def _get_partial_pressures(pin, fO2_shift, ddict):
    """Partial pressures [bar] of all 11 species from the 4 primaries.

    Internal helper that takes ``fO2_shift`` as an explicit argument
    instead of reading ``ddict['fO2_shift_IW']``. The user-facing
    ``get_partial_pressures`` is a thin wrapper that reads fO2_shift
    from ddict. Both share the same physics; the explicit form supports
    the authoritative-O solver where fO2_shift is an unknown rather
    than a config input.

    `pin` provides H2O, CO2, N2, S2; the other 7 species are derived
    via equilibrium constants and the fO2 buffer.
    """

    p_d = {s: 0.0 for s in volatile_species}

    p_d['H2O'] = pin['H2O']

    if is_included('H2', ddict):
        gamma = ModifiedKeq('janaf_H2')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['H2'] = gamma * pin['H2O']

    p_d['CO2'] = pin['CO2']

    if is_included('CO', ddict):
        gamma = ModifiedKeq('janaf_CO')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['CO'] = gamma * pin['CO2']

    if is_included('H2', ddict) and is_included('CH4', ddict):
        gamma = ModifiedKeq('schaefer_CH4')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['CH4'] = gamma * pin['CO2'] * p_d['H2'] ** 2.0

    p_d['N2'] = pin['N2']

    if is_included('NH3', ddict):
        gamma = ModifiedKeq('janaf_NH3')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['NH3'] = (gamma * pin['N2'] * p_d['H2'] ** 3) ** 0.5

    # O2 is set unconditionally from the fO2 buffer regardless of the
    # `O2_included` flag, because every sulfur and carbon equilibrium
    # below reads p_O2 (SO2 reaction, CO2-CO couple, CH4 couple).
    fO2_model = OxygenFugacity()
    p_d['O2'] = 10.0 ** fO2_model(ddict['T_magma'], fO2_shift)

    p_d['S2'] = pin['S2']

    if is_included('SO2', ddict):
        gamma = ModifiedKeq('janaf_SO2')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['SO2'] = (gamma * pin['S2'] * p_d['O2'] ** 2) ** 0.5

    if is_included('H2S', ddict) and is_included('H2', ddict):
        gamma = ModifiedKeq('janaf_H2S')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['H2S'] = (gamma * pin['S2'] * p_d['H2'] ** 2) ** 0.5

    # Silent clip: solver Monte-Carlo restarts can produce negative or
    # non-finite `pin`, which then propagates through the sqrt expressions;
    # the downstream mass tallies require non-negative, finite, real
    # pressures.
    # `max(0.0, p_d[k])` alone is insufficient: builtin max preserves NaN
    # whenever NaN is the second argument (NaN < 0 is False, so max
    # returns the second arg unchanged), and Python complex numbers do
    # not define `<` at all. The explicit checks below pin every output
    # to a non-negative real (clipping NaN, +/-inf, and any complex
    # intermediate produced by sqrt-of-negative paths down to 0).
    for k in p_d.keys():
        v = p_d[k]
        if isinstance(v, complex) or not np.isfinite(v):
            p_d[k] = 0.0
        else:
            p_d[k] = max(0.0, float(v))

    return p_d


def get_partial_pressures(pin, ddict):
    """Partial pressures [bar] of all 11 species from the 4 primaries.

    Thin wrapper around ``_get_partial_pressures`` that reads fO2_shift
    from ``ddict['fO2_shift_IW']``. Preserved bit-for-bit identical to
    the pre-refactor function for all callers.
    """
    return _get_partial_pressures(pin, ddict['fO2_shift_IW'], ddict)


def get_total_pressure(p_d):
    """Sum partial pressures to get total pressure"""
    return sum(p_d.values())


def atmosphere_mean_molar_mass(p_d):
    """Mean molar mass of the atmosphere [g/mol].

    When ``ptot`` collapses to ~0 (every partial pressure clipped to zero
    by the NaN-aware guard in ``_get_partial_pressures``), the division
    would raise ZeroDivisionError. Return a sentinel value of 1.0 g/mol;
    downstream callers in ``_atmosphere_mass`` multiply by ``p_d[k]=0``
    so all element masses come out zero, which propagates a clean
    "atmosphere is empty here" signal to the mass-balance residual.
    """

    ptot = get_total_pressure(p_d)

    if ptot < 1e-30:
        return 1.0

    mu_atm = 0
    for key, value in p_d.items():
        mu_atm += molar_mass[key] * value
    mu_atm /= ptot

    return mu_atm


def _atmosphere_mass(pin, fO2_shift, ddict):
    """Atmospheric mass of volatiles and totals for H, C, N, O, S.

    Internal helper that takes ``fO2_shift`` as an explicit argument.
    The user-facing ``atmosphere_mass`` is a thin wrapper that reads
    fO2_shift from ddict.

    CALLIOPE stores pressures in bar throughout; the only conversion
    to SI Pa happens here (factor 1e5) when computing column mass
    `kg = p_Pa * 4 pi R^2 / g`.
    """

    p_d = _get_partial_pressures(pin, fO2_shift, ddict)
    mu_atm = atmosphere_mean_molar_mass(p_d)

    mass_atm_d = {}
    for key, value in p_d.items():
        mass_atm_d[key] = value * 1.0e5 / ddict['gravity']
        mass_atm_d[key] *= 4.0 * np.pi * ddict['radius'] ** 2.0
        mass_atm_d[key] *= molar_mass[key] / mu_atm

    # Per-element kg accumulators are built in mol-of-element first
    # (atom counts: H2O=2H, CH4=4H, H2S=2H, NH3=3H+1N, SO2=2O+1S, ...)
    # then multiplied by the element molar mass at the end.
    mass_atm_d['H'] = 2 * mass_atm_d['H2O'] / molar_mass['H2O']
    if is_included('H2', ddict):
        mass_atm_d['H'] += 2 * mass_atm_d['H2'] / molar_mass['H2']
    if is_included('CH4', ddict):
        mass_atm_d['H'] += 4 * mass_atm_d['CH4'] / molar_mass['CH4']
    if is_included('H2S', ddict):
        mass_atm_d['H'] += 2 * mass_atm_d['H2S'] / molar_mass['H2S']
    if is_included('NH3', ddict):
        mass_atm_d['H'] += 3 * mass_atm_d['NH3'] / molar_mass['NH3']
    mass_atm_d['H'] *= molar_mass['H']

    mass_atm_d['C'] = mass_atm_d['CO2'] / molar_mass['CO2']
    if is_included('CO', ddict):
        mass_atm_d['C'] += mass_atm_d['CO'] / molar_mass['CO']
    if is_included('CH4', ddict):
        mass_atm_d['C'] += mass_atm_d['CH4'] / molar_mass['CH4']
    mass_atm_d['C'] *= molar_mass['C']

    mass_atm_d['N'] = 2 * mass_atm_d['N2'] / molar_mass['N2']
    if is_included('NH3', ddict):
        mass_atm_d['N'] += mass_atm_d['NH3'] / molar_mass['NH3']
    mass_atm_d['N'] *= molar_mass['N']

    mass_atm_d['O'] = mass_atm_d['H2O'] / molar_mass['H2O']
    mass_atm_d['O'] += 2 * mass_atm_d['O2'] / molar_mass['O2']
    # CO2 is one of the four primary species in `pin`; mass_atm_d['CO2']
    # is populated unconditionally and contributes to the C tally
    # without gating. The O contribution must match for element
    # bookkeeping to be symmetric.
    mass_atm_d['O'] += mass_atm_d['CO2'] / molar_mass['CO2'] * 2.0
    if is_included('CO', ddict):
        mass_atm_d['O'] += mass_atm_d['CO'] / molar_mass['CO']
    if is_included('SO2', ddict):
        mass_atm_d['O'] += mass_atm_d['SO2'] / molar_mass['SO2'] * 2.0
    mass_atm_d['O'] *= molar_mass['O']

    mass_atm_d['S'] = 2 * mass_atm_d['S2'] / molar_mass['S2']
    if is_included('SO2', ddict):
        mass_atm_d['S'] += mass_atm_d['SO2'] / molar_mass['SO2']
    if is_included('H2S', ddict):
        mass_atm_d['S'] += mass_atm_d['H2S'] / molar_mass['H2S']
    mass_atm_d['S'] *= molar_mass['S']

    for e in element_list:
        mass_atm_d[e] = max(0.0, mass_atm_d[e])

    return mass_atm_d


def atmosphere_mass(pin, ddict):
    """Atmospheric mass of volatiles and totals for H, C, N, O, S.

    Thin wrapper around ``_atmosphere_mass`` that reads fO2_shift from
    ``ddict['fO2_shift_IW']``. Preserved bit-for-bit identical to the
    pre-refactor function for all callers.
    """
    return _atmosphere_mass(pin, ddict['fO2_shift_IW'], ddict)


def _dissolved_mass(pin, fO2_shift, ddict):
    """Volatile masses in the (molten) mantle.

    Internal helper that takes ``fO2_shift`` as an explicit argument
    instead of reading ``ddict['fO2_shift_IW']``. Two solubility laws
    (``SolubilityN2('dasgupta')`` and ``SolubilityS2()``) consume fO2_shift
    directly, so the shift must flow through to them too. The user-facing
    ``dissolved_mass`` is a thin wrapper that reads fO2_shift from ddict.
    """

    mass_int_d = {}

    p_d = _get_partial_pressures(pin, fO2_shift, ddict)
    ptot = get_total_pressure(p_d)

    # Henry's-law / power-law solubility laws return ppmw of the
    # species in the melt; prefactor converts ppmw -> kg dissolved.
    prefactor = 1e-6 * ddict['M_mantle'] * ddict['Phi_global']

    sol_H2O = SolubilityH2O()
    ppmw_H2O = sol_H2O(p_d['H2O'])
    mass_int_d['H2O'] = prefactor * ppmw_H2O

    sol_CO2 = SolubilityCO2()
    ppmw_CO2 = sol_CO2(p_d['CO2'], ddict['T_magma'])
    mass_int_d['CO2'] = prefactor * ppmw_CO2

    if is_included('CO', ddict):
        sol_CO = SolubilityCO()
        ppmw_CO = sol_CO(p_d['CO'], ptot)
        mass_int_d['CO'] = prefactor * ppmw_CO
    else:
        mass_int_d['CO'] = 0.0

    if is_included('CH4', ddict):
        sol_CH4 = SolubilityCH4()
        ppmw_CH4 = sol_CH4(p_d['CH4'], ptot)
        mass_int_d['CH4'] = prefactor * ppmw_CH4
    else:
        mass_int_d['CH4'] = 0.0

    # Override class default 'libourel'; dasgupta carries fO2 + p_total
    # dependence which the linear Libourel law does not.
    sol_N2 = SolubilityN2('dasgupta')
    ppmw_N2 = sol_N2(p_d['N2'], ptot, ddict['T_magma'], fO2_shift)
    mass_int_d['N2'] = prefactor * ppmw_N2

    sol_S2 = SolubilityS2()
    ppmw_S2 = sol_S2(p_d['S2'], ddict['T_magma'], fO2_shift)
    mass_int_d['S2'] = prefactor * ppmw_S2

    # No SolubilityH2S / NH3 / SO2 / O2 / H2 in CALLIOPE; these
    # species do not partition into the melt phase and contribute
    # zero to the dissolved-element tallies below.
    mass_int_d['H'] = mass_int_d['H2O'] * 2 / molar_mass['H2O']
    if is_included('CH4', ddict):
        mass_int_d['H'] += mass_int_d['CH4'] * 4 / molar_mass['CH4']
    mass_int_d['H'] *= molar_mass['H']

    mass_int_d['C'] = mass_int_d['CO2'] / molar_mass['CO2']
    if is_included('CO', ddict):
        mass_int_d['C'] += mass_int_d['CO'] / molar_mass['CO']
    if is_included('CH4', ddict):
        mass_int_d['C'] += mass_int_d['CH4'] / molar_mass['CH4']
    mass_int_d['C'] *= molar_mass['C']

    mass_int_d['N'] = mass_int_d['N2']

    mass_int_d['O'] = mass_int_d['H2O'] / molar_mass['H2O']
    mass_int_d['O'] += mass_int_d['CO2'] / molar_mass['CO2'] * 2.0
    if is_included('CO', ddict):
        mass_int_d['O'] += mass_int_d['CO'] / molar_mass['CO']
    mass_int_d['O'] *= molar_mass['O']

    mass_int_d['S'] = mass_int_d['S2']

    for e in element_list:
        mass_int_d[e] = max(0.0, mass_int_d[e])

    return mass_int_d


def dissolved_mass(pin, ddict):
    """Volatile masses in the (molten) mantle.

    Thin wrapper around ``_dissolved_mass`` that reads fO2_shift from
    ``ddict['fO2_shift_IW']``. Preserved bit-for-bit identical to the
    pre-refactor function for all callers.
    """
    return _dissolved_mass(pin, ddict['fO2_shift_IW'], ddict)


def func(pin_arr, ddict, mass_target_d):
    """Mass-balance residual [kg per element] for the four primary partial pressures [bar]."""

    pin_dict = {'H2O': pin_arr[0], 'CO2': pin_arr[1], 'N2': pin_arr[2], 'S2': pin_arr[3]}

    mass_atm_d = atmosphere_mass(pin_dict, ddict)
    mass_int_d = dissolved_mass(pin_dict, ddict)

    res_l = [0.0] * 4
    for i, vol in enumerate(['H', 'C', 'N', 'S']):
        res_l[i] = mass_atm_d[vol] + mass_int_d[vol] - mass_target_d[vol]

    return res_l


def obj(pin_arr, ddict, mass_target_d):
    """Function to compute the residual of the mass balance given the partial pressures [bar]"""

    res_l = func(pin_arr, ddict, mass_target_d)
    return np.dot(res_l, res_l) ** 0.5


def func_authoritative_O(x_arr, ddict, mass_target_d):
    """5-residual vector for the authoritative-O solver mode.

    Mass-balance residual [kg per element] over five unknowns:
    ``[pH2O, pCO2, pN2, pS2, fO2_shift]``. The first four equations are
    the usual H, C, N, S mass balances; the fifth is the O mass balance
    that was implicit in the chemistry (because fO2 was an input) and is
    now an explicit constraint (because fO2 is an unknown).

    Parameters
    ----------
    x_arr : array_like, length 5
        ``[pH2O_bar, pCO2_bar, pN2_bar, pS2_bar, fO2_shift_IW]``. The
        first four are partial pressures in bar; the fifth is the
        IW-buffer offset in log10 units (typical range -6 to +8).
    ddict : dict
        Coupler options dict. Reads everything except ``fO2_shift_IW``,
        which is taken from ``x_arr[4]`` to expose it as an unknown.
    mass_target_d : dict
        Target elemental mass inventories [kg]. MUST include the keys
        ``'H'``, ``'C'``, ``'N'``, ``'S'``, ``'O'`` (all five). Missing
        ``'O'`` raises ``KeyError``.

    Returns
    -------
    list of float, length 5
        Residuals ``(atm_kg + dissolved_kg) - target_kg`` for H, C, N,
        S, O in that order.
    """

    pin_dict = {'H2O': x_arr[0], 'CO2': x_arr[1], 'N2': x_arr[2], 'S2': x_arr[3]}
    fO2_shift = x_arr[4]

    mass_atm_d = _atmosphere_mass(pin_dict, fO2_shift, ddict)
    mass_int_d = _dissolved_mass(pin_dict, fO2_shift, ddict)

    res_l = [0.0] * 5
    for i, elem in enumerate(['H', 'C', 'N', 'S', 'O']):
        res_l[i] = mass_atm_d[elem] + mass_int_d[elem] - mass_target_d[elem]

    return res_l


def obj_authoritative_O(x_arr, ddict, mass_target_d):
    """Scalar objective for trust-constr fallback in the authoritative-O solver."""
    res_l = func_authoritative_O(x_arr, ddict, mass_target_d)
    return np.dot(res_l, res_l) ** 0.5


def get_initial_pressures(target_d, p_max=P_GUESS_MAX_BAR):
    """Cold-start guesses for the four primary partial pressures [bar].

    Log-uniform draw over [P_GUESS_MIN_BAR, ``p_max``] bar, covering ~17
    orders of magnitude from trace-volatile undersaturation up to the
    default ~100 kbar (the upper end of magma-ocean surface-pressure
    regimes). `target_d` is accepted for API stability but not consulted.

    Parameters
    ----------
    target_d : dict
        Accepted for API parity; not consulted.
    p_max : float, default ``P_GUESS_MAX_BAR``
        Upper bound of the log-uniform pressure draw [bar]. Raise it for
        high-pressure (e.g. sub-Neptune) cases whose surface pressure can
        exceed the default; the solver box (``P_CEILING_BAR``) is unchanged.
    """
    hi = np.log10(p_max)
    lo = np.log10(P_GUESS_MIN_BAR)
    pH2O = 10 ** np.random.uniform(low=lo, high=hi)
    pCO2 = 10 ** np.random.uniform(low=lo, high=hi)
    pN2 = 10 ** np.random.uniform(low=lo, high=hi)
    pS2 = 10 ** np.random.uniform(low=lo, high=hi)

    return pH2O, pCO2, pN2, pS2


def get_initial_pressures_with_fO2(
    target_d, fO2_hint, restart=False, rng=None, p_max=P_GUESS_MAX_BAR
):
    """Cold-start guesses for the five unknowns of the authoritative-O solver.

    Returns ``[pH2O, pCO2, pN2, pS2, fO2_shift]``. The four pressures use
    the same log-uniform draw as ``get_initial_pressures`` over
    ``[P_GUESS_MIN_BAR, p_max]`` bar. The fifth element is ``fO2_hint`` on
    the first attempt; on solver restart (``restart=True``) it is redrawn
    from a uniform distribution over ``[FO2_GUESS_MIN, FO2_GUESS_MAX]``,
    which covers the reducing-mantle to highly-oxidized regimes likely to
    be encountered.

    Parameters
    ----------
    target_d : dict
        Target elemental mass inventories. Accepted for API parity with
        ``get_initial_pressures`` but not consulted.
    fO2_hint : float
        Initial guess for the IW-buffer offset (log10 units). Typical
        PROTEUS user values lie in ``[-4, +6]``.
    restart : bool, default False
        When True, redraw fO2_shift from ``Uniform(FO2_GUESS_MIN,
        FO2_GUESS_MAX)``. When False, return ``fO2_hint`` unchanged.
    rng : np.random.Generator or None, default None
        Random number generator for the log-uniform pressure draw and
        (when ``restart=True``) the fO2 redraw. When None, the global
        ``np.random`` state is used. The authoritative-O entry point
        threads a seeded generator through this argument to make solver
        outcomes reproducible across calls.
    p_max : float, default ``P_GUESS_MAX_BAR``
        Upper bound of the log-uniform pressure draw [bar]. Raise it for
        high-pressure (e.g. sub-Neptune) cases; the solver box
        (``P_CEILING_BAR``) is unchanged.

    Returns
    -------
    tuple of 5 floats
        ``(pH2O, pCO2, pN2, pS2, fO2_shift)``.
    """
    if rng is None:
        rng = np.random

    hi = np.log10(p_max)
    lo = np.log10(P_GUESS_MIN_BAR)
    pH2O = 10 ** rng.uniform(low=lo, high=hi)
    pCO2 = 10 ** rng.uniform(low=lo, high=hi)
    pN2 = 10 ** rng.uniform(low=lo, high=hi)
    pS2 = 10 ** rng.uniform(low=lo, high=hi)

    if restart:
        fO2 = rng.uniform(low=FO2_GUESS_MIN, high=FO2_GUESS_MAX)
    else:
        fO2 = float(fO2_hint)

    return pH2O, pCO2, pN2, pS2, fO2


def get_target_from_params(ddict):

    N_ocean_moles = ddict['hydrogen_earth_oceans']
    CH_ratio = ddict['CH_ratio']
    Nitrogen = ddict['nitrogen_ppmw']
    Sulfur = ddict['sulfur_ppmw']

    H_kg = N_ocean_moles * ocean_moles * molar_mass['H2']
    C_kg = CH_ratio * H_kg
    N_kg = Nitrogen * 1.0e-6 * ddict['M_mantle']
    S_kg = Sulfur * 1.0e-6 * ddict['M_mantle']
    target_d = {'H': H_kg, 'C': C_kg, 'N': N_kg, 'S': S_kg}
    return target_d


def get_target_from_pressures(ddict):

    target_d = {}

    pin_dict = {}
    for vol in volatile_species:
        if is_included(vol, ddict):
            pin_dict[vol] = ddict[vol + '_initial_bar']

    p_tot = np.sum(list(pin_dict.values()))
    if p_tot < 1.0e-3:
        raise Exception('Initial surface pressure too low! (%.2e bar)' % p_tot)

    # Per-element short-circuits: when no S-, C-, or N-bearing primary
    # is initially present, pin the corresponding target to zero
    # rather than dragging the dissolved-mass solubility laws through
    # divide-by-near-zero paths (1e-20 bar threshold is well below any
    # solver-relevant pressure).
    ptot_S = pin_dict['S2']
    if is_included('SO2', ddict):
        ptot_S += pin_dict['SO2']
    if is_included('H2S', ddict):
        ptot_S += pin_dict['H2S']
    if ptot_S < 1.0e-20:
        target_d['S'] = 0.0

    ptot_C = pin_dict['CO2']
    if is_included('CO', ddict):
        ptot_C += pin_dict['CO']
    if ptot_C < 1.0e-20:
        target_d['C'] = 0.0

    ptot_N = pin_dict['N2']
    if is_included('NH3', ddict):
        ptot_N += pin_dict['NH3']
    if ptot_N < 1.0e-20:
        target_d['N'] = 0.0

    mass_atm_d = atmosphere_mass(pin_dict, ddict)
    mass_int_d = dissolved_mass(pin_dict, ddict)

    for vol in ['H', 'C', 'N', 'S']:
        if vol in target_d.keys():
            continue
        target_d[vol] = mass_atm_d[vol] + mass_int_d[vol]

    return target_d


def equilibrium_atmosphere(
    target_d,
    ddict,
    hide_warnings=True,
    rtol=1e-5,
    atol=1e10,
    xtol=1e-8,
    p_guess=None,
    nsolve=1500,
    nguess=7500,
    print_result=True,
    opt_solver=True,
    p_max=P_GUESS_MAX_BAR,
):
    """Solve for surface partial pressures assuming melt-vapour equilibrium.

    Parameters
    ----------
    target_d : dict
        Target elemental mass inventories [kg], with keys 'H', 'C', 'N', 'S'.
    ddict : dict
        Dictionary of coupler options variables (planet, magma state, inclusion flags).
    hide_warnings : bool, default True
        Hide floating point runtime warnings raised by `scipy` for poor guesses.
    rtol : float, default 1e-5
        Relative tolerance for mass conservation.
    atol : float, default 1e10
        Absolute tolerance for mass conservation [kg].
    xtol : float, default 1e-8
        Relative tolerance for fsolve.
    p_guess : dict or None, default None
        Dictionary of initial guess for primary-species partial pressures [bar].
        Must contain the keys 'H2O', 'CO2', 'N2', 'S2', each mapping to a
        finite real number. If None, an internal Monte-Carlo log-uniform
        draw is used. A non-dict value raises TypeError; missing keys or
        non-finite values raise ValueError.
    nsolve : int, default 1500
        Maximum number of inner-solver iterations per attempt.
    nguess : int, default 7500
        Maximum number of Monte-Carlo restarts before giving up.
    print_result : bool, default True
        If True, log final outgassed partial pressures at INFO level.
    opt_solver : bool, default True
        If True, alternate between fsolve and trust-constr on each restart.
    p_max : float, default ``P_GUESS_MAX_BAR``
        Upper bound [bar] of the Monte-Carlo cold-start pressure draw. Raise it
        for high-pressure (e.g. sub-Neptune) cases whose surface pressure can
        exceed the default; the solver box (``P_CEILING_BAR``) is unchanged.

    Returns
    -------
    partial_pressures : dict
        Volatile partial pressures [bar] keyed `<species>_bar`, plus per-species
        reservoir masses [kg], elemental totals, residuals, and atmospheric
        diagnostics (P_surf, M_atm, atm_kg_per_mol, ratios).
    """

    if print_result:
        log.info('Solving for equilibrium partial pressures at surface')
    log.debug('    target masses: %s' % str(target_d))

    # Hard ub = P_CEILING_BAR is well above any realistic magma-ocean
    # surface pressure; it prevents trust-constr from exploring
    # non-physical regions during a poor cold start.
    lb = [0.0] * 4
    ub = [P_CEILING_BAR] * 4

    if p_guess is None:
        x0 = get_initial_pressures(target_d, p_max=p_max)
    else:
        # Validate up front so a missing key surfaces as a clear ValueError
        # rather than a bare KeyError from the tuple construction below.
        # The isinstance check handles cases where a caller passes a list,
        # a pandas Series, or accidentally a non-dict object; without it,
        # the membership test below would raise an opaque TypeError.
        if not isinstance(p_guess, dict):
            raise TypeError(f'p_guess must be a dict or None, got {type(p_guess).__name__}.')
        required = ('H2O', 'CO2', 'N2', 'S2')
        missing = [k for k in required if k not in p_guess]
        if missing:
            raise ValueError(
                f'p_guess is missing required keys: {missing}. '
                f'Expected all of {list(required)}.'
            )
        x0 = (p_guess['H2O'], p_guess['CO2'], p_guess['N2'], p_guess['S2'])

        # Reject non-finite values up front. NaN > 1e-10 evaluates False,
        # which would silently collapse ub to 1.0 and propagate NaN through
        # the solver to produce garbage output with no error signal.
        for k, v in zip(required, x0):
            if not np.isfinite(v):
                raise ValueError(f'p_guess[{k!r}] must be a finite real number, got {v!r}.')

        # Zero or near-zero guess collapses ub from 1e7 to 1.0 to keep
        # trust-constr from wandering inside a degenerate slot.
        for i in range(4):
            ub[i] = ub[i] if (x0[i] > 1e-10) else 1.0

    bounds = opt.Bounds(lb=lb, ub=ub)

    tolerance = np.amax(list(target_d.values())) * rtol + atol + TRUNC_MASS
    log.debug('Required tolerance: %g' % tolerance)

    with warnings.catch_warnings():
        # Solver Monte-Carlo restarts produce poor guesses that trip
        # RuntimeWarning / UserWarning; the bad attempts are discarded
        # so the warnings should not reach the caller.
        if hide_warnings:
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            warnings.filterwarnings('ignore', category=UserWarning)

        solver: int = 0
        for count in range(nguess):
            if solver == 0:
                sol, _, ier, _ = opt.fsolve(
                    func, x0, args=(ddict, target_d), maxfev=nsolve, xtol=xtol, full_output=True
                )
                success = bool(ier == 1)
            else:
                result = opt.minimize(
                    obj,
                    x0,
                    args=(ddict, target_d),
                    method='trust-constr',
                    bounds=bounds,
                    options={'maxiter': nsolve, 'xtol': xtol},
                )
                success = result.success
                sol = result.x

            # Reject solver-claimed successes that miss the residual
            # tolerance: fsolve and trust-constr converge by their own
            # criteria, which are not the kg-mass-balance criterion.
            this_resid = func(sol, ddict, target_d)
            loss = np.amax(np.abs(this_resid))
            if loss > tolerance:
                if success:
                    log.debug('Solution rejected by residual')
                    log.debug('    d(i=%d) = %.2e kg' % (np.argmax(this_resid), loss))
                success = False

            if success:
                break

            x0 = get_initial_pressures(target_d, p_max=p_max)

            # Alternate fsolve <-> trust-constr on each restart so a
            # basin one solver cannot escape gets a chance from the
            # other. Disable via opt_solver=False to pin to fsolve.
            if opt_solver:
                solver = 1 - solver

    if not success:
        raise RuntimeError(
            'Could not find solution for volatile abundances (max attempts, %d)' % nguess
        )

    log.debug('    Initial guess attempt number = %d' % count)

    res_l = func(sol, ddict, target_d)
    log.debug('    Residuals: %s' % res_l)

    sol_dict = {'H2O': sol[0], 'CO2': sol[1], 'N2': sol[2], 'S2': sol[3]}
    p_d = get_partial_pressures(sol_dict, ddict)

    # Pass the 4-key primary dict (sol_dict), not the expanded p_d:
    # atmosphere_mass and dissolved_mass each call get_partial_pressures
    # internally, so feeding them p_d would re-run gas-phase chemistry
    # on a dict that already contains the secondaries (read but ignored).
    mass_atm_d = atmosphere_mass(sol_dict, ddict)
    mass_int_d = dissolved_mass(sol_dict, ddict)

    # CALLIOPE has no solid-mantle reservoir; every `_kg_solid` field
    # is 0.0 and exists only for schema parity with the downstream
    # PROTEUS hf_row which carries solid/melt/atmosphere splits.
    outdict = {'M_atm': 0.0, 'P_surf': 0.0}
    for s in volatile_species:
        outdict[s + '_bar'] = 0.0
        outdict[s + '_kg_atm'] = 0.0
        outdict[s + '_kg_liquid'] = 0.0
        outdict[s + '_kg_solid'] = 0.0
        outdict[s + '_kg_total'] = 0.0

        if s in p_d.keys():
            outdict[s + '_bar'] = p_d[s]
            outdict['P_surf'] += outdict[s + '_bar']

    P_surf = outdict['P_surf']
    for s in volatile_species:
        outdict[s + '_vmr'] = (
            (outdict[s + '_bar'] / P_surf) if P_surf > P_SURF_FLOOR_BAR else 0.0
        )

        if print_result:
            log.info(
                '    %-6s : %-8.2f bar (%.2e VMR)'
                % (s, outdict[s + '_bar'], outdict[s + '_vmr'])
            )

    all = [s for s in volatile_species]
    all.extend(['H', 'C', 'N', 'S', 'O'])
    for s in all:
        tot_kg = 0.0

        if s in mass_atm_d.keys():
            outdict[s + '_kg_atm'] = mass_atm_d[s]
            tot_kg += mass_atm_d[s]

        if s in mass_int_d.keys():
            outdict[s + '_kg_liquid'] = mass_int_d[s]
            outdict[s + '_kg_solid'] = 0.0
            tot_kg += mass_int_d[s]

        outdict[s + '_kg_total'] = tot_kg

    for s in volatile_species:
        outdict['M_atm'] += outdict[s + '_kg_atm']

    outdict['atm_kg_per_mol'] = 0.0
    for s in volatile_species:
        outdict[s + '_mol_atm'] = outdict[s + '_kg_atm'] / molar_mass[s]
        outdict[s + '_mol_solid'] = outdict[s + '_kg_solid'] / molar_mass[s]
        outdict[s + '_mol_liquid'] = outdict[s + '_kg_liquid'] / molar_mass[s]
        outdict[s + '_mol_total'] = (
            outdict[s + '_mol_atm'] + outdict[s + '_mol_solid'] + outdict[s + '_mol_liquid']
        )

        outdict['atm_kg_per_mol'] += outdict[s + '_vmr'] * molar_mass[s]

    for e1 in element_list:
        for e2 in element_list:
            if e1 == e2:
                continue
            em1 = outdict[e1 + '_kg_atm']
            em2 = outdict[e2 + '_kg_atm']
            if em2 == 0:
                continue
            outdict['%s/%s_atm' % (e1, e2)] = em1 / em2

    outdict['H_res'] = res_l[0]
    outdict['C_res'] = res_l[1]
    outdict['N_res'] = res_l[2]
    outdict['S_res'] = res_l[3]

    return outdict


def equilibrium_atmosphere_authoritative_O(
    target_d,
    ddict,
    fO2_hint=4.0,
    hide_warnings=True,
    rtol=1e-5,
    atol=1e10,
    xtol=1e-8,
    p_guess=None,
    nsolve=1500,
    nguess=7500,
    print_result=True,
    opt_solver=True,
    random_seed=None,
    p_max=P_GUESS_MAX_BAR,
):
    """Solve for partial pressures AND fO2 given total elemental masses including O.

    Authoritative-oxygen solver mode. Unlike ``equilibrium_atmosphere``
    (which takes fO2 as a config input via ``ddict['fO2_shift_IW']``),
    this entry point treats fO2 as a fifth unknown and solves a 5x5
    nonlinear mass-balance system. The user supplies a target O mass
    alongside H/C/N/S, and the solver returns the partial pressures
    plus the IW-buffer offset (``fO2_shift_derived``) that produces
    that equilibrium.

    Use this when the science model declares atmospheric+dissolved O
    as a budget (e.g. mantle FeO inventory, or whole-planet O accounting
    where atmospheric escape and ingassing debit the same O reservoir).
    For the legacy mode where fO2 is buffered to a user-specified IW
    offset and O is derived, use ``equilibrium_atmosphere`` instead.

    Parameters
    ----------
    target_d : dict
        Target elemental mass inventories [kg]. MUST contain the keys
        ``'H'``, ``'C'``, ``'N'``, ``'S'``, ``'O'``. Missing ``'O'``
        raises ``KeyError``.
    ddict : dict
        Coupler options dict (planet, magma state, inclusion flags).
        ``ddict['fO2_shift_IW']`` is IGNORED by this entry point; the
        value is treated as a solver unknown initialised from
        ``fO2_hint`` instead.
    fO2_hint : float, default 4.0
        Initial guess for the IW-buffer offset (log10 units). Provide
        a value close to the expected solution to speed convergence;
        typical PROTEUS values lie in [-4, +6]. The solver's Monte-Carlo
        restarts redraw from Uniform(-6, +8) if the hint does not lead
        to convergence.
    hide_warnings : bool, default True
        Hide floating point runtime warnings raised by `scipy` for poor guesses.
    rtol : float, default 1e-5
        Relative tolerance for mass conservation.
    atol : float, default 1e10
        Absolute tolerance for mass conservation [kg].
    xtol : float, default 1e-8
        Relative tolerance for fsolve.
    p_guess : dict or None, default None
        Initial guess for primary-species partial pressures [bar]. Keys
        must include ``'H2O'``, ``'CO2'``, ``'N2'``, ``'S2'``; the
        optional key ``'fO2_shift_IW'`` overrides ``fO2_hint`` for the
        starting guess. Non-dict raises TypeError; missing required keys
        or non-finite values raise ValueError.
    nsolve : int, default 1500
        Maximum number of inner-solver iterations per attempt.
    nguess : int, default 7500
        Maximum number of Monte-Carlo restarts before giving up.
    print_result : bool, default True
        If True, log final outgassed partial pressures and derived
        fO2 at INFO level.
    opt_solver : bool, default True
        If True, alternate between fsolve and trust-constr on each
        restart so a basin one solver cannot escape gets a chance from
        the other.
    random_seed : int or None, default None
        Seed for the Monte-Carlo restart RNG. ``None`` uses the global
        ``np.random`` state (non-deterministic). An integer seed makes
        solver outcomes reproducible across calls, which is required
        for regression testing and for diffing two runs.
    p_max : float, default ``P_GUESS_MAX_BAR``
        Upper bound [bar] of the Monte-Carlo cold-start pressure draw. Raise it
        for high-pressure (e.g. sub-Neptune) cases whose surface pressure can
        exceed the default; the solver box (``P_CEILING_BAR``) is unchanged.

    Returns
    -------
    partial_pressures : dict
        Volatile partial pressures [bar] keyed ``<species>_bar``, plus
        per-species reservoir masses [kg], elemental totals, residuals,
        and atmospheric diagnostics. Two additions relative to
        ``equilibrium_atmosphere``:

        - ``fO2_shift_derived`` : float
            The IW-buffer offset the solver converged to. Equals
            ``fO2_hint`` only if the hint happened to be the
            self-consistent value.
        - ``O_res`` : float
            5th residual (O mass-balance), in kg. Pairs with the
            existing ``H_res``/``C_res``/``N_res``/``S_res`` keys.

    Raises
    ------
    KeyError
        If ``target_d`` is missing the ``'O'`` key.
    TypeError, ValueError
        If ``p_guess`` fails validation (same contract as
        ``equilibrium_atmosphere``).
    RuntimeError
        If the solver fails to converge after ``nguess`` Monte-Carlo
        restarts. The error message includes the final pressures and
        fO2_shift attempt for diagnosis.

    Notes
    -----
    Mathematical model. The four existing mass-balance equations for
    H/C/N/S are extended with a fifth for O. The four primary
    partial pressures (H2O, CO2, N2, S2) are joined by fO2_shift as a
    fifth unknown. All seven derived partial pressures (H2, CO, CH4,
    NH3, O2, SO2, H2S) and the two fO2-coupled solubility laws (N2
    Dasgupta, S2 Gaillard) consume fO2_shift through the same
    physics functions ``_get_partial_pressures``, ``_atmosphere_mass``,
    ``_dissolved_mass`` that ``equilibrium_atmosphere`` uses.

    Examples
    --------
    Reproducing an equilibrium_atmosphere result through the new mode:

    >>> # First, run the legacy mode at fO2_shift_IW = +4
    >>> ddict = {..., 'fO2_shift_IW': 4.0}
    >>> out_legacy = equilibrium_atmosphere(target_d_HCNS, ddict)
    >>> # Then, run the new mode with the implied O budget
    >>> target_d_HCNSO = dict(target_d_HCNS,
    ...                      O=out_legacy['O_kg_total'])
    >>> out_new = equilibrium_atmosphere_authoritative_O(
    ...     target_d_HCNSO, ddict, fO2_hint=4.0)
    >>> abs(out_new['fO2_shift_derived'] - 4.0) < 0.01  # round-trip
    True
    """

    required_elements = ('H', 'C', 'N', 'S', 'O')

    # Contract check: every required element key must be present, finite,
    # and non-negative. A missing key raises KeyError (matching the legacy
    # "target_d must include 'O'" message). Non-finite or negative values
    # are user-error and raise ValueError before the solver wastes effort.
    missing = [e for e in required_elements if e not in target_d]
    if missing:
        raise KeyError(
            'target_d is missing required element keys: %s. '
            'Authoritative-O mode requires all of %s. Got keys: %s'
            % (missing, list(required_elements), sorted(target_d.keys()))
        )
    for e in required_elements:
        v = target_d[e]
        if not np.isfinite(v):
            raise ValueError('target_d[%r] must be a finite real number [kg], got %r.' % (e, v))
        if v < 0:
            raise ValueError('target_d[%r] must be non-negative [kg], got %r.' % (e, v))

    # Validate fO2_hint and the planet/state parameters consumed from
    # ddict by the residual chain. The solver evaluates the residual at
    # arbitrarily wild trial points, so it cannot recover from a bad
    # entry value; failing fast here gives a useful error rather than
    # a ZeroDivisionError from deep inside the chemistry path.
    if not np.isfinite(fO2_hint):
        raise ValueError(
            'fO2_hint must be a finite real number (log10 IW offset), got %r.' % fO2_hint
        )
    if not (FO2_HARD_MIN <= fO2_hint <= FO2_HARD_MAX):
        raise ValueError(
            'fO2_hint=%.3f is outside the solver bounds [%+g, %+g]. '
            'Pick a value in [%+g, %+g] for physically realistic mantle '
            'redox states.'
            % (fO2_hint, FO2_HARD_MIN, FO2_HARD_MAX, FO2_GUESS_MIN, FO2_GUESS_MAX)
        )

    for required_ddict_key in ('M_mantle', 'Phi_global', 'T_magma', 'gravity', 'radius'):
        if required_ddict_key not in ddict:
            raise KeyError(
                'ddict is missing required key %r. Authoritative-O mode '
                'requires M_mantle, Phi_global, T_magma, gravity, radius.' % required_ddict_key
            )

    M_mantle = ddict['M_mantle']
    if not (np.isfinite(M_mantle) and M_mantle > 0):
        raise ValueError("ddict['M_mantle']=%r must be a positive finite mass [kg]." % M_mantle)

    Phi_global = ddict['Phi_global']
    if not (np.isfinite(Phi_global) and 0.0 <= Phi_global <= 1.0):
        raise ValueError(
            "ddict['Phi_global']=%r must lie in [0, 1] (melt mass fraction)." % Phi_global
        )

    T_magma = ddict['T_magma']
    if not (np.isfinite(T_magma) and T_magma > 0):
        raise ValueError(
            "ddict['T_magma']=%r must be a positive finite temperature [K]." % T_magma
        )

    if nguess < 1:
        raise ValueError('nguess must be >= 1, got %d.' % nguess)
    if nsolve < 1:
        raise ValueError('nsolve must be >= 1, got %d.' % nsolve)

    # Seeded RNG for the Monte-Carlo restart draws. random_seed=None
    # falls back to the global np.random state to preserve historical
    # non-deterministic behaviour for callers that do not opt in to
    # reproducibility.
    rng = np.random.default_rng(random_seed) if random_seed is not None else np.random

    if print_result:
        log.info(
            'Solving for equilibrium partial pressures + fO2_shift '
            '(authoritative-O mode, fO2_hint=%.2f)',
            fO2_hint,
        )
    log.debug('    target masses: %s', target_d)

    # Bounds. Pressures: [0, P_CEILING_BAR] bar (same as equilibrium_atmosphere).
    # fO2_shift: [FO2_HARD_MIN, FO2_HARD_MAX] log10 units. The physically
    # meaningful range is roughly [FO2_GUESS_MIN, FO2_GUESS_MAX] (mantle
    # reducing to highly oxidized); the wider hard box gives trust-constr room
    # to explore on poor cold starts without escaping to non-physical territory.
    lb = [0.0, 0.0, 0.0, 0.0, FO2_HARD_MIN]
    ub = [P_CEILING_BAR, P_CEILING_BAR, P_CEILING_BAR, P_CEILING_BAR, FO2_HARD_MAX]

    # Physical pressure ceiling for the acceptance gate, captured before
    # the per-guess ub-collapse below mutates ub. The gate tests against
    # this documented 1e7 bar bound, not a collapsed conditioning value,
    # so a tiny-guess slot whose root legitimately lands above its
    # collapsed 1.0 bar bound is not falsely rejected.
    p_ceiling = ub[0]

    if p_guess is None:
        x0 = get_initial_pressures_with_fO2(target_d, fO2_hint, rng=rng, p_max=p_max)
    else:
        if not isinstance(p_guess, dict):
            raise TypeError(f'p_guess must be a dict or None, got {type(p_guess).__name__}.')
        required = ('H2O', 'CO2', 'N2', 'S2')
        missing = [k for k in required if k not in p_guess]
        if missing:
            raise ValueError(
                f'p_guess is missing required keys: {missing}. '
                f'Expected all of {list(required)}.'
            )
        # fO2_shift_IW in p_guess overrides fO2_hint; absent means use fO2_hint.
        fO2_seed = p_guess.get('fO2_shift_IW', fO2_hint)
        x0 = (p_guess['H2O'], p_guess['CO2'], p_guess['N2'], p_guess['S2'], fO2_seed)

        # Reject non-finite values.
        for k, v in zip(required + ('fO2_shift_IW',), x0):
            if not np.isfinite(v):
                raise ValueError(f'p_guess[{k!r}] must be a finite real number, got {v!r}.')

        # Match the legacy ub-collapse for tiny pressure guesses so
        # trust-constr does not wander in a degenerate slot. fO2 bound
        # stays at the wide range.
        ub_collapsed = list(ub)
        for i in range(4):
            ub_collapsed[i] = ub_collapsed[i] if (x0[i] > 1e-10) else 1.0
        ub = ub_collapsed

    bounds = opt.Bounds(lb=lb, ub=ub)

    # Per-element tolerance: each residual must satisfy
    # ``|res_i| <= max(target_i * rtol, TRUNC_MASS)``. The gate is
    # relative per element, so mass closure is judged against each
    # element's own budget rather than the largest. The absolute floor is
    # the small TRUNC_MASS noise level (10 kg) so a near-zero target does
    # not demand an exactly-zero residual. The floor is deliberately not
    # tied to ``atol`` (the planetary negligible-mass threshold, ~1e16
    # kg), which would dominate the relative gate on small-budget
    # elements such as N and let multi-percent closure errors pass there.
    target_vec = np.array([target_d[e] for e in required_elements])
    elem_tolerance = np.maximum(target_vec * rtol, TRUNC_MASS)
    log.debug('Per-element tolerance: %s kg', elem_tolerance.tolist())

    # `sol` initialised to a sentinel so the post-loop RuntimeError path
    # can format the final attempt even if every iteration crashed
    # before `sol` was assigned. `success` starts False so an early
    # break from a 0-iteration loop (already rejected by the nguess>=1
    # validator above, but a defensive belt) raises RuntimeError rather
    # than UnboundLocalError.
    sol = np.array(x0, dtype=float)
    success = False
    count = 0

    with warnings.catch_warnings():
        if hide_warnings:
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            warnings.filterwarnings('ignore', category=UserWarning)

        solver: int = 0
        for count in range(nguess):
            # Wrap each solver call: fsolve and trust-constr can both
            # raise ZeroDivisionError or FloatingPointError if the
            # residual blows up in a numerically unrecoverable way at a
            # trial point. Catch and treat as a failed attempt so the
            # restart loop continues instead of propagating a crash to
            # the caller.
            try:
                if solver == 0:
                    sol, _, ier, _ = opt.fsolve(
                        func_authoritative_O,
                        x0,
                        args=(ddict, target_d),
                        maxfev=nsolve,
                        xtol=xtol,
                        full_output=True,
                    )
                    success = bool(ier == 1)
                else:
                    result = opt.minimize(
                        obj_authoritative_O,
                        x0,
                        args=(ddict, target_d),
                        method='trust-constr',
                        bounds=bounds,
                        options={'maxiter': nsolve, 'xtol': xtol},
                    )
                    success = result.success
                    sol = result.x
            except (ZeroDivisionError, FloatingPointError, ValueError) as exc:
                log.debug(
                    'Solver attempt %d (method=%s) raised %s: %s; restarting',
                    count,
                    'fsolve' if solver == 0 else 'trust-constr',
                    type(exc).__name__,
                    exc,
                )
                success = False

            # Per-element residual gate. Compute defensively so a
            # post-solver evaluation crash also routes to restart.
            if success:
                try:
                    this_resid = func_authoritative_O(sol, ddict, target_d)
                    resid_abs = np.abs(np.asarray(this_resid))
                    if np.any(resid_abs > elem_tolerance):
                        worst = int(np.argmax(resid_abs - elem_tolerance))
                        log.debug(
                            'Solution rejected by per-element residual: '
                            'element=%s, |res|=%.2e, tol=%.2e',
                            required_elements[worst],
                            resid_abs[worst],
                            elem_tolerance[worst],
                        )
                        success = False
                except (ZeroDivisionError, FloatingPointError, ValueError) as exc:
                    log.debug(
                        'Post-solver residual evaluation raised %s: %s; rejecting attempt %d',
                        type(exc).__name__,
                        exc,
                        count,
                    )
                    success = False

            # Reject a converged-but-non-physical root. The production
            # path runs only the unbounded fsolve (opt_solver=False), so a
            # root can satisfy mass balance yet sit outside the physical
            # box: a derived fO2_shift beyond [-12, +12] (the target O is
            # unreachable at this H/C/N/S/T_magma), a negative partial
            # pressure, or a partial pressure above the 1e7 bar ceiling.
            # trust-constr enforces `bounds`; fsolve does not, so the full
            # box is enforced here before the solution is accepted.
            if success:
                sol_p = np.asarray(sol[:4], dtype=float)
                if not (lb[4] <= sol[4] <= ub[4]):
                    log.debug(
                        'Solution rejected: derived fO2_shift=%.3f outside [%.1f, %.1f]',
                        sol[4],
                        lb[4],
                        ub[4],
                    )
                    success = False
                elif np.any(sol_p < -1.0e-6):
                    log.debug('Solution rejected: negative partial pressure %s bar', sol[:4])
                    success = False
                elif np.any(sol_p > p_ceiling * (1.0 + 1.0e-6)):
                    log.debug(
                        'Solution rejected: partial pressure above %.1e bar ceiling: %s bar',
                        p_ceiling,
                        sol[:4],
                    )
                    success = False

            if success:
                break

            # Restart. Redraw pressures from log-uniform; redraw fO2
            # from Uniform(-6, +8) to give it a chance from a different
            # basin if the hint led to a non-converging region.
            x0 = get_initial_pressures_with_fO2(
                target_d, fO2_hint, restart=True, rng=rng, p_max=p_max
            )

            if opt_solver:
                solver = 1 - solver

    if not success:
        raise RuntimeError(
            'Could not find solution for volatile abundances + fO2 under '
            'authoritative-O mode (max attempts: %d). '
            'Final attempt: pH2O=%.3e bar, pCO2=%.3e bar, pN2=%.3e bar, '
            'pS2=%.3e bar, fO2_shift=%.3f. '
            'Either the target O budget is outside the physically '
            'reachable range at this (H, C, N, S, T_magma), or the '
            'chemistry has a non-monotonic region the solver could not '
            'escape. Consider adjusting fO2_hint or the target masses.'
            % (nguess, sol[0], sol[1], sol[2], sol[3], sol[4])
        )

    log.debug('    Initial guess attempt number = %d', count)

    res_l = func_authoritative_O(sol, ddict, target_d)
    log.debug('    Residuals: %s', res_l)
    log.debug('    Derived fO2_shift: %.4f', sol[4])

    sol_dict = {'H2O': sol[0], 'CO2': sol[1], 'N2': sol[2], 'S2': sol[3]}
    fO2_derived = sol[4]
    p_d = _get_partial_pressures(sol_dict, fO2_derived, ddict)

    mass_atm_d = _atmosphere_mass(sol_dict, fO2_derived, ddict)
    mass_int_d = _dissolved_mass(sol_dict, fO2_derived, ddict)

    # Output dict structure matches equilibrium_atmosphere bit-for-bit
    # plus the two new keys (fO2_shift_derived, O_res). The PROTEUS-side
    # wrapper consumes the same fields regardless of which solver mode
    # ran.
    outdict = {'M_atm': 0.0, 'P_surf': 0.0}
    for s in volatile_species:
        outdict[s + '_bar'] = 0.0
        outdict[s + '_kg_atm'] = 0.0
        outdict[s + '_kg_liquid'] = 0.0
        outdict[s + '_kg_solid'] = 0.0
        outdict[s + '_kg_total'] = 0.0

        if s in p_d.keys():
            outdict[s + '_bar'] = p_d[s]
            outdict['P_surf'] += outdict[s + '_bar']

    P_surf = outdict['P_surf']
    for s in volatile_species:
        outdict[s + '_vmr'] = (
            (outdict[s + '_bar'] / P_surf) if P_surf > P_SURF_FLOOR_BAR else 0.0
        )

        if print_result:
            log.info(
                '    %-6s : %-8.2f bar (%.2e VMR)',
                s,
                outdict[s + '_bar'],
                outdict[s + '_vmr'],
            )

    all_keys = [s for s in volatile_species]
    all_keys.extend(['H', 'C', 'N', 'S', 'O'])
    for s in all_keys:
        tot_kg = 0.0

        if s in mass_atm_d.keys():
            outdict[s + '_kg_atm'] = mass_atm_d[s]
            tot_kg += mass_atm_d[s]

        if s in mass_int_d.keys():
            outdict[s + '_kg_liquid'] = mass_int_d[s]
            outdict[s + '_kg_solid'] = 0.0
            tot_kg += mass_int_d[s]

        outdict[s + '_kg_total'] = tot_kg

    for s in volatile_species:
        outdict['M_atm'] += outdict[s + '_kg_atm']

    outdict['atm_kg_per_mol'] = 0.0
    for s in volatile_species:
        outdict[s + '_mol_atm'] = outdict[s + '_kg_atm'] / molar_mass[s]
        outdict[s + '_mol_solid'] = outdict[s + '_kg_solid'] / molar_mass[s]
        outdict[s + '_mol_liquid'] = outdict[s + '_kg_liquid'] / molar_mass[s]
        outdict[s + '_mol_total'] = (
            outdict[s + '_mol_atm'] + outdict[s + '_mol_solid'] + outdict[s + '_mol_liquid']
        )

        outdict['atm_kg_per_mol'] += outdict[s + '_vmr'] * molar_mass[s]

    for e1 in element_list:
        for e2 in element_list:
            if e1 == e2:
                continue
            em1 = outdict[e1 + '_kg_atm']
            em2 = outdict[e2 + '_kg_atm']
            if em2 == 0:
                continue
            outdict['%s/%s_atm' % (e1, e2)] = em1 / em2

    outdict['H_res'] = res_l[0]
    outdict['C_res'] = res_l[1]
    outdict['N_res'] = res_l[2]
    outdict['S_res'] = res_l[3]
    outdict['O_res'] = res_l[4]
    outdict['fO2_shift_derived'] = fO2_derived

    if print_result:
        log.info('    Derived fO2_shift = %.4f (hint was %.4f)', fO2_derived, fO2_hint)

    return outdict
