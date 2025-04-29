from __future__ import annotations

import logging
import warnings

import numpy as np
from scipy.optimize import fsolve

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

log = logging.getLogger("fwl."+__name__)

# Solve partial pressures functions
# Originally formulated by Dan Bower
# See the related issue on the PROTEUS GitHub page:-
# https://github.com/FormingWorlds/PROTEUS/issues/42
# Paper to cite:-
# https://www.sciencedirect.com/science/article/pii/S0012821X22005301

# Solve for the equilibrium chemistry of a magma ocean atmosphere
# for a given set of solubility and redox relations

TRUNC_MASS = 1e1

def is_included(gas, ddict):
    return bool(ddict[gas+"_included"]>0)

def get_partial_pressures(pin, ddict):
    """Partial pressure of all considered species from oxidised species"""

    # we only need to know pH2O, pCO2, and pN2, since reduced species
    # can be directly determined from equilibrium chemistry

    fO2_shift = ddict["fO2_shift_IW"]

    # return results in dict, to be explicit about which pressure
    # corresponds to which volatile
    p_d = {}

    # H2O
    p_d['H2O'] = pin['H2O']

    # H2
    if is_included("H2", ddict):
        gamma = ModifiedKeq('janaf_H2')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['H2'] = gamma*pin["H2O"]

    # CO2
    p_d['CO2'] = pin['CO2']

    # CO
    if is_included("CO", ddict):
        gamma = ModifiedKeq('janaf_CO')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['CO'] = gamma*pin["CO2"]

    # CH4
    if is_included("H2",ddict) and is_included("CH4",ddict):
        gamma = ModifiedKeq('schaefer_CH4')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['CH4'] = gamma*pin["CO2"]*p_d['H2']**2.0

    # N2
    p_d['N2']  = pin['N2']

    # NH3
    if is_included("NH3", ddict):
        gamma = ModifiedKeq('janaf_NH3')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['NH3']  = (gamma*pin['N2']*p_d['H2']**3)**0.5

    # O2
    fO2_model = OxygenFugacity()
    p_d['O2'] = 10.0**fO2_model(ddict['T_magma'], fO2_shift)

    # S2
    p_d['S2'] = pin['S2']

    # SO2
    if is_included("SO2", ddict):
        gamma = ModifiedKeq('janaf_SO2')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['SO2']  = (gamma*pin['S2']*p_d['O2']**2)**0.5

    # H2S
    if is_included("H2S", ddict):
        gamma = ModifiedKeq('janaf_H2S')
        gamma = gamma(ddict['T_magma'], fO2_shift)
        p_d['H2S']  = (gamma*pin['S2']*p_d['H2']**2)**0.5

    return p_d



def get_total_pressure(p_d):
    """Sum partial pressures to get total pressure"""
    return sum(p_d.values())

def atmosphere_mean_molar_mass(p_d):
    """Mean molar mass of the atmosphere"""

    ptot = get_total_pressure(p_d)

    mu_atm = 0
    for key, value in p_d.items():
        mu_atm += molar_mass[key]*value
    mu_atm /= ptot

    return mu_atm



def atmosphere_mass(pin, ddict):
    """Atmospheric mass of volatiles and totals for H, C, and N"""

    p_d = get_partial_pressures(pin, ddict)
    mu_atm = atmosphere_mean_molar_mass(p_d)

    mass_atm_d = {}
    for key, value in p_d.items():
        # 1.0E5 because pressures are in bar
        mass_atm_d[key] = value*1.0E5/ddict['gravity']
        mass_atm_d[key] *= 4.0*np.pi*ddict['radius']**2.0
        mass_atm_d[key] *= molar_mass[key]/mu_atm

    # total mass of H
    mass_atm_d['H'] = 2*mass_atm_d['H2O'] / molar_mass['H2O']
    if is_included('H2', ddict):
        mass_atm_d['H'] += 2*mass_atm_d['H2'] / molar_mass['H2']
    if is_included('CH4', ddict):
        mass_atm_d['H'] += 4*mass_atm_d['CH4'] / molar_mass['CH4']    # note factor 4 to account for stoichiometry
    if is_included("H2S", ddict):
        mass_atm_d['H'] += 2*mass_atm_d['H2S'] / molar_mass['H2S']
    if is_included("NH3", ddict):
        mass_atm_d['H'] += 3*mass_atm_d['NH3'] / molar_mass['NH3']
    # below converts moles of H2 to mass of H
    mass_atm_d['H'] *= molar_mass['H']

    # total mass of C
    mass_atm_d['C'] = mass_atm_d['CO2'] / molar_mass['CO2']
    if is_included("CO", ddict):
        mass_atm_d['C'] += mass_atm_d['CO'] / molar_mass['CO']
    if is_included("CH4", ddict):
        mass_atm_d['C'] += mass_atm_d['CH4'] / molar_mass['CH4']
    # below converts moles of C to mass of C
    mass_atm_d['C'] *= molar_mass['C']

    # total mass of N
    mass_atm_d['N'] = 2*mass_atm_d['N2']/molar_mass['N2']
    if is_included("NH3", ddict):
        mass_atm_d['N'] += 3*mass_atm_d['NH3'] / molar_mass['NH3']
    mass_atm_d['N'] *= molar_mass['N']

    # total mass of O
    mass_atm_d['O'] = mass_atm_d['H2O'] / molar_mass['H2O']
    if is_included("CO", ddict):
        mass_atm_d['O'] += mass_atm_d['CO'] / molar_mass['CO']
    if is_included("CO2", ddict):
        mass_atm_d['O'] += mass_atm_d['CO2'] / molar_mass['CO2'] * 2.0
    if is_included("SO2", ddict):
        mass_atm_d['O'] += mass_atm_d['SO2'] / molar_mass['SO2'] * 2.0
    # below converts moles of O to mass of O
    mass_atm_d['O'] *= molar_mass['O']

    # total mass of S
    mass_atm_d['S'] = mass_atm_d['S2'] / molar_mass['S2']
    if is_included("SO2", ddict):
        mass_atm_d['S'] += mass_atm_d['SO2'] / molar_mass['SO2']
    if is_included("H2S", ddict):
        mass_atm_d['S'] += mass_atm_d['H2S'] / molar_mass['H2S']
    mass_atm_d['S'] *= molar_mass['S']

    for e in element_list:
        mass_atm_d[e] = max(0.0, mass_atm_d[e])

    return mass_atm_d



def dissolved_mass(pin, ddict):
    """Volatile masses in the (molten) mantle"""

    mass_int_d = {}

    p_d = get_partial_pressures(pin, ddict)
    ptot = get_total_pressure(p_d)

    prefactor = 1E-6*ddict['M_mantle']*ddict['Phi_global']

    # H2O
    sol_H2O = SolubilityH2O() # gets the default solubility model
    ppmw_H2O = sol_H2O(p_d['H2O'])
    mass_int_d['H2O'] = prefactor*ppmw_H2O

    # CO2
    sol_CO2 = SolubilityCO2() # gets the default solubility model
    ppmw_CO2 = sol_CO2(p_d['CO2'], ddict['T_magma'])
    mass_int_d['CO2'] = prefactor*ppmw_CO2

    # CO
    if is_included("CO", ddict):
        sol_CO = SolubilityCO() # gets the default solubility model
        ppmw_CO = sol_CO(p_d["CO"], ptot)
        mass_int_d['CO'] = prefactor*ppmw_CO

    # CH4
    if is_included("CH4", ddict):
        sol_CH4 = SolubilityCH4() # gets the default solubility model
        ppmw_CH4 = sol_CH4(p_d["CH4"], ptot)
        mass_int_d['CH4'] = prefactor*ppmw_CH4

    # N2
    sol_N2 = SolubilityN2("dasgupta") # calculate fO2-dependent solubility
    ppmw_N2 = sol_N2(p_d['N2'], ptot,  ddict['T_magma'],  ddict['fO2_shift_IW'])
    mass_int_d['N2'] = prefactor*ppmw_N2

    # S2
    sol_S2 = SolubilityS2()
    ppmw_S2 = sol_S2(p_d["S2"], ddict["T_magma"], ddict['fO2_shift_IW'])
    mass_int_d['S2'] = prefactor*ppmw_S2

    # now get totals of H, C, N, O, S
    mass_int_d['H'] = mass_int_d['H2O']*2/molar_mass['H2O']
    if is_included("CH4", ddict):
        mass_int_d['H'] += mass_int_d['CH4']*4/molar_mass["CH4"]
    mass_int_d['H'] *= molar_mass['H']

    mass_int_d['C'] = mass_int_d['CO2']/molar_mass['CO2']
    if is_included("CO", ddict):
        mass_int_d['C'] += mass_int_d['CO']/molar_mass['CO']
    if is_included("CH4", ddict):
        mass_int_d['C'] += mass_int_d['CH4']/molar_mass['CH4']
    mass_int_d['C'] *= molar_mass['C']

    mass_int_d['N'] = mass_int_d['N2']

    mass_int_d['O'] = mass_int_d['H2O'] / molar_mass['H2O']
    mass_int_d['O'] += mass_int_d['CO2'] / molar_mass['CO2'] * 2.0
    if is_included("CO", ddict):
        mass_int_d['O'] += mass_int_d['CO'] / molar_mass['CO']
    mass_int_d['O'] *= molar_mass['O']

    mass_int_d['S'] = mass_int_d['S2']

    for e in element_list:
        mass_int_d[e] = max(0.0, mass_int_d[e])

    return mass_int_d

def func(pin_arr, ddict, mass_target_d):
    """Function to compute the residual of the mass balance given the partial pressures [bar]"""

    pin_dict = {
        "H2O" : pin_arr[0],
        "CO2" : pin_arr[1],
        "N2"  : pin_arr[2],
        "S2" : pin_arr[3]
    }

    # get atmospheric masses
    mass_atm_d = atmosphere_mass(pin_dict, ddict)

    # get (molten) mantle masses
    mass_int_d = dissolved_mass(pin_dict, ddict)

    # compute residuals
    res_l = []
    for vol in ['H','C','N','S']:
        # absolute residual
        res_l.append(mass_atm_d[vol] + mass_int_d[vol] - mass_target_d[vol])

    # Debug
    # H_kg = (2*mass_atm_d["H2O"]/molar_mass["H2O"] + 2*mass_atm_d["H2"]/molar_mass["H2"] + 4*mass_atm_d["CH4"]/molar_mass["CH4"]) *molar_mass['H']
    # C_kg = (mass_atm_d["CO2"]/molar_mass["CO2"] + mass_atm_d["CO"]/molar_mass["CO"] + mass_atm_d["CH4"]/molar_mass["CH4"]) * molar_mass['C']
    # N_kg = mass_atm_d["N2"]
    # print("Post:", H_kg, C_kg, N_kg)

    return res_l

def get_log_rand(rng):
    r = np.random.uniform(low=rng[0], high=rng[1])
    return 10.0**r

def get_initial_pressures(target_d):
    """Get initial guesses of partial pressures"""

    # all in bar
    cH2O = [-12, +5]  # range in log10 units
    cCO2 = [-12, +5]
    cN2  = [-12, +5]
    cS2  = [-12, +5]

    pH2O = get_log_rand(cH2O)
    pCO2 = get_log_rand(cCO2)
    pN2  = get_log_rand(cN2 )
    pS2  = get_log_rand(cS2)

    if target_d['H'] < TRUNC_MASS:
        pH2O = 0.0
    if target_d['C'] < TRUNC_MASS:
        pCO2 = 0.0
    if target_d['N'] < TRUNC_MASS:
        pN2  = 0.0
    if target_d['S'] < TRUNC_MASS:
        pS2 = 0.0

    return pH2O, pCO2, pN2, pS2


def get_target_from_params(ddict):

    N_ocean_moles = ddict['hydrogen_earth_oceans']
    CH_ratio =      ddict['CH_ratio']
    Nitrogen =      ddict['nitrogen_ppmw']
    Sulfur =        ddict['sulfur_ppmw']

    H_kg = N_ocean_moles * ocean_moles * molar_mass['H2']
    C_kg = CH_ratio * H_kg
    N_kg = Nitrogen * 1.0E-6 * ddict["M_mantle"]
    S_kg = Sulfur * 1.0E-6 * ddict["M_mantle"]
    target_d = {'H': H_kg, 'C': C_kg, 'N': N_kg, 'S': S_kg}
    return target_d

def get_target_from_pressures(ddict):

    target_d = {}

    # store partial pressures for included gases
    pin_dict = {}
    for vol in volatile_species:
        if is_included(vol, ddict):
            pin_dict[vol] = ddict[vol+"_initial_bar"]

    p_tot = np.sum(list(pin_dict.values()))
    if p_tot < 1.0e-3:
        raise Exception("Initial surface pressure too low! (%.2e bar)"%p_tot)

    # Check if no sulfur is present
    ptot_S = pin_dict["S2"]
    if is_included("SO2", ddict):
        ptot_S += pin_dict["SO2"]
    if is_included("H2S", ddict):
        ptot_S += pin_dict["H2S"]
    if ptot_S < 1.0e-20:
        target_d['S'] = 0.0

    # Check if no carbon is present
    ptot_C = pin_dict["CO2"]
    if is_included("CO", ddict):
        ptot_C += pin_dict["CO"]
    if ptot_C < 1.0e-20:
        target_d['C'] = 0.0

    # Check if no nitrogen is present
    ptot_N = pin_dict["N2"]
    if is_included("NH3", ddict):
        ptot_N += pin_dict["NH3"]
    if ptot_N < 1.0e-20:
        target_d['N'] = 0.0

    # get dissolved+atmosphere masses from partial pressures
    mass_atm_d = atmosphere_mass(pin_dict, ddict)
    mass_int_d = dissolved_mass(pin_dict, ddict)

    for vol in ['H','C','N','S']:
        if vol in target_d.keys():
            continue
        target_d[vol] = mass_atm_d[vol] + mass_int_d[vol]

    return target_d

def equilibrium_atmosphere(target_d, ddict, hide_warnings=True,
                            rtol=1e-5, atol=1e10, xtol=1e-9,
                            p_guess=None, nsolve=1500, nguess=7500):
    """Solves for surface partial pressures assuming melt-vapour eqm


    Parameters
    ----------
        target_d : dict
            Target elemental mass inventories [kg]
        ddict : dict
            Dictionary of coupler options variables

        hide_warnings : bool
            Hide floating point runtime warnings
        rtol : float
            Relative tolerance for mass conservation
        atol : float
            Absolute tolerance for mass conservation
        xtol : float
            Relative tolerance for fsolve
        p_guess : dict
            Dictionary of initial guess for partial pressures [bar]
        nsolve : int
            Maximum number of iterations allowed by fsolve
        nguess : int
            Maximum number of guesses before giving up

    Returns
    ----------
        partial_pressures : dict
            Dictionary of: volatile partial pressures [Pa], and corresponding reservoir masses [kg]
    """


    log.info("Solving for equilibrium partial pressures at surface")
    log.debug("    target masses: %s"%str(target_d))

    # solver parameters
    count = 0
    ier = 0

    # initial guess
    if p_guess is None:
        x0 = get_initial_pressures(target_d)
    else:
        x0 = (p_guess["H2O"], p_guess["CO2"], p_guess["N2"], p_guess["S2"])

    # do calculation
    success = False
    with warnings.catch_warnings():
        # Suppress warnings from solver, since they are triggered when
        # the model makes a poor guess for the composition. These are then discarded,
        # so the warning should not propagate anywhere. Errors are still printed.
        if hide_warnings:
            warnings.filterwarnings("ignore", category=RuntimeWarning)

        # could in principle result in an infinite loop, if randomising
        # the ic never finds the physical solution (but in practice,
        # this doesn't seem to happen)
        for count in range(nguess):

            # for non-dimensionalising within the solver
            scalars = [1.0, 1.0, 1.0, 1.0]
            for i,x in enumerate(x0):
                if x < TRUNC_MASS:
                    scalars[i] = 1e-6
                else:
                    scalars[i] = x*0.5

            # call solver
            sol, info, ier, msg = fsolve(func, x0, args=(ddict, target_d),
                                            full_output=True,
                                            epsfcn=1e-1, diag=scalars,
                                            xtol=xtol, maxfev=nsolve)

            # solver converged?
            success = bool(ier == 1)

            # if any negative pressures, report ier!=1
            if any(sol<0):
                # sometimes, a solution exists with negative pressures, which is clearly non-physical.
                # Here, assert we must have positive pressures.
                success = False

            # check residuals
            this_resid = func(sol, ddict, target_d)
            tolerance = np.amax(list(target_d.values())) * rtol + atol + TRUNC_MASS
            loss = np.amax(np.abs(this_resid))
            if loss > tolerance:
                if success:
                    log.debug("Rejected by residual, d(i=%d) = %.2e kg"%(np.argmax(this_resid), loss))
                success = False

            # break if success!
            if success:
                break

            # new initial guess for solver
            x0 = get_initial_pressures(target_d)

    if not success:
        raise RuntimeError("Could not find solution for volatile abundances (max attempts, %d)" % nguess)

    log.debug("    Initial guess attempt number = %d" % count)

    sol_dict = {
        "H2O" : sol[0],
        "CO2" : sol[1],
        "N2"  : sol[2],
        "S2" : sol[3]
    }

    # Final partial pressures [bar]
    p_d        = get_partial_pressures(sol_dict, ddict)

    # Final masses [kg]
    mass_atm_d = atmosphere_mass(p_d, ddict)
    mass_int_d = dissolved_mass(p_d, ddict)

    # Residuals [relative]
    res_l      = func(sol, ddict, target_d)
    log.debug("    Residuals: %s"%res_l)

    # Output dict
    outdict = {"M_atm":0.0}

    # Initialise and store partial pressures
    outdict["P_surf"] = 0.0
    for s in volatile_species:

        # Defaults
        outdict[s+"_bar"]   = 0.0      # surface partial pressure [bar]
        outdict[s+"_kg_atm"]    = 0.0      # kg in atmosphere
        outdict[s+"_kg_liquid"] = 0.0      # kg in liquid
        outdict[s+"_kg_solid"]  = 0.0      # kg in solid (not handled here)
        outdict[s+"_kg_total"]  = 0.0      # kg (total)

        # Store partial pressures
        if s in p_d.keys():
            outdict[s+"_bar"] = p_d[s]  # store as bar
            outdict["P_surf"] += outdict[s+"_bar"]

    # Store VMRs (=mole fractions) and total atmosphere
    for s in volatile_species:
        outdict[s+"_vmr"] = outdict[s+"_bar"]/outdict["P_surf"]

        log.info("    %-6s : %-8.2f bar (%.2e VMR)" % (s,outdict[s+"_bar"], outdict[s+"_vmr"]))

    # Store masses of both gases and elements
    all = [s for s in volatile_species]
    all.extend(["H","C","N","S","O"])
    for s in all:
        tot_kg = 0.0

        if s in mass_atm_d.keys():
            outdict[s+"_kg_atm"] = mass_atm_d[s]
            tot_kg += mass_atm_d[s]

        if s in mass_int_d.keys():
            outdict[s+"_kg_liquid"] = mass_int_d[s]
            outdict[s+"_kg_solid"] = 0.0
            tot_kg += mass_int_d[s]

        outdict[s+"_kg_total"] = tot_kg

    # Total atmosphere mass
    for s in volatile_species:
        outdict["M_atm"] += outdict[s+"_kg_atm"]

    # Store moles of gases and atmosphere total mmw
    outdict["atm_kg_per_mol"] = 0.0
    for s in volatile_species:
        outdict[s+"_mol_atm"]    = outdict[s+"_kg_atm"] / molar_mass[s]
        outdict[s+"_mol_solid"]  = outdict[s+"_kg_solid"] / molar_mass[s]
        outdict[s+"_mol_liquid"] = outdict[s+"_kg_liquid"] / molar_mass[s]
        outdict[s+"_mol_total"]  = outdict[s+"_mol_atm"] + outdict[s+"_mol_solid"] + outdict[s+"_mol_liquid"]

        outdict["atm_kg_per_mol"] += outdict[s+"_vmr"] * molar_mass[s]

    # Calculate elemental ratios (by mass)
    for e1 in element_list:
        for e2 in element_list:
            if e1==e2:
                continue
            em1 = outdict[e1+"_kg_atm"]
            em2 = outdict[e2+"_kg_atm"]
            if em2 == 0:
                continue  # avoid division by zero
            outdict["%s/%s_atm"%(e1,e2)] = em1/em2

    # Store residuals
    outdict["H_res"] = res_l[0]
    outdict["C_res"] = res_l[1]
    outdict["N_res"] = res_l[2]
    outdict["S_res"] = res_l[3]

    return outdict
