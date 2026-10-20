"""Cross-code thermodynamic and gas-only equilibrium consistency tests.

Verifies:
1. CALLIOPE equilibrium constants (ModifiedKeq) and algebraic secondary gas
   speciation against the 9-coefficient NASA Glenn reference values across
   T in {1500, 2000, 2500, 3000} K within +-0.002 dex of the measured per-temperature
   reaction-data offset (K_OFFSET_EXPECTED_DEX). For H2S, CALLIOPE's janaf_H2S
   (tools/fit_H2S.ipynb) fits g_H2S - g_S2 - g_H2 (unit S2 stoichiometric
   coefficient over 200-4000 K instead of 0.5 S2), producing a +0.29 to -0.39 dex
   per-mole K_p offset across 1500-3000 K; H2S, S2 on 54 cases, and any case-species key
   whose offset shifts by > half its assert margin under a 0.5 S2 refit are held: these keys
   depend on the janaf_H2S constant.
2. Live 72-case gas-only equilibrium cross-comparison against atmodeller (60
   buffered cases with matched fO2 plus 12 authoritative-O cases) on the shared
   volatile species (H2O, H2, O2, CO2, CO, CH4, N2, NH3, S2, SO2, H2S). Elemental
   inventories are defined in moles (n_O = 0.5 n_H + 2 n_C + 2 n_S at the
   stoichiometric H2O + CO2 + SO2 endpoint) and converted to mass per code with
   that code's atomic masses (calliope.constants.molar_mass and molmass.Formula),
   and surface gravity is matched to atmodeller.Planet.surface_gravity.
   Across all 12 authoritative-O cases (1500 to 3000 K), solved log10(fO2) (O2 partial
   pressure) agrees within expected +-0.005 dex (measured -0.0007 to +0.0026).
   H2 follows DeltaK with H2O and O2, CO and CH4 with CO2 and O2, NH3 with N2 and H2;
   the relation holds for the pinned values.
"""

from __future__ import annotations

import importlib
import logging
import math
import os
from pathlib import Path
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

# Measured log10(K_CALLIOPE / K_NASA_Glenn) per mole of product at each temperature,
# asserted within +-K_OFFSET_MARGIN_DEX (0.002 dex); H2S is held:
# these keys depend on the janaf_H2S constant.
K_OFFSET_MARGIN_DEX: float = 0.002
K_OFFSET_HELD_FOR_H2S_FIX: tuple[str, ...] = ('H2S',)
K_OFFSET_EXPECTED_DEX: dict[float, dict[str, float]] = {
    1500.0: {
        'H2': -0.005993,
        'CO': 0.013892,
        'CH4': 0.059526,
        'SO2': 0.007863,
        'NH3': -0.010623,
    },
    2000.0: {
        'H2': 0.002763,
        'CO': -0.006148,
        'CH4': 0.100441,
        'SO2': 0.007004,
        'NH3': 0.033241,
    },
    2500.0: {
        'H2': 0.003928,
        'CO': -0.004536,
        'CH4': 0.098591,
        'SO2': 0.006492,
        'NH3': 0.053011,
    },
    3000.0: {'H2': 0.003666, 'CO': 0.005555, 'CH4': 0.072052, 'SO2': 0.007657, 'NH3': 0.059389},
}

# Keys held: these keys depend on the janaf_H2S constant:
# buffered H2S, S2 on all 60 buffered cases plus the 152 buffered (case_id, species)
# keys whose offset shifts by > 0.002 dex under the S2-coefficient-0.5 janaf_H2S refit.
BUFFERED_SPECIES_HELD_FOR_H2S_FIX: tuple[str, ...] = ('H2S',)
AUTHO_HELD_FOR_H2S_FIX: tuple[str, ...] = ('H2S',)
AUTHO_CASE_KEYS_HELD_FOR_H2S_FIX: dict[str, tuple[str, ...]] = {
    'authO_default_T3000': ('S2',),
    'authO_high_T3000': ('S2',),
    'authO_low_T3000': ('S2',),
}
BUFFERED_CASE_KEYS_HELD_FOR_H2S_FIX: dict[str, tuple[str, ...]] = {
    'buffered_default_T1500_IW+0': ('SO2', 'S2'),
    'buffered_default_T1500_IW+2': ('SO2', 'S2'),
    'buffered_default_T1500_IW+4': (
        'P_surf',
        'H2O',
        'N2',
        'H2',
        'NH3',
        'CH4',
        'CO2',
        'CO',
        'SO2',
        'S2',
    ),
    'buffered_default_T1500_IW-2': ('SO2', 'S2'),
    'buffered_default_T1500_IW-4': ('SO2', 'S2'),
    'buffered_default_T2000_IW+0': ('SO2', 'S2'),
    'buffered_default_T2000_IW+2': ('SO2', 'S2'),
    'buffered_default_T2000_IW-2': ('SO2', 'S2'),
    'buffered_default_T2000_IW-4': ('SO2', 'S2'),
    'buffered_default_T2500_IW+0': (
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_default_T2500_IW+2': (
        'P_surf',
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_default_T2500_IW-2': ('SO2', 'S2'),
    'buffered_default_T2500_IW-4': ('SO2', 'S2'),
    'buffered_default_T3000_IW+0': (
        'P_surf',
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_default_T3000_IW+2': ('SO2', 'S2'),
    'buffered_default_T3000_IW-2': (
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_default_T3000_IW-4': (
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_high_T1500_IW+0': ('SO2', 'S2'),
    'buffered_high_T1500_IW+2': ('SO2', 'S2'),
    'buffered_high_T1500_IW+4': (
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_high_T1500_IW-2': ('SO2', 'S2'),
    'buffered_high_T1500_IW-4': ('SO2', 'S2'),
    'buffered_high_T2000_IW+0': ('SO2', 'S2'),
    'buffered_high_T2000_IW+2': ('SO2', 'S2'),
    'buffered_high_T2000_IW+4': ('SO2', 'S2'),
    'buffered_high_T2000_IW-2': ('SO2', 'S2'),
    'buffered_high_T2000_IW-4': ('SO2', 'S2'),
    'buffered_high_T2500_IW+0': ('SO2', 'S2'),
    'buffered_high_T2500_IW+2': (
        'P_surf',
        'H2O',
        'N2',
        'H2',
        'NH3',
        'CH4',
        'CO2',
        'CO',
        'SO2',
        'S2',
    ),
    'buffered_high_T2500_IW-2': ('SO2', 'S2'),
    'buffered_high_T2500_IW-4': ('SO2', 'S2'),
    'buffered_high_T3000_IW+0': (
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_high_T3000_IW+2': (
        'P_surf',
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_high_T3000_IW-2': ('SO2', 'S2'),
    'buffered_high_T3000_IW-4': ('SO2', 'S2'),
    'buffered_low_T1500_IW+0': ('SO2', 'S2'),
    'buffered_low_T1500_IW+2': ('SO2', 'S2'),
    'buffered_low_T1500_IW+4': (
        'P_surf',
        'H2O',
        'N2',
        'H2',
        'NH3',
        'CH4',
        'CO2',
        'CO',
        'SO2',
        'S2',
    ),
    'buffered_low_T1500_IW-2': ('SO2', 'S2'),
    'buffered_low_T1500_IW-4': ('SO2', 'S2'),
    'buffered_low_T2000_IW+0': ('SO2', 'S2'),
    'buffered_low_T2000_IW+2': ('SO2', 'S2'),
    'buffered_low_T2000_IW-2': ('SO2', 'S2'),
    'buffered_low_T2000_IW-4': ('SO2', 'S2'),
    'buffered_low_T2500_IW+0': (
        'P_surf',
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_low_T2500_IW+2': ('SO2', 'S2'),
    'buffered_low_T2500_IW-2': (
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_low_T2500_IW-4': (
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_low_T3000_IW+0': (
        'P_surf',
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_low_T3000_IW-2': (
        'H2O',
        'H2',
        'NH3',
        'CH4',
        'SO2',
        'S2',
    ),
    'buffered_low_T3000_IW-4': (
        'H2O',
        'N2',
        'H2',
        'NH3',
        'CH4',
        'CO2',
        'CO',
        'SO2',
        'S2',
    ),
}

# Measured per-case log10(P_CALLIOPE / P_atmodeller) offsets at matched moles
# and g_surf = 9.930522 m s^-2
# for all H2S-independent keys (shift <= 0.002 dex), asserted within +-0.005 dex.
GRID_OFFSET_MARGIN_DEX: float = 0.005
AUTHO_MAJOR_MARGIN_DEX: float = 0.001
BUFFERED_EXPECTED_DEX: dict[str, dict[str, float]] = {
    'buffered_default_T1500_IW+0': {
        'P_surf': 0.000059,
        'H2O': 0.003495,
        'CO2': -0.013405,
        'N2': 0.001060,
        'H2': -0.002498,
        'CO': 0.000487,
        'CH4': 0.041125,
        'NH3': -0.013840,
    },
    'buffered_default_T1500_IW+2': {
        'P_surf': 0.000101,
        'H2O': 0.000892,
        'CO2': -0.003977,
        'N2': 0.000164,
        'H2': -0.005101,
        'CO': 0.009915,
        'CH4': 0.045347,
        'NH3': -0.018192,
    },
    'buffered_default_T1500_IW+4': {},
    'buffered_default_T1500_IW-2': {
        'P_surf': -0.002583,
        'H2O': 0.002076,
        'CO2': -0.029318,
        'N2': 0.000579,
        'H2': -0.003917,
        'CO': -0.015427,
        'CH4': 0.022373,
        'NH3': -0.016210,
    },
    'buffered_default_T1500_IW-4': {
        'P_surf': -0.001930,
        'H2O': 0.003201,
        'CO2': -0.047030,
        'N2': 0.000268,
        'H2': -0.002791,
        'CO': -0.033138,
        'CH4': 0.006914,
        'NH3': -0.014676,
    },
    'buffered_default_T2000_IW+0': {
        'P_surf': -0.000438,
        'H2O': -0.001864,
        'CO2': 0.004618,
        'N2': -0.000573,
        'H2': 0.000900,
        'CO': -0.001531,
        'CH4': 0.106858,
        'NH3': 0.034303,
    },
    'buffered_default_T2000_IW+2': {
        'P_surf': 0.000770,
        'H2O': 0.000770,
        'CO2': 0.002408,
        'N2': 0.000360,
        'H2': 0.003533,
        'CO': -0.003740,
        'CH4': 0.109916,
        'NH3': 0.038721,
    },
    'buffered_default_T2000_IW+4': {
        'P_surf': 0.000050,
        'H2O': 0.000025,
        'CO2': 0.000320,
        'N2': 0.000030,
        'SO2': 0.000410,
        'H2': 0.002788,
        'CO': -0.005829,
        'CH4': 0.106337,
        'NH3': 0.037438,
        'S2': -0.013189,
    },
    'buffered_default_T2000_IW-2': {
        'P_surf': -0.000474,
        'H2O': -0.003074,
        'CO2': 0.004897,
        'N2': -0.000419,
        'H2': -0.000310,
        'CO': -0.001251,
        'CH4': 0.104718,
        'NH3': 0.032566,
    },
    'buffered_default_T2000_IW-4': {
        'P_surf': -0.001929,
        'H2O': -0.005236,
        'CO2': -0.001259,
        'N2': -0.000143,
        'H2': -0.002472,
        'CO': -0.007407,
        'CH4': 0.094238,
        'NH3': 0.029461,
    },
    'buffered_default_T2500_IW+0': {
        'P_surf': 0.000809,
        'CO2': 0.003510,
        'N2': -0.000426,
        'CO': -0.001025,
    },
    'buffered_default_T2500_IW+2': {
        'CO2': 0.003023,
        'N2': 0.001426,
        'CO': -0.001513,
    },
    'buffered_default_T2500_IW+4': {
        'P_surf': 0.000010,
        'H2O': -0.000024,
        'CO2': 0.000238,
        'N2': 0.000007,
        'SO2': 0.000056,
        'H2': 0.003904,
        'CO': -0.004298,
        'CH4': 0.106637,
        'NH3': 0.058870,
        'S2': -0.012873,
    },
    'buffered_default_T2500_IW-2': {
        'P_surf': -0.000238,
        'H2O': -0.002902,
        'CO2': 0.003624,
        'N2': -0.000951,
        'H2': 0.001026,
        'CO': -0.000911,
        'CH4': 0.104267,
        'NH3': 0.054074,
    },
    'buffered_default_T2500_IW-4': {
        'P_surf': -0.000137,
        'H2O': -0.003232,
        'CO2': 0.003493,
        'N2': -0.000697,
        'H2': 0.000696,
        'CO': -0.001043,
        'CH4': 0.103476,
        'NH3': 0.053706,
    },
    'buffered_default_T3000_IW+0': {
        'CO2': -0.003524,
        'N2': 0.001099,
        'CO': 0.002030,
    },
    'buffered_default_T3000_IW+2': {
        'P_surf': 0.000328,
        'H2O': 0.000143,
        'CO2': -0.001947,
        'N2': 0.000043,
        'H2': 0.003810,
        'CO': 0.003608,
        'CH4': 0.077725,
        'NH3': 0.065125,
    },
    'buffered_default_T3000_IW+4': {
        'P_surf': -0.000029,
        'H2O': -0.000059,
        'CO2': -0.000327,
        'N2': -0.000032,
        'SO2': -0.000025,
        'H2': 0.003608,
        'CO': 0.005228,
        'CH4': 0.078941,
        'NH3': 0.064785,
        'S2': -0.015363,
    },
    'buffered_default_T3000_IW-2': {
        'P_surf': 0.000528,
        'CO2': -0.007445,
        'N2': -0.002127,
        'CO': -0.001891,
    },
    'buffered_default_T3000_IW-4': {
        'P_surf': 0.000057,
        'CO2': -0.007972,
        'N2': -0.002513,
        'CO': -0.002417,
    },
    'buffered_high_T1500_IW+0': {
        'P_surf': -0.001623,
        'H2O': 0.001206,
        'CO2': -0.041057,
        'N2': 0.002076,
        'H2': -0.004786,
        'CO': -0.027165,
        'CH4': 0.008896,
        'NH3': -0.016765,
    },
    'buffered_high_T1500_IW+2': {
        'P_surf': -0.000353,
        'H2O': 0.000371,
        'CO2': -0.005483,
        'N2': 0.000431,
        'H2': -0.005622,
        'CO': 0.008408,
        'CH4': 0.042798,
        'NH3': -0.018841,
    },
    'buffered_high_T1500_IW+4': {
        'P_surf': -0.001024,
        'H2O': -0.001189,
        'CO2': -0.001059,
        'N2': -0.000483,
        'H2': -0.007182,
        'CO': 0.012833,
    },
    'buffered_high_T1500_IW-2': {
        'P_surf': 0.000063,
        'H2O': 0.005610,
        'CO2': -0.057314,
        'N2': 0.001431,
        'H2': -0.000383,
        'CO': -0.043422,
        'CH4': 0.001447,
        'NH3': -0.010482,
    },
    'buffered_high_T1500_IW-4': {
        'P_surf': 0.000058,
        'H2O': 0.006079,
        'CO2': -0.059556,
        'N2': 0.000991,
        'H2': 0.000086,
        'CO': -0.045664,
        'CH4': 0.000142,
        'NH3': -0.009999,
    },
    'buffered_high_T2000_IW+0': {
        'P_surf': -0.002222,
        'H2O': -0.004157,
        'CO2': -0.001036,
        'N2': -0.001936,
        'H2': -0.001393,
        'CO': -0.007184,
        'CH4': 0.096619,
        'NH3': 0.030183,
    },
    'buffered_high_T2000_IW+2': {
        'P_surf': 0.000141,
        'H2O': -0.000084,
        'CO2': 0.002094,
        'N2': -0.000098,
        'H2': 0.002679,
        'CO': -0.004055,
        'CH4': 0.107893,
        'NH3': 0.037210,
    },
    'buffered_high_T2000_IW+4': {
        'P_surf': 0.000336,
        'H2O': 0.000348,
        'CO2': 0.000465,
        'N2': 0.000168,
        'H2': 0.003111,
        'CO': -0.005683,
        'CH4': 0.107129,
        'NH3': 0.037992,
    },
    'buffered_high_T2000_IW-2': {
        'P_surf': -0.007045,
        'H2O': -0.011838,
        'CO2': -0.028697,
        'N2': -0.001367,
        'H2': -0.009074,
        'CO': -0.034846,
        'CH4': 0.053596,
        'NH3': 0.018946,
    },
    'buffered_high_T2000_IW-4': {
        'P_surf': -0.005106,
        'H2O': -0.010000,
        'CO2': -0.068058,
        'N2': -0.000695,
        'H2': -0.007237,
        'CO': -0.074207,
        'CH4': 0.017909,
        'NH3': 0.022038,
    },
    'buffered_high_T2500_IW+0': {
        'P_surf': -0.000621,
        'H2O': -0.002292,
        'CO2': 0.002856,
        'N2': -0.001717,
        'H2': 0.001636,
        'CO': -0.001680,
        'CH4': 0.104718,
        'NH3': 0.054605,
    },
    'buffered_high_T2500_IW+2': {},
    'buffered_high_T2500_IW+4': {
        'P_surf': 0.000052,
        'H2O': 0.000025,
        'CO2': 0.000260,
        'N2': 0.000024,
        'SO2': 0.000509,
        'H2': 0.003952,
        'CO': -0.004276,
        'CH4': 0.106755,
        'NH3': 0.058951,
        'S2': -0.011966,
    },
    'buffered_high_T2500_IW-2': {
        'P_surf': -0.001747,
        'H2O': -0.005618,
        'CO2': -0.000760,
        'N2': -0.002027,
        'H2': -0.001691,
        'CO': -0.005296,
        'CH4': 0.094449,
        'NH3': 0.049461,
    },
    'buffered_high_T2500_IW-4': {
        'P_surf': -0.005461,
        'H2O': -0.011019,
        'CO2': -0.020763,
        'N2': -0.001161,
        'H2': -0.007092,
        'CO': -0.025299,
        'CH4': 0.063644,
        'NH3': 0.041792,
    },
    'buffered_high_T3000_IW+0': {
        'P_surf': 0.001230,
        'CO2': -0.004896,
        'N2': -0.000937,
        'CO': 0.000659,
    },
    'buffered_high_T3000_IW+2': {
        'CO2': -0.000046,
        'N2': 0.001867,
        'CO': 0.005508,
    },
    'buffered_high_T3000_IW+4': {
        'P_surf': -0.000029,
        'H2O': -0.000056,
        'CO2': -0.000328,
        'N2': -0.000035,
        'SO2': 0.000034,
        'H2': 0.003610,
        'CO': 0.005227,
        'CH4': 0.078945,
        'NH3': 0.064787,
        'S2': -0.015246,
    },
    'buffered_high_T3000_IW-2': {
        'P_surf': -0.000419,
        'H2O': -0.002850,
        'CO2': -0.006795,
        'N2': -0.002063,
        'H2': 0.000816,
        'CO': -0.001240,
        'CH4': 0.066889,
        'NH3': 0.059582,
    },
    'buffered_high_T3000_IW-4': {
        'P_surf': -0.001193,
        'H2O': -0.004302,
        'CO2': -0.010580,
        'N2': -0.001725,
        'H2': -0.000636,
        'CO': -0.005025,
        'CH4': 0.060201,
        'NH3': 0.057573,
    },
    'buffered_low_T1500_IW+0': {
        'P_surf': 0.000976,
        'H2O': 0.004620,
        'CO2': -0.010242,
        'N2': 0.001015,
        'H2': -0.001373,
        'CO': 0.003650,
        'CH4': 0.046538,
        'NH3': -0.012175,
    },
    'buffered_low_T1500_IW+2': {
        'P_surf': -0.000176,
        'H2O': 0.000411,
        'CO2': -0.004067,
        'N2': 0.000021,
        'H2': -0.005582,
        'CO': 0.009825,
        'CH4': 0.044295,
        'NH3': -0.018985,
    },
    'buffered_low_T1500_IW+4': {},
    'buffered_low_T1500_IW-2': {
        'P_surf': 0.000309,
        'H2O': 0.005897,
        'CO2': -0.013585,
        'N2': 0.000462,
        'H2': -0.000096,
        'CO': 0.000307,
        'CH4': 0.045749,
        'NH3': -0.010536,
    },
    'buffered_low_T1500_IW-4': {
        'P_surf': -0.000855,
        'H2O': 0.004807,
        'CO2': -0.017428,
        'N2': 0.000136,
        'H2': -0.001186,
        'CO': -0.003537,
        'CH4': 0.039726,
        'NH3': -0.012334,
    },
    'buffered_low_T2000_IW+0': {
        'P_surf': -0.000272,
        'H2O': -0.001511,
        'CO2': 0.004701,
        'N2': -0.000424,
        'H2': 0.001253,
        'CO': -0.001447,
        'CH4': 0.107647,
        'NH3': 0.034908,
    },
    'buffered_low_T2000_IW+2': {
        'P_surf': 0.001022,
        'H2O': 0.000961,
        'CO2': 0.002554,
        'N2': 0.000523,
        'H2': 0.003725,
        'CO': -0.003595,
        'CH4': 0.110445,
        'NH3': 0.039089,
    },
    'buffered_low_T2000_IW+4': {
        'P_surf': 0.000016,
        'H2O': -0.000013,
        'CO2': 0.000303,
        'N2': 0.000014,
        'SO2': 0.000052,
        'H2': 0.002750,
        'CO': -0.005846,
        'CH4': 0.106245,
        'NH3': 0.037373,
        'S2': -0.013904,
    },
    'buffered_low_T2000_IW-2': {
        'P_surf': -0.000197,
        'H2O': -0.002610,
        'CO2': 0.005748,
        'N2': -0.000289,
        'H2': 0.000153,
        'CO': -0.000400,
        'CH4': 0.106496,
        'NH3': 0.033326,
    },
    'buffered_low_T2000_IW-4': {
        'P_surf': -0.000007,
        'H2O': -0.002639,
        'CO2': 0.005993,
        'N2': -0.000075,
        'H2': 0.000125,
        'CO': -0.000155,
        'CH4': 0.106684,
        'NH3': 0.033390,
    },
    'buffered_low_T2500_IW+0': {
        'CO2': 0.004133,
        'N2': 0.000297,
        'CO': -0.000403,
    },
    'buffered_low_T2500_IW+2': {
        'P_surf': 0.000355,
        'H2O': 0.000092,
        'CO2': 0.001739,
        'N2': 0.000155,
        'H2': 0.004019,
        'CO': -0.002797,
        'CH4': 0.108369,
        'NH3': 0.059117,
    },
    'buffered_low_T2500_IW+4': {
        'P_surf': 0.000005,
        'H2O': -0.000029,
        'CO2': 0.000236,
        'N2': 0.000005,
        'SO2': 0.000010,
        'H2': 0.003899,
        'CO': -0.004300,
        'CH4': 0.106625,
        'NH3': 0.058862,
        'S2': -0.012964,
    },
    'buffered_low_T2500_IW-2': {
        'P_surf': 0.000252,
        'CO2': 0.002772,
        'N2': -0.001700,
        'CO': -0.001764,
    },
    'buffered_low_T2500_IW-4': {
        'P_surf': 0.000038,
        'CO2': 0.002656,
        'N2': -0.001886,
        'CO': -0.001880,
    },
    'buffered_low_T3000_IW+0': {
        'CO2': -0.003236,
        'N2': 0.001470,
        'CO': 0.002319,
    },
    'buffered_low_T3000_IW+2': {
        'P_surf': -0.000188,
        'H2O': -0.000440,
        'CO2': -0.002216,
        'N2': -0.000219,
        'SO2': 0.000368,
        'H2': 0.003226,
        'CO': 0.003339,
        'CH4': 0.076289,
        'NH3': 0.064119,
        'S2': -0.014578,
    },
    'buffered_low_T3000_IW+4': {
        'P_surf': -0.000010,
        'H2O': -0.000044,
        'CO2': -0.000312,
        'N2': -0.000016,
        'SO2': -0.000016,
        'H2': 0.003622,
        'CO': 0.005243,
        'CH4': 0.078985,
        'NH3': 0.064814,
        'S2': -0.015345,
    },
    'buffered_low_T3000_IW-2': {
        'P_surf': 0.000772,
        'CO2': -0.007732,
        'N2': -0.002289,
        'CO': -0.002177,
    },
    'buffered_low_T3000_IW-4': {
        'P_surf': 0.000072,
    },
}

AUTHO_EXPECTED_DEX: dict[str, dict[str, float]] = {
    'authO_default_T1500': {
        'P_surf': 0.000003,
        'H2O': 0.000004,
        'CO2': 0.000003,
        'N2': 0.000003,
        'SO2': 0.000003,
        'O2': 0.000098,
        'H2': -0.006038,
        'CO': 0.013845,
        'CH4': 0.047354,
        'NH3': -0.019679,
        'S2': -0.015914,
    },
    'authO_default_T2000': {
        'P_surf': 0.000003,
        'H2O': 0.000002,
        'CO2': 0.000020,
        'N2': 0.000004,
        'SO2': 0.000004,
        'O2': -0.000723,
        'H2': 0.003127,
        'CO': -0.005767,
        'CH4': 0.107438,
        'NH3': 0.037933,
        'S2': -0.012555,
    },
    'authO_default_T2500': {
        'P_surf': 0.000003,
        'H2O': -0.000014,
        'CO2': 0.000125,
        'N2': 0.000004,
        'SO2': 0.000011,
        'O2': -0.000270,
        'H2': 0.004048,
        'CO': -0.004276,
        'CH4': 0.107082,
        'NH3': 0.059085,
        'S2': -0.012422,
    },
    'authO_default_T3000': {
        'P_surf': 0.000003,
        'H2O': -0.000074,
        'CO2': -0.000527,
        'N2': -0.000039,
        'SO2': 0.000043,
        'O2': 0.002620,
        'H2': 0.002282,
        'CO': 0.003718,
        'CH4': 0.073470,
        'NH3': 0.062793,
    },
    'authO_high_T1500': {
        'P_surf': 0.000003,
        'H2O': 0.000003,
        'CO2': 0.000003,
        'N2': 0.000003,
        'SO2': 0.000003,
        'O2': 0.000111,
        'H2': -0.006045,
        'CO': 0.013840,
        'CH4': 0.047329,
        'NH3': -0.019688,
        'S2': -0.015939,
    },
    'authO_high_T2000': {
        'P_surf': 0.000003,
        'H2O': 0.000003,
        'CO2': 0.000011,
        'N2': 0.000004,
        'SO2': 0.000004,
        'O2': -0.000732,
        'H2': 0.003132,
        'CO': -0.005771,
        'CH4': 0.107449,
        'NH3': 0.037941,
        'S2': -0.012537,
    },
    'authO_high_T2500': {
        'P_surf': 0.000003,
        'H2O': -0.000005,
        'CO2': 0.000060,
        'N2': 0.000003,
        'SO2': 0.000011,
        'O2': -0.000406,
        'H2': 0.004126,
        'CO': -0.004273,
        'CH4': 0.107309,
        'NH3': 0.059201,
        'S2': -0.012150,
    },
    'authO_high_T3000': {
        'P_surf': 0.000003,
        'H2O': -0.000032,
        'CO2': -0.000260,
        'N2': -0.000020,
        'SO2': 0.000062,
        'O2': 0.002460,
        'H2': 0.002405,
        'CO': 0.004065,
        'CH4': 0.074142,
        'NH3': 0.062986,
    },
    'authO_low_T1500': {
        'P_surf': 0.000003,
        'H2O': 0.000004,
        'CO2': 0.000002,
        'N2': 0.000003,
        'SO2': 0.000003,
        'O2': 0.000091,
        'H2': -0.006035,
        'CO': 0.013848,
        'CH4': 0.047367,
        'NH3': -0.019674,
        'S2': -0.015901,
    },
    'authO_low_T2000': {
        'P_surf': 0.000003,
        'H2O': -0.000000,
        'CO2': 0.000040,
        'N2': 0.000004,
        'SO2': 0.000004,
        'O2': -0.000711,
        'H2': 0.003119,
        'CO': -0.005753,
        'CH4': 0.107429,
        'NH3': 0.037920,
        'S2': -0.012580,
    },
    'authO_low_T2500': {
        'P_surf': 0.000003,
        'H2O': -0.000035,
        'CO2': 0.000263,
        'N2': 0.000004,
        'SO2': 0.000011,
        'O2': -0.000145,
        'H2': 0.003965,
        'CO': -0.004201,
        'CH4': 0.106929,
        'NH3': 0.058960,
        'S2': -0.012671,
    },
    'authO_low_T3000': {
        'P_surf': 0.000003,
        'H2O': -0.000161,
        'CO2': -0.001038,
        'N2': -0.000076,
        'SO2': 0.000011,
        'O2': 0.002529,
        'H2': 0.002241,
        'CO': 0.003252,
        'CH4': 0.072967,
        'NH3': 0.062713,
    },
}

GAS_SPECIES_STOICH: dict[str, dict[str, int]] = {
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
}


def _import_sibling(module_name: str) -> Any:
    """Import sibling equilibrium package, failing when CHEM_CROSS_CODE_REQUIRED=1 is set."""
    if os.environ.get('CHEM_CROSS_CODE_REQUIRED') == '1':
        return importlib.import_module(module_name)
    return pytest.importorskip(module_name)


def _assert_near_dex(
    actual_dex: float,
    expected_dex: float,
    margin_dex: float,
    label: str,
) -> None:
    """Assert that a signed log10 ratio [dex] lies within expected_dex +- margin_dex."""
    lo_dex = expected_dex - margin_dex
    hi_dex = expected_dex + margin_dex
    assert lo_dex <= actual_dex <= hi_dex, (
        f'{label}: log10 offset {actual_dex:+.6f} dex outside '
        f'[{lo_dex:+.6f}, {hi_dex:+.6f}] dex (expected {expected_dex:+.6f} +- {margin_dex:.4f})'
    )


def _reconstruct_gas_element_masses(
    cal_out: dict[str, float],
    r_int: float,
    g_surf: float,
    elements: tuple[str, ...],
) -> dict[str, float]:
    """Rebuild gas-phase element totals [kg] from solved species partial pressures."""
    p_tot = sum(float(cal_out[f'{sp}_bar']) for sp in GAS_SPECIES_STOICH)
    mu_atm = (
        sum(float(cal_out[f'{sp}_bar']) * float(molar_mass[sp]) for sp in GAS_SPECIES_STOICH)
        / p_tot
    )
    p_to_m = 1.0e5 * 4.0 * math.pi * (r_int**2) / g_surf
    totals = {el: 0.0 for el in elements}
    for sp, stoich in GAS_SPECIES_STOICH.items():
        m_sp = float(molar_mass[sp])
        sp_mass = float(cal_out[f'{sp}_bar']) * p_to_m * (m_sp / mu_atm)
        for el, count in stoich.items():
            if el in totals:
                totals[el] += sp_mass * (count * float(molar_mass[el]) / m_sp)
    return totals


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


def _guarded_log10_ratio(p1: float, p2: float) -> float:
    """Return log10(p1 / p2), or nan if either pressure is non-positive."""
    if p1 <= 0.0 or p2 <= 0.0:
        return float('nan')
    return math.log10(p1 / p2)


def _solve_calliope_gas_only(
    target: dict[str, float],
    ddict: dict[str, float | int],
    autho: bool,
) -> tuple[dict[str, float], bool]:
    """Solve CALLIOPE gas-only equilibrium with nguess=1000, tracking whether preseed is used."""
    np.random.seed(42)
    kw: dict[str, Any] = {
        'rtol': 1.0e-5,
        'atol': 1.0e10,
        'xtol': 1.0e-8,
        'nguess': 1000,
        'opt_solver': False,
        'print_result': False,
    }
    fo2 = float(ddict['fO2_shift_IW'])
    try:
        if autho:
            return (
                equilibrium_atmosphere_authoritative_O(
                    target, ddict, fO2_hint=fo2, random_seed=42, **kw
                ),
                False,
            )
        return equilibrium_atmosphere(target, ddict, **kw), False
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
        return (
            equilibrium_atmosphere_authoritative_O(
                target, ddict, fO2_hint=float(res.x[4]), p_guess=p_g, random_seed=42, **kw
            ),
            True,
        )
    return equilibrium_atmosphere(target, ddict, p_guess=p_g, **kw), True


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

        # Solve gas-only equilibrium and verify reaction quotients in log10 space.
        ddict = _build_gas_only_ddict(
            m_mantle=4.0311e24,
            g_surf=9.930522397640205,
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
            nguess=1000,
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

        product_quotients = {
            'H2': math.log10(q_h2) - ref_dict['janaf_H2'],
            'CO': math.log10(q_co) - ref_dict['janaf_CO'],
            'CH4': math.log10(q_ch4) - ref_dict['schaefer_CH4'],
            'SO2': math.log10(q_so2) - 0.5 * ref_dict['janaf_SO2'],
            'NH3': math.log10(q_nh3) - 0.5 * ref_dict['janaf_NH3'],
            'H2S': math.log10(q_h2s) - 0.5 * ref_dict['janaf_H2S'],
        }
        exp_k = K_OFFSET_EXPECTED_DEX[t_k]
        for sp, d_val in product_quotients.items():
            if sp in K_OFFSET_HELD_FOR_H2S_FIX:
                h2s_report_rows.append(f'{t_label}: d_log10_K_{sp}={d_val:+.6f} dex')
                if int(t_k) == 1500 and sp == 'H2S':
                    assert abs(d_val - 0.287432) <= 0.02, (
                        'CALLIOPE janaf_H2S constant offset has changed; '
                        'H2S-dependent holds must be revisited'
                    )
            else:
                _assert_near_dex(
                    d_val,
                    exp_k[sp],
                    K_OFFSET_MARGIN_DEX,
                    f'{t_label} {sp}',
                )

    h2s_summary = '; '.join(h2s_report_rows)
    record_property('h2s_k_offsets_dex', h2s_summary)
    print(f'H2S equilibrium-constant offsets vs NASA Glenn: {h2s_summary}')


def test_calliope_vs_atmodeller_72_case_gas_only_grid(record_property, tmp_path):
    """Verify CALLIOPE vs atmodeller cross-code agreement over all 72 gas-only grid cases."""
    atmodeller = _import_sibling('atmodeller')
    from atmodeller.containers import FixedFugacityConstraint, SolverParameters
    from molmass import Formula

    mass_ocean = float(ocean_moles * molar_mass['H2'])
    m_earth = 5.972e24
    r_earth = 6.335439e6
    cmf = 0.325
    m_mantle = (1.0 - cmf) * m_earth
    atom_mass_atm = {el: Formula(el).mass * 1.0e-3 for el in ('H', 'C', 'N', 'S', 'O')}

    # Define inventories in moles (with exact 0.5*n_H + 2*n_C + 2*n_S stoichiometric O)
    # and convert to kg per code with that code's atomic masses.
    moles_inventories: dict[str, dict[str, float]] = {}
    for label, scale in (('low', 0.1), ('default', 1.0), ('high', 10.0)):
        h_ref = scale * mass_ocean
        n_h = h_ref / float(molar_mass['H'])
        n_c = h_ref / float(molar_mass['C'])
        n_n = (0.5 * h_ref) / float(molar_mass['N'])
        n_s = h_ref / float(molar_mass['S'])
        n_o = 0.5 * n_h + 2.0 * n_c + 2.0 * n_s
        moles_inventories[label] = {'H': n_h, 'C': n_c, 'N': n_n, 'S': n_s, 'O': n_o}

    buffered_cases: list[dict[str, Any]] = []
    for inv_label in ('low', 'default', 'high'):
        m_inv = moles_inventories[inv_label]
        for t_k in (1500.0, 2000.0, 2500.0, 3000.0):
            for fo2 in (-4.0, -2.0, 0.0, 2.0, 4.0):
                buffered_cases.append(
                    {
                        'case_id': f'buffered_{inv_label}_T{int(t_k)}_IW{fo2:+.0f}',
                        'T_K': t_k,
                        'is_1500_K': int(t_k) == 1500,
                        'fO2_shift_IW': fo2,
                        'cal_targets': {
                            el: m_inv[el] * float(molar_mass[el]) for el in ('H', 'C', 'N', 'S')
                        },
                        'atm_targets': {
                            el: m_inv[el] * atom_mass_atm[el] for el in ('H', 'C', 'N', 'S')
                        },
                    }
                )
    autho_cases: list[dict[str, Any]] = []
    for inv_label in ('low', 'default', 'high'):
        m_inv = moles_inventories[inv_label]
        for t_k in (1500.0, 2000.0, 2500.0, 3000.0):
            autho_cases.append(
                {
                    'case_id': f'authO_{inv_label}_T{int(t_k)}',
                    'T_K': t_k,
                    'is_1500_K': int(t_k) == 1500,
                    'fO2_shift_IW': 2.0,
                    'cal_targets': {
                        el: m_inv[el] * float(molar_mass[el])
                        for el in ('H', 'C', 'N', 'S', 'O')
                    },
                    'atm_targets': {
                        el: m_inv[el] * atom_mass_atm[el] for el in ('H', 'C', 'N', 'S', 'O')
                    },
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

    # Solve 60 buffered cases in atmodeller and take g_surf from Planet.surface_gravity.
    n_buf = len(buffered_cases)
    planet_buf = atmodeller.Planet(
        planet_mass=np.full(n_buf, m_earth),
        core_mass_fraction=np.full(n_buf, cmf),
        mantle_melt_fraction=np.zeros(n_buf),
        surface_radius=np.full(n_buf, r_earth),
        temperature=np.array([c['T_K'] for c in buffered_cases]),
        pressure=np.full(n_buf, np.nan),
    )
    g_surf = float(np.squeeze(planet_buf.surface_gravity)[0])
    f_o2_arr = np.array(
        [float(10.0 ** fo2_fischer(c['T_K'], c['fO2_shift_IW'])) for c in buffered_cases]
    )
    mdl_buf = atmodeller.EquilibriumModel(atmodeller.SpeciesNetwork(gas_species))
    mdl_buf.solve(
        state=planet_buf,
        fugacity_constraints={'O2_g': FixedFugacityConstraint(f_o2_arr)},
        mass_constraints={
            el: np.array([c['atm_targets'][el] for c in buffered_cases])
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
            el: np.array([c['atm_targets'][el] for c in autho_cases])
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

    report_env = os.environ.get('CHEM_SULFUR_REPORT_PATH')
    report_path = Path(report_env) if report_env else (tmp_path / 'sulfur_72_case_report.tsv')
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text('case_id\td_H2S_dex\td_S2_dex\td_SO2_dex\n', encoding='utf-8')

    sulfur_report_rows: list[str] = []
    held_for_h2s_rows: list[str] = []
    preseed_count = 0

    for idx, case in enumerate(buffered_cases):
        cid = case['case_id']
        ddict = _build_gas_only_ddict(
            m_mantle=m_mantle,
            g_surf=g_surf,
            r_int=r_earth,
            t_magma=case['T_K'],
            fo2_shift=case['fO2_shift_IW'],
        )
        cal_out, used_preseed = _solve_calliope_gas_only(
            case['cal_targets'], ddict, autho=False
        )
        if used_preseed:
            preseed_count += 1

        d_h2s_case = _guarded_log10_ratio(
            float(cal_out['H2S_bar']), float(asdict_buf['H2S_g']['pressure'][idx])
        )
        d_s2_case = _guarded_log10_ratio(
            float(cal_out['S2_bar']), float(asdict_buf['S2_g']['pressure'][idx])
        )
        d_so2_case = _guarded_log10_ratio(
            float(cal_out['SO2_bar']), float(asdict_buf['O2S_g']['pressure'][idx])
        )
        sulfur_report_rows.append(
            f'{cid}: d_H2S={d_h2s_case:+.4f}, d_S2={d_s2_case:+.4f}, d_SO2={d_so2_case:+.4f}'
        )
        with report_path.open('a', encoding='utf-8') as f:
            f.write(f'{cid}\t{d_h2s_case:+.6f}\t{d_s2_case:+.6f}\t{d_so2_case:+.6f}\n')

        # Rebuild gas element totals from solved partial pressures and check against target.
        rec_masses = _reconstruct_gas_element_masses(
            cal_out, r_int=r_earth, g_surf=g_surf, elements=('H', 'C', 'N', 'S')
        )
        hcns_max = max(float(case['cal_targets'][k]) for k in ('H', 'C', 'N', 'S'))
        cal_tol = hcns_max * 1.0e-5 + 1.0e10 + float(TRUNC_MASS)
        for el in ('H', 'C', 'N', 'S'):
            assert abs(rec_masses[el] - float(case['cal_targets'][el])) <= cal_tol, (
                f'{cid} reconstructed {el} mass mismatch'
            )

        exp_case = BUFFERED_EXPECTED_DEX[cid]
        held_case = BUFFERED_CASE_KEYS_HELD_FOR_H2S_FIX.get(cid, ())
        d_psurf = math.log10(float(cal_out['P_surf']) / float(ptot_buf[idx]))
        if 'P_surf' in held_case:
            held_for_h2s_rows.append(f'{cid}:P_surf={d_psurf:+.6f}')
        else:
            _assert_near_dex(
                d_psurf,
                exp_case['P_surf'],
                GRID_OFFSET_MARGIN_DEX,
                f'{cid} P_surf',
            )
        for sp_cal, sp_atm in (
            ('H2O', 'H2O_g'),
            ('N2', 'N2_g'),
            ('H2', 'H2_g'),
            ('NH3', 'H3N_g'),
            ('CH4', 'CH4_g'),
            ('CO2', 'CO2_g'),
            ('CO', 'CO_g'),
            ('SO2', 'O2S_g'),
            ('H2S', 'H2S_g'),
            ('S2', 'S2_g'),
        ):
            d_sp = math.log10(
                float(cal_out[f'{sp_cal}_bar']) / float(asdict_buf[sp_atm]['pressure'][idx])
            )
            if sp_cal in BUFFERED_SPECIES_HELD_FOR_H2S_FIX or sp_cal in held_case:
                held_for_h2s_rows.append(f'{cid}:{sp_cal}={d_sp:+.6f}')
            else:
                _assert_near_dex(
                    d_sp,
                    exp_case[sp_cal],
                    GRID_OFFSET_MARGIN_DEX,
                    f'{cid} {sp_cal}',
                )

        q_so2_cal = float(cal_out['SO2_bar']) / (
            (float(cal_out['S2_bar']) ** 0.5) * float(cal_out['O2_bar'])
        )
        q_so2_atm = float(asdict_buf['O2S_g']['pressure'][idx]) / (
            (float(asdict_buf['S2_g']['pressure'][idx]) ** 0.5)
            * float(asdict_buf['O2_g']['pressure'][idx])
        )
        _assert_near_dex(
            math.log10(q_so2_cal / q_so2_atm),
            K_OFFSET_EXPECTED_DEX[case['T_K']]['SO2'],
            K_OFFSET_MARGIN_DEX,
            f'{cid} Q_SO2',
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
        cal_out, used_preseed = _solve_calliope_gas_only(
            case['cal_targets'],
            ddict,
            autho=True,
        )
        if used_preseed:
            preseed_count += 1

        d_h2s_case = _guarded_log10_ratio(
            float(cal_out['H2S_bar']), float(asdict_auth['H2S_g']['pressure'][idx])
        )
        d_s2_case = _guarded_log10_ratio(
            float(cal_out['S2_bar']), float(asdict_auth['S2_g']['pressure'][idx])
        )
        d_so2_case = _guarded_log10_ratio(
            float(cal_out['SO2_bar']), float(asdict_auth['O2S_g']['pressure'][idx])
        )
        sulfur_report_rows.append(
            f'{cid}: d_H2S={d_h2s_case:+.4f}, d_S2={d_s2_case:+.4f}, d_SO2={d_so2_case:+.4f}'
        )
        with report_path.open('a', encoding='utf-8') as f:
            f.write(f'{cid}\t{d_h2s_case:+.6f}\t{d_s2_case:+.6f}\t{d_so2_case:+.6f}\n')

        rec_masses = _reconstruct_gas_element_masses(
            cal_out, r_int=r_earth, g_surf=g_surf, elements=('H', 'C', 'N', 'S', 'O')
        )
        for el in ('H', 'C', 'N', 'S', 'O'):
            cal_el_tol = max(float(case['cal_targets'][el]) * 1.0e-5, float(TRUNC_MASS))
            assert abs(rec_masses[el] - float(case['cal_targets'][el])) <= cal_el_tol, (
                f'{cid} reconstructed {el} mass mismatch'
            )

        if case['is_1500_K']:
            # A 1e-9 relative O error shifts log10 fO2 by at most 2.7e-5 dex
            # (largest at 1500 K).
            cal_o_tgt = float(case['cal_targets']['O'])
            cal_o_rel = abs(rec_masses['O'] - cal_o_tgt) / cal_o_tgt
            assert cal_o_rel <= 1.0e-9, (
                f'{cid} CALLIOPE relative O residual {cal_o_rel:.4e} > 1.0e-9'
            )

            stoich_o = {'H2O_g': 1, 'O2_g': 2, 'CO2_g': 2, 'CO_g': 1, 'O2S_g': 2}
            atm_o_tgt = float(case['atm_targets']['O'])
            atm_o_rec = sum(
                float(asdict_auth[sp_g]['gas_mass'][idx])
                * (
                    stoich_o[sp_g]
                    * float(atom_mass_atm['O'])
                    / float(asdict_auth[sp_g]['molar_mass'][idx])
                )
                for sp_g in stoich_o
            )
            atm_o_rel = abs(atm_o_rec - atm_o_tgt) / atm_o_tgt
            assert atm_o_rel <= 1.0e-9, (
                f'{cid} atmodeller relative O residual {atm_o_rel:.4e} > 1.0e-9'
            )

        exp_case = AUTHO_EXPECTED_DEX[cid]
        _assert_near_dex(
            math.log10(float(cal_out['P_surf']) / float(ptot_auth[idx])),
            exp_case['P_surf'],
            AUTHO_MAJOR_MARGIN_DEX,
            f'{cid} P_surf',
        )
        held_auth_case = AUTHO_CASE_KEYS_HELD_FOR_H2S_FIX.get(cid, ())
        for sp_cal, sp_atm in (
            ('H2O', 'H2O_g'),
            ('CO2', 'CO2_g'),
            ('N2', 'N2_g'),
            ('SO2', 'O2S_g'),
            ('O2', 'O2_g'),
            ('H2', 'H2_g'),
            ('CO', 'CO_g'),
            ('CH4', 'CH4_g'),
            ('NH3', 'H3N_g'),
            ('H2S', 'H2S_g'),
            ('S2', 'S2_g'),
        ):
            if sp_cal in AUTHO_HELD_FOR_H2S_FIX or sp_cal in held_auth_case:
                continue
            margin = (
                AUTHO_MAJOR_MARGIN_DEX
                if sp_cal in ('H2O', 'CO2', 'N2', 'SO2', 'S2')
                else GRID_OFFSET_MARGIN_DEX
            )
            _assert_near_dex(
                math.log10(
                    float(cal_out[f'{sp_cal}_bar'])
                    / float(asdict_auth[sp_atm]['pressure'][idx])
                ),
                exp_case[sp_cal],
                margin,
                f'{cid} {sp_cal}',
            )

        q_so2_cal = float(cal_out['SO2_bar']) / (
            (float(cal_out['S2_bar']) ** 0.5) * float(cal_out['O2_bar'])
        )
        q_so2_atm = float(asdict_auth['O2S_g']['pressure'][idx]) / (
            (float(asdict_auth['S2_g']['pressure'][idx]) ** 0.5)
            * float(asdict_auth['O2_g']['pressure'][idx])
        )
        _assert_near_dex(
            math.log10(q_so2_cal / q_so2_atm),
            K_OFFSET_EXPECTED_DEX[case['T_K']]['SO2'],
            K_OFFSET_MARGIN_DEX,
            f'{cid} Q_SO2',
        )

    record_property('calliope_preseed_count', preseed_count)
    assert preseed_count == 0, (
        f'Expected 0 preseeded solves with nguess=1000, got {preseed_count}'
    )

    record_property('buffered_held_for_h2s_fix_count', len(held_for_h2s_rows))
    record_property('buffered_held_for_h2s_fix_dex', '; '.join(held_for_h2s_rows))

    record_property('sulfur_72_case_report_dex', '\n'.join(sulfur_report_rows))
    print('72-case sulfur species log10(CALLIOPE/atmodeller) report [dex]:')
    for row in sulfur_report_rows:
        print(f'  {row}')
    assert report_path.is_file() and report_path.stat().st_size > 0
