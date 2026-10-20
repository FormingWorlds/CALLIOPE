"""Cross-code thermodynamic and gas-only equilibrium consistency tests.

Scope: gas-phase cross-code consistency only (Phi_global = 0.0; the melt is
checked for element balance only in test_element_balance_grid.py, and no test
compares melt solubility across codes).

Verifies:
1. CALLIOPE equilibrium constants (ModifiedKeq) and algebraic secondary gas
   speciation against the 9-coefficient NASA Glenn reference values across
   T in {1500, 2000, 2500, 3000} K within +-0.002 dex of the measured
   per-temperature reaction-data offset (K_OFFSET_EXPECTED_DEX), for all six reactions.
2. Live 72-case gas-only equilibrium cross-comparison against atmodeller (60
   buffered cases with matched fO2 plus 12 authoritative-O cases) on the shared
   volatile species (H2O, H2, O2, CO2, CO, CH4, N2, NH3, S2, SO2, H2S). This
   72-case grid is a regression snapshot of thermodynamic data that differ by
   construction: the secondary-species pins follow from the K offsets and the
   primary pins to within 1.5e-6 dex, and the CH4 and NH3 K differences remain
   open. In the 60 buffered cases both codes receive CALLIOPE's
   OxygenFugacity('fischer') fO2 (Fischer et al. (2011) calculate the 1-bar
   Fe-FeO curve from Chase (1998) (Fig. 6) and add the pressure integral of
   their Eq. 2; their supplement (Table S2) gives a P-T polynomial. The source
   of this 2-term form is not stated in the paper body, and no independent
   numerical pin of the Fischer 1-bar fit exists in the repository test data).
   Buffered CALLIOPE solves draw initial guesses from a local default_rng(42)
   copy of get_initial_pressures so multi-start draws do not mutate global NumPy
   RNG state while keeping the pinned roots unchanged. Elemental
   inventories are defined in moles (n_O = 0.5 n_H + 2 n_C + 2 n_S at the
   stoichiometric H2O + CO2 + SO2 endpoint) and converted to mass per code with
   that code's atomic masses (calliope.constants.molar_mass and molmass.Formula),
   and surface gravity is matched to atmodeller.Planet.surface_gravity. Across
   all 12 authoritative-O cases (1500 to 3000 K), solved log10(fO2) (O2 partial
   pressure) is pinned to its per-case measured offset (-0.000732 to +0.002620
   dex) +-0.005 dex (-0.005732 to +0.007620 dex total). H2 follows DeltaK with
   H2O and O2, CO with CO2 and O2, CH4 with CO2, O2 and H2 (d_CH4 = DeltaK +
   d_CO2 - d_O2 + 2 d_H2), NH3 with N2 and H2; the relation holds for the
   pinned values.
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

import calliope.solve as solve_mod
from calliope.chemistry import ModifiedKeq
from calliope.constants import molar_mass, noble_gases, ocean_moles, volatile_species
from calliope.oxygen_fugacity import OxygenFugacity
from calliope.solve import (
    TRUNC_MASS,
    atmosphere_mass,
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
)

logging.getLogger('fwl.calliope').setLevel(logging.WARNING)


def _check_cross_code_required_env(module_name: str) -> None:
    """Validate CHEM_CROSS_CODE_REQUIRED and enforce module presence when set to '1'."""
    req = os.environ.get('CHEM_CROSS_CODE_REQUIRED')
    if req not in (None, '', '0', '1'):
        pytest.fail(f'Invalid CHEM_CROSS_CODE_REQUIRED={req!r}; expected unset, "0", or "1".')
    if req == '1':
        importlib.import_module(module_name)


_check_cross_code_required_env('atmodeller')
_check_cross_code_required_env('molmass')
atmodeller = pytest.importorskip('atmodeller')
molmass = pytest.importorskip('molmass')
from atmodeller.containers import FixedFugacityConstraint, SolverParameters  # noqa: E402
from atmodeller.output import OutputSolution  # noqa: E402
from molmass import Formula  # noqa: E402

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
# asserted within +-K_OFFSET_MARGIN_DEX (0.002 dex).
K_OFFSET_MARGIN_DEX: float = 0.002
K_OFFSET_EXPECTED_DEX: dict[float, dict[str, float]] = {
    1500.0: {
        'H2': -0.005993,
        'CO': 0.013892,
        'CH4': 0.059525,
        'SO2': 0.007862,
        'NH3': 0.012940,
        'H2S': -0.014997,
    },
    2000.0: {
        'H2': 0.002763,
        'CO': -0.006148,
        'CH4': 0.100441,
        'SO2': 0.007005,
        'NH3': 0.030070,
        'H2S': -0.013703,
    },
    2500.0: {
        'H2': 0.003927,
        'CO': -0.004536,
        'CH4': 0.098591,
        'SO2': 0.006492,
        'NH3': 0.033799,
        'H2S': -0.014854,
    },
    3000.0: {
        'H2': 0.003666,
        'CO': 0.005555,
        'CH4': 0.072052,
        'SO2': 0.007657,
        'NH3': 0.029483,
        'H2S': -0.017231,
    },
}

# Measured log10(P_CALLIOPE / P_atmodeller) offsets at matched moles.
GRID_OFFSET_MARGIN_DEX: float = 0.005
AUTHO_MAJOR_MARGIN_DEX: float = 0.001
BUFFERED_EXPECTED_DEX: dict[str, dict[str, float]] = {
    'buffered_default_T1500_IW+0': {
        'P_surf': 0.000048,
        'H2O': 0.003468,
        'CO2': -0.013395,
        'N2': 0.000743,
        'H2': -0.002525,
        'CO': 0.000497,
        'CH4': 0.041081,
        'NH3': 0.009524,
        'SO2': 0.026244,
        'H2S': 0.000859,
        'S2': 0.036763,
    },
    'buffered_default_T1500_IW+2': {
        'P_surf': 0.000135,
        'H2O': 0.000950,
        'CO2': -0.003962,
        'N2': 0.000124,
        'H2': -0.005043,
        'CO': 0.009929,
        'CH4': 0.045478,
        'NH3': 0.005439,
        'SO2': 0.027950,
        'H2S': 0.000048,
        'S2': 0.040175,
    },
    'buffered_default_T1500_IW+4': {
        'P_surf': 0.000624,
        'H2O': 0.000836,
        'N2': 0.000321,
        'H2': -0.005157,
        'NH3': 0.005366,
        'CH4': 0.048979,
        'CO2': -0.000234,
        'CO': 0.013658,
        'SO2': 0.018509,
        'H2S': -0.009507,
        'S2': 0.021292,
    },
    'buffered_default_T1500_IW-2': {
        'P_surf': -0.002583,
        'H2O': 0.002056,
        'CO2': -0.029284,
        'N2': 0.000264,
        'H2': -0.003937,
        'CO': -0.015392,
        'CH4': 0.022369,
        'NH3': 0.007167,
        'SO2': 0.027154,
        'H2S': 0.000358,
        'S2': 0.038584,
    },
    'buffered_default_T1500_IW-4': {
        'P_surf': -0.001929,
        'H2O': 0.003185,
        'CO2': -0.046986,
        'N2': 0.000033,
        'H2': -0.002807,
        'CO': -0.033094,
        'CH4': 0.006925,
        'NH3': 0.008746,
        'SO2': 0.025788,
        'H2S': 0.000121,
        'S2': 0.035852,
    },
    'buffered_default_T2000_IW+0': {
        'P_surf': -0.000453,
        'H2O': -0.001897,
        'CO2': 0.004615,
        'N2': -0.000562,
        'H2': 0.000867,
        'CO': -0.001533,
        'CH4': 0.106790,
        'NH3': 0.031089,
        'SO2': 0.019237,
        'H2S': -0.000604,
        'S2': 0.024465,
    },
    'buffered_default_T2000_IW+2': {
        'P_surf': 0.000268,
        'H2O': 0.000068,
        'CO2': 0.002171,
        'N2': 0.000125,
        'H2': 0.002831,
        'CO': -0.003977,
        'CH4': 0.108275,
        'NH3': 0.034379,
        'SO2': 0.014226,
        'H2S': -0.003650,
        'S2': 0.014444,
    },
    'buffered_default_T2000_IW+4': {
        'P_surf': 0.000026,
        'H2O': -0.000002,
        'CO2': 0.000308,
        'N2': 0.000019,
        'SO2': 0.000163,
        'H2': 0.002762,
        'CO': -0.005840,
        'CH4': 0.106273,
        'NH3': 0.034222,
        'S2': -0.013683,
        'H2S': -0.017783,
    },
    'buffered_default_T2000_IW-2': {
        'P_surf': -0.000475,
        'H2O': -0.003084,
        'CO2': 0.004901,
        'N2': -0.000395,
        'H2': -0.000321,
        'CO': -0.001247,
        'CH4': 0.104701,
        'NH3': 0.029391,
        'SO2': 0.020721,
        'H2S': -0.000308,
        'S2': 0.027433,
    },
    'buffered_default_T2000_IW-4': {
        'P_surf': -0.001929,
        'H2O': -0.005245,
        'CO2': -0.001252,
        'N2': -0.000119,
        'H2': -0.002481,
        'CO': -0.007401,
        'CH4': 0.094227,
        'NH3': 0.026289,
        'SO2': 0.023118,
        'H2S': -0.000071,
        'S2': 0.032226,
    },
    'buffered_default_T2500_IW+0': {
        'P_surf': -0.000629,
        'CO2': 0.003150,
        'N2': -0.000739,
        'CO': -0.001386,
        'H2O': -0.002407,
        'H2': 0.001520,
        'NH3': 0.035710,
        'CH4': 0.104781,
        'SO2': 0.017491,
        'H2S': -0.002335,
        'S2': 0.021998,
    },
    'buffered_default_T2500_IW+2': {
        'CO2': 0.001667,
        'N2': 0.000075,
        'CO': -0.002869,
        'P_surf': 0.000208,
        'H2O': -0.000083,
        'H2': 0.003845,
        'NH3': 0.039604,
        'CH4': 0.107948,
        'SO2': 0.003374,
        'H2S': -0.014127,
        'S2': -0.006236,
    },
    'buffered_default_T2500_IW+4': {
        'P_surf': 0.000005,
        'H2O': -0.000029,
        'CO2': 0.000236,
        'N2': 0.000005,
        'SO2': 0.000010,
        'H2': 0.003899,
        'CO': -0.004300,
        'CH4': 0.106625,
        'NH3': 0.039650,
        'S2': -0.012965,
        'H2S': -0.017437,
    },
    'buffered_default_T2500_IW-2': {
        'P_surf': -0.000426,
        'H2O': -0.003929,
        'CO2': 0.003974,
        'N2': -0.000532,
        'H2': -0.000001,
        'CO': -0.000562,
        'CH4': 0.102563,
        'NH3': 0.033531,
        'SO2': 0.020120,
        'H2S': -0.001227,
        'S2': 0.027255,
    },
    'buffered_default_T2500_IW-4': {
        'P_surf': -0.000153,
        'H2O': -0.004051,
        'CO2': 0.003996,
        'N2': -0.000133,
        'H2': -0.000123,
        'CO': -0.000540,
        'CH4': 0.102340,
        'NH3': 0.033547,
        'SO2': 0.020667,
        'H2S': -0.000802,
        'S2': 0.028351,
    },
    'buffered_default_T3000_IW+0': {
        'CO2': -0.005388,
        'N2': -0.000711,
        'CO': 0.000166,
        'P_surf': -0.000521,
        'H2O': -0.001956,
        'H2': 0.001711,
        'NH3': 0.031694,
        'CH4': 0.070085,
        'SO2': 0.016341,
        'H2S': -0.006835,
        'S2': 0.017369,
    },
    'buffered_default_T3000_IW+2': {
        'P_surf': -0.000203,
        'H2O': -0.000456,
        'CO2': -0.002224,
        'N2': -0.000230,
        'H2': 0.003210,
        'CO': 0.003331,
        'CH4': 0.076249,
        'NH3': 0.034183,
        'SO2': 0.000290,
        'H2S': -0.021387,
        'S2': -0.014734,
    },
    'buffered_default_T3000_IW+4': {
        'P_surf': -0.000029,
        'H2O': -0.000059,
        'CO2': -0.000327,
        'N2': -0.000032,
        'SO2': -0.000031,
        'H2': 0.003607,
        'CO': 0.005228,
        'CH4': 0.078939,
        'NH3': 0.034878,
        'S2': -0.015376,
        'H2S': -0.021312,
    },
    'buffered_default_T3000_IW-2': {
        'P_surf': -0.000450,
        'CO2': -0.005982,
        'N2': -0.000580,
        'CO': -0.000427,
        'H2O': -0.003578,
        'H2': 0.000088,
        'NH3': 0.029325,
        'CH4': 0.066246,
        'SO2': 0.021141,
        'H2S': -0.003659,
        'S2': 0.026968,
    },
    'buffered_default_T3000_IW-4': {
        'P_surf': -0.000042,
        'CO2': -0.005692,
        'N2': -0.000156,
        'CO': -0.000137,
        'H2O': -0.003558,
        'H2': 0.000109,
        'NH3': 0.029568,
        'CH4': 0.066578,
        'SO2': 0.021686,
        'H2S': -0.003093,
        'S2': 0.028058,
    },
    'buffered_high_T1500_IW+0': {
        'P_surf': -0.001691,
        'H2O': 0.000999,
        'CO2': -0.040687,
        'N2': 0.000229,
        'H2': -0.004993,
        'CO': -0.026796,
        'CH4': 0.008852,
        'NH3': 0.005565,
        'SO2': 0.028506,
        'H2S': 0.000653,
        'S2': 0.041287,
    },
    'buffered_high_T1500_IW+2': {
        'P_surf': -0.000380,
        'H2O': 0.000318,
        'CO2': -0.005479,
        'N2': -0.000092,
        'H2': -0.005674,
        'CO': 0.008413,
        'CH4': 0.042699,
        'NH3': 0.004383,
        'SO2': 0.028533,
        'H2S': -0.000001,
        'S2': 0.041341,
    },
    'buffered_high_T1500_IW+4': {
        'P_surf': 0.000163,
        'H2O': 0.000290,
        'CO2': -0.000465,
        'N2': 0.000086,
        'H2': -0.005702,
        'CO': 0.013427,
        'NH3': 0.004430,
        'CH4': 0.047656,
        'SO2': 0.026947,
        'H2S': -0.001614,
        'S2': 0.038170,
    },
    'buffered_high_T1500_IW-2': {
        'P_surf': 0.000040,
        'H2O': 0.005398,
        'CO2': -0.056763,
        'N2': -0.000591,
        'H2': -0.000595,
        'CO': -0.042871,
        'CH4': 0.001573,
        'NH3': 0.011753,
        'SO2': 0.023975,
        'H2S': 0.000520,
        'S2': 0.032224,
    },
    'buffered_high_T1500_IW-4': {
        'P_surf': 0.000055,
        'H2O': 0.005883,
        'CO2': -0.059005,
        'N2': -0.001003,
        'H2': -0.000110,
        'CO': -0.045113,
        'CH4': 0.000302,
        'NH3': 0.012275,
        'SO2': 0.023158,
        'H2S': 0.000189,
        'S2': 0.030591,
    },
    'buffered_high_T2000_IW+0': {
        'P_surf': -0.002218,
        'H2O': -0.004148,
        'CO2': -0.001039,
        'N2': -0.001809,
        'H2': -0.001384,
        'CO': -0.007187,
        'CH4': 0.096634,
        'NH3': 0.027089,
        'SO2': 0.021385,
        'H2S': -0.000707,
        'S2': 0.028761,
    },
    'buffered_high_T2000_IW+2': {
        'P_surf': 0.000047,
        'H2O': -0.000219,
        'CO2': 0.002050,
        'N2': -0.000125,
        'H2': 0.002544,
        'CO': -0.004098,
        'CH4': 0.107580,
        'NH3': 0.033824,
        'SO2': 0.017634,
        'H2S': -0.000530,
        'S2': 0.021259,
    },
    'buffered_high_T2000_IW+4': {
        'P_surf': 0.000133,
        'H2O': 0.000118,
        'CO2': 0.000363,
        'N2': 0.000067,
        'H2': 0.002881,
        'CO': -0.005785,
        'CH4': 0.106567,
        'NH3': 0.034425,
        'SO2': 0.001386,
        'H2S': -0.016441,
        'S2': -0.011237,
    },
    'buffered_high_T2000_IW-2': {
        'P_surf': -0.007045,
        'H2O': -0.011829,
        'CO2': -0.028712,
        'N2': -0.001219,
        'H2': -0.009066,
        'CO': -0.034860,
        'CH4': 0.053598,
        'NH3': 0.015862,
        'SO2': 0.029338,
        'H2S': -0.000435,
        'S2': 0.044667,
    },
    'buffered_high_T2000_IW-4': {
        'P_surf': -0.005106,
        'H2O': -0.009994,
        'CO2': -0.068076,
        'N2': -0.000583,
        'H2': -0.007230,
        'CO': -0.074224,
        'CH4': 0.017905,
        'NH3': 0.018933,
        'SO2': 0.028031,
        'H2S': 0.000093,
        'S2': 0.042054,
    },
    'buffered_high_T2500_IW+0': {
        'P_surf': -0.000820,
        'H2O': -0.002707,
        'CO2': 0.002803,
        'N2': -0.001375,
        'H2': 0.001221,
        'CO': -0.001733,
        'CH4': 0.103835,
        'NH3': 0.034943,
        'SO2': 0.019225,
        'H2S': -0.000901,
        'S2': 0.025466,
    },
    'buffered_high_T2500_IW+2': {
        'P_surf': 0.000348,
        'H2O': 0.000099,
        'N2': 0.000083,
        'H2': 0.004026,
        'NH3': 0.039880,
        'CH4': 0.108388,
        'CO2': 0.001745,
        'CO': -0.002791,
        'SO2': 0.011084,
        'H2S': -0.006235,
        'S2': 0.009185,
    },
    'buffered_high_T2500_IW+4': {
        'P_surf': 0.000009,
        'H2O': -0.000024,
        'CO2': 0.000238,
        'N2': 0.000004,
        'SO2': 0.000051,
        'H2': 0.003903,
        'CO': -0.004298,
        'CH4': 0.106635,
        'NH3': 0.039655,
        'S2': -0.012882,
        'H2S': -0.017392,
    },
    'buffered_high_T2500_IW-2': {
        'P_surf': -0.001759,
        'H2O': -0.005707,
        'CO2': -0.000732,
        'N2': -0.001381,
        'H2': -0.001779,
        'CO': -0.005268,
        'CH4': 0.094301,
        'NH3': 0.030440,
        'SO2': 0.022602,
        'H2S': -0.000523,
        'S2': 0.032220,
    },
    'buffered_high_T2500_IW-4': {
        'P_surf': -0.005451,
        'H2O': -0.011112,
        'CO2': -0.020658,
        'N2': -0.000561,
        'H2': -0.007185,
        'CO': -0.025193,
        'CH4': 0.063564,
        'NH3': 0.022741,
        'SO2': 0.028440,
        'H2S': -0.000091,
        'S2': 0.043895,
    },
    'buffered_high_T3000_IW+0': {
        'P_surf': -0.000766,
        'CO2': -0.005518,
        'N2': -0.001132,
        'CO': 0.000037,
        'H2O': -0.002328,
        'H2': 0.001338,
        'NH3': 0.030925,
        'CH4': 0.069211,
        'SO2': 0.021613,
        'H2S': -0.001937,
        'S2': 0.027911,
    },
    'buffered_high_T3000_IW+2': {
        'CO2': -0.002075,
        'N2': -0.000118,
        'CO': 0.003479,
        'P_surf': 0.000075,
        'H2O': -0.000149,
        'H2': 0.003518,
        'NH3': 0.034701,
        'CH4': 0.077012,
        'SO2': 0.003886,
        'H2S': -0.017484,
        'S2': -0.007541,
    },
    'buffered_high_T3000_IW+4': {
        'P_surf': -0.000035,
        'H2O': -0.000063,
        'CO2': -0.000331,
        'N2': -0.000037,
        'SO2': -0.000030,
        'H2': 0.003604,
        'CO': 0.005224,
        'CH4': 0.078929,
        'NH3': 0.034870,
        'S2': -0.015373,
        'H2S': -0.021314,
    },
    'buffered_high_T3000_IW-2': {
        'P_surf': -0.000630,
        'H2O': -0.003926,
        'CO2': -0.006455,
        'N2': -0.001052,
        'H2': -0.000260,
        'CO': -0.000900,
        'CH4': 0.065077,
        'NH3': 0.028567,
        'SO2': 0.024260,
        'H2S': -0.000887,
        'S2': 0.033207,
    },
    'buffered_high_T3000_IW-4': {
        'P_surf': -0.001182,
        'H2O': -0.005154,
        'CO2': -0.009942,
        'N2': -0.000541,
        'H2': -0.001487,
        'CO': -0.004388,
        'CH4': 0.059135,
        'NH3': 0.026982,
        'SO2': 0.025927,
        'H2S': -0.000448,
        'S2': 0.036539,
    },
    'buffered_low_T1500_IW+0': {
        'P_surf': 0.000986,
        'H2O': 0.004648,
        'CO2': -0.010242,
        'N2': 0.000981,
        'H2': -0.001345,
        'CO': 0.003650,
        'CH4': 0.046594,
        'NH3': 0.011414,
        'SO2': 0.025154,
        'H2S': 0.000950,
        'S2': 0.034584,
    },
    'buffered_low_T1500_IW+2': {
        'P_surf': 0.000178,
        'H2O': 0.001021,
        'CO2': -0.003932,
        'N2': 0.000151,
        'H2': -0.004971,
        'CO': 0.009960,
        'CH4': 0.045652,
        'NH3': 0.005559,
        'SO2': 0.027144,
        'H2S': -0.000687,
        'S2': 0.038563,
    },
    'buffered_low_T1500_IW+4': {
        'P_surf': 0.000442,
        'H2O': 0.000586,
        'N2': 0.000229,
        'H2': -0.005406,
        'NH3': 0.004946,
        'CH4': 0.048387,
        'CO2': -0.000327,
        'CO': 0.013565,
        'SO2': 0.005775,
        'H2S': -0.022490,
        'S2': -0.004174,
    },
    'buffered_low_T1500_IW-2': {
        'P_surf': 0.000310,
        'H2O': 0.005908,
        'CO2': -0.013589,
        'N2': 0.000415,
        'H2': -0.000085,
        'CO': 0.000303,
        'CH4': 0.045767,
        'NH3': 0.013021,
        'SO2': 0.023356,
        'H2S': 0.000412,
        'S2': 0.030987,
    },
    'buffered_low_T1500_IW-4': {
        'P_surf': -0.000855,
        'H2O': 0.004817,
        'CO2': -0.017436,
        'N2': 0.000089,
        'H2': -0.001176,
        'CO': -0.003544,
        'CH4': 0.039739,
        'NH3': 0.011222,
        'SO2': 0.024115,
        'H2S': 0.000079,
        'S2': 0.032504,
    },
    'buffered_low_T2000_IW+0': {
        'P_surf': -0.000390,
        'H2O': -0.001776,
        'CO2': 0.004684,
        'N2': -0.000439,
        'H2': 0.000987,
        'CO': -0.001465,
        'CH4': 0.107099,
        'NH3': 0.031331,
        'SO2': 0.017854,
        'H2S': -0.001867,
        'S2': 0.021698,
    },
    'buffered_low_T2000_IW+2': {
        'P_surf': 0.000401,
        'H2O': 0.000179,
        'CO2': 0.002244,
        'N2': 0.000213,
        'H2': 0.002942,
        'CO': -0.003905,
        'CH4': 0.108570,
        'NH3': 0.034590,
        'SO2': 0.007047,
        'H2S': -0.010718,
        'S2': 0.000085,
    },
    'buffered_low_T2000_IW+4': {
        'P_surf': 0.000014,
        'H2O': -0.000016,
        'CO2': 0.000303,
        'N2': 0.000013,
        'SO2': 0.000028,
        'H2': 0.002748,
        'CO': -0.005847,
        'CH4': 0.106238,
        'NH3': 0.034198,
        'S2': -0.013954,
        'H2S': -0.017933,
    },
    'buffered_low_T2000_IW-2': {
        'P_surf': -0.000212,
        'H2O': -0.002712,
        'CO2': 0.005787,
        'N2': -0.000248,
        'H2': 0.000051,
        'CO': -0.000361,
        'CH4': 0.106331,
        'NH3': 0.030023,
        'SO2': 0.019605,
        'H2S': -0.001052,
        'S2': 0.025200,
    },
    'buffered_low_T2000_IW-4': {
        'P_surf': -0.000007,
        'H2O': -0.002723,
        'CO2': 0.006045,
        'N2': -0.000021,
        'H2': 0.000040,
        'CO': -0.000104,
        'CH4': 0.106567,
        'NH3': 0.030120,
        'SO2': 0.019869,
        'H2S': -0.000799,
        'S2': 0.025729,
    },
    'buffered_low_T2500_IW+0': {
        'CO2': 0.003211,
        'N2': -0.000619,
        'CO': -0.001325,
        'P_surf': -0.000474,
        'H2O': -0.002160,
        'H2': 0.001768,
        'NH3': 0.036141,
        'CH4': 0.105338,
        'SO2': 0.012259,
        'H2S': -0.007319,
        'S2': 0.011534,
    },
    'buffered_low_T2500_IW+2': {
        'P_surf': -0.000020,
        'H2O': -0.000333,
        'CO2': 0.001543,
        'N2': -0.000041,
        'H2': 0.003594,
        'CO': -0.002993,
        'CH4': 0.107323,
        'NH3': 0.039170,
        'SO2': 0.000404,
        'H2S': -0.017348,
        'S2': -0.012176,
    },
    'buffered_low_T2500_IW+4': {
        'P_surf': 0.000005,
        'H2O': -0.000029,
        'CO2': 0.000236,
        'N2': 0.000005,
        'SO2': 0.000005,
        'H2': 0.003898,
        'CO': -0.004300,
        'CH4': 0.106624,
        'NH3': 0.039649,
        'S2': -0.012973,
        'H2S': -0.017442,
    },
    'buffered_low_T2500_IW-2': {
        'P_surf': -0.000398,
        'CO2': 0.003957,
        'N2': -0.000506,
        'CO': -0.000579,
        'H2O': -0.003772,
        'H2': 0.000155,
        'NH3': 0.033779,
        'CH4': 0.102858,
        'SO2': 0.016053,
        'H2S': -0.005138,
        'S2': 0.019122,
    },
    'buffered_low_T2500_IW-4': {
        'P_surf': -0.000027,
        'CO2': 0.004401,
        'N2': -0.000132,
        'CO': -0.000135,
        'H2O': -0.003762,
        'H2': 0.000166,
        'NH3': 0.033982,
        'CH4': 0.103324,
        'SO2': 0.016524,
        'H2S': -0.004656,
        'S2': 0.020064,
    },
    'buffered_low_T3000_IW+0': {
        'CO2': -0.005249,
        'N2': -0.000538,
        'CO': 0.000306,
        'P_surf': -0.000336,
        'H2O': -0.001829,
        'H2': 0.001837,
        'NH3': 0.031970,
        'CH4': 0.070477,
        'SO2': 0.009067,
        'H2S': -0.013983,
        'S2': 0.002821,
    },
    'buffered_low_T3000_IW+2': {
        'P_surf': -0.000242,
        'H2O': -0.000502,
        'CO2': -0.002244,
        'N2': -0.000247,
        'SO2': -0.000193,
        'H2': 0.003164,
        'CO': 0.003310,
        'CH4': 0.076137,
        'NH3': 0.034106,
        'S2': -0.015700,
        'H2S': -0.021916,
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
        'NH3': 0.034909,
        'S2': -0.015346,
        'H2S': -0.021282,
    },
    'buffered_low_T3000_IW-2': {
        'P_surf': -0.000442,
        'CO2': -0.006049,
        'N2': -0.000598,
        'CO': -0.000494,
        'H2O': -0.003456,
        'H2': 0.000210,
        'NH3': 0.029500,
        'CH4': 0.066423,
        'SO2': 0.014223,
        'H2S': -0.010454,
        'S2': 0.013132,
    },
    'buffered_low_T3000_IW-4': {
        'P_surf': -0.000054,
        'H2O': -0.003444,
        'N2': -0.000209,
        'H2': 0.000223,
        'NH3': 0.029712,
        'CH4': 0.066749,
        'CO2': -0.005748,
        'CO': -0.000194,
        'SO2': 0.014727,
        'H2S': -0.009938,
        'S2': 0.014140,
    },
}

AUTHO_EXPECTED_DEX: dict[str, dict[str, float]] = {
    'authO_default_T1500': {
        'P_surf': 0.000003,
        'H2O': 0.000004,
        'CO2': 0.000003,
        'N2': 0.000003,
        'SO2': 0.000003,
        'O2': 0.000085,
        'H2': -0.006032,
        'CO': 0.013852,
        'CH4': 0.047380,
        'NH3': 0.003895,
        'S2': -0.015889,
        'H2S': -0.028973,
    },
    'authO_default_T2000': {
        'P_surf': 0.000003,
        'H2O': 0.000002,
        'CO2': 0.000020,
        'N2': 0.000004,
        'SO2': 0.000004,
        'O2': -0.000718,
        'H2': 0.003124,
        'CO': -0.005769,
        'CH4': 0.107427,
        'NH3': 0.034758,
        'S2': -0.012566,
        'H2S': -0.016863,
    },
    'authO_default_T2500': {
        'P_surf': 0.000003,
        'H2O': -0.000014,
        'CO2': 0.000125,
        'N2': 0.000004,
        'SO2': 0.000005,
        'O2': -0.000181,
        'H2': 0.004004,
        'CO': -0.004319,
        'CH4': 0.106906,
        'NH3': 0.039806,
        'S2': -0.012612,
        'H2S': -0.017156,
    },
    'authO_default_T3000': {
        'P_surf': 0.000003,
        'H2O': -0.000077,
        'CO2': -0.000512,
        'N2': -0.000039,
        'SO2': -0.000030,
        'O2': 0.002875,
        'H2': 0.002152,
        'CO': 0.003605,
        'CH4': 0.072969,
        'NH3': 0.032692,
        'H2S': -0.025641,
        'S2': -0.021125,
    },
    'authO_high_T1500': {
        'P_surf': 0.000003,
        'H2O': 0.000003,
        'CO2': 0.000003,
        'N2': 0.000003,
        'SO2': 0.000003,
        'O2': 0.000085,
        'H2': -0.006031,
        'CO': 0.013853,
        'CH4': 0.047381,
        'NH3': 0.003895,
        'S2': -0.015887,
        'H2S': -0.028972,
    },
    'authO_high_T2000': {
        'P_surf': 0.000003,
        'H2O': 0.000003,
        'CO2': 0.000011,
        'N2': 0.000004,
        'SO2': 0.000004,
        'O2': -0.000721,
        'H2': 0.003126,
        'CO': -0.005777,
        'CH4': 0.107426,
        'NH3': 0.034761,
        'S2': -0.012560,
        'H2S': -0.016857,
    },
    'authO_high_T2500': {
        'P_surf': 0.000003,
        'H2O': -0.000005,
        'CO2': 0.000060,
        'N2': 0.000003,
        'SO2': 0.000004,
        'O2': -0.000221,
        'H2': 0.004033,
        'CO': -0.004364,
        'CH4': 0.106939,
        'NH3': 0.039850,
        'S2': -0.012533,
        'H2S': -0.017088,
    },
    'authO_high_T3000': {
        'P_surf': 0.000003,
        'H2O': -0.000034,
        'CO2': -0.000246,
        'N2': -0.000019,
        'SO2': -0.000010,
        'O2': 0.002971,
        'H2': 0.002146,
        'CO': 0.003823,
        'CH4': 0.073128,
        'NH3': 0.032693,
        'H2S': -0.025722,
        'S2': -0.021276,
    },
    'authO_low_T1500': {
        'P_surf': 0.000003,
        'H2O': 0.000004,
        'CO2': 0.000002,
        'N2': 0.000003,
        'SO2': 0.000003,
        'O2': 0.000085,
        'H2': -0.006032,
        'CO': 0.013851,
        'CH4': 0.047379,
        'NH3': 0.003895,
        'S2': -0.015889,
        'H2S': -0.028973,
    },
    'authO_low_T2000': {
        'P_surf': 0.000003,
        'H2O': -0.000000,
        'CO2': 0.000040,
        'N2': 0.000004,
        'SO2': 0.000004,
        'O2': -0.000708,
        'H2': 0.003117,
        'CO': -0.005753,
        'CH4': 0.107424,
        'NH3': 0.034748,
        'S2': -0.012585,
        'H2S': -0.016879,
    },
    'authO_low_T2500': {
        'P_surf': 0.000003,
        'H2O': -0.000035,
        'CO2': 0.000264,
        'N2': 0.000004,
        'SO2': 0.000005,
        'O2': -0.000103,
        'H2': 0.003944,
        'CO': -0.004221,
        'CH4': 0.106845,
        'NH3': 0.039716,
        'S2': -0.012769,
        'H2S': -0.017295,
    },
    'authO_low_T3000': {
        'P_surf': 0.000003,
        'H2O': -0.000163,
        'CO2': -0.001023,
        'N2': -0.000076,
        'SO2': -0.000067,
        'O2': 0.002663,
        'H2': 0.002172,
        'CO': 0.003200,
        'CH4': 0.072709,
        'NH3': 0.032702,
        'H2S': -0.025446,
        'S2': -0.020774,
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
        f'[{lo_dex:+.6f}, {hi_dex:+.6f}] dex '
        f'(expected {expected_dex:+.6f} +- {margin_dex:.4f})'
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
) -> dict[str, float]:
    """Solve CALLIOPE gas-only equilibrium with nguess=1000 using a local RNG."""
    kw: dict[str, Any] = {
        'rtol': 1.0e-5,
        'atol': 1.0e10,
        'xtol': 1.0e-8,
        'nguess': 1000,
        'opt_solver': False,
        'print_result': False,
    }
    fo2 = float(ddict['fO2_shift_IW'])
    if autho:
        return equilibrium_atmosphere_authoritative_O(
            target, ddict, fO2_hint=fo2, random_seed=42, **kw
        )
    rng = np.random.default_rng(42)
    orig_get_init = solve_mod.get_initial_pressures

    def _seeded_get_init(
        target_d_arg: dict[str, float],
        p_guess_max: float = solve_mod.P_GUESS_MAX_BAR,
        n_extra: int = 0,
    ) -> tuple[float, ...]:
        hi = np.log10(p_guess_max)
        lo = np.log10(solve_mod.P_GUESS_MIN_BAR)
        return tuple(
            float(10.0 ** rng.uniform(low=lo, high=hi)) for _ in range(4 + int(n_extra))
        )

    solve_mod.get_initial_pressures = _seeded_get_init
    try:
        return equilibrium_atmosphere(target, ddict, **kw)
    finally:
        solve_mod.get_initial_pressures = orig_get_init


def test_required_env_guards(monkeypatch):
    """CHEM_CROSS_CODE_REQUIRED accepts unset, '0' and '1' and enforces the module on '1'."""
    for bad_val in ('true', 'yes', '2', 'false'):
        monkeypatch.setenv('CHEM_CROSS_CODE_REQUIRED', bad_val)
        with pytest.raises(pytest.fail.Exception, match='Invalid CHEM_CROSS_CODE_REQUIRED'):
            _check_cross_code_required_env('atmodeller')
    for ok_val in ('0', '1', ''):
        monkeypatch.setenv('CHEM_CROSS_CODE_REQUIRED', ok_val)
        _check_cross_code_required_env('math')
    monkeypatch.setenv('CHEM_CROSS_CODE_REQUIRED', '1')
    with pytest.raises(ModuleNotFoundError):
        _check_cross_code_required_env('_nonexistent_cross_code_module_r62')
    monkeypatch.setenv('CHEM_CROSS_CODE_REQUIRED', '0')
    _check_cross_code_required_env('_nonexistent_cross_code_module_r62')


def test_calliope_equilibrium_constants_and_reaction_quotients_vs_nasa_glenn(record_property):
    """Verify CALLIOPE ModifiedKeq and secondary gas speciation.

    Checked against NASA Glenn reference data.
    """
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
            _assert_near_dex(
                d_val,
                exp_k[sp],
                K_OFFSET_MARGIN_DEX,
                f'{t_label} {sp}',
            )


def test_calliope_vs_atmodeller_72_case_gas_only_grid(record_property, tmp_path):
    """Verify CALLIOPE vs atmodeller cross-code agreement over all 72 gas-only grid cases."""
    report_env = os.environ.get('CHEM_SULFUR_REPORT_PATH')
    report_path = Path(report_env) if report_env else (tmp_path / 'sulfur_72_case_report.tsv')
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text('case_id\td_H2S_dex\td_S2_dex\td_SO2_dex\n', encoding='utf-8')

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
                        'fO2_shift_IW': fo2,
                        'cal_targets': {
                            el: m_inv[el] * float(molar_mass[el])
                            for el in (
                                'H',
                                'C',
                                'N',
                                'S',
                            )
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
        ch4_g_combo = -(g['CH4_g'] + g['O2_g'] - g['CO2_g'] - 2.0 * g['H2_g']) / ln10
        assert ch4_g_combo == pytest.approx(ref_dict['schaefer_CH4'], rel=0.0, abs=1.0e-6)
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
    assert isinstance(mdl_buf.output, OutputSolution)
    buf_success = np.asarray(mdl_buf.output.multi_attempt_solution.solver_success, dtype=bool)
    buf_failed = [buffered_cases[i]['case_id'] for i in np.where(~buf_success)[0]]
    assert np.all(buf_success), f'atmodeller buffered solve failed on {buf_failed}'
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
    assert isinstance(mdl_auth.output, OutputSolution)
    sol_auth = mdl_auth.output.multi_attempt_solution
    auth_success = np.asarray(sol_auth.solver_success, dtype=bool)
    auth_failed = [autho_cases[i]['case_id'] for i in np.where(~auth_success)[0]]
    assert np.all(auth_success), f'atmodeller authO solve failed on {auth_failed}'
    asdict_auth = mdl_auth.output.asdict()
    ptot_auth = np.asarray(mdl_auth.output.total_pressure())

    sulfur_report_rows: list[str] = []
    expected_sulfur: dict[str, tuple[float, float, float]] = {}

    for idx, case in enumerate(buffered_cases):
        cid = case['case_id']
        ddict = _build_gas_only_ddict(
            m_mantle=m_mantle,
            g_surf=g_surf,
            r_int=r_earth,
            t_magma=case['T_K'],
            fo2_shift=case['fO2_shift_IW'],
        )
        cal_out = _solve_calliope_gas_only(case['cal_targets'], ddict, autho=False)

        d_h2s_case = _guarded_log10_ratio(
            float(cal_out['H2S_bar']), float(asdict_buf['H2S_g']['pressure'][idx])
        )
        d_s2_case = _guarded_log10_ratio(
            float(cal_out['S2_bar']), float(asdict_buf['S2_g']['pressure'][idx])
        )
        d_so2_case = _guarded_log10_ratio(
            float(cal_out['SO2_bar']), float(asdict_buf['O2S_g']['pressure'][idx])
        )
        expected_sulfur[cid] = (d_h2s_case, d_s2_case, d_so2_case)
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
        d_psurf = _guarded_log10_ratio(float(cal_out['P_surf']), float(ptot_buf[idx]))
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
            d_sp = _guarded_log10_ratio(
                float(cal_out[f'{sp_cal}_bar']), float(asdict_buf[sp_atm]['pressure'][idx])
            )
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
        cal_out = _solve_calliope_gas_only(
            case['cal_targets'],
            ddict,
            autho=True,
        )

        d_h2s_case = _guarded_log10_ratio(
            float(cal_out['H2S_bar']), float(asdict_auth['H2S_g']['pressure'][idx])
        )
        d_s2_case = _guarded_log10_ratio(
            float(cal_out['S2_bar']), float(asdict_auth['S2_g']['pressure'][idx])
        )
        d_so2_case = _guarded_log10_ratio(
            float(cal_out['SO2_bar']), float(asdict_auth['O2S_g']['pressure'][idx])
        )
        expected_sulfur[cid] = (d_h2s_case, d_s2_case, d_so2_case)
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

        # Measured relative O residuals (12 authO cases): CALLIOPE <= 1.17e-10 (rtol=1e-5),
        # atmodeller <= 4.52e-11 (atol=1e-6); 1e-9 gives 8.5x/22x margins (<= 2.7e-5 dex fO2).
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
            _guarded_log10_ratio(float(cal_out['P_surf']), float(ptot_auth[idx])),
            exp_case['P_surf'],
            AUTHO_MAJOR_MARGIN_DEX,
            f'{cid} P_surf',
        )
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
            d_sp = _guarded_log10_ratio(
                float(cal_out[f'{sp_cal}_bar']), float(asdict_auth[sp_atm]['pressure'][idx])
            )
            margin = (
                AUTHO_MAJOR_MARGIN_DEX
                if sp_cal in ('H2O', 'CO2', 'N2', 'SO2', 'S2')
                else GRID_OFFSET_MARGIN_DEX
            )
            _assert_near_dex(
                d_sp,
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

    record_property('sulfur_72_case_report_dex', '\n'.join(sulfur_report_rows))
    print('72-case sulfur species log10(CALLIOPE/atmodeller) report [dex]:')
    for row in sulfur_report_rows:
        print(f'  {row}')
    assert report_path.is_file()
    report_lines = [
        line.strip()
        for line in report_path.read_text(encoding='utf-8').splitlines()
        if line.strip()
    ]
    n_rows = len(report_lines)
    assert n_rows == 73, f'Expected 1 header + 72 data rows, got {n_rows}'
    expected_header = 'case_id\td_H2S_dex\td_S2_dex\td_SO2_dex'
    assert report_lines[0] == expected_header, (
        f'Header mismatch: got {report_lines[0]!r}, expected {expected_header!r}'
    )
    expected_case_ids = [c['case_id'] for c in buffered_cases] + [
        c['case_id'] for c in autho_cases
    ]
    assert len(set(expected_case_ids)) == 72
    row_case_ids = [line.split('\t')[0] for line in report_lines[1:]]
    assert set(row_case_ids) == set(expected_case_ids), 'Report case ID set mismatch'
    assert len(row_case_ids) == 72, f'Expected 72 rows, got {len(row_case_ids)}'
    s_targets_by_id = {
        c['case_id']: float(c['cal_targets']['S']) for c in buffered_cases + autho_cases
    }
    nan_pair_skips = 0
    for line in report_lines[1:]:
        parts = line.split('\t')
        cid = parts[0]
        vals = tuple(float(x) for x in parts[1:])
        exp_vals = expected_sulfur[cid]
        case_nan_skips = 0
        for v_read, v_exp, sp in zip(vals, exp_vals, ('H2S', 'S2', 'SO2'), strict=True):
            if math.isnan(v_read) and math.isnan(v_exp):
                case_nan_skips += 1
                nan_pair_skips += 1
                continue
            assert abs(v_read - v_exp) <= 1.0e-6, (
                f'{cid} {sp} read-back offset {v_read} != {v_exp}'
            )
        if s_targets_by_id[cid] > 0.0:
            assert case_nan_skips == 0, (
                f'{cid}: expected 0 NaN-pair skips for S target > 0, got {case_nan_skips}'
            )
    # Every grid case has an S target > 0, so the report has zero NaN-pair skips.
    assert nan_pair_skips == 0, f'Expected 0 total NaN-pair skips, got {nan_pair_skips}'
