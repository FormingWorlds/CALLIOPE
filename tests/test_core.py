from __future__ import annotations

import math

import pytest

from calliope.chemistry import ModifiedKeq
from calliope.constants import (
    M_earth,
    R_earth,
    element_list,
    molar_mass,
    ocean_moles,
    volatile_species,
)
from calliope.oxygen_fugacity import OxygenFugacity
from calliope.solubility import SolubilityH2O
from calliope.structure import calculate_mantle_mass

# ---------- Oxygen fugacity tests ----------


def test_oxygen_fugacity_oneill_value_2000K():
    of = OxygenFugacity('oneill')
    val = of(2000.0)
    # Expected from formula: 2*(-244118+115.559*T-8.474*T*ln(T))/(ln(10)*8.31441*T)
    # For T=2000 K this is -7.407823842131363
    assert val == pytest.approx(-7.407823842131363, rel=1e-3, abs=5e-3)


def test_oxygen_fugacity_fischer_value_2000K():
    of = OxygenFugacity('fischer')
    val = of(2000.0)
    # 6.94059 - (28.1808e3)/T -> ~ -7.14981 at 2000 K
    assert val == pytest.approx(-7.1498, rel=1e-4, abs=1e-3)


def test_oxygen_fugacity_shift_is_additive():
    of = OxygenFugacity('oneill')
    base = of(1800.0, 0.0)
    shifted = of(1800.0, 0.75)
    assert shifted == pytest.approx(base + 0.75, abs=1e-10)


def test_oxygen_fugacity_monotonic_in_T():
    of = OxygenFugacity('oneill')
    low = of(1500.0)
    high = of(3000.0)
    # Becomes less negative (increases) with temperature for this model
    assert high > low


# ---------- Modified equilibrium constant tests ----------


def test_modified_keq_janaf_H2_numeric_and_shift_dependence():
    T = 2000.0
    mk = ModifiedKeq('janaf_H2', fO2_model='oneill')
    g0 = mk(T, fO2_shift=0.0)
    # Precomputed expectation ~1.469 at 2000 K with oneill
    assert g0 == pytest.approx(1.469, rel=1e-2, abs=2e-2)

    # Increasing fO2 (more oxidizing) should decrease Geq for positive fO2 stoichiometry
    g_shift = mk(T, fO2_shift=+1.0)
    assert g_shift < g0

    # Different fO2 model should change result
    mk2 = ModifiedKeq('janaf_H2', fO2_model='fischer')
    g_fischer = mk2(T, fO2_shift=0.0)
    assert g_fischer != pytest.approx(g0, rel=1e-6)


def test_modified_keq_janaf_CO_positive():
    T = 1800.0
    mk = ModifiedKeq('janaf_CO', fO2_model='oneill')
    g = mk(T, fO2_shift=0.0)
    assert g > 0.0
    assert math.isfinite(g)


def test_modified_keq_schaefer_models_return_finite():
    T = 2200.0
    mk_h = ModifiedKeq('schaefer_H', fO2_model='oneill')
    mk_c = ModifiedKeq('schaefer_C', fO2_model='oneill')
    mk_ch4 = ModifiedKeq('schaefer_CH4', fO2_model='oneill')
    assert math.isfinite(mk_h(T, 0.0))
    assert math.isfinite(mk_c(T, 0.0))
    assert math.isfinite(mk_ch4(T, 0.0))


# ---------- Structure tests ----------


def test_calculate_mantle_mass_earth_like():
    # Using Earth values and the internal earth_fr, earth_fm in model:
    # mantle_mass ≈ (1 - 0.325) * M_earth
    mantle = calculate_mantle_mass(R_earth, M_earth, corefrac=0.55)
    assert mantle == pytest.approx((1.0 - 0.325) * M_earth, rel=1e-6)


def test_calculate_mantle_mass_negative_raises():
    # Choose a mass smaller than the implied core mass to trigger exception
    with pytest.raises(Exception):
        calculate_mantle_mass(R_earth, 0.1 * M_earth, corefrac=0.55)


def test_calculate_mantle_mass_decreases_with_corefrac():
    m1 = calculate_mantle_mass(R_earth, M_earth, corefrac=0.4)
    m2 = calculate_mantle_mass(R_earth, M_earth, corefrac=0.6)
    assert m2 < m1


# ---------- Solubility tests (H2O-only; other species are incomplete) ----------


def test_solubility_h2o_default_peridotite_sqrt_law():
    s = SolubilityH2O()  # default peridotite
    # peridotite: 524 * p^0.5
    assert s(0.0) == 0.0
    assert s(100.0) == pytest.approx(524.0 * 10.0, rel=1e-12)


def test_solubility_h2o_other_parameterizations():
    s = SolubilityH2O('basalt_dixon')
    assert s(100.0) == pytest.approx(965.0 * 10.0, rel=1e-12)

    s = SolubilityH2O('basalt_wilson')
    # 215 * p^0.7, with p=100 -> 215 * 10^1.4
    expected = 215.0 * (100.0**0.7)
    assert s(100.0) == pytest.approx(expected, rel=1e-12)

    s = SolubilityH2O('anorthite_diopside')
    assert s(100.0) == pytest.approx(727.0 * 10.0, rel=1e-12)

    s = SolubilityH2O('lunar_glass')
    assert s(100.0) == pytest.approx(683.0 * 10.0, rel=1e-12)


# ---------- Constants and metadata tests ----------


def test_molar_mass_contains_expected_species():
    required = {'H2O', 'CO2', 'H2', 'CH4', 'CO', 'N2', 'O2', 'SO2', 'H2S', 'S2', 'NH3'}
    assert required.issubset(set(molar_mass.keys()))
    # All molar masses should be positive
    assert all(m > 0 for m in molar_mass.values())


def test_volatile_species_and_elements_defined():
    assert isinstance(volatile_species, list) and len(volatile_species) > 0
    assert {'H', 'O', 'C', 'N', 'S'}.issubset(set(element_list))


def test_ocean_moles_positive():
    assert ocean_moles > 0.0
