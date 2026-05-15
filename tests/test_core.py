from __future__ import annotations

import math

import pytest

from calliope.chemistry import ModifiedKeq
from calliope.constants import (
    element_list,
    molar_mass,
    ocean_moles,
    volatile_species,
)

pytestmark = pytest.mark.unit

# Oxygen-fugacity tests live in tests/test_oxygen_fugacity.py per the
# 1:1 source-to-test mirroring rule (see
# .github/.claude/rules/calliope-tests.md section 12).


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


# Structure tests live in tests/test_structure.py per the 1:1
# source-to-test mirroring rule (see .github/.claude/rules/calliope-tests.md
# section 12).


# Solubility tests live in tests/test_solubility.py per the 1:1
# source-to-test mirroring rule (see .github/.claude/rules/calliope-tests.md
# section 12).


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
