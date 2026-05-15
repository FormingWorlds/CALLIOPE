from __future__ import annotations

import pytest

from calliope.constants import (
    element_list,
    molar_mass,
    ocean_moles,
    volatile_species,
)

pytestmark = pytest.mark.unit

# Per-source tests live in their 1:1-mirrored files (see
# .github/.claude/rules/calliope-tests.md section 12):
#   chemistry.py        -> tests/test_chemistry.py
#   oxygen_fugacity.py  -> tests/test_oxygen_fugacity.py
#   solubility.py       -> tests/test_solubility.py
#   solve.py            -> tests/test_solve.py
#   structure.py        -> tests/test_structure.py
# This file retains the metadata tests on the utility module
# `constants.py`, which is exempt from the 1:1 mirroring requirement.


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
