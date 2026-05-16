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
    """The `molar_mass` table includes every volatile species CALLIOPE
    tracks (H2O, CO2, H2, CH4, CO, N2, O2, SO2, H2S, S2, NH3) with a
    positive value for each."""
    required = {'H2O', 'CO2', 'H2', 'CH4', 'CO', 'N2', 'O2', 'SO2', 'H2S', 'S2', 'NH3'}
    assert required.issubset(set(molar_mass.keys()))
    # All molar masses should be positive
    assert all(m > 0 for m in molar_mass.values())


def test_volatile_species_and_elements_defined():
    """`volatile_species` is a non-empty list and `element_list` contains
    the five CHNOS atoms PROTEUS budgets are written in."""
    assert isinstance(volatile_species, list) and len(volatile_species) > 0
    assert {'H', 'O', 'C', 'N', 'S'}.issubset(set(element_list))


def test_ocean_moles_positive():
    """One Earth ocean is roughly 7.7e22 moles of H2O. The constant must
    be positive, finite, and in a physically plausible order of magnitude."""
    assert ocean_moles > 0.0
    # Discrimination guard: a stub that returned 1.0 (or any small constant)
    # would pass the bare positivity check. The constant must be in the
    # 1e22 - 1e23 mole range that matches the Earth-ocean reference.
    assert 1e22 < ocean_moles < 1e23
