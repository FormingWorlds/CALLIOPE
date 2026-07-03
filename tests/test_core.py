from __future__ import annotations

import pytest

from calliope.constants import (
    element_list,
    element_list_chnos,
    molar_mass,
    noble_gases,
    ocean_moles,
    volatile_species,
)

pytestmark = pytest.mark.unit

# Per-source tests live in their 1:1-mirrored files:
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
    """`volatile_species` is a non-empty list and `element_list` is the
    complete element registry: the five CHNOS atoms PROTEUS budgets are
    written in, plus the noble gases."""
    assert isinstance(volatile_species, list) and len(volatile_species) > 0
    assert {'H', 'O', 'C', 'N', 'S'}.issubset(set(element_list))
    # element_list is CHNOS plus the noble gases, and the two subsets are
    # disjoint (a noble gas is never a reacting CHNOS element).
    assert set(noble_gases).issubset(set(element_list))
    assert set(element_list_chnos).isdisjoint(set(noble_gases))
    assert set(element_list) == set(element_list_chnos) | set(noble_gases)


def test_ocean_moles_positive():
    """One Earth ocean is roughly 7.7e22 moles of H2O. The constant must
    be positive, finite, and in a physically plausible order of magnitude."""
    assert ocean_moles > 0.0
    # Discrimination guard: a stub that returned 1.0 (or any small constant)
    # would pass the bare positivity check. The constant must be in the
    # 1e22 - 1e23 mole range that matches the Earth-ocean reference.
    assert 1e22 < ocean_moles < 1e23
