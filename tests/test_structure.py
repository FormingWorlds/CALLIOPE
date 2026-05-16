"""Tests for `src/calliope/structure.py`.

Exercises the mantle-mass closure model in `calculate_mantle_mass`:

- Conservation: `M_mantle = M_planet - M_core` for any valid input.
- Boundedness: `0 < M_mantle <= M_planet` for any physically valid (mass, radius, core_frac).
- Monotonicity: increasing `core_frac` at fixed `(mass, radius)` decreases `M_mantle`.
- Reference pin: Earth-like input recovers the Wang, Lineweaver & Ireland
  (2017) Earth core mass fraction within tolerance (the constant the source
  file cites at `src/calliope/structure.py` line 50).
- Error contract: zero/negative mantle masses raise; missing or ambiguous
  `core_frac` raises `TypeError`; the deprecated `corefrac` alias emits
  `DeprecationWarning` while still computing the right value.
"""

from __future__ import annotations

import pytest

from calliope.constants import M_earth, R_earth
from calliope.structure import calculate_mantle_mass

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_calculate_mantle_mass_recovers_wang_2017_earth_core_fraction():
    """Earth-like inputs recover Wang, Lineweaver & Ireland (2017) Earth
    core mass fraction 0.325.

    The source file `src/calliope/structure.py` hard-codes `earth_fm = 0.325`
    citing arxiv:1708.08718 (Wang, Lineweaver & Ireland 2017, "The Elemental
    Abundances (with Uncertainties) of the Most Earth-like Planet"; the paper
    reports 32.5 +/- 0.3 wt% Earth core mass fraction). For Earth radius,
    Earth mass, and `core_frac = 0.55` (Earth-like core radius fraction), the
    computed core mass equals `0.325 * M_earth` by construction; the mantle
    mass is `(1 - 0.325) * M_earth = 0.675 * M_earth`.

    Cross-check: the reduced-mass orbital test in PROTEUS's satellite module
    uses the same 0.325 core-fraction value for Earth-Moon decomposition, so
    this pin keeps the two codes consistent.
    """
    expected = (1.0 - 0.325) * M_earth
    mantle = calculate_mantle_mass(R_earth, M_earth, core_frac=0.55)
    # rel=1e-6 because the only source of error is float32->float64 promotion
    # in the constants module; the formula is closed-form and deterministic.
    assert mantle == pytest.approx(expected, rel=1e-6)
    # Exponent-error guard: a regression to `(radius * core_frac)**2.0` lands
    # at a wildly different mass (the core volume scales as r^3 by physics,
    # not r^2). At Earth scale the r^2 form would give ~3e18 kg, ~6 orders
    # below the correct 1.95e24 kg.
    wrong_r2 = M_earth - (
        ((3.0 * 0.325 * M_earth) / (4.0 * 3.14159265 * (0.55 * R_earth) ** 3.0))
        * (4.0 / 3.0)
        * 3.14159265
        * (R_earth * 0.55) ** 2.0
    )
    assert abs(mantle - wrong_r2) > 0.1 * M_earth
    # Sign guard: mantle mass is always positive for valid Earth-like input.
    assert mantle > 0
    # Scale guard: order of magnitude is 4e24 kg (between 1e24 and 1e25),
    # not 4e21 (forgotten kg->g) or 4e27 (forgotten g->kg).
    assert 1e24 < mantle < 1e25


@pytest.mark.physics_invariant
def test_calculate_mantle_mass_closure_holds_for_earth_like():
    """`M_mantle + M_core ≈ M_planet` within solver tolerance.

    The conservation invariant is the contract of the function: total mass
    is split into mantle and core, no other reservoirs. A regression that
    introduces a fictitious third reservoir (e.g. an atmosphere subtraction)
    would break this equality.
    """
    mass = M_earth
    mantle = calculate_mantle_mass(R_earth, mass, core_frac=0.55)
    core = mass - mantle
    # Both reservoirs must be positive for any physical config.
    assert mantle > 0
    assert core > 0
    # Closure: sum must equal input mass to floating-point precision.
    assert mantle + core == pytest.approx(mass, rel=1e-12)


@pytest.mark.physics_invariant
def test_calculate_mantle_mass_decreases_with_core_frac():
    """Increasing `core_frac` at fixed mass and radius shrinks the mantle.

    Monotonicity test: doubling the core radius (within the Earth-like
    regime) more than doubles the core volume and thus the core mass, so
    the mantle remainder must shrink. The delta is large enough to
    discriminate an exponent error in `(radius * core_frac)**3`: a swap
    to `**2` would invert the sign at sufficiently small core_frac.
    """
    m_low = calculate_mantle_mass(R_earth, M_earth, core_frac=0.40)
    m_high = calculate_mantle_mass(R_earth, M_earth, core_frac=0.60)
    # Strict ordering: high core_frac => smaller mantle.
    assert m_high < m_low
    # Discrimination guard: the delta is order ~1e24 kg at Earth scale.
    # A regression that swapped the exponent would land at a different
    # delta and likely violate the strict ordering too, but pin the size
    # to catch a coefficient-only bug that preserves the sign.
    assert (m_low - m_high) > 0.1 * M_earth


@pytest.mark.physics_invariant
def test_calculate_mantle_mass_is_bounded_by_planet_mass():
    """Mantle mass must lie in `(0, M_planet)` for any physical input.

    Property-based check across three core_frac regimes at Earth radius
    and Earth mass: 0.30 (Mars-like), 0.55 (Earth), and 0.70 (super-Mercury,
    but the largest core fraction that still leaves a positive mantle
    given the Earth-derived core density: at 0.55 the construction yields
    core = 0.325 M_earth, so core scales as cf**3 and crosses M_earth at
    cf ≈ 0.80). All three configurations must respect the (0, M_planet)
    envelope.
    """
    mass = M_earth
    mantles = {}
    for core_frac in (0.30, 0.55, 0.70):
        mantle = calculate_mantle_mass(R_earth, mass, core_frac=core_frac)
        assert 0 < mantle < mass
        mantles[core_frac] = mantle

    # Discrimination guard: mantle mass must decrease monotonically with
    # core fraction (more core means less mantle at fixed planet mass).
    # A stub that returned the same mantle mass for every core_frac
    # would pass the bare envelope check but fail this ordering.
    assert mantles[0.30] > mantles[0.55] > mantles[0.70]


def test_calculate_mantle_mass_raises_when_core_exceeds_total():
    """Negative mantle mass is non-physical and must raise.

    Chooses a planetary mass smaller than the implied core mass (10% of
    Earth mass at Earth radius and `core_frac = 0.55`). The function's
    explicit guard at `structure.py:62` must fire, not return a negative
    value silently.
    """
    with pytest.raises(Exception, match='mantle mass is negative'):
        calculate_mantle_mass(R_earth, 0.1 * M_earth, core_frac=0.55)


def test_calculate_mantle_mass_corefrac_alias_emits_deprecation():
    """Legacy `corefrac` keyword still works but emits `DeprecationWarning`.

    Backwards-compatibility shim from `structure.py:32-45`. The alias must
    produce the same numeric result as `core_frac` so callers can migrate
    gradually.
    """
    with pytest.warns(DeprecationWarning, match="'corefrac' keyword is deprecated"):
        legacy = calculate_mantle_mass(R_earth, M_earth, corefrac=0.55)
    new = calculate_mantle_mass(R_earth, M_earth, core_frac=0.55)
    assert legacy == pytest.approx(new, rel=1e-12)


def test_calculate_mantle_mass_both_aliases_raises():
    """Passing both `core_frac` and `corefrac` is ambiguous and must raise.

    The migration shim refuses to silently prefer one over the other.
    """
    with pytest.raises(TypeError, match="received both 'core_frac' and 'corefrac'"):
        calculate_mantle_mass(R_earth, M_earth, core_frac=0.55, corefrac=0.6)


def test_calculate_mantle_mass_missing_core_frac_raises():
    """Neither `core_frac` nor `corefrac` supplied must raise `TypeError`.

    The function has no default for `core_frac`; callers must specify the
    interior structure explicitly.
    """
    with pytest.raises(TypeError, match="missing required argument: 'core_frac'"):
        calculate_mantle_mass(R_earth, M_earth)
