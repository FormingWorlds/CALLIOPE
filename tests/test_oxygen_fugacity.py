"""Tests for `src/calliope/oxygen_fugacity.py`.

Exercises the `OxygenFugacity` IW-buffer dispatcher and its two
underlying fits:

- Reference pin: Fischer et al. (2011) IW value at T = 2000 K against
  the closed-form `6.94059 - 28.1808e3 / T`, with a discrimination
  guard against the O'Neill & Eggins (2002) IW at the same T. The
  guard catches a regression that silently dispatches to the wrong
  buffer (the canonical buffer-flip trap from
  `.github/.claude/rules/calliope-tests.md` Section 16).
- Monotonicity: `log10(fO2)` is monotonic in T along each buffer
  over the 1500-3000 K range.
- Symmetry: `fO2_shift` is strictly additive: `of(T, dIW) = of(T, 0) + dIW`.
- Boundedness: `log10(fO2)` is finite for any valid (T > 0, dIW finite).
- Error contract: T <= 0 raises `ValueError` mentioning the divergence
  in the underlying formulae; unknown buffer names raise
  `AttributeError` at construction.

See `.github/.claude/rules/calliope-tests.md` sections 1-3 for the
anti-happy-path, discrimination-guard, and physics-invariant rules.
"""

from __future__ import annotations

import math

import pytest

from calliope.oxygen_fugacity import OxygenFugacity

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_oxygen_fugacity_fischer_value_at_2000K_matches_published_fit():
    """Fischer 2011 IW at T = 2000 K equals 6.94059 - 28.1808e3 / 2000.

    The source file `src/calliope/oxygen_fugacity.py` line 28 implements
    Fischer et al. (2011, EPSL 304, 496) Eq. 2 as `6.94059 - 28.1808e3 / T`.
    At T = 2000 K this gives `log10(fO2) = -7.14981`.

    Cross-buffer discrimination guard: O'Neill & Eggins (2002) at the
    same T gives ~-7.4078 (a 0.26 dex offset). A regression that
    silently dispatches to 'oneill' instead of 'fischer' would land
    outside the tolerance, catching the buffer-flip trap from
    `calliope-tests.md` Section 16.
    """
    of = OxygenFugacity('fischer')
    val = of(2000.0)
    expected = 6.94059 - 28.1808e3 / 2000.0  # -7.14981
    # rel=1e-4 matches the published precision of the Fischer 2011 fit.
    assert val == pytest.approx(expected, rel=1e-4)
    # Wrong-buffer guard: O'Neill at 2000 K gives -7.4078; a regression
    # to that buffer would land 0.26 dex away.
    wrong_oneill = -7.4078
    assert abs(val - wrong_oneill) > 0.2
    # Sign guard: log10(fO2) at the IW buffer is always negative under
    # standard conditions (T < ~10000 K).
    assert val < 0
    # Scale guard: order of magnitude is -7, not -70 (forgotten log10)
    # or -0.7 (factor-10 unit slip on the temperature coefficient).
    assert -10 < val < -3


@pytest.mark.physics_invariant
def test_oxygen_fugacity_oneill_value_at_2000K_matches_published_fit():
    """O'Neill & Eggins (2002) IW at T = 2000 K matches the closed form.

    Implements Eq. 11 of O'Neill & Eggins (2002, J. Chem. Thermodyn. 34, 1311):
    `2 * (-244118 + 115.559*T - 8.474*T*ln(T)) / (ln(10) * 8.31441 * T)`.
    At T = 2000 K this evaluates to ~-7.4078. The discrimination guard
    against Fischer is the mirror of the test above.
    """
    of = OxygenFugacity('oneill')
    val = of(2000.0)
    expected = -7.407823842131363
    assert val == pytest.approx(expected, rel=1e-3, abs=5e-3)
    # Wrong-buffer guard: Fischer at 2000 K is -7.14981.
    wrong_fischer = -7.14981
    assert abs(val - wrong_fischer) > 0.2
    # Sign and scale guards.
    assert val < 0
    assert -10 < val < -3


@pytest.mark.physics_invariant
def test_oxygen_fugacity_shift_is_strictly_additive():
    """`fO2_shift` adds linearly to the buffer value for both buffers.

    The dispatcher at `__call__` returns `callmodel(T) + fO2_shift`. A
    regression that multiplies instead of adds, or that applies the
    shift twice, would break this identity.
    """
    for buffer_name in ('fischer', 'oneill'):
        of = OxygenFugacity(buffer_name)
        base = of(1800.0, fO2_shift=0.0)
        plus_half = of(1800.0, fO2_shift=0.5)
        plus_three = of(1800.0, fO2_shift=3.0)
        # Strict additivity to floating-point precision.
        assert plus_half == pytest.approx(base + 0.5, abs=1e-12)
        assert plus_three == pytest.approx(base + 3.0, abs=1e-12)
        # Negative shifts work too (reduced fO2).
        minus_two = of(1800.0, fO2_shift=-2.0)
        assert minus_two == pytest.approx(base - 2.0, abs=1e-12)


@pytest.mark.physics_invariant
def test_oxygen_fugacity_monotonic_in_T():
    """Both buffers produce `log10(fO2)` strictly increasing with T over 1500-3000 K.

    A hot magma ocean is more reducing on an absolute scale than a cool
    one at the IW buffer, but the IW buffer itself is defined by the
    Fe-FeO equilibrium and its log10(fO2) becomes less negative (closer
    to zero) as T increases. The chosen window spans realistic surface
    temperatures from solidification (~1500 K) to early Earth magma
    ocean (~3000 K), giving a delta large enough to resolve a regression
    that flipped the sign of the slope.
    """
    for buffer_name in ('fischer', 'oneill'):
        of = OxygenFugacity(buffer_name)
        low = of(1500.0)
        mid = of(2250.0)
        high = of(3000.0)
        # Strict ordering: monotonic, not just non-decreasing.
        assert low < mid < high
        # Discrimination: the 1500 K -> 3000 K delta is ~2 dex for
        # Fischer (28180.8/1500 - 28180.8/3000 = 9.39). A regression
        # that flipped the slope sign would invert the ordering above
        # but pin the magnitude to catch a coefficient-only bug too.
        assert (high - low) > 1.0


@pytest.mark.physics_invariant
def test_oxygen_fugacity_finite_over_realistic_temperature_range():
    """`log10(fO2)` is finite for any T in the realistic 800-5000 K range.

    Boundedness check: no nan, no inf, no complex intermediate. Catches a
    regression that introduces a `log(negative)` or `0/0` along an
    unexpected code path.
    """
    for buffer_name in ('fischer', 'oneill'):
        of = OxygenFugacity(buffer_name)
        for T in (800.0, 1500.0, 2500.0, 3500.0, 5000.0):
            for dIW in (-3.0, 0.0, 1.5):
                val = of(T, fO2_shift=dIW)
                assert math.isfinite(val)
                # Bounded order-of-magnitude: log10(fO2) on the IW buffer
                # stays in [-40, +5] over the 800-5000 K window for any
                # dIW in [-3, +1.5]. The lower edge is set by the
                # T = 800 K / dIW = -3 corner (~-31). Pin the envelope so
                # a unit-conversion bug (e.g. log10 -> ln, factor 2.3)
                # surfaces.
                assert -40 < val < 5


@pytest.mark.parametrize('bad_T', [0.0, -1.0, -300.0])
def test_oxygen_fugacity_nonpositive_T_raises(bad_T):
    """`T <= 0` raises `ValueError` for both buffers.

    Without the guard, T = 0 silently propagates `nan` through every
    downstream equilibrium constant via the `1/T` and `T * log(T)` terms.
    The error message must name the temperature so users debugging an IC
    file can find the offending input.
    """
    for buffer_name in ('fischer', 'oneill'):
        of = OxygenFugacity(buffer_name)
        with pytest.raises(ValueError, match='Temperature must be positive'):
            of(bad_T)


def test_oxygen_fugacity_unknown_buffer_name_raises():
    """Constructing with an unknown buffer name raises `AttributeError`.

    The dispatcher uses `getattr(self, model)` so a typo lands as an
    `AttributeError` at construction (eager), not at first call (lazy).
    Eager failure is preferable: it surfaces config typos before any
    chemistry step runs.
    """
    with pytest.raises(AttributeError):
        OxygenFugacity('hirschmann')  # not implemented in CALLIOPE
    with pytest.raises(AttributeError):
        OxygenFugacity('typo_fisher')
