"""Property-based exploration tests for CALLIOPE solubility laws.

These tests use the ``hypothesis`` library to fuzz the input space of
the Dasgupta N2 and Gaillard S2 solubility laws and the O'Neill IW
buffer, asserting that the invariants pinned in test_invariants.py
hold over the full physically-relevant parameter space (not just the
fixed parametric points).

The tests are marked ``slow`` because each hypothesis-driven test runs
many iterations (default 200) and the cumulative wall-time is several
seconds per test. They are not part of the PR gate; the nightly run
exercises them.
"""

from __future__ import annotations

import math

import pytest

# Gate the `hypothesis` import behind ``pytest.importorskip``: the
# Docker-based PR image installs with ``pip install --no-deps`` and
# would otherwise fail collection on an unconditional ``import
# hypothesis`` at module top.
pytest.importorskip('hypothesis')

from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

from calliope.chemistry import ModifiedKeq  # noqa: E402
from calliope.oxygen_fugacity import OxygenFugacity  # noqa: E402
from calliope.solubility import SolubilityN2, SolubilityS2  # noqa: E402

# Fixed-seed (derandomized) profile so the property-based exploration
# replays identically across hypothesis versions. A failure surfaced
# here must be reproducible from the same input sequence rather than
# depending on the per-run seed strategy.
settings.register_profile('calliope_deterministic', derandomize=True)
settings.load_profile('calliope_deterministic')

pytestmark = [pytest.mark.slow, pytest.mark.timeout(3600)]


# Physical bounds for the fuzzing strategies. These cover the
# magma-ocean regime CALLIOPE is targeted at, intentionally wider than
# the calibration footprints of the individual laws so the tests also
# exercise extrapolation behaviour.
T_RANGE = (1500.0, 2500.0)
DIW_RANGE = (-6.0, +6.0)
P_RANGE = (1e-6, 1e4)
P_TOT_RANGE = (1e-6, 1e5)


# ===========================================================================
# Dasgupta N2: monotonicity + finiteness across the physical space
# ===========================================================================


@given(
    p_N2=st.floats(min_value=P_RANGE[0], max_value=P_RANGE[1]),
    p_tot=st.floats(min_value=P_TOT_RANGE[0], max_value=P_TOT_RANGE[1]),
    T=st.floats(min_value=T_RANGE[0], max_value=T_RANGE[1]),
    dIW_low=st.floats(min_value=DIW_RANGE[0], max_value=-0.1),
    dIW_high=st.floats(min_value=0.1, max_value=DIW_RANGE[1]),
)
@settings(deadline=None, max_examples=200)
def test_dasgupta_monotonic_with_oxidation(p_N2, p_tot, T, dIW_low, dIW_high):
    """For any (p_N2, p_tot, T), the Dasgupta solubility at a more
    reducing fO2 is >= the value at a more oxidising fO2."""
    N2 = SolubilityN2('dasgupta')
    val_low = N2.dasgupta(p_N2, p_tot, T, dIW_low)
    val_high = N2.dasgupta(p_N2, p_tot, T, dIW_high)
    assert val_low >= val_high, (
        f'Dasgupta at dIW={dIW_low:.2f} ({val_low:.4e}) is less than '
        f'at dIW={dIW_high:.2f} ({val_high:.4e}) '
        f'(p_N2={p_N2:.4e}, p_tot={p_tot:.4e}, T={T:.1f})'
    )

    # Discrimination guard: monotonicity alone is satisfied by a constant
    # function. When the dIW separation is large enough (>= 4 dex), the
    # reduced-N branch (which carries an exp(-1.6 dIW) factor) gives a
    # meaningful gap between the two values.
    if (dIW_high - dIW_low) >= 4.0:
        assert val_low > val_high or val_low == val_high == 0.0, (
            f'Dasgupta values at dIW_low={dIW_low} and dIW_high={dIW_high} '
            f'are equal but nonzero ({val_low:.4e}); expected strict '
            f'inequality for a 4-dex dIW span'
        )


@given(
    p_N2=st.floats(min_value=P_RANGE[0], max_value=P_RANGE[1]),
    p_tot=st.floats(min_value=P_TOT_RANGE[0], max_value=P_TOT_RANGE[1]),
    T=st.floats(min_value=T_RANGE[0], max_value=T_RANGE[1]),
    dIW=st.floats(min_value=DIW_RANGE[0], max_value=DIW_RANGE[1]),
)
@settings(deadline=None, max_examples=200)
def test_dasgupta_finite_nonnegative(p_N2, p_tot, T, dIW):
    """Dasgupta produces finite, non-negative values everywhere in
    the physically-relevant input space."""
    N2 = SolubilityN2('dasgupta')
    val = N2.dasgupta(p_N2, p_tot, T, dIW)
    assert math.isfinite(val), (
        f'Dasgupta NaN/Inf at p_N2={p_N2:.4e}, p_tot={p_tot:.4e}, T={T:.1f}, dIW={dIW:.2f}'
    )
    assert val >= 0.0


# ===========================================================================
# Gaillard S2: monotonicity + finiteness across the physical space
# ===========================================================================


@given(
    p_S2=st.floats(min_value=1e-19, max_value=P_RANGE[1]),
    T=st.floats(min_value=T_RANGE[0], max_value=T_RANGE[1]),
    dIW_low=st.floats(min_value=DIW_RANGE[0], max_value=-0.1),
    dIW_high=st.floats(min_value=0.1, max_value=DIW_RANGE[1]),
)
@settings(deadline=None, max_examples=200)
def test_gaillard_monotonic_with_oxidation(p_S2, T, dIW_low, dIW_high):
    """Gaillard at a more reducing fO2 is >= value at more oxidising."""
    S2 = SolubilityS2('gaillard')
    val_low = S2.gaillard(p_S2, T, dIW_low)
    val_high = S2.gaillard(p_S2, T, dIW_high)
    assert val_low >= val_high, (
        f'Gaillard at dIW={dIW_low:.2f} ({val_low:.4e}) less than '
        f'at dIW={dIW_high:.2f} ({val_high:.4e}) '
        f'(p_S2={p_S2:.4e}, T={T:.1f})'
    )

    # Discrimination guard: when both endpoints are above the p_S2 < 1e-20
    # floor and the dIW separation is >= 4 dex, the +0.5 ln(p_S2/fO2) term
    # gives a strict inequality, not just >=.
    if val_low > 0.0 and (dIW_high - dIW_low) >= 4.0:
        assert val_low > val_high, (
            f'Gaillard at dIW_low={dIW_low}, dIW_high={dIW_high} should '
            f'differ strictly when both endpoints are above the floor'
        )


@given(
    p_S2=st.floats(min_value=1e-19, max_value=P_RANGE[1]),
    T=st.floats(min_value=T_RANGE[0], max_value=T_RANGE[1]),
    dIW=st.floats(min_value=DIW_RANGE[0], max_value=DIW_RANGE[1]),
)
@settings(deadline=None, max_examples=200)
def test_gaillard_finite_nonnegative(p_S2, T, dIW):
    """Hypothesis fuzz: Gaillard S2 solubility is finite and non-negative
    everywhere in the (p_S2, T, dIW) physical input space."""
    S2 = SolubilityS2('gaillard')
    val = S2.gaillard(p_S2, T, dIW)
    assert math.isfinite(val)
    assert val >= 0.0


# ===========================================================================
# OxygenFugacity: continuity + agreement between O'Neill and Fischer
# in the calibration overlap
# ===========================================================================


@given(
    T=st.floats(min_value=T_RANGE[0], max_value=T_RANGE[1]),
    dIW=st.floats(min_value=DIW_RANGE[0], max_value=DIW_RANGE[1]),
)
@settings(deadline=None, max_examples=200)
def test_oneill_finite_over_bounds(T, dIW):
    """O'Neill & Eggins (2002) IW buffer + shift is finite and in a
    physically plausible range everywhere in the physical T-dIW space."""
    val = OxygenFugacity('oneill')(T, dIW)
    assert math.isfinite(val)

    # Discrimination guard: log10(fO2) for IW + dIW over the physical
    # T-dIW window must be in roughly [-30, +5]. A stub returning 1.0
    # for any input would pass the finite check; this band excludes that.
    assert -30.0 < val < 5.0, (
        f"O'Neill log10(fO2) at T={T:.1f}, dIW={dIW:.2f} = {val:.3f} "
        f'outside the physically plausible window [-30, +5]'
    )


@given(
    T=st.floats(min_value=1800.0, max_value=2200.0),
)
@settings(deadline=None, max_examples=100)
def test_oneill_fischer_disagree_by_less_than_dex(T):
    """Within the 1800-2200 K window (overlap of the two
    calibrations) the O'Neill and Fischer buffers agree to within
    1 log10 unit, consistent with the doc claim."""
    oneill_val = OxygenFugacity('oneill')(T, 0.0)
    fischer_val = OxygenFugacity('fischer')(T, 0.0)
    diff = abs(oneill_val - fischer_val)
    assert diff < 1.0

    # Discrimination guard: confirm each buffer returns a physically
    # plausible value independently. A stub that returned the same constant
    # for both buffers would give diff = 0 and pass the diff < 1 check
    # trivially. Both buffers at T in [1800, 2200] K with dIW=0 should
    # return log10(fO2) in roughly [-12, -5] (the upper edge widens as
    # T approaches 2200 K).
    assert -12.0 < oneill_val < -5.0, (
        f"O'Neill log10(fO2) at T={T:.1f} = {oneill_val:.3f} "
        f'outside the expected [-12, -5] window for IW'
    )
    assert -12.0 < fischer_val < -5.0, (
        f'Fischer log10(fO2) at T={T:.1f} = {fischer_val:.3f} '
        f'outside the expected [-12, -5] window for IW'
    )


# ===========================================================================
# ModifiedKeq: G_eq is finite and well-behaved over the input space
# ===========================================================================


@given(
    T=st.floats(min_value=T_RANGE[0], max_value=T_RANGE[1]),
    dIW=st.floats(min_value=DIW_RANGE[0], max_value=DIW_RANGE[1]),
)
@settings(deadline=None, max_examples=200)
def test_modified_keq_finite_positive(T, dIW):
    """Every modified equilibrium constant is finite and strictly
    positive across the physical input space (since K_eq = 10**(...))"""
    for method in (
        'janaf_H2',
        'janaf_CO',
        'schaefer_CH4',
        'janaf_SO2',
        'janaf_H2S',
        'janaf_NH3',
    ):
        Keq = ModifiedKeq(method)
        val = Keq(T, dIW)
        assert math.isfinite(val), f'{method}(T={T}, dIW={dIW}) = {val} is not finite'
        assert val > 0.0
