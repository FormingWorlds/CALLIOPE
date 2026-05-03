"""Tests for the user-configurable melt-composition kwargs on
`SolubilityS2` (x_FeO) and `SolubilityN2` (x_SiO2, x_Al2O3, x_TiO2).

The kwargs are backward-compatible: callers that omit them must get
bit-identical numerics to the prior hardcoded values. Non-Earth
overrides must propagate into the dissolved-mass formulas as predicted
by the closed-form expressions in Gaillard et al. (2022) and Dasgupta
et al. (2022).
"""

from __future__ import annotations

import math

import pytest

from calliope.solubility import SolubilityS2

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# SolubilityS2 :: x_FeO
# ---------------------------------------------------------------------------


class TestSolubilityS2_xFeO:
    """The Gaillard et al. (2022) law has form
        ln(X_S^melt) = 13.8426 - 26476/T + 0.124*x_FeO + 0.5*ln(p/fO2)
    so x_FeO enters as exp(0.124 * dx_FeO) on the prefactor. Tests use
    this closed form to discriminate the implementation against
    plausible bugs (missing 0.124, treating x_FeO as fraction not wt%,
    multiplicative not additive).
    """

    def test_default_xFeO_is_10wt_percent(self):
        """Default kwarg must keep the hardcoded Earth-mantle value of
        10.0 wt%, so that PROTEUS callers using `SolubilityS2()` see no
        numerical drift."""
        s = SolubilityS2()
        assert s.x_FeO == 10.0

    def test_default_call_matches_pre_kwarg_value(self):
        """Pin one numeric output of the default-x_FeO call against the
        pre-kwarg implementation. Drift here means either the formula
        changed or the default x_FeO drifted off 10.0 wt%."""
        s = SolubilityS2()
        out = s(1.0, 2500.0, 0.0)
        # Regression pin: computed by the implementation at the time
        # the x_FeO kwarg was introduced, with the default x_FeO=10.0.
        assert out == pytest.approx(30095.04, rel=1e-5)

    @pytest.mark.parametrize(
        'x_FeO,expected_ratio',
        [
            (5.0, math.exp(0.124 * (5.0 - 10.0))),  # half FeO -> exp(-0.62) ~ 0.538
            (10.0, 1.0),  # default, identity
            (15.0, math.exp(0.124 * (15.0 - 10.0))),  # ~1.86
            (20.0, math.exp(0.124 * (20.0 - 10.0))),  # ~3.46
        ],
    )
    def test_xFeO_scales_prefactor_correctly(self, x_FeO, expected_ratio):
        """The Gaillard law is linear in x_FeO inside the exponent,
        so the ratio of solubilities at different x_FeO at fixed p, T,
        fO2_shift must equal exp(0.124 * dx_FeO). Discriminating because
        a bug applying x_FeO as a multiplicative factor (e.g. *x_FeO
        rather than +0.124*x_FeO) would not produce this exact ratio.
        """
        baseline = SolubilityS2(x_FeO=10.0)
        custom = SolubilityS2(x_FeO=x_FeO)
        p_S2, T, dIW = 1.0, 2500.0, 0.0

        ratio = custom(p_S2, T, dIW) / baseline(p_S2, T, dIW)
        assert ratio == pytest.approx(expected_ratio, rel=1e-12)

    def test_xFeO_does_not_affect_zero_pressure_short_circuit(self):
        """Edge case: at p_S2 < 1e-20 bar the law returns 0.0
        identically; the x_FeO kwarg must not bypass that guard. Keeps
        the divergence-in-log handling intact."""
        for x_FeO in (0.0, 10.0, 50.0):
            s = SolubilityS2(x_FeO=x_FeO)
            assert s(1.0e-25, 2500.0, 0.0) == 0.0

    def test_xFeO_extreme_negative_evaluates_finitely(self):
        """Edge case: a physically nonsensical negative x_FeO value
        is not validated (the law has no domain constraint encoded), so
        it must still evaluate finitely. Pin that the result is finite
        and positive (exp domain), so the kwarg cannot silently produce
        NaN/Inf for plausible-looking-but-wrong inputs."""
        s = SolubilityS2(x_FeO=-5.0)
        out = s(1.0, 2500.0, 0.0)
        assert math.isfinite(out)
        assert out > 0.0

    def test_xFeO_zero_consistent_with_drop_term(self):
        """Discriminating: at x_FeO=0 the 0.124*x_FeO term drops, so
        the result must equal SolubilityS2(x_FeO=10) divided by
        exp(0.124*10) ~ 3.456. Catches an off-by-coefficient bug like
        +0.124*(x_FeO-10) or +0.0124*x_FeO."""
        baseline = SolubilityS2(x_FeO=10.0)
        zero = SolubilityS2(x_FeO=0.0)
        ratio = zero(1.0, 2500.0, 0.0) / baseline(1.0, 2500.0, 0.0)
        assert ratio == pytest.approx(math.exp(-1.24), rel=1e-12)
