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
import warnings

import pytest

from calliope.solubility import SolubilityN2, SolubilityS2

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
        changed or the default x_FeO drifted off 10.0 wt%.

        Hidden coupling: this pin depends on OxygenFugacity('oneill')
        evaluated at T=2500 K, fO2_shift=0. Any audit of the IW buffer
        coefficients in oxygen_fugacity.py will require regenerating
        this number — the failure mode here surfaces in test_solubility,
        not test_oxygen_fugacity.
        """
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
        it must still evaluate finitely. isfinite catches NaN and Inf;
        the prior `out > 0.0` check was redundant with that since
        np.exp of any finite real is positive."""
        s = SolubilityS2(x_FeO=-5.0)
        out = s(1.0, 2500.0, 0.0)
        assert math.isfinite(out)

    def test_xFeO_zero_consistent_with_drop_term(self):
        """Discriminating: at x_FeO=0 the 0.124*x_FeO term drops, so
        the result must equal SolubilityS2(x_FeO=10) divided by
        exp(0.124*10) ~ 3.456. Catches an off-by-coefficient bug like
        +0.124*(x_FeO-10) or +0.0124*x_FeO."""
        baseline = SolubilityS2(x_FeO=10.0)
        zero = SolubilityS2(x_FeO=0.0)
        ratio = zero(1.0, 2500.0, 0.0) / baseline(1.0, 2500.0, 0.0)
        assert ratio == pytest.approx(math.exp(-1.24), rel=1e-12)


# ---------------------------------------------------------------------------
# SolubilityN2 :: x_SiO2 / x_Al2O3 / x_TiO2
# ---------------------------------------------------------------------------


class TestSolubilityN2_meltComposition:
    """The Dasgupta et al. (2022) law adds a molecular term whose
    prefactor c_melt = exp(4.67 + 7.11 x_SiO2 - 13.06 x_Al2O3 - 120.67 x_TiO2)
    is precomputed in __init__ and stored as `dasfac_2`. Tests pin the
    default values and verify the closed-form scaling with each kwarg
    independently, in line with the published expression."""

    DEFAULT_SIO2 = 0.56
    DEFAULT_AL2O3 = 0.11
    DEFAULT_TIO2 = 0.01
    DEFAULT_DASFAC = math.exp(4.67 + 7.11 * 0.56 - 13.06 * 0.11 - 120.67 * 0.01)  # ~ 406.79

    def test_default_kwargs_give_pre_existing_dasfac(self):
        """No-arg call must reproduce the previously hardcoded
        dasfac_2 = exp(4.67 + 7.11*0.56 - 13.06*0.11 - 120.67*0.01) ~ 406.79.
        Pin the closed-form value so any coefficient drift surfaces."""
        s = SolubilityN2('dasgupta')
        assert s.dasfac_2 == pytest.approx(self.DEFAULT_DASFAC, rel=1e-12)
        # Numeric pin against the actual computed value at kwarg
        # introduction time (catches off-by-1 in the constants too).
        assert s.dasfac_2 == pytest.approx(406.791, rel=1e-4)

    def test_default_kwargs_attribute_values(self):
        """Pin x_SiO2/Al2O3/TiO2 attributes to their documented
        defaults; a future maintainer changing these would have to
        update this test, signalling the doc-text needs the same
        update."""
        s = SolubilityN2('dasgupta')
        assert s.x_SiO2 == self.DEFAULT_SIO2
        assert s.x_Al2O3 == self.DEFAULT_AL2O3
        assert s.x_TiO2 == self.DEFAULT_TIO2

    @pytest.mark.parametrize(
        'kwarg,delta,coef',
        [
            ('x_SiO2', 0.10, 7.11),  # +0.10 SiO2 -> exp(0.711) ~ 2.04
            ('x_Al2O3', 0.05, -13.06),  # +0.05 Al2O3 -> exp(-0.653) ~ 0.520
            ('x_TiO2', 0.01, -120.67),  # +0.01 TiO2 -> exp(-1.2067) ~ 0.299
        ],
    )
    def test_each_kwarg_scales_dasfac_independently(self, kwarg, delta, coef):
        """Each composition kwarg enters c_melt with its own
        coefficient. Bumping one kwarg by `delta` must scale dasfac_2
        by exp(coef*delta), independent of the other two. Catches
        copy-paste swaps among the three coefficients (e.g. using
        7.11 for Al2O3 instead of SiO2)."""
        kwargs_default = {
            'x_SiO2': self.DEFAULT_SIO2,
            'x_Al2O3': self.DEFAULT_AL2O3,
            'x_TiO2': self.DEFAULT_TIO2,
        }
        kwargs_modified = dict(kwargs_default)
        kwargs_modified[kwarg] += delta

        s_default = SolubilityN2('dasgupta', **kwargs_default)
        s_modified = SolubilityN2('dasgupta', **kwargs_modified)

        ratio = s_modified.dasfac_2 / s_default.dasfac_2
        assert ratio == pytest.approx(math.exp(coef * delta), rel=1e-12)

    def test_dasgupta_call_uses_new_dasfac(self):
        """Discriminating: the Dasgupta call adds `pb_N2 * dasfac_2` on
        top of an exponential redox-dependent term. Pick conditions
        where the dasfac term carries appreciable weight (high p,
        oxidising) and verify changing dasfac changes ppmw in the
        predicted direction."""
        p_N2 = 1000.0  # bar
        p_total = 5000.0  # bar
        T = 1800.0
        dIW = +5.0  # strongly oxidising; suppresses the redox term

        s_default = SolubilityN2('dasgupta')
        s_high_SiO2 = SolubilityN2('dasgupta', x_SiO2=0.66)  # +0.10

        ppmw_default = s_default(p_N2, p_total, T, dIW)
        ppmw_high = s_high_SiO2(p_N2, p_total, T, dIW)

        # The dasfac term scales by exp(7.11 * 0.10) = ~ 2.036; the
        # redox term is unchanged. At dIW=+5 the redox term is heavily
        # suppressed, so dasfac dominates and the ratio approaches 2.04.
        ratio = ppmw_high / ppmw_default
        assert 1.5 < ratio < 2.1

    def test_libourel_unaffected_by_composition_kwargs(self):
        """The Libourel law has no melt-composition dependence; its
        output must be unchanged regardless of x_SiO2/Al2O3/TiO2.
        Discriminating: a refactor that accidentally threaded the new
        kwargs into power_law would be caught here."""
        p_N2 = 100.0
        baseline = SolubilityN2('libourel')(p_N2)
        custom = SolubilityN2('libourel', x_SiO2=0.99, x_Al2O3=0.99, x_TiO2=0.99)(p_N2)
        assert custom == pytest.approx(baseline, rel=1e-12)

    def test_zero_composition_yields_pure_4_67_prefactor(self):
        """Edge: setting all three to zero collapses dasfac_2 to
        exp(4.67) ~ 106.7. Catches a bug where the constant 4.67 was
        accidentally mixed into a kwarg coefficient."""
        s = SolubilityN2('dasgupta', x_SiO2=0.0, x_Al2O3=0.0, x_TiO2=0.0)
        assert s.dasfac_2 == pytest.approx(math.exp(4.67), rel=1e-12)

    def test_libourel_call_path_does_not_break_with_extreme_dasfac_inputs(self):
        """Edge: even if a user supplies pathological composition
        values that would overflow the dasgupta exponential (e.g.
        enormous SiO2 of 10.0, well outside any physical mole fraction),
        the libourel law must still evaluate to the Henry-law output.
        Discriminating: dasfac_2 is gated on composition='dasgupta',
        so libourel callers must NOT pay the exp() precompute cost AND
        must NOT see any RuntimeWarning emitted at construction time.
        """
        with warnings.catch_warnings():
            warnings.simplefilter('error')  # any warning becomes an error
            s = SolubilityN2('libourel', x_SiO2=10.0)
        # libourel path stores no precomputed dasfac_2
        assert s.dasfac_2 is None
        out = s(50.0)
        assert math.isfinite(out)
        assert out == pytest.approx(0.0611 * 50.0, rel=1e-12)


# ---------------------------------------------------------------------------
# Backward-compatibility regression tests (PROTEUS callers)
# ---------------------------------------------------------------------------


class TestBackwardCompatibility:
    """Pin numeric outputs of the no-arg SolubilityS2() and
    SolubilityN2() constructors used in `solve.dissolved_mass` so any
    silent default-value change surfaces here."""

    def test_solve_dissolved_mass_callsites_use_defaults(self):
        """Discriminating: `solve.dissolved_mass` instantiates
        `SolubilityS2()` and `SolubilityN2('dasgupta')` with no
        composition kwargs. Pin both default-instantiated objects to
        their pre-kwarg numeric outputs at a fixed (p, T, fO2_shift)
        triple. A drift here means a future change to the default
        kwarg values broke PROTEUS-side runs.

        Hidden coupling: the S2 pin below couples to
        OxygenFugacity('oneill'). An IW-buffer coefficient audit will
        require regenerating these numbers.
        """
        # SolubilityS2 default at (p_S2, T, dIW)
        s2 = SolubilityS2()
        ppmw_S2 = s2(1.0, 2500.0, 0.0)
        assert ppmw_S2 == pytest.approx(30095.04, rel=1e-5)

        # SolubilityN2('dasgupta') default at (p_N2, p_tot, T, dIW)
        n2 = SolubilityN2('dasgupta')
        ppmw_N2 = n2(50.0, 200.0, 2000.0, 0.0)
        # Regression value frozen from the pre-kwarg implementation.
        assert ppmw_N2 == pytest.approx(2.14133, rel=1e-5)
