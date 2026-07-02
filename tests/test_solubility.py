"""Tests for `src/calliope/solubility.py`.

Exercises the Henry's-law solubility models for H2O (Sossi 2023 default,
Dixon 1995, Hamilton 1964 / Wilson and Head 1981, Newcombe 2017), S2
(Gaillard 2022), and N2 (Dasgupta 2022, Libourel).

- Reference pins: H2O peridotite default against the Sossi et al. (2023)
  `524 * p^0.5` constant; S2 default against the Gaillard et al. (2022)
  Earth-mantle numerics; N2 dasgupta prefactor against the published
  composition coefficients.
- Conservation: zero partial pressure returns identically zero; Henry's
  identity in the linear regime.
- Monotonicity: dissolved ppmw increases with partial pressure for every
  H2O parameterization.
- Closed-form scaling: composition kwargs (`x_FeO`, `x_SiO2`, `x_Al2O3`,
  `x_TiO2`) enter the exponential prefactors as `exp(coef * delta)`;
  the ratio of solubilities at two compositions matches the closed-form
  factor to floating-point precision.
- Edge cases: zero-pressure short-circuit on S2, negative composition
  values evaluate finitely, libourel path unaffected by composition kwargs.
"""

from __future__ import annotations

import math
import warnings

import pytest

from calliope.constants import molar_mass, noble_gases
from calliope.solubility import (
    JAMBON86_STP_HENRY,
    STP_MOLAR_VOLUME_CM3,
    SolubilityH2O,
    SolubilityN2,
    SolubilityNobleGas,
    SolubilityS2,
    jambon86_ppmw_per_bar,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


# ---------------------------------------------------------------------------
# SolubilityH2O :: H2O parameterizations
# ---------------------------------------------------------------------------


class TestSolubilityH2O:
    """Power-law H2O solubilities. Each parameterization carries a published
    `const` and `exponent`; the test pins both via the closed-form
    `const * p^exponent` identity at chosen pressures."""

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    def test_peridotite_default_matches_sossi_2023_fit(self):
        """Sossi et al. (2023) peridotite H2O fit: `ppmw = 524 * p^0.5`.

        `SolubilityH2O('peridotite')` is the package-wide default; the
        constant 524 ppmw/bar^0.5 and the 0.5 exponent come from Sossi
        et al. (2023).

        Discrimination guards: a regression that swapped to the Dixon
        1995 basalt constant (965 ppmw/bar^0.5) at the same exponent
        would change the result by factor 1.84; an exponent flip
        (1.0 vs 0.5) would change it by `sqrt(100) = 10x` at p = 100 bar.
        """
        s = SolubilityH2O()  # peridotite default
        # Zero-pressure boundary: identically zero by power-law definition.
        assert s(0.0) == 0.0
        # At p = 100 bar: 524 * sqrt(100) = 5240 ppmw.
        val = s(100.0)
        expected = 524.0 * 10.0
        assert val == pytest.approx(expected, rel=1e-12)
        # Wrong-law guard: basalt_dixon would give 9650, off by ~84%.
        assert abs(val - 965.0 * 10.0) > 1000.0
        # Wrong-exponent guard: p^1.0 instead of p^0.5 would give 52400.
        assert abs(val - 524.0 * 100.0) > 1000.0
        # Sign + scale guards: positive ppmw, order 1e3 to 1e4 at 100 bar.
        assert val > 0
        assert 1e3 < val < 1e4

    @pytest.mark.physics_invariant
    def test_basalt_dixon_matches_published_constant(self):
        """Dixon et al. (1995) basalt H2O: `ppmw = 965 * p^0.5`."""
        s = SolubilityH2O('basalt_dixon')
        val = s(100.0)
        expected = 965.0 * 10.0
        assert val == pytest.approx(expected, rel=1e-12)
        # Wrong-law guard against peridotite default (524) at same p.
        assert abs(val - 524.0 * 10.0) > 1000.0

    @pytest.mark.physics_invariant
    def test_basalt_wilson_uses_non_half_exponent(self):
        """Hamilton (1964) / Wilson and Head (1981) basalt: `ppmw = 215 * p^0.7`.

        Non-square-root exponent: discriminating, because a regression
        that defaults all H2O laws to the 0.5 exponent would land at
        2150 ppmw at p = 100 bar instead of the correct ~5403 ppmw.
        """
        s = SolubilityH2O('basalt_wilson')
        val = s(100.0)
        expected = 215.0 * (100.0**0.7)  # ~5403 ppmw
        assert val == pytest.approx(expected, rel=1e-12)
        # Wrong-exponent guard: p^0.5 would give 2150.
        assert abs(val - 215.0 * 10.0) > 1000.0

    @pytest.mark.physics_invariant
    def test_anorthite_diopside_and_lunar_glass_match_newcombe_2017(self):
        """Newcombe et al. (2017): anorthite-diopside 727 ppmw/bar^0.5,
        lunar glass 683 ppmw/bar^0.5. Both follow the sqrt law."""
        s_ad = SolubilityH2O('anorthite_diopside')
        s_lg = SolubilityH2O('lunar_glass')
        val_ad = s_ad(100.0)
        val_lg = s_lg(100.0)
        assert val_ad == pytest.approx(727.0 * 10.0, rel=1e-12)
        assert val_lg == pytest.approx(683.0 * 10.0, rel=1e-12)
        # The two are close (~6% apart), so a regression that swapped
        # them would not be caught by a single-law check; pin both.
        assert val_ad != pytest.approx(val_lg, rel=1e-3)

    @pytest.mark.physics_invariant
    def test_h2o_monotonic_in_pressure_for_every_parameterization(self):
        """All five parameterizations are strictly increasing in p.

        Henry's law sign convention: higher partial pressure -> more
        dissolved ppmw. A regression that flipped the sign of the
        exponent (very unlikely but possible if power_law gained a
        negative-exponent default) would invert this ordering.
        """
        for name in (
            'peridotite',
            'basalt_dixon',
            'basalt_wilson',
            'anorthite_diopside',
            'lunar_glass',
        ):
            s = SolubilityH2O(name)
            low = s(10.0)
            mid = s(100.0)
            high = s(1000.0)
            assert low < mid < high
            # Discrimination: high - low spans ~2 orders for sqrt laws,
            # ~4 orders for p^0.7. Pin a positive delta.
            assert (high - low) > 0


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
        assert s.x_FeO == pytest.approx(10.0)

        # Discrimination guard: confirm the default x_FeO actually flows
        # through to the call path. A constructor that stored 10 on the
        # attribute but used a different value internally would pass the
        # bare attribute pin.
        s_explicit = SolubilityS2(x_FeO=10.0)
        assert s(1.0, 2500.0, 0.0) == pytest.approx(s_explicit(1.0, 2500.0, 0.0), rel=1e-12)

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    def test_default_call_matches_gaillard_2022_earth_mantle_value(self):
        """Pin S2 ppmw under Gaillard et al. (2022) Earth-mantle defaults.

        At p_S2 = 1 bar, T = 2500 K, fO2_shift = 0, x_FeO = 10 wt%, and
        the default Fischer 2011 IW buffer:

        ln(X) = 13.8426 - 26476/2500 + 0.124*10 + 0.5*ln(1/fO2_bar)
              = 13.8426 - 10.5904 + 1.24 + 0.5*ln(10^-IW(2500))
              ~ 9.479
        ppmw = exp(9.479) ~ 13086

        Hidden coupling: the pin depends on the default IW buffer
        (Fischer 2011 since 2026-05). Any change to the IW buffer
        coefficients in oxygen_fugacity.py will require regenerating
        this number; the failure surfaces here, not in
        test_oxygen_fugacity.
        """
        s = SolubilityS2()
        out = s(1.0, 2500.0, 0.0)
        # Regression pin: computed at default x_FeO=10.0 with Fischer 2011.
        assert out == pytest.approx(13085.87, rel=1e-5)
        # Sign and scale guards: positive ppmw, order 1e4 at these conditions.
        assert out > 0
        assert 1e3 < out < 1e5

    @pytest.mark.physics_invariant
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

        # Discrimination guard: the wrong formula (multiplicative *x_FeO/10
        # instead of additive +0.124*x_FeO) would give ratio = x_FeO/10.
        # At x_FeO=5 the wrong ratio is 0.5 (correct exp(-0.62)~0.538);
        # at x_FeO=15 the wrong ratio is 1.5 (correct exp(0.62)~1.86);
        # at x_FeO=20 the wrong ratio is 2.0 (correct exp(1.24)~3.46).
        # The identity case x_FeO=10 is excluded since both formulas coincide.
        if x_FeO != 10.0:
            wrong_ratio_multiplicative = x_FeO / 10.0
            assert abs(ratio - wrong_ratio_multiplicative) > 0.03

    def test_xFeO_does_not_affect_zero_pressure_short_circuit(self):
        """Edge case: at p_S2 < 1e-20 bar the law returns 0.0
        identically; the x_FeO kwarg must not bypass that guard. Keeps
        the divergence-in-log handling intact."""
        for x_FeO in (0.0, 10.0, 50.0):
            s = SolubilityS2(x_FeO=x_FeO)
            assert s(1.0e-25, 2500.0, 0.0) == pytest.approx(0.0, abs=1e-30)

        # Discrimination guard: the floor branch at low p_S2 must be distinct
        # from the main path. At p_S2 = 1.0 bar each of the three x_FeO values
        # gives a different, nonzero result; a stub that hard-coded 0.0 would
        # fail here.
        for x_FeO in (0.0, 10.0, 50.0):
            s = SolubilityS2(x_FeO=x_FeO)
            assert s(1.0, 2500.0, 0.0) > 0.0

    def test_xFeO_extreme_negative_evaluates_finitely(self):
        """Edge case: a physically nonsensical negative x_FeO value
        is not validated (the law has no domain constraint encoded), so
        it must still evaluate finitely and positively. The output must
        also follow the published linear-in-x_FeO law: negative x_FeO
        gives less solubility than x_FeO = 0."""
        s = SolubilityS2(x_FeO=-5.0)
        out = s(1.0, 2500.0, 0.0)
        assert math.isfinite(out)

        # Discrimination guard: the linear-in-x_FeO law predicts that
        # x_FeO = -5 gives exp(-0.62) ~ 0.538 times the x_FeO = 0 value.
        # A clamped-at-zero implementation (treating negative x_FeO as 0)
        # would give the same result as x_FeO = 0; a sign-flip bug would
        # give exp(+0.62) ~ 1.86 times.
        zero = SolubilityS2(x_FeO=0.0)(1.0, 2500.0, 0.0)
        assert out == pytest.approx(zero * math.exp(0.124 * -5.0), rel=1e-12)

    @pytest.mark.physics_invariant
    def test_xFeO_zero_consistent_with_drop_term(self):
        """Discriminating: at x_FeO=0 the 0.124*x_FeO term drops, so
        the result must equal SolubilityS2(x_FeO=10) divided by
        exp(0.124*10) ~ 3.456. Catches an off-by-coefficient bug like
        +0.124*(x_FeO-10) or +0.0124*x_FeO."""
        baseline = SolubilityS2(x_FeO=10.0)
        zero = SolubilityS2(x_FeO=0.0)
        ratio = zero(1.0, 2500.0, 0.0) / baseline(1.0, 2500.0, 0.0)
        assert ratio == pytest.approx(math.exp(-1.24), rel=1e-12)

        # Discrimination guard: an off-by-10 in the coefficient
        # (+0.0124 * x_FeO instead of +0.124 * x_FeO) would give a ratio
        # of exp(-0.124) ~ 0.883. The correct ratio is exp(-1.24) ~ 0.289;
        # the gap is 0.59, well outside any tolerance.
        wrong_coef_ratio = math.exp(-0.124)
        assert abs(ratio - wrong_coef_ratio) > 0.1


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

    @pytest.mark.physics_invariant
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

        # Discrimination guard: each kwarg's coefficient is distinct, so a
        # copy-paste bug that used a sibling coefficient would yield a
        # different ratio. The three coefficients are +7.11 (SiO2), -13.06
        # (Al2O3), -120.67 (TiO2); confirm the measured ratio does not
        # match either of the other two predictions.
        sibling_coefs = [c for c in (7.11, -13.06, -120.67) if c != coef]
        for wrong_coef in sibling_coefs:
            wrong_ratio = math.exp(wrong_coef * delta)
            assert abs(ratio - wrong_ratio) > 1e-3 * abs(ratio), (
                f'Ratio {ratio:.4e} matches sibling coefficient {wrong_coef} (expected {coef})'
            )

    @pytest.mark.physics_invariant
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

        # Discrimination guard: an implementation that did NOT thread the
        # new x_SiO2 kwarg into dasfac_2 would leave ratio ~ 1.0; a wrong
        # sign on the coefficient would give exp(-7.11 * 0.10) ~ 0.49.
        # Both failure modes are excluded by the [1.5, 2.1] band, but
        # tighten by confirming the gap from 1.0 is meaningful.
        assert abs(ratio - 1.0) > 0.3, (
            f'Ratio {ratio:.4f} is too close to 1.0; the x_SiO2 kwarg '
            f'may not be threading through to the call path'
        )

    def test_libourel_unaffected_by_composition_kwargs(self):
        """The Libourel law has no melt-composition dependence; its
        output must be unchanged regardless of x_SiO2/Al2O3/TiO2.
        Discriminating: a refactor that accidentally threaded the new
        kwargs into power_law would be caught here."""
        p_N2 = 100.0
        baseline = SolubilityN2('libourel')(p_N2)
        custom = SolubilityN2('libourel', x_SiO2=0.99, x_Al2O3=0.99, x_TiO2=0.99)(p_N2)
        assert custom == pytest.approx(baseline, rel=1e-12)

        # Discrimination guard: a stub that returned 0 for both calls would
        # pass the equality check trivially. Confirm the libourel output is
        # physically meaningful (positive, finite, in a reasonable ppmw range
        # for p_N2 = 100 bar).
        assert baseline > 0.0
        assert math.isfinite(baseline)

    def test_zero_composition_yields_pure_4_67_prefactor(self):
        """Edge: setting all three to zero collapses dasfac_2 to
        exp(4.67) ~ 106.7. Catches a bug where the constant 4.67 was
        accidentally mixed into a kwarg coefficient."""
        s = SolubilityN2('dasgupta', x_SiO2=0.0, x_Al2O3=0.0, x_TiO2=0.0)
        assert s.dasfac_2 == pytest.approx(math.exp(4.67), rel=1e-12)

        # Discrimination guard: nearby constants (4.6, 5.0) would give
        # prefactors of exp(4.6) ~ 99.5 or exp(5.0) ~ 148.4, distinguishable
        # from exp(4.67) ~ 106.7 at the 7-40% level. Confirm the value is
        # in the narrow band around exp(4.67).
        assert s.dasfac_2 < math.exp(4.75)
        assert s.dasfac_2 > math.exp(4.60)

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
        their numeric outputs at a fixed (p, T, fO2_shift) triple. A
        drift here means a future change to the default kwarg values
        broke PROTEUS-side runs.

        Hidden coupling: the S2 pin below couples to the default IW
        buffer (Fischer et al. 2011). A change to the IW-buffer
        coefficients will require regenerating these numbers.
        """
        # SolubilityS2 default at (p_S2, T, dIW)
        s2 = SolubilityS2()
        ppmw_S2 = s2(1.0, 2500.0, 0.0)
        assert ppmw_S2 == pytest.approx(13085.87, rel=1e-5)

        # SolubilityN2('dasgupta') default at (p_N2, p_tot, T, dIW)
        n2 = SolubilityN2('dasgupta')
        ppmw_N2 = n2(50.0, 200.0, 2000.0, 0.0)
        # Regression value: N2 solubility under dasgupta has no fO2
        # dependence at this dIW=0 evaluation, so this is buffer-agnostic.
        assert ppmw_N2 == pytest.approx(2.14133, rel=1e-5)


# ---------------------------------------------------------------------------
# SolubilityNobleGas :: Jambon, Weill & Braun (1986) Henry's law
# ---------------------------------------------------------------------------


class TestSolubilityNobleGas:
    """Noble gases dissolve by a single linear Henry's law,
    `ppmw = const * p`, with `const` derived from the Jambon et al. (1986)
    STP-volume solubility constant, the molar mass, and the STP molar
    volume. Tests pin the derived constant, the linearity, the physical
    ordering set by the tabulated data, and the melt-independence, and
    cross-check the constant against atmodeller."""

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    def test_every_gas_constant_pinned_to_independent_literal(self):
        """Pin each noble gas Henry constant [ppmw/bar] to an independently
        hand-computed literal, so a fat-fingered entry in `JAMBON86_STP_HENRY`
        or `molar_mass` for ANY of the five gases fails here even on a CI
        image without atmodeller installed.

        Each literal is `(k_STP / 2.24e4) * M[g/mol] * 1e6` evaluated by hand
        from the Jambon et al. (1986) constants and the molar masses:
        He 0.100065, Ne 0.225220, Ar 0.105220, Kr 0.112229, Xe 0.099642.
        """
        expected = {
            'He': 0.100065,
            'Ne': 0.225220,
            'Ar': 0.105220,
            'Kr': 0.112229,
            'Xe': 0.099642,
        }
        for gas in noble_gases:
            assert jambon86_ppmw_per_bar(gas) == pytest.approx(expected[gas], rel=1e-4)

        # Discrimination guard: the constants are close for several pairs
        # (He/Ar/Xe all near 0.1), so pin that each gas differs from every
        # other by more than the tolerance. A doubled Ar constant (0.105 ->
        # 0.210) would move it next to Ne and fail its own pin; this loop
        # catches a value copied from the wrong gas.
        for a in noble_gases:
            for b in noble_gases:
                if a == b:
                    continue
                if abs(expected[a] - expected[b]) > 1e-3:
                    assert jambon86_ppmw_per_bar(a) != pytest.approx(expected[b], rel=1e-4), (
                        f'{a} constant collides with {b}'
                    )

        # Unit-slip guard: forgetting the kg->g factor would make every
        # constant 1000x too small (order 1e-4). All five are order 1e-1.
        for gas in noble_gases:
            assert 1e-2 < jambon86_ppmw_per_bar(gas) < 1e0

    def test_conversion_matches_closed_form_for_every_gas(self):
        """`jambon86_ppmw_per_bar` reproduces the closed-form chain
        `(k_STP / V_STP) * M[g/mol] * 1e6`, confirming the helper wires the
        module constants together as documented. This is a construction
        check; the independent numeric anchor is the literal-pin test above.
        """
        for gas in noble_gases:
            expected = (
                (JAMBON86_STP_HENRY[gas] / STP_MOLAR_VOLUME_CM3)
                * (molar_mass[gas] * 1.0e3)
                * 1.0e6
            )
            assert jambon86_ppmw_per_bar(gas) == pytest.approx(expected, rel=1e-12)

        # Discrimination guard: a unit slip that forgot the kg->g factor
        # (molar_mass in kg/mol, not g/mol) would make every constant 1000x
        # too small. Pin the He magnitude against that failure mode.
        assert jambon86_ppmw_per_bar('He') > 1.0e-3

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    def test_helium_constant_matches_jambon_1986(self):
        """Pin the He Henry constant against the Jambon et al. (1986)
        tholeiitic-basalt abstract value of 56e-5 cm3 STP/g/bar, which is the
        value atmodeller also uses.

        const = (56e-5 / 2.24e4) * 4.0026 * 1e6 ~ 0.100065 ppmw/bar.

        Discrimination guards below rule out the plausible bug classes:
        wrong molar mass, wrong STP volume, wrong power, sign, and scale.
        """
        const = jambon86_ppmw_per_bar('He')
        assert const == pytest.approx(0.100065, rel=1e-4)
        # Wrong-molar-mass guard: using Ne's molar mass (20.1797) instead of
        # He's (4.0026) would give ~0.5045, a factor of ~5 larger.
        wrong_mass = (56e-5 / STP_MOLAR_VOLUME_CM3) * 20.1797 * 1.0e6
        assert abs(const - wrong_mass) > 0.1
        # Sign guard: solubility constants are strictly positive.
        assert const > 0
        # Scale guard: order 1e-1 ppmw/bar, not 1e2 (missing STP-volume
        # division) or 1e-4 (kg/g slip).
        assert 1e-2 < const < 1e0

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    def test_noble_constants_cross_check_against_atmodeller(self):
        """Cross-implementation check: CALLIOPE and atmodeller derive the
        Jambon et al. (1986) Henry constants from the same primitives for all
        five gases, so the two backends must agree to floating point. This is
        what makes the coupled cross-backend comparison a true parity test
        rather than two independent calibrations.
        """
        pytest.importorskip('atmodeller')
        from atmodeller.solubility import get_solubility_models

        models = get_solubility_models()
        for gas in noble_gases:
            cal = jambon86_ppmw_per_bar(gas)
            model = models[f'{gas}_basalt_jambon86']
            # atmodeller's power law is linear in fugacity, so the ppmw at
            # 1 bar is the Henry constant itself.
            try:
                atm = float(model.concentration(1.0, temperature=1500.0))
            except TypeError:
                atm = float(model.concentration(1.0))
            assert cal == pytest.approx(atm, rel=1e-10), f'{gas}: {cal} vs {atm}'

        # Discrimination guard: the five constants are distinct, so a
        # transcription that reused He's value for every gas would fail the
        # per-gas comparison. Confirm they span a real range.
        vals = [jambon86_ppmw_per_bar(g) for g in noble_gases]
        assert max(vals) / min(vals) > 2.0

    @pytest.mark.physics_invariant
    def test_henry_law_is_linear_in_pressure(self):
        """Noble gas solubility is linear (exponent 1): doubling the partial
        pressure doubles the dissolved concentration exactly, whereas a
        square-root law such as the H2O or S2 parameterization would scale by
        `sqrt(2)`. The linear form is pinned here so a stray exponent is
        caught.
        """
        s = SolubilityNobleGas('Ar')
        low = s(10.0)
        high = s(20.0)
        # Linear: exactly a factor of 2.
        assert high == pytest.approx(2.0 * low, rel=1e-12)
        # Wrong-exponent guard: a p^0.5 law would give a factor of sqrt(2).
        assert abs(high / low - 2.0**0.5) > 0.2
        # Zero-pressure boundary: identically zero dissolved gas.
        assert s(0.0) == 0.0

    @pytest.mark.physics_invariant
    def test_dissolved_concentration_monotonic_and_positive(self):
        """Every noble gas is strictly increasing and positive in pressure,
        the Henry's-law sign convention. Discriminating because a sign flip
        on the constant would invert the ordering."""
        for gas in noble_gases:
            s = SolubilityNobleGas(gas)
            a, b, c = s(1.0), s(100.0), s(1.0e4)
            assert 0.0 < a < b < c
            # Linearity pins the spacing: c/a == 1e4 exactly.
            assert c == pytest.approx(1.0e4 * a, rel=1e-12)

    def test_melt_composition_independence(self):
        """The Jambon law carries no melt-composition, temperature, or
        redox dependence: the class takes only the gas name and the call
        takes only pressure. Verify the mapping from gas to constant is the
        published ordering (Ne most soluble by weight, Xe least of the pair
        with He), which a swapped-constant bug would violate."""
        consts = {g: SolubilityNobleGas(g)(1.0) for g in noble_gases}
        # Ne has the largest ppmw/bar of the five (25e-5 STP constant times
        # its molar mass dominates); Xe the smallest.
        assert consts['Ne'] == max(consts.values())
        assert consts['Xe'] == min(consts.values())
        # Discrimination: He and Ar are close (0.1001 vs 0.1052) but not
        # equal; a bug that collapsed all gases to one constant would fail.
        assert consts['He'] != pytest.approx(consts['Ar'], rel=1e-3)

    def test_rejects_non_noble_gas(self):
        """The constructor validates its argument: a non-noble species
        (e.g. a CHNOS gas passed by mistake) raises ValueError rather than
        silently building a broken model."""
        with pytest.raises(ValueError, match='not a noble gas'):
            SolubilityNobleGas('H2O')
        with pytest.raises(ValueError, match='not a noble gas'):
            SolubilityNobleGas('CO2')
        # A valid gas must still build and evaluate.
        assert SolubilityNobleGas('He')(1.0) > 0.0
