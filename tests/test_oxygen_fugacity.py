"""Tests for `src/calliope/oxygen_fugacity.py`.

Exercises the `OxygenFugacity` IW-buffer dispatcher and its three
underlying fits:

- Cross-calibration and regression check: coded Fischer IW value at
  T = 2000 K against the closed-form `6.94059 - 28.1808e3 / T`, with a
  discrimination guard against the O'Neill & Eggins (2002) IW at the same
  T. Fischer et al. (2011) calculate the 1-bar Fe-FeO curve from Chase
  (1998) (Fig. 6) and add the pressure integral of their Eq. 2; their
  supplement (Table S2) gives a P-T polynomial. The source of this 2-term
  form is not stated in the paper body, and no independent numerical pin
  of the Fischer 1-bar fit exists in the repository test data. The guard
  catches a regression that silently dispatches to the wrong buffer; a
  default change in `oxygen_fugacity.py` would otherwise silently shift
  every PROTEUS-side number pinned to the previous default.
- Reference pin: Hirschmann (2021) IW at 1 bar and T = 1473 K against
  the IW value plotted in their Fig. 5, with guards against both the
  Fischer and O'Neill buffers.
- Analytical limit: the Hirschmann (2021) low-pressure slope matches
  `2 dV / (R T ln 10)` from the Fe and FeO molar volumes.
- Continuity: the fcc/bcc and hcp branches of Hirschmann (2021) meet at
  the fcc-hcp iron boundary within the published fit mismatch.
- Monotonicity: `log10(fO2)` is monotonic in T along each buffer
  over the 1500-3000 K range, and in P along the Hirschmann buffer.
- Symmetry: `fO2_shift` is strictly additive: `of(T, dIW) = of(T, 0) + dIW`.
- Boundedness: `log10(fO2)` is finite for any valid (T > 0, dIW finite).
- Error contract: T <= 0 raises `ValueError` mentioning the divergence
  in the underlying formulae; negative pressure raises `ValueError`;
  unknown buffer names, `None` and non-model attribute names raise
  `ValueError` at construction.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from calliope.oxygen_fugacity import OxygenFugacity

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_oxygen_fugacity_fischer_value_at_2000K_cross_check_and_regression(validation_plot):
    """Coded Fischer IW at T = 2000 K, cross-checked against O'Neill & Eggins (2002).

    Fischer et al. (2011) calculate the 1-bar Fe-FeO curve from Chase (1998)
    (Fig. 6) and add the pressure integral of their Eq. 2; their supplement
    (Table S2) gives a P-T polynomial. The source of this 2-term form is not
    stated in the paper body, and no independent numerical pin of the Fischer
    1-bar fit exists in the repository test data. At T = 2000 K the coded
    `'fischer'` (-7.14981) and `'oneill'` (-7.40782) expressions differ by
    0.258 dex; pinning both the 0.258 dex offset and the closed-form value
    provides a regression check and a wrong-buffer discrimination guard. Both
    values cluster with other published IW fits (Frost 1991 ~ -7.04; O'Neill
    1988 ~ -7.55) within the ~0.5 dex scatter of the literature.
    """
    T = 2000.0
    fischer = OxygenFugacity('fischer')(T)
    oneill = OxygenFugacity('oneill')(T)
    # Regression check on the offset between the two coded buffer expressions.
    assert (fischer - oneill) == pytest.approx(0.258, abs=0.02)
    # Direct regression check on the coded Fischer 1-bar expression at 2000 K.
    coded = 6.94059 - 28.1808e3 / T
    assert fischer == pytest.approx(coded, rel=1e-4)
    # Wrong-buffer guard: a silent dispatch to 'oneill' lands 0.26 dex
    # away from the Fischer value.
    assert abs(fischer - (-7.4078)) > 0.2
    # Sign guard: log10(fO2) at the IW buffer is always negative under
    # standard conditions (T < ~10000 K).
    assert fischer < 0
    # Scale guard: order of magnitude is -7, not -70 (forgotten log10)
    # or -0.7 (factor-10 unit slip on the temperature coefficient).
    assert -10 < fischer < -3

    def figure():
        import matplotlib.pyplot as plt

        grid = np.linspace(1900.0, 2100.0, 21)
        fig, ax = plt.subplots()
        for model, ls in (('fischer', '-'), ('oneill', '--'), ('hirschmann', ':')):
            ax.plot(grid, [OxygenFugacity(model)(t) for t in grid], ls=ls, label=model)
        ax.plot([T], [coded], 's', ms=10, mfc='none', label='coded form')
        ax.plot([T], [fischer], 'o', label='code, fischer')
        ax.plot([T], [oneill], 'x', label='code, oneill')
        ax.set(xlabel='T [K]', ylabel=r'$\log_{10} f\mathrm{O_2}$ [bar] at IW, 1 bar')
        ax.legend()
        fig.tight_layout()
        return fig

    validation_plot(figure, 'oxygen_fugacity_fischer_2000k')


@pytest.mark.physics_invariant
def test_oxygen_fugacity_oneill_value_at_2000K_matches_published_fit():
    """O'Neill & Eggins (2002) IW at T = 2000 K matches the closed form.

    Builds log10 fO2 from the free energy of Fe + 0.5 O2 = FeO(liq) of
    O'Neill & Eggins (2002, Chem. Geol. 186, 151), Table 6 on p. 163:
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
    for buffer_name in ('fischer', 'oneill', 'hirschmann'):
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
    """All buffers produce `log10(fO2)` strictly increasing with T over 1500-3000 K.

    A hot magma ocean is more reducing on an absolute scale than a cool
    one at the IW buffer, but the IW buffer itself is defined by the
    Fe-FeO equilibrium and its log10(fO2) becomes less negative (closer
    to zero) as T increases. The chosen window spans realistic surface
    temperatures from solidification (~1500 K) to early Earth magma
    ocean (~3000 K), giving a delta large enough to resolve a regression
    that flipped the sign of the slope.
    """
    for buffer_name in ('fischer', 'oneill', 'hirschmann'):
        of = OxygenFugacity(buffer_name)
        low = of(1500.0)
        mid = of(2250.0)
        high = of(3000.0)
        # Strict ordering: monotonic, not just non-decreasing.
        assert low < mid < high
        # Discrimination: 1500 K -> 3000 K delta exceeds 1 dex for both fits.
        assert (high - low) > 1.0


@pytest.mark.physics_invariant
def test_oxygen_fugacity_finite_over_realistic_temperature_range():
    """`log10(fO2)` is finite for any T in the realistic 800-5000 K range.

    Boundedness check: no nan, no inf, no complex intermediate. Catches a
    regression that introduces a `log(negative)` or `0/0` along an
    unexpected code path.
    """
    for buffer_name in ('fischer', 'oneill', 'hirschmann'):
        of = OxygenFugacity(buffer_name)
        for T in (800.0, 1500.0, 2500.0, 3500.0, 5000.0):
            for dIW in (-3.0, 0.0, 1.5):
                val = of(T, fO2_shift=dIW)
                assert math.isfinite(val)
                # log10(fO2) stays in [-40, +5] over 800-5000 K for dIW in [-3, +1.5];
                # catches a log10 -> ln unit-conversion slip (factor 2.3).
                assert -40 < val < 5


@pytest.mark.parametrize('bad_T', [0.0, -1.0, -300.0])
def test_oxygen_fugacity_nonpositive_T_raises(bad_T):
    """`T <= 0` raises `ValueError` for every buffer.

    Without the guard, T = 0 silently propagates `nan` through every
    downstream equilibrium constant via the `1/T` and `T * log(T)` terms.
    The error message must name the temperature so users debugging an IC
    file can find the offending input.
    """
    for buffer_name in ('fischer', 'oneill', 'hirschmann'):
        of = OxygenFugacity(buffer_name)
        with pytest.raises(ValueError, match='Temperature must be positive'):
            of(bad_T)


def test_oxygen_fugacity_unknown_buffer_name_raises():
    """Constructing with a name outside FO2_MODELS raises `ValueError`.

    The check runs at construction (eager), not at first call (lazy), so a
    config typo surfaces before any chemistry step runs. `None` and the
    names of non-model attributes (`__call__`, which a bare attribute
    lookup would resolve) are refused the same way; every name in
    FO2_MODELS still constructs and evaluates to a finite fO2.
    """
    from calliope.oxygen_fugacity import FO2_MODELS

    for bad in ('frost1991', 'typo_fisher', None, '__call__', 'Hirschmann'):
        with pytest.raises(ValueError, match='Unknown IW buffer'):
            OxygenFugacity(bad)
    for name in FO2_MODELS:
        assert np.isfinite(OxygenFugacity(name)(2000.0))


@pytest.mark.physics_invariant
def test_default_fo2_model_is_shared_by_oxygen_fugacity_and_chemistry():
    """OxygenFugacity and chemistry.ModifiedKeq default to the same single
    DEFAULT_FO2_MODEL constant, so a future change to the default cannot be
    applied to one but not the other.

    The behavioural check pins that the shared default actually dispatches to
    the Hirschmann (2021) buffer: the bare-default OxygenFugacity matches the
    Hirschmann value at 1473 K, and a regression that flipped the default to
    Fischer or O'Neill would land more than 0.2 dex away, outside the
    tolerance. 1473 K is chosen because Hirschmann and Fischer nearly cross
    near 2500 K, where they could not be told apart."""
    from calliope.chemistry import ModifiedKeq
    from calliope.oxygen_fugacity import DEFAULT_FO2_MODEL, OxygenFugacity

    assert DEFAULT_FO2_MODEL == 'hirschmann'

    # Both classes carry the same default in their signature: changing one
    # without the other is the trap this guards against.
    assert OxygenFugacity.__init__.__defaults__ == (DEFAULT_FO2_MODEL,)
    assert ModifiedKeq.__init__.__defaults__ == (DEFAULT_FO2_MODEL,)

    # The bare default dispatches to Hirschmann at the 1 bar reference: pin
    # the value and guard against both wrong-buffer regressions.
    T = 1473.0
    default_val = OxygenFugacity()(T)
    assert default_val == pytest.approx(-11.957, abs=2e-3)
    assert default_val == pytest.approx(OxygenFugacity('hirschmann').hirschmann(T), abs=1e-12)
    fischer_val = 6.94059 - 28.1808e3 / T
    assert abs(default_val - fischer_val) > 0.2  # default is not Fischer
    oneill_val = OxygenFugacity('oneill')(T)
    assert abs(default_val - oneill_val) > 0.2  # default is not O'Neill


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_oxygen_fugacity_hirschmann_value_at_1473K_matches_published_figure(validation_plot):
    """Hirschmann (2021) IW at 1 bar and T = 1473 K matches their Fig. 5.

    Fig. 5 of Hirschmann (2021, GCA 313, 74) plots the 1473 K IW value of
    O'Neill & Pownceby (1993) at log10 fO2 ~ -11.94 (read off the figure,
    so the tolerance is 0.05 dex), and their Fig. 4 shows the Table 1 fit
    sits within ~0.02 dex of that calibration at this T. The coded fit
    gives -11.957. Fischer (-12.19) and O'Neill & Eggins (-11.70) both
    land more than 0.2 dex away, so a silent dispatch to either fails.
    """
    of = OxygenFugacity('hirschmann')
    T, fig5, fig5_tol = 1473.0, -11.94, 0.05
    val = of(T)
    assert val == pytest.approx(fig5, abs=fig5_tol)
    # Wrong-buffer guards.
    fischer = OxygenFugacity('fischer')(T)
    oneill = OxygenFugacity('oneill')(T)
    assert abs(val - fischer) > 0.2
    assert abs(val - oneill) > 0.2
    # The dispatcher evaluates the 1 bar reference, not some other pressure:
    # 1 kbar already moves the buffer by ~0.04 dex at this T.
    assert val == pytest.approx(of.hirschmann(T, P_bar=1.0), abs=1e-12)
    kbar = of.hirschmann(T, P_bar=1.0e3)
    assert abs(val - kbar) > 0.03
    # Scale guard against a dropped log10 or a T-unit slip.
    assert -14 < val < -10

    def figure():
        import matplotlib.pyplot as plt

        grid = np.linspace(T - 30.0, T + 30.0, 25)
        fig, ax = plt.subplots(figsize=(8.4, 4.4))
        for model, ls in (('hirschmann', '-'), ('fischer', '--'), ('oneill', ':')):
            ax.plot(grid, [OxygenFugacity(model)(t) for t in grid], ls=ls, label=model)
        ax.errorbar(
            [T],
            [fig5],
            yerr=fig5_tol,
            fmt='s',
            mfc='none',
            capsize=4,
            label=r'Fig. 5, $\pm$0.05 dex',
        )
        ax.plot([T], [val], 'o', label='code, hirschmann')
        ax.plot([T], [fischer], 'x', label='code, fischer')
        ax.plot([T], [oneill], '+', ms=10, label='code, oneill')
        ax.plot([T], [kbar], '^', label='code, hirschmann at 1 kbar')
        ax.set(xlabel='T [K]', ylabel=r'$\log_{10} f\mathrm{O_2}$ [bar] at IW, 1 bar')
        ax.legend(fontsize='small', loc='upper left', bbox_to_anchor=(1.02, 1.0))
        fig.tight_layout()
        return fig

    validation_plot(figure, 'oxygen_fugacity_hirschmann_1473k')


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_oxygen_fugacity_hirschmann_low_pressure_slope_matches_volume_change(validation_plot):
    """Low-pressure slope of the Hirschmann (2021) IW buffer follows dV.

    For Fe + 1/2 O2 = FeO, d(log10 fO2)/dP = 2 dV / (R T ln 10), with
    dV = V(FeO) - V(Fe). V(FeO) = 12.256 cm3/mol is the stoichiometric
    intercept in Hirschmann (2021) Appendix A; V(Fe) = 55.845 / 7.874
    cm3/mol from the density of alpha iron. Room-temperature volumes
    ignore thermal expansion and stoichiometry change, so the fit slope
    over 0-1 GPa is expected within 15%, not exactly (it runs 7.5, 8.2 and
    10.4% high at 1200, 2000 and 3000 K).
    A bar/GPa slip in the pressure conversion changes the slope by 1e4.
    """
    of = OxygenFugacity('hirschmann')
    dV = (12.256 - 55.845 / 7.874) * 1e-6  # m3/mol
    R = 8.314462618
    temps, ratios, tol = (1200.0, 2000.0, 3000.0), [], 0.15
    for T in temps:
        slope = of.hirschmann(T, P_bar=1.0e4) - of.hirschmann(T, P_bar=1.0)  # dex/GPa
        expected = 2 * dV * 1e9 / (R * T * math.log(10))
        assert slope == pytest.approx(expected, rel=tol)
        # Sign guard: compressing the assemblage favours the denser Fe
        # side, so the buffer must move to higher fO2.
        assert slope > 0
        ratios.append(slope / expected)

    def figure():
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.axhspan(1.0 - tol, 1.0 + tol, color='0.85', label='tolerance, 15 %')
        ax.axhline(1.0, color='0.5', ls=':')
        ax.plot(temps, ratios, 'o', label='fit slope / volume estimate')
        ax.set(xlabel='T [K]', ylabel='slope ratio, 1 bar to 1 GPa', ylim=(0.7, 1.3))
        ax.legend()
        fig.tight_layout()
        return fig

    validation_plot(figure, 'oxygen_fugacity_hirschmann_pressure_slope')


@pytest.mark.physics_invariant
def test_oxygen_fugacity_hirschmann_branches_continuous_at_fcc_hcp_boundary():
    """The fcc/bcc and hcp fits agree at the fcc-hcp iron boundary.

    Hirschmann (2021) fit the two iron polymorph fields separately, with a
    reported maximum mismatch against the thermodynamic model of
    0.028 dex. Both branches fit the same underlying model, so the jump at
    the Eq. 18 boundary must stay within twice that mismatch; the fit's own
    jump is at most 0.0185 dex at these temperatures. A sign error in any
    nonzero Table 1 coefficient raises the largest jump to at least 0.069 dex
    (the square-root term of fcc/bcc c); the next smallest is 0.38 dex.
    """
    of = OxygenFugacity('hirschmann')
    for T in (1000.0, 1500.0, 2000.0, 2500.0, 3000.0):
        P_gpa = -18.64 + 0.04359 * T - 5.069e-6 * T**2
        # The boundary lies inside the 0-100 GPa calibration window.
        assert 0 < P_gpa < 100
        below = of.hirschmann(T, P_bar=P_gpa * 1e4 * (1 - 1e-6))
        above = of.hirschmann(T, P_bar=P_gpa * 1e4 * (1 + 1e-6))
        assert abs(above - below) < 0.056
        assert math.isfinite(below) and math.isfinite(above)


@pytest.mark.physics_invariant
def test_oxygen_fugacity_hirschmann_monotonic_in_pressure_and_vectorised():
    """Hirschmann (2021) IW rises strictly with P at fixed T, across the
    fcc-hcp switch, and the array path matches the scalar path.

    The 1 bar to 100 GPa sweep at 2000 K crosses the fcc-hcp boundary
    (~48 GPa), so a branch-selection bug that picked the wrong polymorph
    on one side would show up as a non-monotonic step or as an array /
    scalar mismatch.
    """
    of = OxygenFugacity('hirschmann')
    P_bar = np.array([1.0, 1e3, 1e4, 1e5, 4.8e5, 4.9e5, 7e5, 1e6])
    vals = of.hirschmann(2000.0, P_bar=P_bar)
    assert vals.shape == P_bar.shape
    assert np.all(np.diff(vals) > 0)
    # Scalar and array evaluations agree element by element.
    scalars = [of.hirschmann(2000.0, P_bar=p) for p in P_bar]
    np.testing.assert_allclose(vals, scalars, rtol=0, atol=1e-12)
    # Edge case: zero pressure is allowed (the P**0.5 term is defined).
    assert math.isfinite(of.hirschmann(2000.0, P_bar=0.0))
    # Error contract: negative pressure has no physical meaning and would
    # give nan through the P**0.5 term.
    with pytest.raises(ValueError, match='Pressure must be non-negative'):
        of.hirschmann(2000.0, P_bar=-1.0)


@pytest.mark.physics_invariant
def test_oxygen_fugacity_hirschmann_hcp_field_and_boundary_placement():
    """Hirschmann (2021) at 2000 K, where Eq. 18 puts the fcc-hcp boundary at
    48.26 GPa. Deep in the hcp field, at 90 GPa, the buffer is the hcp row of
    Table 1 (13.897), 0.66 dex above the fcc/bcc row (13.234), so a fit that
    never selects hcp fails. Each side of the boundary is pinned about 6.5 GPa
    away: 55 GPa is hcp (6.654, fcc 6.561) and 42 GPa is fcc/bcc (3.732, hcp
    3.676), so a boundary moved by 10 GPa either way puts one of them on the
    wrong branch. The pins are hand evaluations of the Table 1 rows, not
    values read from the paper.
    """
    of = OxygenFugacity('hirschmann')
    deep = of.hirschmann(2000.0, P_bar=90.0e4)
    assert deep == pytest.approx(13.897, abs=1e-3)
    assert abs(deep - 13.234) > 0.6  # fcc/bcc row at the same point

    above = of.hirschmann(2000.0, P_bar=55.0e4)
    below = of.hirschmann(2000.0, P_bar=42.0e4)
    assert above == pytest.approx(6.654, abs=1e-3)
    assert abs(above - 6.561) > 0.08  # fcc/bcc row
    assert below == pytest.approx(3.732, abs=1e-3)
    assert abs(below - 3.676) > 0.05  # hcp row
    # Boundedness: the buffer rises with pressure through the transition.
    assert below < above < deep


@pytest.mark.physics_invariant
def test_oxygen_fugacity_hirschmann_dispatcher_matches_the_full_fit_at_1bar():
    """The dispatcher evaluates Hirschmann (2021) at 1 bar from precomputed
    fcc/bcc coefficients, skipping the pressure polynomials and the fcc-hcp
    branch test. It must give exactly ``hirschmann(T, P_bar=1.0)`` across
    300-9000 K, including both ends of the held boundary range (1000 K, ~4300 K)
    and the ~451 K and ~8148 K points where the unheld parabola reaches 1 bar,
    and the shift must stay additive. At 10 GPa the full fit differs, so the
    fast path is a 1 bar shortcut only."""
    of = OxygenFugacity('hirschmann')
    temperatures = np.concatenate(
        [np.geomspace(300.0, 9000.0, 200), [451.0, 1000.0, 4300.0, 8148.0]]
    )
    for T in temperatures:
        assert of(T) == pytest.approx(of.hirschmann(T, P_bar=1.0), rel=1e-15, abs=1e-12), T
        assert of(T, 3.0) == pytest.approx(of.hirschmann(T) + 3.0, abs=1e-12), T
    assert abs(of.hirschmann(2000.0, P_bar=1.0e5) - of(2000.0)) > 0.5
    with pytest.raises(ValueError, match='Temperature must be positive'):
        of(-1.0)


@pytest.mark.physics_invariant
def test_oxygen_fugacity_hirschmann_1bar_stays_on_fcc_branch_outside_boundary_fit():
    """The Eq. 18 parabola drops below 1 bar under ~451 K and above ~8148 K.
    With the boundary held at its 1000 K and vertex (~4300 K) values there,
    the 1 bar buffer stays on the fcc/bcc row: -85.418 at 300 K (hcp row
    -85.200) and 2.317 at 8150 K (hcp row -40.125, a 42 dex fall). Across
    8148 K the 1 bar buffer is continuous and keeps rising with T. The full fit
    is evaluated, as the 1 bar dispatcher skips the boundary.
    """
    of = OxygenFugacity('hirschmann')
    iw = of.hirschmann
    assert iw(300.0) == pytest.approx(-85.418, abs=1e-3)
    assert abs(iw(300.0) - (-85.200)) > 0.2  # hcp row
    assert iw(8150.0) == pytest.approx(2.317, abs=1e-3)
    assert abs(iw(8140.0) - iw(8150.0)) < 0.01
    assert iw(8000.0) < iw(8150.0) < iw(9000.0)
    with pytest.raises(ValueError, match='Temperature must be positive'):
        of(0.0)
