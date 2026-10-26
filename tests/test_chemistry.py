"""Tests for `src/calliope/chemistry.py`.

Exercises the `ModifiedKeq` equilibrium-constant dispatcher and its
two reaction families:

- JANAF tables: `janaf_H2`, `janaf_CO`, `janaf_SO2`, `janaf_H2S`,
  `janaf_NH3`. Coefficients are fits to the NIST JANAF Thermochemical
  Tables over 1500-3000 K (`janaf_H2`, `janaf_CO`) or 700-5000 K
  (`janaf_H2S`, `janaf_NH3`, checked against every JANAF row in that range).
- Schaefer & Fegley series: `schaefer_H`, `schaefer_C`, `schaefer_CH4`.

Reference pin: `janaf_H2` Geq at T = 2000 K under the legacy O'Neill
2002 IW buffer matches the closed-form `10^(Keq - 0.5 * log10_fO2)`
within published precision. Includes a wrong-model discrimination
guard against `schaefer_H` at the same conditions.

Physics invariants:
- Keq > 0 for every model over the 1500-3000 K range (Geq is `10^x`).
- Monotonic decrease in Geq with positive `fO2_shift` for models with
  positive `fO2_stoich` (more oxidising -> less reduced product).
- All JANAF and Schaefer call returns are finite for the documented T range.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pytest

from calliope import chemistry
from calliope.chemistry import ModifiedKeq

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_modified_keq_janaf_H2_matches_closed_form_at_2000K_with_oneill(validation_plot):
    """`janaf_H2` Geq at T = 2000 K under O'Neill 2002 IW matches the
    closed-form `10^(Keq - 0.5 * log10_fO2)`.

    The `janaf_H2` coefficients (`a = -13152.48`, `b = 3.038586`, stoich
    coefficient 0.5) are fits to the NIST JANAF Thermochemical Tables
    for the reaction `H2O = H2 + 0.5 O2` over 1500-3000 K.

    At T = 2000 K with `fO2_model='oneill'` and `fO2_shift = 0`:

    ```
    Keq = 10^(-13152.48 / 2000 + 3.038586) = 10^-3.5378 ~ 2.90e-4
    log10(fO2)_oneill(2000) ~ -7.4078
    Geq = 10^(Keq - 0.5 * fO2) = 10^(-3.5378 - 0.5 * -7.4078)
        = 10^(0.1661) ~ 1.467
    ```

    Discrimination guard: `schaefer_H` at the same conditions has
    coefficients `(-12794 / T + 2.7768, 0.5)` -> Keq = 10^-3.6202 and
    Geq = 10^(0.0837) ~ 1.213. A regression that silently dispatched
    to the wrong reaction would land 0.25 units away.
    """
    T = 2000.0
    mk = ModifiedKeq('janaf_H2', fO2_model='oneill')
    g = mk(T, fO2_shift=0.0)
    # Closed-form expectation: 1.467.
    expected = 10 ** (
        (-13152.477779978302 / T + 3.038586383273608)
        - 0.5
        * (2 * (-244118 + 115.559 * T - 8.474 * T * math.log(T)) / (math.log(10) * 8.31441 * T))
    )
    assert g == pytest.approx(expected, rel=1e-6)
    # Wrong-model discrimination: schaefer_H at the same T gives ~1.213.
    mk_schaefer = ModifiedKeq('schaefer_H', fO2_model='oneill')
    wrong = mk_schaefer(T, fO2_shift=0.0)
    assert abs(g - wrong) > 0.2
    # Sign + scale guards: Geq is positive (10^x) and of order unity.
    assert g > 0
    assert 0.5 < g < 5.0

    def figure():
        import matplotlib.pyplot as plt

        grid = np.linspace(1500.0, 3000.0, 61)
        fig, ax = plt.subplots()
        ax.plot(grid, [math.log10(mk(t, fO2_shift=0.0)) for t in grid], label='janaf_H2')
        ax.plot(
            grid,
            [math.log10(mk_schaefer(t, fO2_shift=0.0)) for t in grid],
            ls='--',
            label='schaefer_H',
        )
        ax.plot([T], [math.log10(expected)], 's', ms=10, mfc='none', label='closed form')
        ax.plot([T], [math.log10(g)], 'o', label='code, janaf_H2')
        ax.plot([T], [math.log10(wrong)], 'x', label='code, schaefer_H')
        ax.set(xlabel='T [K]', ylabel=r'$\log_{10} G_\mathrm{eq}$ (H$_2$/H$_2$O)')
        ax.legend()
        fig.tight_layout()
        return fig

    validation_plot(figure, 'chemistry_janaf_h2_2000k')


@pytest.mark.physics_invariant
def test_modified_keq_janaf_H2_decreases_with_oxidising_shift():
    """Positive `fO2_shift` (more oxidising) decreases Geq for any
    reaction with positive `fO2_stoich`.

    The dispatcher returns `Geq = 10^(Keq - fO2_stoich * log10_fO2)`.
    For `janaf_H2` (stoich = 0.5), increasing `fO2_shift` raises
    `log10_fO2`, lowering `Keq - 0.5 * fO2`, so `Geq` drops.
    Discrimination: a sign flip on `fO2_stoich` would invert this
    ordering.
    """
    T = 2000.0
    mk = ModifiedKeq('janaf_H2', fO2_model='oneill')
    g_reducing = mk(T, fO2_shift=-2.0)
    g_neutral = mk(T, fO2_shift=0.0)
    g_oxidising = mk(T, fO2_shift=+2.0)
    # Strict ordering: more oxidising -> lower Geq for stoich > 0.
    assert g_reducing > g_neutral > g_oxidising
    # Discrimination: the delta from -2 to +2 spans roughly 2 dex on
    # log10(Geq) for stoich = 0.5, so the ratio g_reducing / g_oxidising
    # should be ~100. Pin a wide envelope to catch coefficient-only bugs.
    assert g_reducing / g_oxidising > 10


@pytest.mark.physics_invariant
def test_modified_keq_fO2_model_choice_changes_result():
    """Switching the underlying IW buffer (Fischer vs O'Neill) changes
    Geq at the same (T, fO2_shift). Buffer-flip propagation guard.

    A regression that silently dispatched to the wrong fO2 model would
    leave this test passing for a wrong reason if the test only checked
    one buffer; instead, pin both and assert the difference exceeds the
    rel=1e-6 tolerance.
    """
    T = 2000.0
    g_oneill = ModifiedKeq('janaf_H2', fO2_model='oneill')(T, fO2_shift=0.0)
    g_fischer = ModifiedKeq('janaf_H2', fO2_model='fischer')(T, fO2_shift=0.0)
    # The two buffers differ by ~0.26 dex at 2000 K; for stoich = 0.5
    # the Geq ratio differs by ~10^0.13 ~ 1.35.
    assert g_oneill != pytest.approx(g_fischer, rel=1e-6)
    # Both must be positive and finite.
    assert g_oneill > 0 and g_fischer > 0
    assert math.isfinite(g_oneill) and math.isfinite(g_fischer)


@pytest.mark.physics_invariant
def test_modified_keq_janaf_CO_positive_and_finite():
    """`janaf_CO` Geq is positive and finite at T = 1800 K."""
    T = 1800.0
    mk = ModifiedKeq('janaf_CO', fO2_model='oneill')
    g = mk(T, fO2_shift=0.0)
    # Geq = 10^x is always > 0 for finite x; this catches a regression
    # that introduced a `log()` of a non-positive intermediate.
    assert g > 0.0
    assert math.isfinite(g)
    # Discrimination: at T = 1800 K and oneill, the closed form gives
    # log10(Keq) = -14467.51/1800 + 4.348 = -3.689, log10(fO2) ~ -8.34,
    # log10(Geq) = -3.689 + 0.5 * 8.34 = 0.481, Geq ~ 3.0.
    assert 1.0 < g < 10.0


@pytest.mark.physics_invariant
def test_modified_keq_schaefer_models_finite_over_calliope_use_range():
    """All three Schaefer models return finite, positive Geq for T in
    [1500, 3000] K (the documented CALLIOPE use range)."""
    for T in (1500.0, 2200.0, 3000.0):
        for model in ('schaefer_H', 'schaefer_C', 'schaefer_CH4'):
            mk = ModifiedKeq(model, fO2_model='oneill')
            g = mk(T, fO2_shift=0.0)
            assert math.isfinite(g)
            assert g > 0


def test_modified_keq_unknown_model_raises():
    """Constructing with an unknown reaction model raises `AttributeError`.

    Dispatcher uses `getattr(self, Keq_model)` so a typo lands at
    construction, not at first call. Eager failure: chemistry-step
    typos surface at config parse, not mid-simulation.
    """
    with pytest.raises(AttributeError):
        ModifiedKeq('typo_janaf')
    with pytest.raises(AttributeError):
        ModifiedKeq('janaf_O3')  # not an implemented reaction


# JANAF fit range from the default PROTEUS outgassing floor (700 K) to 5000 K.
# Tolerance = max fit residual plus margin [dex].
JANAF_FIT_RANGE_K = (700.0, 5000.0)
JANAF_FIT_TOL_DEX = {'H2S': 0.025, 'NH3': 0.04}
# Wrong constants (doubled form) the row test must reject: H2S with the S2 coefficient 1
# instead of 1/2, and an NH3 fit that misses the rows by up to 0.158 dex (H2S 1.853 dex).
DISCRIMINATING_FIT = {'H2S': (13462.03094, 7.24546062), 'NH3': (5328.031246, 11.98476092)}


def _per_mole_log_k(product, T):
    """Per-mole log10 K from the doubled-form `janaf_<product>` method."""
    log_k2, stoich = getattr(ModifiedKeq(f'janaf_{product}'), f'janaf_{product}')(T)
    assert stoich == 0
    return 0.5 * log_k2


def _rows_in_fit_range(janaf_log_k, product):
    """JANAF rows (T, per-mole log10 K) inside the fit range."""
    T, y = janaf_log_k(product)
    keep = (T >= JANAF_FIT_RANGE_K[0]) & (T <= JANAF_FIT_RANGE_K[1])
    return T[keep], y[keep]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
@pytest.mark.parametrize('product', ['H2S', 'NH3'], ids=['hydrogen_sulfide', 'ammonia'])
def test_janaf_constant_matches_janaf_rows_over_fit_range(
    janaf_log_k, product, validation_plot
):
    """The doubled-form constant matches every NIST-JANAF row from 700 to 5000 K.

    Reference: JANAF 4th ed. log Kf columns in `tools/` (H2S: log Kf(H2S) -
    0.5 log Kf(S2); NH3: log Kf(NH3)). Both range ends are table rows, so the edges are tested.
    """
    T, janaf = _rows_in_fit_range(janaf_log_k, product)
    assert T[0] == pytest.approx(JANAF_FIT_RANGE_K[0])
    assert T[-1] == pytest.approx(JANAF_FIT_RANGE_K[1])
    fit = _per_mole_log_k(product, T)
    assert np.max(np.abs(fit - janaf)) < JANAF_FIT_TOL_DEX[product]
    # Discrimination: the wrong fits miss the rows by far more.
    a, b = DISCRIMINATING_FIT[product]
    assert np.max(np.abs(0.5 * (a / T - b) - janaf)) > 3 * JANAF_FIT_TOL_DEX[product]
    # Sign: the formation reactions are exothermic, so log K falls with T.
    assert np.all(np.diff(fit) < 0)

    def figure():
        import matplotlib.pyplot as plt

        tol = JANAF_FIT_TOL_DEX[product]
        fig, (top, bottom) = plt.subplots(2, 1, sharex=True, height_ratios=[2, 1])
        top.plot(T, janaf, 'o', ms=3, label='JANAF rows')
        top.plot(T, fit, label=f'janaf_{product}')
        top.plot(T, 0.5 * (a / T - b), ls='--', label='wrong fit')
        top.set(ylabel=r'$\log_{10} K$ per mole')
        top.legend()
        bottom.plot(T, fit - janaf, 'o', ms=3)
        bottom.axhline(tol, color='0.5', ls=':')
        bottom.axhline(-tol, color='0.5', ls=':')
        bottom.set(xlabel='T [K]', ylabel='code - JANAF [dex]')
        fig.tight_layout()
        return fig

    validation_plot(figure, f'chemistry_janaf_{product.lower()}_rows')


@pytest.mark.reference_pinned
@pytest.mark.parametrize('product', ['H2S', 'NH3'], ids=['hydrogen_sulfide', 'ammonia'])
def test_janaf_coefficients_are_least_squares_fit_of_janaf_rows(
    janaf_log_k, product, validation_plot
):
    """The coefficients equal a least-squares a/T - b fit of the JANAF rows in range.

    A changed fit range, form or coefficient moves the fit away from the code.
    """
    T, janaf = _rows_in_fit_range(janaf_log_k, product)
    slope, intercept = np.polyfit(1.0 / T, 2.0 * janaf, 1)
    fit = 2.0 * _per_mole_log_k(product, T)
    # atol because the doubled H2S log K crosses zero near 1830 K.
    np.testing.assert_allclose(fit, slope / T + intercept, rtol=1e-6, atol=1e-6)
    assert fit[0] > fit[-1]

    def figure():
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(1e4 / T, slope / T + intercept, label='least-squares fit of the JANAF rows')
        ax.plot(1e4 / T, fit, 'o', ms=3, label=f'janaf_{product}')
        ax.set(xlabel=r'$10^4 / T$ [K$^{-1}$]', ylabel=r'$2 \log_{10} K$')
        ax.legend()
        fig.tight_layout()
        return fig

    validation_plot(figure, f'chemistry_janaf_{product.lower()}_fit')


@pytest.fixture
def fit_range_warnings(monkeypatch, caplog):
    """Fresh once-per-process warning state; returns the log capture."""
    monkeypatch.setattr(chemistry, '_WARNED_BELOW_FIT', set())
    caplog.set_level(logging.WARNING, logger='fwl.calliope.chemistry')
    return caplog


@pytest.mark.parametrize('product', ['H2S', 'NH3'], ids=['hydrogen_sulfide', 'ammonia'])
def test_janaf_fit_warns_once_below_fit_range(fit_range_warnings, product):
    """Below 700 K the fit extrapolates: one warning per process and constant, naming the lowest T."""
    method = getattr(ModifiedKeq(f'janaf_{product}'), f'janaf_{product}')
    method(JANAF_FIT_RANGE_K[0])
    assert fit_range_warnings.records == []
    for T in (np.array([800.0, 650.0]), 699.0, 300.0):
        method(T)
    assert len(fit_range_warnings.records) == 1
    record = fit_range_warnings.records[0]
    assert record.levelno == logging.WARNING
    assert f'janaf_{product}' in record.getMessage()
    assert 'T = 650.0 K' in record.getMessage()
    assert f'fit range starts at {JANAF_FIT_RANGE_K[0]:g} K' in record.getMessage()
    other = 'NH3' if product == 'H2S' else 'H2S'
    getattr(ModifiedKeq(f'janaf_{other}'), f'janaf_{other}')(650.0)
    assert len(fit_range_warnings.records) == 2


def test_janaf_fit_warning_edges_and_other_constants(fit_range_warnings):
    """No warning from 700 K up, at inf, for NaN, for an empty array, or for the other constants."""
    for name in 'janaf_H2 janaf_CO janaf_SO2 schaefer_CH4 schaefer_C schaefer_H'.split():
        getattr(ModifiedKeq(name), name)(500.0)
    for name in ('janaf_H2S', 'janaf_NH3'):
        method = getattr(ModifiedKeq(name), name)
        method(np.array([700.0, 1500.0, 5000.0, np.inf]))
        method(2000)
        method(np.array([700, 5000]))
        method(float('nan'))
        assert method(np.array([]))[0].size == 0
    assert fit_range_warnings.records == []
    assert chemistry._WARNED_BELOW_FIT == set()
    # A NaN entry is skipped; the largest float below 700 K still warns.
    just_below = math.nextafter(JANAF_FIT_RANGE_K[0], 0.0)
    ModifiedKeq('janaf_H2S').janaf_H2S(np.array([np.nan, just_below]))
    assert len(fit_range_warnings.records) == 1
    assert f'T = {just_below} K' in fit_range_warnings.records[0].getMessage()
    assert chemistry._WARNED_BELOW_FIT == {'janaf_H2S'}


@pytest.mark.parametrize(
    'T',
    [
        650.0,
        np.float64(650.0),
        np.float32(650.0),
        650,
        np.int64(650),
        np.array(650.0),
        np.array([900.0, 650.0, 800.0]),
        np.array([[900, 650], [800, 720]]),
    ],
    ids=['float', 'float64', 'float32', 'int', 'int64', 'zero_d', 'minimum_inside', 'int_2d'],
)
def test_janaf_fit_warning_for_every_numeric_form_of_T(fit_range_warnings, T):
    """Each numeric form of T gives a finite constant and one warning that names its lowest entry."""
    log_k, stoich = ModifiedKeq('janaf_H2S').janaf_H2S(T)
    assert np.all(np.isfinite(log_k))
    assert stoich == 0
    assert len(fit_range_warnings.records) == 1
    assert 'T = 650.0 K' in fit_range_warnings.records[0].getMessage()
