"""Tests for `src/calliope/chemistry.py`.

Exercises the `ModifiedKeq` equilibrium-constant dispatcher and its
two reaction families:

- JANAF tables: `janaf_H2`, `janaf_CO`, `janaf_SO2`, `janaf_H2S`,
  `janaf_NH3`. Coefficients are fits to the NIST JANAF Thermochemical
  Tables over 1500-3000 K (`janaf_H2`, `janaf_CO`) or 700-5000 K
  (`janaf_H2S`, checked against every JANAF row in that range).
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

import math

import numpy as np
import pytest

from calliope.chemistry import ModifiedKeq

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_modified_keq_janaf_H2_matches_closed_form_at_2000K_with_oneill():
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


# JANAF fit range: the PROTEUS outgassing floor (700 K) to above the hottest
# magma-ocean surfaces in coupled runs. Tolerance = max fit residual plus margin [dex].
JANAF_FIT_RANGE_K = (700.0, 5000.0)
JANAF_FIT_TOL_DEX = {'H2S': 0.025}
# Fit with the S2 coefficient 1 instead of 1/2 over 200-4000 K (doubled form).
WRONG_S2_COEFFICIENT_FIT = {'H2S': (13462.03094, 7.24546062)}


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
@pytest.mark.parametrize('product', ['H2S'], ids=['hydrogen_sulfide'])
def test_janaf_constant_matches_janaf_rows_over_fit_range(janaf_log_k, product):
    """The doubled-form constant matches every NIST-JANAF row from 700 to 5000 K.

    Reference: JANAF 4th ed. log Kf columns in `tools/` (H2S: log Kf(H2S) -
    0.5 log Kf(S2)). Both range ends are table rows, so the edges are tested.
    """
    T, janaf = _rows_in_fit_range(janaf_log_k, product)
    assert T[0] == pytest.approx(JANAF_FIT_RANGE_K[0])
    assert T[-1] == pytest.approx(JANAF_FIT_RANGE_K[1])
    resid = _per_mole_log_k(product, T) - janaf
    assert np.max(np.abs(resid)) < JANAF_FIT_TOL_DEX[product]
    # Discrimination: the S2-coefficient-1 fit misses the rows by far more.
    a, b = WRONG_S2_COEFFICIENT_FIT[product]
    assert np.max(np.abs(0.5 * (a / T - b) - janaf)) > 4 * JANAF_FIT_TOL_DEX[product]
    # Sign: the formation reactions are exothermic, so log K falls with T.
    assert np.all(np.diff(_per_mole_log_k(product, T)) < 0)
    # Scale: per-mole log10 K stays within a few dex (ln instead of log10 is x2.3).
    assert np.all(np.abs(_per_mole_log_k(product, T)) < 10)


@pytest.mark.reference_pinned
@pytest.mark.parametrize('product', ['H2S'], ids=['hydrogen_sulfide'])
def test_janaf_coefficients_are_least_squares_fit_of_janaf_rows(janaf_log_k, product):
    """The coefficients equal a least-squares a/T - b fit of the JANAF rows in range.

    A changed fit range, form or coefficient moves the fit away from the code.
    """
    T, janaf = _rows_in_fit_range(janaf_log_k, product)
    slope, intercept = np.polyfit(1.0 / T, 2.0 * janaf, 1)
    T1, T2 = JANAF_FIT_RANGE_K
    k1, k2 = 2.0 * _per_mole_log_k(product, T1), 2.0 * _per_mole_log_k(product, T2)
    a = (k1 - k2) / (1.0 / T1 - 1.0 / T2)
    assert a == pytest.approx(slope, rel=1e-6)
    assert a / T1 - k1 == pytest.approx(-intercept, rel=1e-6)
    # Sign and scale guards: exothermic slope; a in K is thousands, not tens.
    assert a > 0
    assert 1e3 < a < 1e5
