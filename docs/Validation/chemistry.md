# Validation: `src/calliope/chemistry.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.chemistry` against published sources.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_chemistry.py::test_modified_keq_janaf_H2_matches_closed_form_at_2000K_with_oneill` | NIST JANAF Thermochemical Tables (4th ed.), fits for `H2O = H2 + 0.5 O2` over 1500-3000 K | `src/calliope/chemistry.py:41-43` (`janaf_H2`) | Pins `Geq(janaf_H2)` at T = 2000 K under the O'Neill 2002 IW buffer against the closed-form `10^(Keq - 0.5 * log10_fO2)`; includes a wrong-reaction discrimination guard against `schaefer_H` at the same conditions. |

## Re-derivation note

`ModifiedKeq` returns `Geq = 10^(Keq(T) - fO2_stoich * log10_fO2(T, dIW))`.
For `janaf_H2` the coefficients `(a = -13152.48, b = 3.038586, stoich = 0.5)`
fit the JANAF tables for `H2O = H2 + 0.5 O2` over the documented
1500-3000 K CALLIOPE-use range.

At T = 2000 K, `fO2_model='oneill'`, `fO2_shift = 0`:

```
Keq = 10^(-13152.48 / 2000 + 3.038586)
    = 10^(-6.5762 + 3.0386)
    = 10^-3.5376
    ~ 2.901e-4
```

The O'Neill 2002 IW buffer at the same T:

```
log10(fO2) = 2 * (-244118 + 115.559 * 2000 - 8.474 * 2000 * ln(2000))
             / (ln(10) * 8.31441 * 2000)
           ~ -7.4078
```

Combining:

```
log10(Geq) = -3.5376 - 0.5 * (-7.4078) = 0.1663
Geq = 10^0.1663 ~ 1.467
```

The test reconstructs the closed form symbolically (without intermediate
rounding) and pins to `rel=1e-6` precision.

Wrong-reaction discrimination guard: `schaefer_H` at the same T:

```
Keq_schaefer = 10^(-12794 / 2000 + 2.7768)
             = 10^-3.6232 ~ 2.380e-4
Geq_schaefer = 10^(-3.6232 + 3.7039) = 10^0.0807 ~ 1.205
```

The 0.26 difference between `Geq_janaf_H2 ~ 1.467` and `Geq_schaefer_H ~ 1.205`
exceeds the test's discrimination guard threshold (`abs(g - wrong) > 0.2`).

## Anchor type

Published benchmark (JANAF tables) plus cross-reaction discrimination
(janaf vs schaefer) plus a separate buffer-flip-propagation test
(`test_modified_keq_fO2_model_choice_changes_result`) that asserts
switching from O'Neill to Fischer changes Geq beyond rel=1e-6.

## Cross-references

- `src/calliope/chemistry.py:11-22`: `ModifiedKeq` constructor and `__call__`.
- `src/calliope/chemistry.py:24-62`: per-reaction coefficient methods.
- `docs/Explanations/equilibrium_chemistry.md`: user-facing concept page.
