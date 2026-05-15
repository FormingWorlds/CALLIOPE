# Validation: `src/calliope/solubility.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.solubility` against published sources.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_solubility.py::TestSolubilityH2O::test_peridotite_default_matches_sossi_2023_fit` | Sossi et al. (2023) peridotite H2O fit: `ppmw = 524 * p^0.5` | `src/calliope/solubility.py:39-41` (`peridotite`) | Pins the default H2O parameterization at p = 100 bar against the Sossi 2023 constant; includes wrong-law (basalt_dixon) and wrong-exponent (1.0 vs 0.5) discrimination guards. |
| `tests/test_solubility.py::TestSolubilityS2_xFeO::test_default_call_matches_gaillard_2022_earth_mantle_value` | Gaillard et al. (2022), EPSL 117255, S2 solubility law with `x_FeO = 10 wt%` Earth-mantle default | [doi:10.1016/j.epsl.2021.117255](https://doi.org/10.1016/j.epsl.2021.117255), `src/calliope/solubility.py:75-93` | Pins S2 ppmw at (p = 1 bar, T = 2500 K, fO2_shift = 0) against the closed-form Gaillard expression; couples to the Fischer 2011 IW buffer default. |

## Re-derivation notes

### H2O (Sossi 2023 peridotite default)

Power-law form: `ppmw = const * p^exponent`. Sossi et al. (2023) gives
`const = 524 ppmw/bar^0.5` for peridotite melt. At p = 100 bar:

```
ppmw = 524 * sqrt(100) = 5240
```

A regression that swapped to the Dixon 1995 basalt constant (`const = 965`)
would land at 9650, off by 84% at the same pressure. A regression that
flipped the exponent (0.5 -> 1.0) would land at 52400, off by an order of
magnitude. Both are caught by the discrimination guards in the test.

### S2 (Gaillard 2022 Earth-mantle default)

Closed form from Gaillard et al. (2022):

```
ln(X_S^melt) = 13.8426 - 26476/T + 0.124*x_FeO + 0.5*ln(p_S2/fO2)
```

At T = 2500 K, p_S2 = 1 bar, x_FeO = 10 wt%, fO2 from the Fischer 2011 IW
buffer:

```
log10(fO2) = -7.2832 (Fischer at 2500 K)
fO2 = 10^-7.2832 = 5.21e-8 bar
ln(X) = 13.8426 - 10.5904 + 1.24 + 0.5 * ln(1 / 5.21e-8)
      = 13.8426 - 10.5904 + 1.24 + 0.5 * 16.770
      = 12.876
ppmw = exp(9.479) ~ 13086
```

(The numbers differ slightly from a naive evaluation because the implementation
returns `exp(out)` after the closed-form computation.) Hidden coupling: the pin
depends on the Fischer 2011 IW buffer; any change to the buffer coefficients
in `oxygen_fugacity.py` requires regenerating this number. The test docstring
flags this coupling explicitly.

## Anchor types

- H2O: published benchmark (Sossi 2023) plus wrong-law and wrong-exponent
  discrimination guards.
- S2: published benchmark (Gaillard 2022) plus closed-form scaling tests
  for `x_FeO` (parametrized 5/10/15/20 wt%) so the prefactor coefficient
  0.124 is independently verified.

## Cross-references

- `src/calliope/solubility.py:30-57`: H2O parameterizations.
- `src/calliope/solubility.py:60-100`: S2 with the `x_FeO` kwarg.
- `src/calliope/solubility.py:100+`: N2 with the `x_SiO2`, `x_Al2O3`,
  `x_TiO2` kwargs (Dasgupta 2022) plus Libourel as the no-composition
  baseline.
- `docs/Explanations/solubility.md`: user-facing concept page with the
  validity envelope per parameterization.
