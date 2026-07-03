# Validation: `src/calliope/solubility.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.solubility` against published sources.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_solubility.py::TestSolubilityH2O::test_peridotite_default_matches_sossi_2023_fit` | Sossi et al. (2023) [^cite-sossi2023] peridotite H2O fit: `ppmw = 524 * p^0.5` | `src/calliope/solubility.py:39-41` (`peridotite`) | Pins the default H2O parameterization at p = 100 bar against the Sossi 2023 constant; includes wrong-law (Dixon et al. 1995 [^cite-dixon1995] basalt) and wrong-exponent (1.0 vs 0.5) discrimination guards. |
| `tests/test_solubility.py::TestSolubilityS2_xFeO::test_default_call_matches_gaillard_2022_earth_mantle_value` | Gaillard et al. (2022) [^cite-gaillard2022], EPSL 117255, S2 solubility law with `x_FeO = 10 wt%` Earth-mantle default | [doi:10.1016/j.epsl.2021.117255](https://doi.org/10.1016/j.epsl.2021.117255), `src/calliope/solubility.py:75-93` | Pins S2 ppmw at (p = 1 bar, T = 2500 K, fO2_shift = 0) against the closed-form Gaillard expression; couples to the Fischer 2011 [^cite-fischer2011] IW buffer default. |
| `tests/test_solubility.py::TestSolubilityNobleGas::test_every_gas_constant_pinned_to_independent_literal` | Jambon, Weill & Braun (1986) [^cite-jambon1986] noble gas Henry constants in tholeiitic basalt | [doi:10.1016/0016-7037(86)90193-6](https://doi.org/10.1016/0016-7037(86)90193-6), `src/calliope/solubility.py` (`jambon86_ppmw_per_bar`, `SolubilityNobleGas`) | Pins each of the five noble gas Henry constants [ppmw/bar] to a hand-computed literal with cross-gas discrimination guards, so a mistyped constant fails on a CI image without atmodeller. |
| `tests/test_solubility.py::TestSolubilityNobleGas::test_noble_constants_cross_check_against_atmodeller` | Cross-implementation cross-check against atmodeller | `src/calliope/solubility.py` | Confirms CALLIOPE and atmodeller derive identical Jambon et al. (1986) constants from the same primitives (skips when atmodeller is not installed). |

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
log10(fO2) = -4.3317 (Fischer at T=2500 K, dIW=0)
fO2 = 10^-4.3317 = 4.659e-5 bar
ln(p_S2/fO2) = ln(1 / 4.659e-5) = 9.9742
ln(X) = 13.8426 - 10.5904 + 1.24 + 0.5 * 9.9742
      = 3.2522 + 1.24 + 4.9871
      = 9.4793
ppmw = exp(9.4793) ~ 13086
```

Hidden coupling: the pin depends on the Fischer 2011 IW buffer; any change
to the buffer coefficients in `oxygen_fugacity.py` requires regenerating this
number. The test docstring flags this coupling explicitly. The reference
value in the test (13085.87) matches the closed-form ppmw above to the
4-significant-figure precision of the Fischer coefficients.

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
  `x_TiO2` kwargs (Dasgupta 2022 [^cite-dasgupta2022]) plus Libourel 2003 [^cite-libourel2003] as the no-composition baseline.
- `docs/Explanations/solubility.md`: user-facing concept page with the
  validity envelope per parameterization.

## References

 [^cite-dasgupta2022]: R. Dasgupta, E. Falksen, A. Pal, C. Sun, *[The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions](https://doi.org/10.1016/j.gca.2022.09.012)*, Geochimica et Cosmochimica Acta, 336, 291-307, 2022. [SciX](https://scixplorer.org/abs/2022GeCoA.336..291D/abstract).
 [^cite-dixon1995]: J. E. Dixon, E. M. Stolper, J. R. Holloway, *[An experimental study of water and carbon dioxide solubilities in mid-ocean ridge basaltic liquids. Part I: Calibration and solubility models](https://doi.org/10.1093/oxfordjournals.petrology.a037267)*, Journal of Petrology, 36(6), 1607-1631, 1995. [SciX](https://scixplorer.org/abs/1995JPet...36.1607D/abstract).
 [^cite-fischer2011]: R. A. Fischer, A. J. Campbell, G. A. Shofner, O. T. Lord, P. Dera, V. B. Prakapenka, *[Equation of state and phase diagram of FeO](https://doi.org/10.1016/j.epsl.2011.02.025)*, Earth and Planetary Science Letters, 304, 496-502, 2011. [SciX](https://scixplorer.org/abs/2011E%26PSL.304..496F/abstract).
 [^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
 [^cite-jambon1986]: A. Jambon, H. Weill, J. Braun, *[Solubility of He, Ne, Ar, Kr and Xe in a basalt melt in the range 1250-1600 C. Geochemical implications](https://doi.org/10.1016/0016-7037(86)90193-6)*, Geochimica et Cosmochimica Acta, 50(3), 401-408, 1986. [SciX](https://scixplorer.org/abs/1986GeCoA..50..401J/abstract).
 [^cite-libourel2003]: G. Libourel, B. Marty, F. Humbert, *[Nitrogen solubility in basaltic melt. Part I. Effect of oxygen fugacity](https://doi.org/10.1016/S0016-7037(03)00259-X)*, Geochimica et Cosmochimica Acta, 67(21), 4123-4135, 2003. [SciX](https://scixplorer.org/abs/2003GeCoA..67.4123L/abstract).
 [^cite-sossi2023]: P. A. Sossi, P. M. E. Tollan, J. Badro, D. J. Bower, *[Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets](https://doi.org/10.1016/j.epsl.2022.117894)*, Earth and Planetary Science Letters, 601, 117894, 2023. [SciX](https://scixplorer.org/abs/2023E%26PSL.60117894S/abstract).
