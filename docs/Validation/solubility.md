# Validation: `src/calliope/solubility.py`

## What this page checks

This page checks the solubility laws of CALLIOPE.
A solubility law gives the mass fraction of a gas that dissolves in the silicate melt, in ppmw (parts per million by weight), as a function of the partial pressure of that gas \[bar\] and, for some gases, of the temperature and the oxygen fugacity.
The oxygen fugacity $f\mathrm{O_2}$ \[bar\] measures how oxidising the melt is; CALLIOPE gives it as a shift $\Delta\mathrm{IW}$, in log10 units, from the iron-wüstite (IW) buffer (see the [oxygen fugacity page](oxygen_fugacity.md)).
Validation means that a test runs the code on a case where another source gives the right answer, and checks the result within a stated tolerance.
If a law is wrong, the melt holds the wrong amount of a volatile, so the atmosphere above it gets too much or too little of it.
Each case below says why the check exists, what it tests, how we know the right answer and why it matters, and a box "Test details" ends it.

**Terms used on this page.**

- **Regression check**: a test that compares the code with a stored value, a result the code gave before. It fails when the result changes, but it does not show that the result is right.
- **Cross-code check**: a test that compares CALLIOPE with a second, independent code.
- **atmodeller**: an outgassing code that PROTEUS can use in place of CALLIOPE.

## Water in peridotite melt

**Why this check exists.** `SolubilityH2O()` with no argument uses the peridotite law, the default of the package. A wrong constant, the constant of another melt composition, or a wrong exponent changes the dissolved water at every pressure.

**What it tests.** The law is a power law, $\mathrm{ppmw} = c\,p^{n}$, with $c$ = 524 ppmw bar$^{-0.5}$ and $n$ = 0.5. At $p$ = 100 bar the code must give 5240 ppmw, to a relative 1e-12, and at $p$ = 0 it must give 0. The test also requires that the result differs by more than 1000 ppmw from the basalt law of Dixon et al. (1995)[^cite-dixon1995] ($c$ = 965, 9650 ppmw) and from an exponent of 1 (52400 ppmw), and that it lies between 1e3 and 1e4 ppmw.

**How we know the right answer.** Sossi et al. (2023)[^cite-sossi2023] fit $c$ = 524 ppmw bar$^{-0.5}$ and $n$ = 0.5 to their experiments on peridotite melt (a published benchmark).

**Why it matters.** The solubility of water decides how much water stays in the magma ocean and how much goes into the atmosphere.

![Dissolved water against pressure on log axes, the code law and the Dixon law, with the value and the 2 wrong values at 100 bar](../assets/validation/solubility_h2o_peridotite.svg#only-light){ width="80%" }
![Dissolved water against pressure on log axes, the code law and the Dixon law, with the value and the 2 wrong values at 100 bar](../assets/validation/solubility_h2o_peridotite_dark.svg#only-dark){ width="80%" }

**How to read the plot.** Both axes are logarithmic. The solid line is the code law and the dashed line is the Dixon basalt law of the code; both give context. The test asserts the code value at 100 bar (circle) against $524 \times 100^{0.5}$ = 5240 ppmw (square), to a relative 1e-12, too small to see. The 2 crosses at 100 bar are the wrong values that the test holds: 9650 ppmw from the Dixon constant and 52400 ppmw from an exponent of 1. A pass means that the circle sits in the square and both crosses lie more than 1000 ppmw from it.

??? info "Test details"
    - Reference: Sossi et al. (2023)[^cite-sossi2023] peridotite H2O fit: `ppmw = 524 * p^0.5`.
    - Tolerance: relative 1e-12 at 100 bar; more than 1000 ppmw from the Dixon et al. (1995) basalt law and from an exponent of 1.
    - Test id: `tests/test_solubility.py::TestSolubilityH2O::test_peridotite_default_matches_sossi_2023_fit`
    - Command: `pytest tests/test_solubility.py -k peridotite_default_matches_sossi`

??? note "Worked numbers"
    At p = 100 bar:

    ```
    ppmw = 524 * sqrt(100) = 5240
    ```

    The Dixon 1995 basalt constant (`const = 965`) would give 9650, 84 % higher at the same pressure. An exponent of 1.0 would give 52400, 10 times higher.

## Sulfur at 2500 K

**Why this check exists.** The S$_2$ law depends on the temperature, the FeO content of the melt and the oxygen fugacity, which comes from the IW buffer. A wrong coefficient, an FeO content in the wrong unit, or a buffer other than the default changes the dissolved sulfur.

**What it tests.** `SolubilityS2()` with its defaults, at $p_\mathrm{S_2}$ = 1 bar, $T$ = 2500 K and $\Delta\mathrm{IW}$ = 0, must give 13248.67 ppmw to a relative 1e-5. The default FeO content is $x_\mathrm{FeO}$ = 10 wt%, and the default buffer is that of Hirschmann (2021)[^cite-hirschmann2021]. The test also requires that the result differs by more than 100 ppmw from the value with the Fischer et al. (2011)[^cite-fischer2011] buffer (13085.87 ppmw), and that it lies between 1e3 and 1e5 ppmw.

**How we know the right answer.** The value 13248.67 is the output of the code, so this case is a regression check. The value agrees with a hand evaluation of the law. Gaillard et al. (2022)[^cite-gaillard2022] give the closed form $\ln X_\mathrm{S} = 13.8426 - 26476/T + 0.124\,x_\mathrm{FeO} + 0.5 \ln(p_\mathrm{S_2}/f\mathrm{O_2})$, where $X_\mathrm{S}$ is the sulfur content of the melt in ppmw. With the Hirschmann buffer at 2500 K and 1 bar, $\log_{10} f\mathrm{O_2}$ = -4.34247, and the closed form gives 13249 ppmw.

**Why it matters.** The dissolved sulfur sets how much S$_2$, SO$_2$ and H$_2$S the atmosphere can hold. The value also depends on the default IW buffer: a change of the default or of its coefficients changes this number, and the test then fails here.

![Dissolved sulfur against temperature with the Hirschmann and the Fischer buffer, with the pinned value at 2500 K](../assets/validation/solubility_s2_2500k.svg#only-light){ width="80%" }
![Dissolved sulfur against temperature with the Hirschmann and the Fischer buffer, with the pinned value at 2500 K](../assets/validation/solubility_s2_2500k_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The lines are the code law with the Hirschmann buffer (solid) and with the Fischer buffer (dashed) and give context. The test asserts the marks at 2500 K: the code value (circle), the pinned value 13248.67 ppmw (square) and the Fischer value 13085.87 ppmw (cross). A pass means that the circle agrees with the square to a relative 1e-5, too small to see, and the cross lies more than 100 ppmw below.

**Known limit.** As a regression check, the case fails when a coefficient, the FeO unit or the buffer changes. It does not test the published law itself.

??? info "Test details"
    - Reference: Gaillard et al. (2022)[^cite-gaillard2022], EPSL 117255, S2 solubility law with `x_FeO = 10 wt%` Earth-mantle default ([doi:10.1016/j.epsl.2021.117255](https://doi.org/10.1016/j.epsl.2021.117255)); Hirschmann (2021)[^cite-hirschmann2021] IW buffer.
    - Tolerance: relative 1e-5; more than 100 ppmw from the Fischer-buffer value.
    - Test id: `tests/test_solubility.py::TestSolubilityS2_xFeO::test_default_call_matches_gaillard_2022_earth_mantle_value`
    - Command: `pytest tests/test_solubility.py -k gaillard_2022_earth_mantle`
    - Companion tests in `TestSolubilityS2_xFeO` check the closed-form scaling with $x_\mathrm{FeO}$ at 5, 10, 15 and 20 wt%, so the coefficient 0.124 is checked on its own.

??? note "Worked numbers"
    At T = 2500 K, p_S2 = 1 bar, x_FeO = 10 wt%, fO2 from the default Hirschmann 2021 IW buffer at 1 bar:

    ```
    log10(fO2) = -4.34247 (Hirschmann at T=2500 K, P=1 bar, dIW=0)
    fO2 = 10^-4.34247 = 4.545e-5 bar
    ln(p_S2/fO2) = ln(1 / 4.545e-5) = 9.9989
    ln(X) = 13.8426 - 10.5904 + 1.24 + 0.5 * 9.9989
          = 3.2522 + 1.24 + 4.9995
          = 9.4917
    ppmw = exp(9.4917) ~ 13249
    ```

    The test value 13248.67 matches this closed form. The Fischer 2011 buffer (log10(fO2) = -4.3317 at 2500 K) gives 13085.87, 1.2 % lower.

## Noble gas Henry constants

**Why this check exists.** A noble gas follows Henry's law: its dissolved amount is proportional to its partial pressure, and the factor is the Henry constant. CALLIOPE turns the Henry constants of Jambon et al. (1986)[^cite-jambon1986], given as a volume of gas at standard temperature and pressure (STP) per gram of melt, into ppmw per bar. A mistyped constant, a wrong molar mass, a constant copied from another gas, or a missing kg-to-g factor gives a wrong solubility for that gas.

**What it tests.** For each of He, Ne, Ar, Kr and Xe, `jambon86_ppmw_per_bar` must equal a literal value to a relative 1e-4: He 0.100065, Ne 0.225220, Ar 0.105220, Kr 0.112229 and Xe 0.099642 ppmw bar$^{-1}$. No gas may match the value of another gas, where the two literals differ by more than 0.001. Every constant must lie between 0.01 and 1, which catches a missing kg-to-g factor (1000 times too small).

**How we know the right answer.** Each literal is $(k_\mathrm{STP}/V_\mathrm{STP})\,M \times 10^6$, computed by hand from the constants of Jambon et al. (1986) and the molar masses. Here $k_\mathrm{STP}$ \[cm$^3$ STP g$^{-1}$ bar$^{-1}$\] is the published constant (for He the value 56e-5 of the abstract, which atmodeller also uses; Table 5 of the paper gives 56.5e-5), $V_\mathrm{STP}$ = 2.24e4 cm$^3$ mol$^{-1}$ is the molar volume of an ideal gas at STP, and $M$ \[g mol$^{-1}$\] is the molar mass. The test does not need atmodeller, so it runs on every CI image.

**Why it matters.** The dissolved fraction of each noble gas decides how much of it reaches the atmosphere.

![Henry constants of the 5 noble gases from the code, with the hand-computed literals](../assets/validation/solubility_noble_constants.svg#only-light){ width="80%" }
![Henry constants of the 5 noble gases from the code, with the hand-computed literals](../assets/validation/solubility_noble_constants_dark.svg#only-dark){ width="80%" }

**How to read the plot.** Each grey column is the constant of the code for one gas, and each mark is the literal that the test asserts. A pass means that every mark sits at the top of its column, to a relative 1e-4.

??? info "Test details"
    - Reference: Jambon, Weill and Braun (1986)[^cite-jambon1986], noble gas Henry constants in tholeiitic basalt ([doi:10.1016/0016-7037(86)90193-6](https://doi.org/10.1016/0016-7037(86)90193-6)); code in `jambon86_ppmw_per_bar` and `SolubilityNobleGas`.
    - Tolerance: relative 1e-4 per gas; no match with another gas; $0.01 < \mathrm{constant} < 1$ ppmw bar$^{-1}$.
    - Test id: `tests/test_solubility.py::TestSolubilityNobleGas::test_every_gas_constant_pinned_to_independent_literal`, `tests/test_solubility.py::TestSolubilityNobleGas::test_helium_constant_matches_jambon_1986`
    - Command: `pytest tests/test_solubility.py -k "every_gas_constant_pinned or helium_constant_matches"`
    - The helium test also requires that the He constant differs by more than 0.1 from the value with the molar mass of Ne (about 0.5045).

## Noble gas constants against atmodeller

**Why this check exists.** CALLIOPE and atmodeller both compute the Jambon et al. (1986) constants. If the two codes differ, a comparison of the two outgassing codes measures a difference in the constants, not in the physics.

**What it tests.** For each noble gas, the CALLIOPE constant must equal the atmodeller concentration at 1 bar, to a relative 1e-10. The largest constant must be more than 2 times the smallest, so a transcription that used one value for every gas fails.

**How we know the right answer.** atmodeller derives its constants from the same primitives, and its law is linear in the fugacity, so its concentration at 1 bar is the constant itself (a cross-code check).

**Why it matters.** The [cross-code comparison](../Explanations/cross_backend_comparison.md) of CALLIOPE and atmodeller relies on the 2 codes using the same constants.

![CALLIOPE Henry constants against the atmodeller concentrations at 1 bar for the 5 noble gases](../assets/validation/solubility_noble_atmodeller.svg#only-light){ width="80%" }
![CALLIOPE Henry constants against the atmodeller concentrations at 1 bar for the 5 noble gases](../assets/validation/solubility_noble_atmodeller_dark.svg#only-dark){ width="80%" }

**How to read the plot.** Both axes are logarithmic. Each mark is one gas: the atmodeller concentration at 1 bar on the x-axis and the CALLIOPE constant on the y-axis. A pass means that every mark lies on the dotted 1:1 line, to a relative 1e-10.

**Known limit.** The test is skipped when atmodeller is not installed, so it runs only where atmodeller is present, such as the nightly cross-code job.

??? info "Test details"
    - Reference: cross-code check against atmodeller.
    - Tolerance: relative 1e-10 per gas; the largest constant more than 2 times the smallest.
    - Test id: `tests/test_solubility.py::TestSolubilityNobleGas::test_noble_constants_cross_check_against_atmodeller`
    - Command: `pytest tests/test_solubility.py -k noble_constants_cross_check` (skips without atmodeller)

## Related code and pages

- `src/calliope/solubility.py` (`SolubilityH2O`): the H2O laws.
- `src/calliope/solubility.py` (`SolubilityS2`): the S2 law with the `x_FeO` argument.
- `src/calliope/solubility.py` (`SolubilityN2`): the N2 law with the `x_SiO2`, `x_Al2O3` and `x_TiO2` arguments (Dasgupta et al. 2022[^cite-dasgupta2022]), and Libourel et al. (2003)[^cite-libourel2003] as the law without composition.
- [Solubility](../Explanations/solubility.md): the explanation page, with the range of validity of each law.

## References

[^cite-dasgupta2022]: R. Dasgupta, E. Falksen, A. Pal, C. Sun, *[The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions](https://doi.org/10.1016/j.gca.2022.09.012)*, Geochimica et Cosmochimica Acta, 336, 291-307, 2022. [SciX](https://scixplorer.org/abs/2022GeCoA.336..291D/abstract).
[^cite-dixon1995]: J. E. Dixon, E. M. Stolper, J. R. Holloway, *[An experimental study of water and carbon dioxide solubilities in mid-ocean ridge basaltic liquids. Part I: Calibration and solubility models](https://doi.org/10.1093/oxfordjournals.petrology.a037267)*, Journal of Petrology, 36(6), 1607-1631, 1995. [SciX](https://scixplorer.org/abs/1995JPet...36.1607D/abstract).
[^cite-fischer2011]: R. A. Fischer, A. J. Campbell, G. A. Shofner, O. T. Lord, P. Dera, V. B. Prakapenka, *[Equation of state and phase diagram of FeO](https://doi.org/10.1016/j.epsl.2011.02.025)*, Earth and Planetary Science Letters, 304, 496-502, 2011. [SciX](https://scixplorer.org/abs/2011E%26PSL.304..496F/abstract).
[^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
[^cite-hirschmann2021]: M. M. Hirschmann, *[Iron-wüstite revisited: A revised calibration accounting for variable stoichiometry and the effects of pressure](https://doi.org/10.1016/j.gca.2021.08.039)*, Geochimica et Cosmochimica Acta, 313, 74-84, 2021.
[^cite-jambon1986]: A. Jambon, H. Weill, J. Braun, *[Solubility of He, Ne, Ar, Kr and Xe in a basalt melt in the range 1250-1600 C. Geochemical implications](https://doi.org/10.1016/0016-7037(86)90193-6)*, Geochimica et Cosmochimica Acta, 50(3), 401-408, 1986. [SciX](https://scixplorer.org/abs/1986GeCoA..50..401J/abstract).
[^cite-libourel2003]: G. Libourel, B. Marty, F. Humbert, *[Nitrogen solubility in basaltic melt. Part I. Effect of oxygen fugacity](https://doi.org/10.1016/S0016-7037(03)00259-X)*, Geochimica et Cosmochimica Acta, 67(21), 4123-4135, 2003. [SciX](https://scixplorer.org/abs/2003GeCoA..67.4123L/abstract).
[^cite-sossi2023]: P. A. Sossi, P. M. E. Tollan, J. Badro, D. J. Bower, *[Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets](https://doi.org/10.1016/j.epsl.2022.117894)*, Earth and Planetary Science Letters, 601, 117894, 2023. [SciX](https://scixplorer.org/abs/2023E%26PSL.60117894S/abstract).
