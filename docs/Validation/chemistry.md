# Validation: `src/calliope/chemistry.py`

## What this page checks

This page checks the equilibrium constants that CALLIOPE uses to split the gas between reduced and oxidised species.
An equilibrium constant $K$ fixes the ratio of the partial pressures of the products and reactants of a reaction at a given temperature.
CALLIOPE folds the oxygen fugacity into it: `ModifiedKeq` returns $G_\mathrm{eq} = 10^{\log_{10} K - \nu \log_{10} f\mathrm{O_2}}$, and for example the H$_2$ pressure of the gas is $G_\mathrm{eq}$ of `janaf_H2` times the H$_2$O pressure (`src/calliope/solve.py`).
Validation means that a test runs the code on a case where another source gives the right answer, and checks the result within a stated tolerance.
If a constant is wrong, the gas has the wrong ratio of H$_2$ to H$_2$O, CO to CO$_2$, or H$_2$S and NH$_3$ to their parents, and the atmosphere that PROTEUS passes to the climate model is wrong.
Each case below says why the check exists, what it tests, how we know the right answer and why it matters, and a box "Test details" ends it.

**Terms used on this page.**

- **Oxygen fugacity** $f\mathrm{O_2}$ \[bar\]: a measure of how oxidising the melt is, in the form of an effective O$_2$ pressure.
- **IW buffer**: the $f\mathrm{O_2}$ at which iron metal and wüstite (FeO) are in equilibrium. $\Delta\mathrm{IW}$ is the shift of $f\mathrm{O_2}$ from it in log10 units. The [oxygen fugacity page](oxygen_fugacity.md) checks the buffers.
- **$\nu$**: the O$_2$ factor that a reaction method returns: 0.5 for `janaf_H2`, `janaf_CO` and `schaefer_H`, 1 for `schaefer_CH4`, and 0 for `janaf_SO2`, whose O$_2$ the solve adds through the O$_2$ pressure, and for `janaf_H2S` and `janaf_NH3`, which contain no O$_2$.
- **dex**: a difference in log10 units; 0.025 dex is a factor of $10^{0.025} \approx 1.06$.
- **JANAF**: the NIST-JANAF Thermochemical Tables[^cite-chase1998], tables of formation data of gases and solids against temperature.
- **$K_f$**: the formation constant of a species, from its elements in their reference states.

## H2 to H2O ratio at 2000 K

**Why this check exists.** `ModifiedKeq` picks a reaction by name and combines its constant with the oxygen fugacity of a buffer. A typo in the `janaf_H2` coefficients, a wrong sign of the $f\mathrm{O_2}$ term, a wrong stoichiometric factor, a wrong buffer, or a call to the wrong reaction changes $G_\mathrm{eq}$ without an error message.

**What it tests.** The test evaluates $G_\mathrm{eq}$ of `janaf_H2`, for the reaction H$_2$O = H$_2$ + 0.5 O$_2$, at $T$ = 2000 K and $\Delta\mathrm{IW}$ = 0. $\Delta\mathrm{IW}$ is the shift of the oxygen fugacity from the iron-wüstite (IW) buffer in log10 units; the test uses the buffer of O'Neill and Eggins (2002)[^cite-oneilleggins2002]. The code gives 1.4664. The test computes the same value in closed form, $10^{\log_{10} K - 0.5 \log_{10} f\mathrm{O_2}}$, from the coefficients and the O'Neill and Eggins buffer, and requires agreement to a relative 1e-6. It also requires that `schaefer_H`, the constant of Schaefer and Fegley (2017)[^cite-schaeferfegley2017] for the same reaction, differs by more than 0.2 (it gives 1.2126), and that $G_\mathrm{eq}$ lies between 0.5 and 5.

**How we know the right answer.** The closed form uses the published O'Neill and Eggins buffer, which the test writes out itself, and the `janaf_H2` coefficients, which are a fit to the NIST-JANAF tables[^cite-chase1998] for this reaction over 1500 to 3000 K. The test does not compare the coefficients with the JANAF tables; it checks how the code combines them with the buffer.

**Why it matters.** The H$_2$/H$_2$O ratio of the outgassed gas comes from this constant, and the H$_2$ pressure also enters CH$_4$, NH$_3$ and H$_2$S in the solve.

??? info "Test details"
    - Reference: NIST-JANAF Thermochemical Tables[^cite-chase1998] (4th ed.), fits for `H2O = H2 + 0.5 O2` over 1500-3000 K; O'Neill and Eggins (2002) IW buffer (closed form).
    - Tolerance: relative 1e-6; `schaefer_H` must differ by more than 0.2; $0.5 < G_\mathrm{eq} < 5$.
    - Test id: `tests/test_chemistry.py::test_modified_keq_janaf_H2_matches_closed_form_at_2000K_with_oneill`
    - Command: `pytest tests/test_chemistry.py -k janaf_H2_matches_closed_form`
    - Companion test: `test_modified_keq_fO2_model_choice_changes_result` checks that a switch from the O'Neill to the Fischer buffer changes $G_\mathrm{eq}$ by more than a relative 1e-6.

??? note "Worked numbers"
    `ModifiedKeq` returns `Geq = 10^(Keq(T) - fO2_stoich * log10_fO2(T, dIW))`.
    For `janaf_H2` the coefficients are `(a = -13152.48, b = 3.038586, stoich = 0.5)`.

    At T = 2000 K, `fO2_model='oneill'`, `fO2_shift = 0`:

    ```
    Keq = 10^(-13152.48 / 2000 + 3.038586)
        = 10^(-6.5762 + 3.0386)
        = 10^-3.5376
        ~ 2.900e-4
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

    The test builds the closed form without intermediate rounding. `schaefer_H` at the same T:

    ```
    Keq_schaefer = 10^(-12794 / 2000 + 2.7768)
                 = 10^-3.6202 ~ 2.398e-4
    Geq_schaefer = 10^(-3.6202 + 3.7039) = 10^0.0837 ~ 1.213
    ```

    The difference of 0.25 between the two is above the 0.2 of the test.

## H2S constant against the JANAF tables

**Why this check exists.** `janaf_H2S` gives $\log_{10} K$ of S$_2$ + 2 H$_2$ = 2 H$_2$S as $a/T - b$, where $a$ \[K\] and $b$ are fit coefficients. A wrong coefficient, a wrong stoichiometric coefficient of S$_2$ (1 per mole of H$_2$S in place of 1/2), or a fit over the wrong temperature range gives the wrong H$_2$S pressure.

**What it tests.** Per mole of H$_2$S, that is half of $a/T - b$, the code's $\log_{10} K$ must lie within 0.025 dex of the JANAF value at each of the 44 table rows from 700 to 5000 K, the 2 ends included. The JANAF value per mole is $\log K_f(\mathrm{H_2S}) - 0.5 \log K_f(\mathrm{S_2})$, where $K_f$ is the formation constant of a species. The largest miss of the code is 0.018 dex. To show that the test can see a wrong fit, it also checks a set of wrong coefficients, those with S$_2$ coefficient 1: they must miss the rows by more than 0.075 dex, and they miss by up to 1.853 dex. Last, $\log_{10} K$ must fall with $T$, as for an exothermic reaction.

**How we know the right answer.** The JANAF rows of H$_2$S and S$_2$ are in `tools/janaf_H2S.dat` and `tools/janaf_S2.dat`, taken from the NIST-JANAF tables[^cite-chase1998]. Both tables use the same sulfur reference state, which cancels in the difference.

**Why it matters.** The solve sets the H$_2$S pressure to $\sqrt{G_\mathrm{eq}\, p_\mathrm{S_2}\, p_\mathrm{H_2}^2}$, so a wrong constant changes the sulfur in a hydrogen-rich gas.

??? info "Test details"
    - Reference: NIST-JANAF Thermochemical Tables[^cite-chase1998] (4th ed.), log Kf of H2S and S2 in `tools/`.
    - Tolerance: 0.025 dex at every row from 700 to 5000 K; the discriminating coefficients must miss by more than 0.075 dex (3 times the tolerance). Last compared 2026-10-20.
    - Test id: `tests/test_chemistry.py::test_janaf_constant_matches_janaf_rows_over_fit_range[hydrogen_sulfide]`
    - Command: `pytest "tests/test_chemistry.py::test_janaf_constant_matches_janaf_rows_over_fit_range[hydrogen_sulfide]"`

## H2S coefficients are a least-squares fit

**Why this check exists.** The case above allows a miss of 0.025 dex, so coefficients from a different fit could still pass it. This case pins the coefficients to the one fit they come from.

**What it tests.** The test fits a straight line of $2 \log_{10} K$ against $1/T$ to the JANAF rows from 700 to 5000 K by least squares. At every row, this line must equal the code's $a/T - b$ to a relative and absolute 1e-6. The absolute part is needed because the H$_2$S value crosses zero near 1830 K.

**How we know the right answer.** The fit of `tools/fit_H2S.ipynb` over the same JANAF rows gives the coefficients; the test repeats it with `numpy.polyfit`.

**Why it matters.** A changed fit range, form or coefficient moves the fit away from the code, and this case then fails even when the change stays inside the 0.025 dex of the row test.

??? info "Test details"
    - Reference: the JANAF rows of the case above; `tools/fit_H2S.ipynb`.
    - Tolerance: relative 1e-6 and absolute 1e-6 at every row. Last compared 2026-10-20.
    - Test id: `tests/test_chemistry.py::test_janaf_coefficients_are_least_squares_fit_of_janaf_rows[hydrogen_sulfide]`
    - Command: `pytest "tests/test_chemistry.py::test_janaf_coefficients_are_least_squares_fit_of_janaf_rows[hydrogen_sulfide]"`

## NH3 constant against the JANAF tables

**Why this check exists.** `janaf_NH3` gives $\log_{10} K$ of N$_2$ + 3 H$_2$ = 2 NH$_3$ as $a/T - b$, with fit coefficients $a$ and $b$. A wrong coefficient or fit range gives the wrong NH$_3$ pressure.

**What it tests.** Per mole of NH$_3$, the code's $\log_{10} K$ must lie within 0.04 dex of $\log K_f(\mathrm{NH_3})$ at every JANAF row from 700 to 5000 K (44 rows). The largest miss of the code is 0.034 dex. A set of wrong coefficients must miss the rows by more than 0.12 dex, which shows that the test can see a wrong fit; they miss by up to 0.158 dex. $\log_{10} K$ must fall with $T$.

**How we know the right answer.** The JANAF rows of NH$_3$ are in `tools/janaf_NH3.dat`, from the NIST-JANAF tables[^cite-chase1998]. N$_2$ and H$_2$ are reference states, so $\log K_f(\mathrm{NH_3})$ is the constant per mole of NH$_3$.

**Why it matters.** The solve sets the NH$_3$ pressure to $\sqrt{G_\mathrm{eq}\, p_\mathrm{N_2}\, p_\mathrm{H_2}^3}$, so a wrong constant changes the nitrogen in a hydrogen-rich gas.

??? info "Test details"
    - Reference: NIST-JANAF Thermochemical Tables[^cite-chase1998] (4th ed.), log Kf of NH3 in `tools/`.
    - Tolerance: 0.04 dex at every row from 700 to 5000 K; the discriminating coefficients must miss by more than 0.12 dex (3 times the tolerance). Last compared 2026-10-20.
    - Test id: `tests/test_chemistry.py::test_janaf_constant_matches_janaf_rows_over_fit_range[ammonia]`
    - Command: `pytest "tests/test_chemistry.py::test_janaf_constant_matches_janaf_rows_over_fit_range[ammonia]"`

## NH3 coefficients are a least-squares fit

**Why this check exists.** As for H$_2$S, the 0.04 dex of the row test leaves room for coefficients from a different fit.

**What it tests.** The test fits a straight line of $2 \log_{10} K$ against $1/T$ to the JANAF rows from 700 to 5000 K by least squares. At every row, this line must equal the code's $a/T - b$ to a relative and absolute 1e-6.

**How we know the right answer.** The fit of `tools/fit_NH3.ipynb` over the same rows gives the coefficients; the test repeats it with `numpy.polyfit`.

**Why it matters.** A changed fit range, form or coefficient moves the fit away from the code, and this case then fails.

??? info "Test details"
    - Reference: the JANAF rows of the case above; `tools/fit_NH3.ipynb`.
    - Tolerance: relative 1e-6 and absolute 1e-6 at every row. Last compared 2026-10-20.
    - Test id: `tests/test_chemistry.py::test_janaf_coefficients_are_least_squares_fit_of_janaf_rows[ammonia]`
    - Command: `pytest "tests/test_chemistry.py::test_janaf_coefficients_are_least_squares_fit_of_janaf_rows[ammonia]"`

## All 6 constants against NASA Glenn data

**Why this check exists.** The cases above check `janaf_H2`, `janaf_H2S` and `janaf_NH3` one by one. `janaf_CO`, `janaf_SO2` and `schaefer_CH4` need a check too, and the solve must use each constant in the right mass-action form.

**What it tests.** At 1500, 2000, 2500 and 3000 K the test checks, for `janaf_H2`, `janaf_CO`, `schaefer_CH4`, `janaf_SO2`, `janaf_NH3` and `janaf_H2S` with the Fischer buffer:

- the O$_2$ factor $\nu$ of each reaction method, exactly;
- $G_\mathrm{eq}$ against $10^{\log_{10} K - \nu \log_{10} f\mathrm{O_2}}$, and the ratio of $G_\mathrm{eq}$ at $\Delta\mathrm{IW} = +2$ and $-2$ against $10^{-4\nu}$, to a relative 1e-12;
- after a gas-only solve, each reaction quotient (for example $p_\mathrm{H_2}\, p_\mathrm{O_2}^{0.5} / p_\mathrm{H_2O}$) against the code's $K$, to a relative 1e-12;
- the offset of each quotient from the NASA Glenn value of $\log_{10} K$, to within 0.002 dex of its stored offset.

**How we know the right answer.** The NASA Glenn values are 9-coefficient fits of thermodynamic data, taken from atmodeller 1.0.2; the 72-case test on the [solve page](solve.md#gas-only-grid-against-atmodeller) checks them against atmodeller's own polynomials. The 2 data sets differ by construction, so the stored offsets are measured values (a regression snapshot): the largest is that of CH$_4$, 0.100 dex at 2000 K. The mass-action checks are exact identities of the solve.

**Why it matters.** This case covers `janaf_CO`, `janaf_SO2` and `schaefer_CH4`, which no other case on this page checks, and it fixes how far CALLIOPE's chemistry lies from the data that atmodeller uses.

**Known limit.** `schaefer_C` has no case, and `schaefer_H` appears only as the wrong reaction of the first case. This case runs only in the nightly cross-code job, because it needs atmodeller.

??? info "Test details"
    - Reference: NASA Glenn 9-coefficient thermodynamic data (cross-code check, values from atmodeller 1.0.2).
    - Tolerance: exact O$_2$ factor; relative 1e-12 for $G_\mathrm{eq}$, the shift ratio and each reaction quotient; offset from NASA Glenn within 0.002 dex of the stored offset.
    - Test id: `tests/test_cross_code_consistency.py::test_calliope_equilibrium_constants_and_reaction_quotients_vs_nasa_glenn`
    - Command: `pytest tests/test_cross_code_consistency.py -k nasa_glenn` (integration tier; skips without atmodeller)

## Related code and pages

- `src/calliope/chemistry.py` (`ModifiedKeq.__init__`, `ModifiedKeq.__call__`): constructor and modified constant.
- `src/calliope/chemistry.py` (`ModifiedKeq.schaefer_*`, `ModifiedKeq.janaf_*`): per-reaction coefficient methods.
- [Equilibrium chemistry](../Explanations/equilibrium_chemistry.md): the explanation page.

## References

 [^cite-chase1998]: M. W. Chase, *[NIST-JANAF Thermochemical Tables, 4th edition](https://janaf.nist.gov/)*, Journal of Physical and Chemical Reference Data Monograph 9, 1998.
 [^cite-oneilleggins2002]: H. St. C. O'Neill, S. M. Eggins, *[The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts](https://doi.org/10.1016/S0009-2541(01)00414-4)*, Chemical Geology, 186, 151-181, 2002. [SciX](https://scixplorer.org/abs/2002ChGeo.186..151O/abstract). The IW expression of the code uses the free energy of Fe + 0.5 O$_2$ = FeO(liquid) in their Table 6, p. 163.
 [^cite-schaeferfegley2017]: L. Schaefer, B. Fegley, *[Redox states of initial atmospheres outgassed on rocky planets and planetesimals](https://doi.org/10.3847/1538-4357/aa784f)*, The Astrophysical Journal, 843(2), 120, 2017. [SciX](https://scixplorer.org/abs/2017ApJ...843..120S/abstract).
