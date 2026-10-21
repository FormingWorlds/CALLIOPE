# Oxygen fugacity

CALLIOPE parameterises the redox state of the magma ocean through the oxygen fugacity $f_{\mathrm{O}_2}$ at the surface, expressed in $\log_{10}$ units relative to the iron-wüstite (IW) mineral buffer:

$$
\Delta\mathrm{IW} \equiv \log_{10} f_{\mathrm{O}_2} - \log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T).
$$

Under the buffered-mode entry point [`equilibrium_atmosphere`](mass_balance.md), the user supplies `fO2_shift_IW = ` $\Delta\mathrm{IW}$ as a scalar input and CALLIOPE computes the absolute $\log_{10} f_{\mathrm{O}_2}$ at $T = T_\mathrm{magma}$ by adding the buffer value to the shift. Under the [authoritative-oxygen mode](authoritative_oxygen.md), $\Delta\mathrm{IW}$ is instead a solver unknown that closes the system against a user-supplied total oxygen mass. Both modes share the parameterisations, conventions, and chemistry channels documented on this page; they differ only in whether $\Delta\mathrm{IW}$ is an input or an output.

## The IW mineral buffer

Iron in equilibrium with wüstite (FeO) in the presence of O$_2$ obeys

$$
\mathrm{Fe} + \tfrac{1}{2}\,\mathrm{O_2} \rightleftharpoons \mathrm{FeO},
$$

which fixes a single curve $\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T)$ in $T$-$f_{\mathrm{O}_2}$ space. CALLIOPE supports three parameterisations of this curve.

### Hirschmann (2021), `hirschmann` (default) [^cite-hirschmann2021]

An empirical fit to a thermodynamic model of iron coexisting with wüstite of equilibrium stoichiometry $\mathrm{Fe}_{1-y}\mathrm{O}$, valid from 1000 to 3000 K and from 1 bar to 100 GPa. Hirschmann (2021) [^cite-hirschmann2021] finds that the IW buffer computed from the JANAF properties of stoichiometric FeO is too reducing, by 1.1 log units at 1000 K and 0.2 log units at 3000 K. That comparison is against JANAF itself, not against the two-parameter `fischer` fit coded below: at 1 bar the Hirschmann fit sits $0.42$ dex above `fischer` at 1000 K and $0.11$ dex below it at 3000 K. The fit (their Table 1) is

$$
\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T, P) = a + b\,T + c\,T \ln T + \frac{d}{T},
\qquad m = m_0 + m_1 P + m_2 P^2 + m_3 P^3 + m_4 P^{1/2},
$$

for each parameter $m \in \{a, b, c, d\}$, with $P$ in GPa. Separate coefficient sets apply in the fcc/bcc and hcp iron fields, divided by $P_\mathrm{GPa} = -18.64 + 0.04359\,T - 5.069 \times 10^{-6}\,T^2$ (their Eq. 18). CALLIOPE evaluates this boundary at $T$ held to 1000-4300 K: the parabola falls below 1 bar under about 451 K and above about 8148 K, which would otherwise put the 1 bar buffer on the hcp branch. The boundary is unchanged between the calibration floor (1000 K) and the parabola's vertex (4300 K, 75 GPa).

Implemented as `OxygenFugacity.hirschmann(T, P_bar=1.0)`. The `OxygenFugacity('hirschmann')` dispatcher evaluates the buffer at the 1 bar reference, like the other two buffers, so that $\Delta\mathrm{IW}$ keeps the same meaning. To get the buffer at a given pressure, call the method directly with `P_bar`. Below 1000 K and above 3000 K the fit is an extrapolation that the paper advises against; magma-ocean runs hotter than 3000 K use the default buffer outside its calibrated range. Above 30 to 50 GPa, the effects of FeO metallisation and low-spin Fe$^{2+}$ are not included, so the buffer there should be treated as a metastable reference.

### Fischer et al. (2011), `fischer` [^cite-fischer2011]

A simpler two-parameter fit of the 1-bar IW buffer. Fischer et al. (2011) [^cite-fischer2011] calculate the 1-bar Fe-FeO curve of their Fig. 6 from Chase (1998) [^cite-chase1998] and give a P-T polynomial in their supplement (Table S2); the source of this 2-term form is not stated in the paper body, and no independent pin of it exists in the test data. Fischer's own high-pressure measurements ($\le$200 GPa) extend the buffer to deep-mantle conditions but are not used by CALLIOPE.

$$
\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T) = 6.94059 - \frac{28180.8}{T}.
$$

Implemented as `OxygenFugacity.fischer(T)`.

### O'Neill & Eggins (2002), `oneill` (legacy) [^cite-oneilleggins2002]

A thermochemically-constrained fit derived from low-temperature equilibrium data, expressed as Bower et al. (2022) [^cite-bower2022] Equation (7):

$$
\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T) = \frac{2\left[-244118 + 115.559\,T - 8.474\,T \ln T\right]}{\ln(10)\, R\, T},
$$

with $R = 8.31441$ J K$^{-1}$ mol$^{-1}$. Bower et al. (2022) [^cite-bower2022] adopted this as the "IW buffer to which $f_{\mathrm{O}_2}$ is referenced". CALLIOPE retains it under the model name `oneill` (function `OxygenFugacity.oneill(T)` in `oxygen_fugacity.py`); keep it as the choice when you need to reproduce results from the older literature line.

!!! note "Choice of buffer"
    CALLIOPE defaults to Hirschmann (2021) because it accounts for the variable stoichiometry of wüstite, which the JANAF-based stoichiometric-FeO fits do not, and because it is the calibration that atmodeller's Hirschmann composite uses above 1000 K. At 1 bar the three buffers differ as follows: Hirschmann is $0.23$ dex more oxidising than Fischer at $T = 1500$ K, $0.10$ dex at 2000 K, within $0.01$ dex near 2500 K, and $0.11$ dex more reducing at 3000 K. Against O'Neill, Hirschmann is $0.22$ dex more reducing at 1500 K, $0.35$ dex more oxidising at 2000 K, and $0.95$ dex more oxidising at 3000 K. The choice therefore matters at the few-tenths-of-a-dex level for inferred partial pressures across most of the magma-ocean range. Fischer and O'Neill are available as `OxygenFugacity('fischer')` and `OxygenFugacity('oneill')`, or as `fO2_model='fischer'` and `fO2_model='oneill'` on a solve. See [Backend comparison](cross_backend_comparison.md) for the cross-backend comparison.

## Selecting the buffer

Both solver entry points take the buffer as an argument:

```python
out = equilibrium_atmosphere(target, ddict, fO2_model='fischer')
```

The same choice can be made with the `ddict['fO2_model']` key, which is also what the lower-level helpers (`get_partial_pressures`, `atmosphere_mass`, `dissolved_mass`) read. The buffer applies to the whole solve: the O$_2$ partial pressure, every modified equilibrium constant and the S$_2$ solubility use it, so a solve never mixes two buffers. The Dasgupta N$_2$ law is the exception: it takes the shift itself as its $\Delta\mathrm{IW}$ without converting it (see [Solubility](solubility.md)), so at a fixed absolute $f_{\mathrm{O}_2}$ the dissolved N depends on the buffer chosen. Without either, the solve uses the default `DEFAULT_FO2_MODEL = 'hirschmann'`. An unknown name, or an argument that contradicts `ddict['fO2_model']`, raises `ValueError` from the solver entry points; a `None` argument or `ddict['fO2_model']` there selects the default. The helpers raise `ValueError` for an unknown name or a `None` in `ddict['fO2_model']`.

## How $\Delta\mathrm{IW}$ enters the chemistry

The shift $\Delta\mathrm{IW}$ feeds the equilibrium chemistry through four channels:

### 1. The free $\mathrm{O_2}$ partial pressure

In `solve.get_partial_pressures`:

```python
p_d['O2'] = 10.0 ** OxygenFugacity(buffer)(ddict['T_magma'], fO2_shift)
```

where `buffer` is the IW buffer of the solve (see [Selecting the buffer](#selecting-the-buffer)).

This is what makes O$_2$ a *derived*, not a solved, quantity. The mantle redox state pins $p_\mathrm{O_2}$ at the surface, and through the speciation tree it pins the equilibrium between every oxidised/reduced couple.

### 2. The modified equilibrium constants

For every reaction $A \to B + n_\mathrm{O_2}\,\mathrm{O_2}$, the equilibrium ratio $G_\mathrm{eq} = p_B / p_A$ depends on $f_{\mathrm{O}_2}$ as

$$
G_\mathrm{eq}(T, f_{\mathrm{O}_2}) = 10^{\,\log_{10} K_\mathrm{eq}(T) - n_\mathrm{O_2}\,\log_{10} f_{\mathrm{O}_2}}.
$$

For the H$_2$O-H$_2$ couple ($n_\mathrm{O_2} = +0.5$), reducing conditions ($f_{\mathrm{O}_2}$ smaller, $\Delta\mathrm{IW}$ more negative) drive $G_\mathrm{eq}$ larger, which in turn drives more H$_2$O to dissociate into H$_2$. This is the redox dependence visible in the redox-sweep tutorial, and it is the mechanism behind the H$_2$-dominated, long-lived magma-ocean atmospheres found in Nicholls et al. (2024) [^cite-nicholls2024] Figure 6 at $\Delta\mathrm{IW} \le -1$.

### 3. The S$_2$ Gaillard solubility

The Gaillard et al. (2022) [^cite-gaillard2022] sulfide solubility carries an explicit $\ln f_{\mathrm{O}_2}$ term, so $\Delta\mathrm{IW}$ enters the dissolved-S inventory directly. The implementation in `solubility.SolubilityS2.gaillard` calls back into `OxygenFugacity` with the buffer of the solve to compute the absolute $f_{\mathrm{O}_2}$.

### 4. The N$_2$ Dasgupta solubility

Similarly, the Dasgupta et al. (2022) [^cite-dasgupta2022] N$_2$ solubility includes a $-1.6\,\Delta\mathrm{IW}$ term in its exponent, so reducing conditions sharply increase the dissolved-N inventory. `solubility.SolubilityN2.dasgupta` takes the solve's shift as this $\Delta\mathrm{IW}$ directly, on whichever buffer the solve uses. Dasgupta et al. (2022) define $\Delta\mathrm{IW}$ by the Fe–FeO activity ratio of the melt (their Eq. 2) and name no IW buffer curve, so there is no buffer to convert to. This is the one channel that the `fO2_model` setting does not reach: at a fixed absolute $f_{\mathrm{O}_2}$, switching buffers moves the shift and with it the dissolved N. This is one mechanism by which planet-scale N partitioning is tied to mantle redox; see Nicholls et al. (2026) [^cite-nicholls2026] for an application to L 98-59 d, where the inferred H$_2$-dominated atmosphere with photochemical SO$_2$ implies $\Delta\mathrm{IW}$ between IW-4 and IW-1.

## Reference values for $\Delta\mathrm{IW}$

| Reservoir | $\Delta\mathrm{IW}$ | Source |
|---|---|---|
| Mercury surface | IW-2.8 to IW-5.4 (Fe-based: IW-2.8 to IW-4.5; sulphur-based: IW-5.4 via Namur et al. 2016) | Cartier & Wood (2019) [^cite-cartierwood2019] |
| Mars upper mantle (shergottite source) | $\approx$ IW (specifically IW-1.0 to IW-0.3 for QUE 94201) | Wadhwa (2001) [^cite-wadhwa2001] |
| Mars shergottite parent melts | IW-1.0 to IW+1.9 (variation from crust assimilation) | Wadhwa (2001) [^cite-wadhwa2001] |
| Iron-wüstite buffer | $0$ | by definition |
| Earth's upper mantle (modern) | $\approx$ IW+3.5 | Sossi et al. (2020) [^cite-sossi2020] |
| Earth upper mantle (range) | FMQ$\,\pm\,2$ ($\approx$ IW+1.5 to IW+5.5) | Frost & McCammon (2008) [^cite-frostmccammon2008] |
| Earth mantle at $\sim 8$ GPa | $\approx$ FMQ$-5$ ($\approx$ IW-1.5) | Frost & McCammon (2008) [^cite-frostmccammon2008] |
| Earth transition zone ($\sim$14-23 GPa) | just below IW | Frost & McCammon (2008) [^cite-frostmccammon2008] |
| Earth lower mantle ($>$23 GPa) | metal-saturated ($\sim$1 wt% Fe$^0$); at or below IW | Frost & McCammon (2008) [^cite-frostmccammon2008] |

CALLIOPE's PROTEUS-side default is `fO2_shift_IW = 4.0`, consistent with a near-surface terrestrial composition. Nicholls et al. (2024) [^cite-nicholls2024] Table 2 explored $\Delta\mathrm{IW} \in \{-5, -3, -1, 0, +1, +3, +5\}$ on a 7-point grid and demonstrated that the resulting atmospheric composition spans the full range from H$_2$-dominated reduced atmospheres (TRAPPIST-1 c-like) to H$_2$O/CO$_2$-dominated oxidised atmospheres (Earth-like).

## Limitations

- **No $f_{\mathrm{O}_2}$ evolution**: $\Delta\mathrm{IW}$ is a constant input, not a state variable. In reality the mantle $f_{\mathrm{O}_2}$ should evolve with degree of crystallisation, fractional crystallisation depth, and atmospheric escape; CALLIOPE does not capture this and the user is responsible for choosing a representative value or sweeping over a grid.
- **No $f_{\mathrm{O}_2}$ depth profile**: the surface $f_{\mathrm{O}_2}$ alone enters the chemistry. Bower et al. (2022) [^cite-bower2022] §2.3 and Sossi et al. (2020) [^cite-sossi2020] discuss why the *interface* fugacity (rather than the deep-mantle value) is the relevant choice; this assumption is consistent with CALLIOPE's ideal-gas, single-temperature treatment but breaks down if a Fe-FeO equilibrium curve in the deep mantle differs by more than a few dex.
- **No solid-FeO buffering**: when $\Phi_\mathrm{global} \to 0$, there is no melt to buffer $f_{\mathrm{O}_2}$ against the user-prescribed value. CALLIOPE keeps using the shift regardless of melt fraction, which is a reasonable bookkeeping choice but should not be over-interpreted physically.

## See also

- [Authoritative-oxygen mode](authoritative_oxygen.md): how $\Delta\mathrm{IW}$ is recovered as a solver unknown when O is supplied as a budget instead.
- [Equilibrium chemistry](equilibrium_chemistry.md): the species-by-species speciation tree
- [Solubility laws](solubility.md): where $f_{\mathrm{O}_2}$ enters the S and N solubility paths
- [API reference for `calliope.oxygen_fugacity`](../Reference/api/calliope.oxygen_fugacity.md)

 [^cite-bower2022]: D. J. Bower, K. Hakim, P. A. Sossi, P. Sanan, *[Retention of water in terrestrial magma oceans and carbon-rich early atmospheres](https://doi.org/10.3847/PSJ/ac5fb1)*, The Planetary Science Journal, 3(4), 93, 2022. [SciX](https://scixplorer.org/abs/2022PSJ.....3...93B/abstract).
 [^cite-cartierwood2019]: C. Cartier, B. J. Wood, *[The role of reducing conditions in building Mercury](https://doi.org/10.2138/gselements.15.1.39)*, Elements, 15(1), 39–45, 2019. [SciX](https://scixplorer.org/abs/2019Eleme..15...39C/abstract).
 [^cite-chase1998]: M. W. Chase, *[NIST-JANAF Thermochemical Tables, 4th edition](https://janaf.nist.gov/)*, Journal of Physical and Chemical Reference Data Monograph 9, 1998.
 [^cite-dasgupta2022]: R. Dasgupta, E. Falksen, A. Pal, C. Sun, *[The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions](https://doi.org/10.1016/j.gca.2022.09.012)*, Geochimica et Cosmochimica Acta, 336, 291–307, 2022. [SciX](https://scixplorer.org/abs/2022GeCoA.336..291D/abstract).
 [^cite-fischer2011]: R. A. Fischer, A. J. Campbell, G. A. Shofner, O. T. Lord, P. Dera, V. B. Prakapenka, *[Equation of state and phase diagram of FeO](https://doi.org/10.1016/j.epsl.2011.02.025)*, Earth and Planetary Science Letters, 304, 496–502, 2011. [SciX](https://scixplorer.org/abs/2011E%26PSL.304..496F/abstract).
 [^cite-frostmccammon2008]: D. J. Frost, C. A. McCammon, *[The redox state of Earth's mantle](https://doi.org/10.1146/annurev.earth.36.031207.124322)*, Annual Review of Earth and Planetary Sciences, 36, 389–420, 2008. [SciX](https://scixplorer.org/abs/2008AREPS..36..389F/abstract).
 [^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
 [^cite-hirschmann2021]: M. M. Hirschmann, *[Iron-wüstite revisited: A revised calibration accounting for variable stoichiometry and the effects of pressure](https://doi.org/10.1016/j.gca.2021.08.039)*, Geochimica et Cosmochimica Acta, 313, 74–84, 2021.
 [^cite-nicholls2024]: H. Nicholls, T. Lichtenberg, D. J. Bower, R. Pierrehumbert, *[Magma ocean evolution at arbitrary redox state](https://doi.org/10.1029/2024JE008576)*, Journal of Geophysical Research: Planets, 129, e2024JE008576, 2024. [SciX](https://scixplorer.org/abs/2024JGRE..12908576N/abstract).
 [^cite-nicholls2026]: H. Nicholls, T. Lichtenberg, R. D. Chatterjee, C. M. Guimond, E. Postolec, R. T. Pierrehumbert, *[Volatile-rich evolution of molten super-Earth L 98-59 d](https://doi.org/10.1038/s41550-026-02815-8)*, Nature Astronomy, 2026. [SciX](https://scixplorer.org/abs/2026NatAs.tmp...61N/abstract). [arXiv](https://arxiv.org/abs/2507.02656).
 [^cite-oneilleggins2002]: H. St. C. O'Neill, S. M. Eggins, *[The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts](https://doi.org/10.1016/S0009-2541(01)00414-4)*, Chemical Geology, 186, 151–181, 2002. [SciX](https://scixplorer.org/abs/2002ChGeo.186..151O/abstract).
 [^cite-sossi2020]: P. A. Sossi, A. D. Burnham, J. Badro, A. Lanzirotti, M. Newville, H. St. C. O'Neill, *[Redox state of Earth's magma ocean and its Venus-like early atmosphere](https://doi.org/10.1126/sciadv.abd1387)*, Science Advances, 6, eabd1387, 2020. [SciX](https://scixplorer.org/abs/2020SciA....6.1387S/abstract).
 [^cite-wadhwa2001]: M. Wadhwa, *[Redox state of Mars' upper mantle and crust from Eu anomalies in shergottite pyroxenes](https://doi.org/10.1126/science.1057594)*, Science, 291, 1527–1530, 2001. [SciX](https://scixplorer.org/abs/2001Sci...291.1527W/abstract).
