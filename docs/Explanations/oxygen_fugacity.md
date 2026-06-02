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

which fixes a single curve $\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T)$ in $T$-$f_{\mathrm{O}_2}$ space. CALLIOPE supports two parameterisations of this curve.

### Fischer et al. (2011), `fischer` (default) [^cite-fischer2011]

A simpler two-parameter fit of the 1-bar IW buffer. The fit reproduces the 1-bar curve in Fischer et al. (2011) [^cite-fischer2011] Fig. 6, which itself derives from Chase (1998) [^cite-chase1998] NIST-JANAF tabulation. Fischer's own high-pressure measurements ($\le$200 GPa) extend the buffer to deep-mantle conditions but are not used by CALLIOPE.

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
    The two parameterisations cross near $T \approx 1800$ K and diverge in opposite directions on either side: at $T = 1500$ K Fischer is about $0.4$ dex *more reducing* than O'Neill; at $T = 3000$ K Fischer is about $1.1$ dex *more oxidising*. The crossover means the difference is small (under $0.05$ dex) near 1800 K but grows to several tenths of a dex by 2400 K and reaches roughly $1$ dex at 3000 K. The choice matters at the few-tenths-of-a-dex level for inferred partial pressures across most of the magma-ocean range. CALLIOPE now defaults to Fischer 2011 because it sits within $\sim$0.2 dex of the Hirschmann composite used by atmodeller across the whole magma-ocean range and so produces cross-backend $\Delta\mathrm{IW}$ values that agree to a few tenths of a dex rather than up to $\sim 1$ dex with the older default. The legacy O'Neill choice remains available for reproducibility of pre-existing results; see [Backend comparison](cross_backend_comparison.md) for the quantitative comparison.

## How $\Delta\mathrm{IW}$ enters the chemistry

The shift $\Delta\mathrm{IW}$ feeds the equilibrium chemistry through four channels:

### 1. The free $\mathrm{O_2}$ partial pressure

In `solve.get_partial_pressures`:

```python
fO2_model = OxygenFugacity()
p_d['O2'] = 10.0 ** fO2_model(ddict['T_magma'], fO2_shift)
```

This is what makes O$_2$ a *derived*, not a solved, quantity. The mantle redox state pins $p_\mathrm{O_2}$ at the surface, and through the speciation tree it pins the equilibrium between every oxidised/reduced couple.

### 2. The modified equilibrium constants

For every reaction $A \to B + n_\mathrm{O_2}\,\mathrm{O_2}$, the equilibrium ratio $G_\mathrm{eq} = p_B / p_A$ depends on $f_{\mathrm{O}_2}$ as

$$
G_\mathrm{eq}(T, f_{\mathrm{O}_2}) = 10^{\,\log_{10} K_\mathrm{eq}(T) - n_\mathrm{O_2}\,\log_{10} f_{\mathrm{O}_2}}.
$$

For the H$_2$O-H$_2$ couple ($n_\mathrm{O_2} = +0.5$), reducing conditions ($f_{\mathrm{O}_2}$ smaller, $\Delta\mathrm{IW}$ more negative) drive $G_\mathrm{eq}$ larger, which in turn drives more H$_2$O to dissociate into H$_2$. This is the redox dependence visible in the redox-sweep tutorial, and it is the mechanism behind the H$_2$-dominated, long-lived magma-ocean atmospheres found in Nicholls et al. (2024) [^cite-nicholls2024] Figure 6 at $\Delta\mathrm{IW} \le -1$.

### 3. The S$_2$ Gaillard solubility

The Gaillard et al. (2022) [^cite-gaillard2022] sulfide solubility carries an explicit $\ln f_{\mathrm{O}_2}$ term, so $\Delta\mathrm{IW}$ enters the dissolved-S inventory directly. The implementation in `solubility.SolubilityS2.gaillard` calls back into `OxygenFugacity()` to compute the absolute $f_{\mathrm{O}_2}$.

### 4. The N$_2$ Dasgupta solubility

Similarly, the Dasgupta et al. (2022) [^cite-dasgupta2022] N$_2$ solubility includes a $-1.6\,\Delta\mathrm{IW}$ term in its exponent, so reducing conditions sharply increase the dissolved-N inventory. This is one mechanism by which planet-scale N partitioning is tied to mantle redox; see Nicholls et al. (2026) [^cite-nicholls2026] for an application to L 98-59 d, where the inferred H$_2$-dominated atmosphere with photochemical SO$_2$ implies $\Delta\mathrm{IW}$ between IW-4 and IW-1.

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
 [^cite-nicholls2024]: H. Nicholls, T. Lichtenberg, D. J. Bower, R. Pierrehumbert, *[Magma ocean evolution at arbitrary redox state](https://doi.org/10.1029/2024JE008576)*, Journal of Geophysical Research: Planets, 129, e2024JE008576, 2024. [SciX](https://scixplorer.org/abs/2024JGRE..12908576N/abstract).
 [^cite-nicholls2026]: H. Nicholls, T. Lichtenberg, R. D. Chatterjee, C. M. Guimond, E. Postolec, R. T. Pierrehumbert, *[Volatile-rich evolution of molten super-Earth L 98-59 d](https://doi.org/10.1038/s41550-026-02815-8)*, Nature Astronomy, 2026. [SciX](https://scixplorer.org/abs/2026NatAs.tmp...61N/abstract). [arXiv](https://arxiv.org/abs/2507.02656).
 [^cite-oneilleggins2002]: H. St. C. O'Neill, S. M. Eggins, *[The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts](https://doi.org/10.1016/S0009-2541(01)00414-4)*, Chemical Geology, 186, 151–181, 2002. [SciX](https://scixplorer.org/abs/2002ChGeo.186..151O/abstract).
 [^cite-sossi2020]: P. A. Sossi, A. D. Burnham, J. Badro, A. Lanzirotti, M. Newville, H. St. C. O'Neill, *[Redox state of Earth's magma ocean and its Venus-like early atmosphere](https://doi.org/10.1126/sciadv.abd1387)*, Science Advances, 6, eabd1387, 2020. [SciX](https://scixplorer.org/abs/2020SciA....6.1387S/abstract).
 [^cite-wadhwa2001]: M. Wadhwa, *[Redox state of Mars' upper mantle and crust from Eu anomalies in shergottite pyroxenes](https://doi.org/10.1126/science.1057594)*, Science, 291, 1527–1530, 2001. [SciX](https://scixplorer.org/abs/2001Sci...291.1527W/abstract).
