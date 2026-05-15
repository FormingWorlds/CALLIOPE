# Solubility laws

CALLIOPE uses gas-melt **solubility** laws to relate the dissolved concentration of each volatile species in the magma ocean to its surface partial pressure. The general form is the modified Henry's law,

$$
X_i^\mathrm{melt} = \alpha_i\, p_i^{\,1/\beta_i},
$$

with $X_i^\mathrm{melt}$ in ppmw, $p_i$ in bar, and species-specific empirical constants $\alpha_i, \beta_i$. Bower et al. (2022)[^cite-bower2022] Equation (1) writes the same relation in terms of fugacity; CALLIOPE assumes ideal-gas behaviour so $f \equiv p$.

This page lists the implemented laws, the experimental sources behind each, and the alternative compositions a user can select. The corresponding code lives in `solubility.py`; the speciation-time call paths are in `solve.dissolved_mass`.

## Per-species inventory

### H$_2$O - `SolubilityH2O(composition=...)`

| Composition | Law | Source | $\alpha$ (ppmw bar$^{-1/\beta}$) | $\beta$ |
|---|---|---|---:|---:|
| `peridotite` (default) | $524\, p^{0.5}$ | Sossi et al. (2023)[^cite-sossi2023] | 524 | 2 |
| `basalt_dixon` | $965\, p^{0.5}$ | Dixon et al. (1995)[^cite-dixon1995], refit by Sossi | 965 | 2 |
| `basalt_wilson` | $215\, p^{0.7}$ | Hamilton (1964)[^cite-hamilton1964]; Wilson & Head (1981)[^cite-wilsonhead1981] | 215 | 1/0.7 |
| `anorthite_diopside` | $727\, p^{0.5}$ | Newcombe et al. (2017)[^cite-newcombe2017] | 727 | 2 |
| `lunar_glass` | $683\, p^{0.5}$ | Newcombe et al. (2017)[^cite-newcombe2017] | 683 | 2 |

!!! note "About the `peridotite` constant 524"
    Sossi et al. (2023)[^cite-sossi2023] report two values for the prefactor depending on the FTIR absorption coefficient used: $\alpha = 524 \pm 16$ ppmw bar$^{-0.5}$ from the basaltic-glass calibration ($\epsilon_{3550} = 6.3$ m$^2$ mol$^{-1}$), and $\alpha = 647$ ppmw bar$^{-0.5}$ from the peridotite-glass calibration. CALLIOPE uses 524 to match the value adopted in Nicholls et al. (2024)[^cite-nicholls2024] and the PROTEUS-side fiducial; the label `peridotite` refers to the *experimental melt composition* (Sossi+2023 used a peridotitic starting composition), not to the spectroscopic-calibration choice. If you want the full peridotite-glass calibration, instantiate `SolubilityH2O('peridotite')` and override the constant manually, or wait for an upstream change.

The choice between peridotite (default) and basalt is one of the larger uncertainties in early magma-ocean modelling. Bower et al. (2022)[^cite-bower2022] Table 1 compares all five compositions across the relevant pressure range; Nicholls et al. (2024)[^cite-nicholls2024] uses peridotite as the fiducial, consistent with CALLIOPE's default.

### CO$_2$ - `SolubilityCO2(composition='basalt_dixon')`

Dixon et al. (1995)[^cite-dixon1995] MORB fit, with an explicit Poynting-like temperature/pressure correction:

$$
X_{\mathrm{CO_2}}^\mathrm{melt}\,[\text{mol fr.}] = 3.8 \times 10^{-7} \cdot p_{\mathrm{CO_2}} \cdot \exp\left(-\frac{23 (p_{\mathrm{CO_2}} - 1)}{83.15\, T_\mathrm{magma}}\right)
$$

then converted from molar to ppmw via the algebraic conversion in Bower et al. (2022)[^cite-bower2022] Equation (3):

$$
X_{\mathrm{CO_2}}^\mathrm{melt}\,[\text{ppmw}] = 10^4 \cdot \frac{4400 X_{\mathrm{CO_2}}^\mathrm{melt}}{36.6 - 44 X_{\mathrm{CO_2}}^\mathrm{melt}}.
$$

This is the only solubility law in CALLIOPE that depends on $T_\mathrm{magma}$; the others ignore the temperature term that experimental data only weakly constrains (Bower et al. 2022[^cite-bower2022] §2.2.1 discussion).

### CO - `SolubilityCO(composition='mafic_armstrong')`

Armstrong et al. (2015)[^cite-armstrong2015] mafic-melt fit:

$$
\log_{10} X_\mathrm{CO}^\mathrm{melt}\,[\text{ppmw}] = -0.738 + 0.876\, \log_{10} p_\mathrm{CO} - 5.44 \times 10^{-5} \cdot p_\mathrm{tot}
$$

The $-5.44 \times 10^{-5} \cdot p_\mathrm{tot}$ term is a total-pressure (Poynting) correction that reduces solubility at high pressures. CO solubility is generally an order of magnitude or more lower than CO$_2$, consistent with the experimental constraints summarised in Yoshioka et al. (2019)[^cite-yoshioka2019].

### CH$_4$ - `SolubilityCH4(composition='basalt_ardia')`

Ardia et al. (2013)[^cite-ardia2013] Fe-free haplobasaltic-melt fit (their Eq. (8) with $\ln K_0 = 4.93$ and $\Delta V = 26.85\,\mathrm{cm}^3$/mol at $T_0 = 1400\,^\circ$C, $P_0 = 1$ bar, plotted as Fig. 11; calibrated over 0.7-3 GPa):

$$
X_\mathrm{CH_4}^\mathrm{melt}\,[\text{ppmw}] = p_\mathrm{CH_4} \cdot \exp\left(4.93 - 1.93\,p_\mathrm{tot}^{[\mathrm{GPa}]}\right)
$$

with both pressures converted from bar to GPa internally before evaluation. The exponential term encodes the Poynting correction: at high total pressures, the molar volume of dissolved CH$_4$ in the melt becomes significantly smaller than that of gas-phase CH$_4$, so solubility *decreases* with increasing $p_\mathrm{tot}$. The factor of $\exp(1.93) \approx 6.9$ between $p_\mathrm{tot} = 1$ bar and $p_\mathrm{tot} = 1$ GPa is verified in `tests/test_stoichiometry.py::TestCH4Solubility`.

### N$_2$ - `SolubilityN2(composition='dasgupta')` (in `dissolved_mass`)

CALLIOPE provides two N$_2$ solubility laws but `dissolved_mass` hard-codes `dasgupta`:

| Composition | Law | Source |
|---|---|---|
| `libourel` | $0.0611\, p_\mathrm{N_2}$ (linear Henry's law) | Libourel et al. (2003)[^cite-libourel2003] |
| `dasgupta` (default) | physical-state-dependent (see below) | Dasgupta et al. (2022)[^cite-dasgupta2022] |

The Dasgupta et al. (2022)[^cite-dasgupta2022] law adds an $f_{\mathrm{O}_2}$-dependent reduced-N contribution on top of the molecular dissolution:

$$
X_\mathrm{N_2}^\mathrm{melt}\,[\text{ppmw}] = \sqrt{p_\mathrm{N_2}^{[\mathrm{GPa}]}} \cdot \exp\left(\frac{5908\sqrt{p_\mathrm{tot}^{[\mathrm{GPa}]}}}{T_\mathrm{magma}} - 1.6\,\Delta\mathrm{IW}\right) + p_\mathrm{N_2}^{[\mathrm{GPa}]} \cdot c_\mathrm{melt}
$$

where the prefactor $c_\mathrm{melt}$ depends on the silicate composition. The mole fractions $x_\mathrm{SiO_2}$, $x_\mathrm{Al_2O_3}$, $x_\mathrm{TiO_2}$ are exposed as constructor kwargs `SolubilityN2(x_SiO2=..., x_Al2O3=..., x_TiO2=...)` and default to the Earth-mantle estimates $x_\mathrm{SiO_2} = 0.56$, $x_\mathrm{Al_2O_3} = 0.11$, $x_\mathrm{TiO_2} = 0.01$, giving $c_\mathrm{melt} = \exp(4.67 + 7.11 \cdot 0.56 - 13.06 \cdot 0.11 - 120.67 \cdot 0.01) \approx 407$. The two terms compete: at low $p_\mathrm{N_2}$ and reducing conditions, the first term (reduced N dissolved as nitride / amide complexes) dominates; at high $p_\mathrm{N_2}$ and oxidising conditions, the second term (molecular N$_2$ Henry's law) dominates. The breakpoint occurs roughly at $p_\mathrm{N_2}^{[\mathrm{GPa}]} \approx 10^{-5}$ for $\Delta\mathrm{IW} = 0$.

### S$_2$ - `SolubilityS2(composition='gaillard')`

Gaillard et al. (2022)[^cite-gaillard2022] sulfide-saturated mafic-melt law:

$$
\log_{e} X_\mathrm{S_2}^\mathrm{melt}\,[\text{ppmw}] = 13.8426 - \frac{26476}{T_\mathrm{magma}} + 0.124\,x_\mathrm{FeO}^{[\text{wt\%}]} + 0.5\,\ln\frac{p_\mathrm{S_2}}{f_{\mathrm{O}_2}}
$$

where $x_\mathrm{FeO}$ is the melt FeO content in wt%, exposed as the constructor kwarg `SolubilityS2(x_FeO=...)` and defaulting to $10$ wt% (Earth-mantle reference). The $0.5\,\ln(p_\mathrm{S_2}/f_{\mathrm{O}_2})$ term ties the solubility to redox state: at fixed $p_\mathrm{S_2}$, more reducing conditions ($f_{\mathrm{O}_2}$ smaller) increase the dissolved sulfur (sulfide regime); more oxidising conditions decrease it.

The implementation refuses to evaluate at $p_\mathrm{S_2} < 10^{-20}$ bar (returns 0.0) to avoid `log(0)` in the rare case where the solver bottoms out at exactly zero S inventory.

## What about the missing laws?

CALLIOPE does not include explicit solubility laws for H$_2$, NH$_3$, SO$_2$, H$_2$S, or O$_2$. Their dissolved masses are computed from the *primary*-species solubilities (H$_2$O for H-bearing species, CO$_2$ for C-bearing species, S$_2$ for S-bearing species, N$_2$ for N-bearing species) via stoichiometric atom-counting in `dissolved_mass()`. This is consistent with Bower et al. (2022)[^cite-bower2022] Section 2.2.3 which sets $\alpha_\mathrm{H_2} = \alpha_\mathrm{CO} = \alpha_\mathrm{CH_4} = 0$ on the grounds that experimentally constrained solubilities are 1-2 dex smaller than those of H$_2$O / CO$_2$ at equivalent fugacities (Hirschmann et al. 2012[^cite-hirschmann2012]; Li et al. 2015[^cite-li2015]; Yoshioka et al. 2019[^cite-yoshioka2019]; Ardia et al. 2013[^cite-ardia2013]); within the CALLIOPE framework the same logic justifies omitting solubility for the three S species and ammonia.

If you need explicit dissolution of reduced species into the melt, the [atmodeller](https://atmodeller.readthedocs.io/) project provides full per-species solubility laws including H$_2$ (Hirschmann et al. 2012[^cite-hirschmann2012]), CO (Yoshioka et al. 2019[^cite-yoshioka2019]), and CH$_4$ (Ardia et al. 2013[^cite-ardia2013]), with non-ideal real-gas activity coefficients.

## Validity envelope

| Species (`composition`) | Source | Calibration $T$ | Calibration $p$ | Calibration $f_{\mathrm{O}_2}$ |
|---|---|---|---|---|
| H$_2$O `peridotite` (default) | Sossi et al. (2023)[^cite-sossi2023] | 2173 K (1900 $^\circ$C) | 1 atm total ($f_\mathrm{H_2O}\le 0.027$ bar) | IW-1.9 to IW+6.0 |
| H$_2$O `basalt_dixon` | Dixon et al. (1995)[^cite-dixon1995] | 1473 K (1200 $^\circ$C) | 176-717 bar $p_\mathrm{H_2O}$ | $\sim$QFM ($\approx$ IW+3.5 to IW+5) |
| H$_2$O `basalt_wilson` | Hamilton et al. (1964)[^cite-hamilton1964] | 1373 K (1100 $^\circ$C) | 1000-6000 bar $p_\mathrm{H_2O}$ | unbuffered for the pressure series (separate 1000-bar buffer series used MH, FMQ, MW buffers) |
| H$_2$O `anorthite_diopside` | Newcombe et al. (2017)[^cite-newcombe2017] | 1623 K (1350 $^\circ$C) | 1 atm | IW-2.3 to IW+4.8 |
| H$_2$O `lunar_glass` | Newcombe et al. (2017)[^cite-newcombe2017] | 1623 K (1350 $^\circ$C) | 1 atm | IW-3.0 to IW+4.8 |
| CO$_2$ `basalt_dixon` (default) | Dixon et al. (1995)[^cite-dixon1995] | 1473 K (1200 $^\circ$C) | $\le$815 bar $p_\mathrm{CO_2}$ | $\sim$QFM ($\approx$ IW+3.5 to IW+5) |
| CO `mafic_armstrong` (default) | Armstrong et al. (2015)[^cite-armstrong2015] | 1673 K (1400 $^\circ$C) | 1.2 GPa $p_\mathrm{tot}$ (1.0-1.2 GPa including Stanley et al. 2014 data co-fit) | IW-3.65 to IW+1.46 |
| CH$_4$ `basalt_ardia` (default) | Ardia et al. (2013)[^cite-ardia2013] | 1673-1723 K (1400-1450 $^\circ$C) | 0.7-3 GPa $p_\mathrm{tot}$ | IW-9.5 to IW-1.4 (IW/Si and IWC/C buffers) |
| N$_2$ `libourel` | Libourel et al. (2003)[^cite-libourel2003] | 1673-1698 K (1400-1425 $^\circ$C) | 1 atm | linear regime $\log_{10}f_{\mathrm{O}_2}\in[-10.7, -0.7]$ ($\approx$ IW-1.3 to IW+9) |
| N$_2$ `dasgupta` (default in `dissolved_mass`) | Dasgupta et al. (2022)[^cite-dasgupta2022] | 1323-2600 K (1050-2327 $^\circ$C) | 1 bar to 8.2 GPa $p_\mathrm{tot}$ | IW-8.3 to IW+8.7 |
| S$_2$ `gaillard` (default) | Gaillard et al. (2022)[^cite-gaillard2022] | not stated by paper (1 atm calibration data) | 1 atm | IW-1 to FMQ+0.1 ($\approx$ IW+3.5) |

The $f_{\mathrm{O}_2}$ column reports the *calibration footprint* (the range of $f_{\mathrm{O}_2}$ over which the experiments that produced the fit were performed), not the law's functional dependence: only the Gaillard S$_2$ and Dasgupta N$_2$ laws carry an explicit $f_{\mathrm{O}_2}$ term, while every other law in the table is an $f_{\mathrm{O}_2}$-independent expression fitted to data from experiments performed at the listed $f_{\mathrm{O}_2}$ range. The choice between H$_2$O laws (peridotite vs basalt) is one of the larger uncertainties in early magma-ocean modelling (the line-26 note expands on the Sossi peridotite vs Dixon basalt prefactor difference); the $f_{\mathrm{O}_2}$ footprint here describes where the data lived, not how robust the fit is across compositions.

The Dasgupta N$_2$ and Gaillard S$_2$ rows quote the ranges that Dasgupta et al. (2022)[^cite-dasgupta2022] Equation 10 (n=137 compiled data) and Gaillard et al. (2022)[^cite-gaillard2022] Equation 10 (refit of O'Neill & Mavrogenes (2002)[^cite-oneillmavrogenes2002] plus 8 other experimental sources, n=369) report directly. The Gaillard refit data are all at 1 atm; the formula is applied at magma-ocean pressures in CALLIOPE without an explicit pressure correction. The Libourel linear regime stops at IW-1.3, below which the same paper documents a sharp transition to chemical (network-bound N$^{3-}$) dissolution with $\sim$5 orders of magnitude higher solubility; the `libourel` law in CALLIOPE captures only the oxidising-end linear regime and underestimates dissolved N below IW-1.3. The CO, CH$_4$ and Dasgupta-N$_2$ rows give the *total*-pressure range, since the Poynting and reduced-N branches depend on $p_\mathrm{tot}$ as well as the species partial pressure.

CALLIOPE deliberately makes no attempt to flag extrapolation: the laws are evaluated formally outside their calibration ranges to keep the solver well-posed. For applications outside the bracket, treat the dissolved masses as upper bounds and check sensitivity by switching solubility laws via the constructor argument.

## See also

- [Equilibrium chemistry](equilibrium_chemistry.md): how the gas-phase speciation feeds into the partial pressures the solubility laws then consume.
- [Mass balance & solver](mass_balance.md): how dissolved + atmospheric masses are summed to close the elemental conservation constraints.
- [Oxygen fugacity](oxygen_fugacity.md): the IW-buffer parameterisations that drive the $f_{\mathrm{O}_2}$ dependence in the Gaillard S$_2$ and Dasgupta N$_2$ laws.
- [Authoritative-oxygen mode](authoritative_oxygen.md): how the Gaillard $\ln(p_\mathrm{S_2}/f_{\mathrm{O}_2})$ and Dasgupta $-1.6\,\Delta\mathrm{IW}$ terms enter the 5-residual mass balance once $\Delta\mathrm{IW}$ becomes a solver unknown.
- [API reference for `calliope.solubility`](../Reference/api/calliope.solubility.md).

[^cite-armstrong2015]: L. S. Armstrong, M. M. Hirschmann, B. D. Stanley, E. G. Falksen, S. D. Jacobsen, *[Speciation and solubility of reduced C-O-H-N volatiles in mafic melt: implications for volcanism, atmospheric evolution, and deep volatile cycles in the terrestrial planets](https://doi.org/10.1016/j.gca.2015.07.007)*, Geochimica et Cosmochimica Acta, 171, 283–302, 2015. [SciX](https://scixplorer.org/abs/2015GeCoA.171..283A/abstract).
[^cite-ardia2013]: P. Ardia, M. M. Hirschmann, A. C. Withers, B. D. Stanley, *[Solubility of CH$_4$ in a synthetic basaltic melt, with applications to atmosphere-magma ocean-core partitioning of volatiles and to the evolution of the Martian atmosphere](https://doi.org/10.1016/j.gca.2013.03.028)*, Geochimica et Cosmochimica Acta, 114, 52–71, 2013. [SciX](https://scixplorer.org/abs/2013GeCoA.114...52A/abstract).
[^cite-bower2022]: D. J. Bower, K. Hakim, P. A. Sossi, P. Sanan, *[Retention of water in terrestrial magma oceans and carbon-rich early atmospheres](https://doi.org/10.3847/PSJ/ac5fb1)*, The Planetary Science Journal, 3(4), 93, 2022. [SciX](https://scixplorer.org/abs/2022PSJ.....3...93B/abstract).
[^cite-dasgupta2022]: R. Dasgupta, E. Falksen, A. Pal, C. Sun, *[The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions](https://doi.org/10.1016/j.gca.2022.09.012)*, Geochimica et Cosmochimica Acta, 336, 291–307, 2022. [SciX](https://scixplorer.org/abs/2022GeCoA.336..291D/abstract).
[^cite-dixon1995]: J. E. Dixon, E. M. Stolper, J. R. Holloway, *[An experimental study of water and carbon dioxide solubilities in mid-ocean ridge basaltic liquids. Part I: Calibration and solubility models](https://doi.org/10.1093/oxfordjournals.petrology.a037267)*, Journal of Petrology, 36(6), 1607–1631, 1995. [SciX](https://scixplorer.org/abs/1995JPet...36.1607D/abstract).
[^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
[^cite-hamilton1964]: D. L. Hamilton, C. W. Burnham, E. F. Osborn, *[The solubility of water and effects of oxygen fugacity and water content on crystallization in mafic magmas](https://doi.org/10.1093/petrology/5.1.21)*, Journal of Petrology, 5(1), 21–39, 1964.
[^cite-hirschmann2012]: M. M. Hirschmann, A. C. Withers, P. Ardia, N. T. Foley, *[Solubility of molecular hydrogen in silicate melts and consequences for volatile evolution of terrestrial planets](https://doi.org/10.1016/j.epsl.2012.06.031)*, Earth and Planetary Science Letters, 345–348, 38–48, 2012. [SciX](https://scixplorer.org/abs/2012E%26PSL.345...38H/abstract).
[^cite-li2015]: Y. Li, R. Dasgupta, K. Tsuno, *[The effects of sulfur, silicon, water, and oxygen fugacity on carbon solubility and partitioning in Fe-rich alloy and silicate melt systems at 3 GPa and 1600 °C: implications for core-mantle differentiation and degassing of magma oceans and reduced planetary mantles](https://doi.org/10.1016/j.epsl.2015.01.017)*, Earth and Planetary Science Letters, 415, 54–66, 2015. [SciX](https://scixplorer.org/abs/2015E%26PSL.415...54L/abstract).
[^cite-libourel2003]: G. Libourel, B. Marty, F. Humbert, *[Nitrogen solubility in basaltic melt. Part I. Effect of oxygen fugacity](https://doi.org/10.1016/S0016-7037(03)00259-X)*, Geochimica et Cosmochimica Acta, 67(21), 4123–4135, 2003. [SciX](https://scixplorer.org/abs/2003GeCoA..67.4123L/abstract).
[^cite-newcombe2017]: M. E. Newcombe, A. Brett, J. R. Beckett, M. B. Baker, S. Newman, Y. Guan, J. M. Eiler, E. M. Stolper, *[Solubility of water in lunar basalt at low pH$_2$O](https://doi.org/10.1016/j.gca.2016.12.026)*, Geochimica et Cosmochimica Acta, 200, 330–352, 2017. [SciX](https://scixplorer.org/abs/2017GeCoA.200..330N/abstract).
[^cite-nicholls2024]: H. Nicholls, T. Lichtenberg, D. J. Bower, R. Pierrehumbert, *[Magma ocean evolution at arbitrary redox state](https://doi.org/10.1029/2024JE008576)*, Journal of Geophysical Research: Planets, 129, e2024JE008576, 2024. [SciX](https://scixplorer.org/abs/2024JGRE..12908576N/abstract).
[^cite-oneillmavrogenes2002]: H. St. C. O'Neill, J. A. Mavrogenes, *[The sulfide capacity and the sulfur content at sulfide saturation of silicate melts at 1400 °C and 1 bar](https://doi.org/10.1093/petrology/43.6.1049)*, Journal of Petrology, 43(6), 1049–1087, 2002. [SciX](https://scixplorer.org/abs/2002JPet...43.1049O/abstract).
[^cite-sossi2023]: P. A. Sossi, P. M. E. Tollan, J. Badro, D. J. Bower, *[Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets](https://doi.org/10.1016/j.epsl.2022.117894)*, Earth and Planetary Science Letters, 601, 117894, 2023. [SciX](https://scixplorer.org/abs/2023E%26PSL.60117894S/abstract).
[^cite-wilsonhead1981]: L. Wilson, J. W. Head, *[Ascent and eruption of basaltic magma on the Earth and Moon](https://doi.org/10.1029/JB086iB04p02971)*, Journal of Geophysical Research, 86(B4), 2971–3001, 1981. [SciX](https://scixplorer.org/abs/1981JGR....86.2971W/abstract).
[^cite-yoshioka2019]: T. Yoshioka, D. Nakashima, T. Nakamura, S. Shcheka, H. Keppler, *[Carbon solubility in silicate melts in equilibrium with a CO-CO$_2$ gas phase and graphite](https://doi.org/10.1016/j.gca.2019.06.007)*, Geochimica et Cosmochimica Acta, 259, 129–143, 2019. [SciX](https://scixplorer.org/abs/2019GeCoA.259..129Y/abstract).
