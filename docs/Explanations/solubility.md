# Solubility laws

CALLIOPE uses gas-melt **solubility** laws to relate the dissolved concentration of each volatile species in the magma ocean to its surface partial pressure. The general form is the modified Henry's law,

$$
X_i^\mathrm{melt} = \alpha_i\, p_i^{\,1/\beta_i},
$$

with $X_i^\mathrm{melt}$ in ppmw, $p_i$ in bar, and species-specific empirical constants $\alpha_i, \beta_i$. Bower et al. (2022) Equation (1) writes the same relation in terms of fugacity; CALLIOPE assumes ideal-gas behaviour so $f \equiv p$.

This page lists the implemented laws, the experimental sources behind each, and the alternative compositions a user can select. The corresponding code lives in `solubility.py`; the speciation-time call paths are in `solve.dissolved_mass`.

## Per-species inventory

### H$_2$O - `SolubilityH2O(composition=...)`

| Composition | Law | Source | $\alpha$ (ppmw bar$^{-1/\beta}$) | $\beta$ |
|---|---|---|---:|---:|
| `peridotite` (default) | $524\, p^{0.5}$ | Sossi et al. (2023) | 524 | 2 |
| `basalt_dixon` | $965\, p^{0.5}$ | Dixon et al. (1995, refit by Sossi) | 965 | 2 |
| `basalt_wilson` | $215\, p^{0.7}$ | Hamilton (1964); Wilson & Head (1981) | 215 | 1/0.7 |
| `anorthite_diopside` | $727\, p^{0.5}$ | Newcombe et al. (2017) | 727 | 2 |
| `lunar_glass` | $683\, p^{0.5}$ | Newcombe et al. (2017) | 683 | 2 |

The choice between peridotite (default) and basalt is one of the larger uncertainties in early magma-ocean modelling. Bower et al. (2022) Table 1 compares all five compositions across the relevant pressure range; Nicholls et al. (2024) uses peridotite as the fiducial, consistent with CALLIOPE's default.

### CO$_2$ - `SolubilityCO2(composition='basalt_dixon')`

Dixon et al. (1995) MORB fit, with an explicit Poynting-like temperature/pressure correction:

$$
X_{\mathrm{CO_2}}^\mathrm{melt}\,[\text{mol fr.}] = 3.8 \times 10^{-7} \cdot p_{\mathrm{CO_2}} \cdot \exp\left(-\frac{23 (p_{\mathrm{CO_2}} - 1)}{83.15\, T_\mathrm{magma}}\right)
$$

then converted from molar to ppmw via the algebraic conversion in Bower et al. (2022) Equation (3):

$$
X_{\mathrm{CO_2}}^\mathrm{melt}\,[\text{ppmw}] = 10^4 \cdot \frac{4400 X_{\mathrm{CO_2}}^\mathrm{melt}}{36.6 - 44 X_{\mathrm{CO_2}}^\mathrm{melt}}.
$$

This is the only solubility law in CALLIOPE that depends on $T_\mathrm{magma}$; the others ignore the temperature term that experimental data only weakly constrains (Bower et al. 2022 §2.2.1 discussion).

### CO - `SolubilityCO(composition='mafic_armstrong')`

Armstrong et al. (2015) mafic-melt fit:

$$
\log_{10} X_\mathrm{CO}^\mathrm{melt}\,[\text{ppmw}] = -0.738 + 0.876\, \log_{10} p_\mathrm{CO} - 5.44 \times 10^{-5} \cdot p_\mathrm{tot}
$$

The $-5.44 \times 10^{-5} \cdot p_\mathrm{tot}$ term is a total-pressure (Poynting) correction that reduces solubility at high pressures. CO solubility is generally an order of magnitude or more lower than CO$_2$, consistent with the experimental constraints summarised in Yoshioka et al. (2019).

### CH$_4$ - `SolubilityCH4(composition='basalt_ardia')`

Ardia et al. (2013) basalt fit (their Fig. 6 best-fit, 0.7-3 GPa):

$$
X_\mathrm{CH_4}^\mathrm{melt}\,[\text{ppmw}] = p_\mathrm{CH_4} \cdot \exp\left(4.93 - 1.93\,p_\mathrm{tot}^{[\mathrm{GPa}]}\right)
$$

with both pressures converted from bar to GPa internally before evaluation. The exponential term encodes the Poynting correction: at high total pressures, the molar volume of dissolved CH$_4$ in the melt becomes significantly smaller than that of gas-phase CH$_4$, so solubility *decreases* with increasing $p_\mathrm{tot}$. The factor of $\exp(1.93) \approx 6.9$ between $p_\mathrm{tot} = 1$ bar and $p_\mathrm{tot} = 1$ GPa is verified in `tests/test_stoichiometry.py::TestCH4Solubility`.

### N$_2$ - `SolubilityN2(composition='dasgupta')` (in `dissolved_mass`)

CALLIOPE provides two N$_2$ solubility laws but `dissolved_mass` hard-codes `dasgupta`:

| Composition | Law | Source |
|---|---|---|
| `libourel` | $0.0611\, p_\mathrm{N_2}$ (linear Henry's law) | Libourel et al. (2003) |
| `dasgupta` (default) | physical-state-dependent (see below) | Dasgupta et al. (2022) |

The Dasgupta et al. (2022) law adds an $f_{\mathrm{O}_2}$-dependent reduced-N contribution on top of the molecular dissolution:

$$
X_\mathrm{N_2}^\mathrm{melt}\,[\text{ppmw}] = \sqrt{p_\mathrm{N_2}^{[\mathrm{GPa}]}} \cdot \exp\left(\frac{5908\sqrt{p_\mathrm{tot}^{[\mathrm{GPa}]}}}{T_\mathrm{magma}} - 1.6\,\Delta\mathrm{IW}\right) + p_\mathrm{N_2}^{[\mathrm{GPa}]} \cdot c_\mathrm{melt}
$$

where the prefactor $c_\mathrm{melt}$ depends on the silicate composition (CALLIOPE uses Earth-mantle estimates: $x_\mathrm{SiO_2} = 0.56$, $x_\mathrm{Al_2O_3} = 0.11$, $x_\mathrm{TiO_2} = 0.01$, giving $c_\mathrm{melt} = \exp(4.67 + 7.11 \cdot 0.56 - 13.06 \cdot 0.11 - 120.67 \cdot 0.01) = e^{6.4...}$). The two terms compete: at low $p_\mathrm{N_2}$ and reducing conditions, the first term (reduced N dissolved as nitride / amide complexes) dominates; at high $p_\mathrm{N_2}$ and oxidising conditions, the second term (molecular N$_2$ Henry's law) dominates. The breakpoint occurs roughly at $p_\mathrm{N_2}^{[\mathrm{GPa}]} \approx 10^{-5}$ for $\Delta\mathrm{IW} = 0$.

### S$_2$ - `SolubilityS2(composition='gaillard')`

Gaillard et al. (2022) sulfide-saturated mafic-melt law:

$$
\log_{e} X_\mathrm{S_2}^\mathrm{melt}\,[\text{ppmw}] = 13.8426 - \frac{26476}{T_\mathrm{magma}} + 0.124\,x_\mathrm{FeO}^{[\text{wt\%}]} + 0.5\,\ln\frac{p_\mathrm{S_2}}{f_{\mathrm{O}_2}}
$$

with $x_\mathrm{FeO} = 10$ wt% hard-coded as an Earth-mantle reference. The $0.5\,\ln(p_\mathrm{S_2}/f_{\mathrm{O}_2})$ term ties the solubility to redox state: at fixed $p_\mathrm{S_2}$, more reducing conditions ($f_{\mathrm{O}_2}$ smaller) increase the dissolved sulfur (sulfide regime); more oxidising conditions decrease it.

The implementation refuses to evaluate at $p_\mathrm{S_2} < 10^{-20}$ bar (returns 0.0) to avoid `log(0)` in the rare case where the solver bottoms out at exactly zero S inventory.

## What about the missing laws?

CALLIOPE does not include explicit solubility laws for H$_2$, NH$_3$, SO$_2$, H$_2$S, or O$_2$. Their dissolved masses are computed from the *primary*-species solubilities (H$_2$O for H-bearing species, CO$_2$ for C-bearing species, S$_2$ for S-bearing species, N$_2$ for N-bearing species) via stoichiometric atom-counting in `dissolved_mass()`. This is consistent with Bower et al. (2022) Section 2.2.3 which sets $\alpha_\mathrm{H_2} = \alpha_\mathrm{CO} = \alpha_\mathrm{CH_4} = 0$ on the grounds that experimentally constrained solubilities are 1-2 dex smaller than those of H$_2$O / CO$_2$ at equivalent fugacities (Hirschmann et al. 2012; Li et al. 2015; Yoshioka et al. 2019; Ardia et al. 2013); within the CALLIOPE framework the same logic justifies omitting solubility for the three S species and ammonia.

If you need explicit dissolution of reduced species into the melt, the [atmodeller](https://atmodeller.readthedocs.io/) project provides full per-species solubility laws including H$_2$ (Hirschmann et al. 2012), CO (Yoshioka et al. 2019), and CH$_4$ (Ardia et al. 2013), with non-ideal real-gas activity coefficients.

## Validity envelope

| Species | Calibration $T$ range | Calibration $p$ range |
|---|---|---|
| H$_2$O peridotite | 2173 K | $\le$1 bar to several kbar |
| H$_2$O basalt (Dixon) | 1473 K | 176-2021 bar |
| CO$_2$ (Dixon) | $\le$2000 K | $\le$815 bar |
| CO (Armstrong) | $\sim$1700 K | $\le$3 GPa |
| CH$_4$ (Ardia) | 1573-1873 K | 0.7-3 GPa |
| N$_2$ (Dasgupta) | 1500-2200 K | 0-3 GPa |
| S$_2$ (Gaillard) | 1473-1773 K, $f_{\mathrm{O}_2} <$ IW+1 | sulfide-saturated regime |

CALLIOPE deliberately makes no attempt to flag extrapolation: the laws are evaluated formally outside their calibration ranges to keep the solver well-posed. For applications outside the bracket, treat the dissolved masses as upper bounds and check sensitivity by switching solubility laws via the constructor argument.

## See also

- [Equilibrium chemistry](equilibrium_chemistry.md): how the gas-phase speciation feeds into the partial pressures the solubility laws then consume.
- [Mass balance & solver](mass_balance.md): how dissolved + atmospheric masses are summed to close the elemental conservation constraints.
- [API reference for `calliope.solubility`](../Reference/api/calliope.solubility.md).
