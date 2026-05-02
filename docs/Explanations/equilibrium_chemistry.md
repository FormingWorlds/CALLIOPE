# Equilibrium chemistry

CALLIOPE assumes the gas phase above the magma ocean is in instantaneous chemical equilibrium at the surface temperature $T = T_\mathrm{magma}$ and at the oxygen fugacity $f_{\mathrm{O}_2}$ set by the magma redox state (see [Oxygen fugacity](oxygen_fugacity.md)). This page documents the six redox couples that produce the seven secondary species from the four primary unknowns.

## The modified equilibrium constant

For a generic gas-phase reaction with stoichiometric form

$$
A \rightleftharpoons B + n_{\mathrm{O}_2}\, \mathrm{O_2},
$$

equilibrium requires

$$
K_\mathrm{eq}(T) = \frac{p_B \cdot f_{\mathrm{O}_2}^{n_{\mathrm{O}_2}}}{p_A}.
$$

CALLIOPE rolls $f_{\mathrm{O}_2}$ into the constant by defining

$$
G_\mathrm{eq}(T, f_{\mathrm{O}_2}) \equiv 10^{\,\log_{10} K_\mathrm{eq}(T) - n_{\mathrm{O}_2}\,\log_{10} f_{\mathrm{O}_2}} = \frac{p_B}{p_A},
$$

so that the secondary partial pressure follows from the primary one as $p_B = G_\mathrm{eq} \cdot p_A$. This is the `ModifiedKeq.__call__` returns in `chemistry.py`:

```python
def __call__(self, T, fO2_shift):
    fO2 = self.fO2(T, fO2_shift)
    Keq, fO2_stoich = self.callmodel(T)
    Geq = 10 ** (Keq - fO2_stoich * fO2)
    return Geq
```

For more complex stoichiometries (CH$_4$, NH$_3$, SO$_2$, H$_2$S, which involve multi-power dependences on the primaries), the same $G_\mathrm{eq}$ form is reused but the substitution into $p_B$ involves additional factors of $p_{\mathrm{H_2}}$, $p_{\mathrm{O_2}}$, etc. - see [Speciation tree](#speciation-tree) below.

## Reaction-by-reaction calibration

The fitted constants and their sources, all in $\log_{10}$ form (CALLIOPE's `Keq, fO2_stoich` tuple):

### H$_2$ from H$_2$O ($\mathrm{H_2O} \rightleftharpoons \mathrm{H_2} + \tfrac{1}{2}\,\mathrm{O_2}$)

$$
\log_{10} K_\mathrm{eq}^{\mathrm{H_2}}(T) = -\frac{13152.4778}{T} + 3.0386
$$

(`janaf_H2` in `chemistry.py`; JANAF fit, valid $1500 \le T \le 3000$ K). Bower et al. (2022) give the equivalent Schaefer & Fegley (2017) fit $-12794/T + 2.7768$ as `schaefer_H`. The two differ by $\sim$1% in the resulting $G_\mathrm{eq}$ across the validation range; CALLIOPE uses the JANAF coefficients by default in `solve.get_partial_pressures`. Stoichiometric coefficient on $f_{\mathrm{O}_2}$ is $+0.5$.

### CO from CO$_2$ ($\mathrm{CO_2} \rightleftharpoons \mathrm{CO} + \tfrac{1}{2}\,\mathrm{O_2}$)

$$
\log_{10} K_\mathrm{eq}^{\mathrm{CO}}(T) = -\frac{14467.5114}{T} + 4.3481
$$

(`janaf_CO`; JANAF fit). Schaefer & Fegley (2017) equivalent in `schaefer_C`: $-14787/T + 4.5472$. Stoichiometric coefficient $+0.5$.

### CH$_4$ from CO$_2$ + 2H$_2$ ($\mathrm{CO_2} + 2\,\mathrm{H_2} \rightleftharpoons \mathrm{CH_4} + \mathrm{O_2}$)

$$
\log_{10} K_\mathrm{eq}^{\mathrm{CH_4}}(T) = -\frac{16276}{T} - 5.4738
$$

(`schaefer_CH4`, IVTHANTHERMO via Schaefer & Fegley 2017). Stoichiometric coefficient $+1.0$. The methane abundance is then $p_{\mathrm{CH_4}} = G_\mathrm{eq} \cdot p_{\mathrm{CO_2}} \cdot p_{\mathrm{H_2}}^2$; this is the only species whose speciation depends on a second primary species (H$_2$ via H$_2$O).

### SO$_2$ from S$_2$ + 2O$_2$ ($\mathrm{S_2} + 2\,\mathrm{O_2} \rightleftharpoons 2\,\mathrm{SO_2}$, doubled form)

$$
\log_{10} K_\mathrm{eq}^{\mathrm{SO_2}}(T) = +\frac{37774}{T} - 7.6128
$$

(`janaf_SO2`; doubled JANAF fit, factor 2 of the formation constant for $\tfrac{1}{2}\mathrm{S_2} + \mathrm{O_2} \to \mathrm{SO_2}$). Stoichiometric coefficient is **0** because the O$_2$ dependence is carried explicitly by the $p_{\mathrm{O_2}}^2$ factor in the recovery formula

$$
p_{\mathrm{SO_2}} = \sqrt{G_\mathrm{eq}^{\mathrm{SO_2}} \cdot p_{\mathrm{S_2}} \cdot p_{\mathrm{O_2}}^2}.
$$

This unusual encoding (zero `fO2_stoich`, but explicit $p_{\mathrm{O_2}}$ in the recovery formula) is documented in `chemistry.janaf_SO2` and is what makes the analytical-vs-code consistency test in `tests/test_stoichiometry.py::TestEquilibriumChemistry::test_SO2_equilibrium` parametric over $T$.

### H$_2$S from S$_2$ + 2H$_2$ ($\mathrm{S_2} + 2\,\mathrm{H_2} \rightleftharpoons 2\,\mathrm{H_2S}$, doubled form)

$$
\log_{10} K_\mathrm{eq}^{\mathrm{H_2S}}(T) = +\frac{13462.0309}{T} - 7.2455
$$

(`janaf_H2S`; doubled JANAF fit, fitted to JANAF table data in `tools/fit_H2S.ipynb`). Stoichiometric coefficient 0; recovery is

$$
p_{\mathrm{H_2S}} = \sqrt{G_\mathrm{eq}^{\mathrm{H_2S}} \cdot p_{\mathrm{S_2}} \cdot p_{\mathrm{H_2}}^2}.
$$

### NH$_3$ from N$_2$ + 3H$_2$ ($\mathrm{N_2} + 3\,\mathrm{H_2} \rightleftharpoons 2\,\mathrm{NH_3}$, doubled form)

$$
\log_{10} K_\mathrm{eq}^{\mathrm{NH_3}}(T) = +\frac{5328.0312}{T} - 11.9848
$$

(`janaf_NH3`; doubled JANAF fit, fitted in `tools/fit_NH3.ipynb`). Stoichiometric coefficient 0; recovery is

$$
p_{\mathrm{NH_3}} = \sqrt{G_\mathrm{eq}^{\mathrm{NH_3}} \cdot p_{\mathrm{N_2}} \cdot p_{\mathrm{H_2}}^3}.
$$

!!! note "On the doubled-form notation"
    The "doubled form" for SO$_2$, H$_2$S, NH$_3$ is a notational trick: rather than parameterising the formation constant $K_f$ for $\tfrac{1}{2}A + B \to AB$, CALLIOPE stores $K_f^2$ (which is $K$ for $A + 2B \to 2AB$). This lets the `Geq * p_S2 * p_O2^2` product in `solve.py` be square-rooted in one go to recover the formation product, which is more numerically stable when partial pressures span 12 dex than computing the formation constant and multiplying by a half-power directly. Whatever it gains in stability it costs in clarity, which is why the test file ties the analytical $K_f$ form back to the code form for every doubled-form species.

## Speciation tree

`get_partial_pressures(pin, ddict)` in `solve.py` walks the eleven species in this order, with conditional gating on the per-species `_included` flag:

| Step | Species | Source | Formula |
|---|---|---|---|
| 1 | H$_2$O | primary input | $p_\mathrm{H_2O} = \mathrm{pin[H_2O]}$ |
| 2 | H$_2$ | reduces H$_2$O | $p_\mathrm{H_2} = G^{\mathrm{H_2}} \cdot p_\mathrm{H_2O}$ |
| 3 | CO$_2$ | primary input | $p_\mathrm{CO_2} = \mathrm{pin[CO_2]}$ |
| 4 | CO | reduces CO$_2$ | $p_\mathrm{CO} = G^{\mathrm{CO}} \cdot p_\mathrm{CO_2}$ |
| 5 | CH$_4$ | requires H$_2$ + CO$_2$ | $p_\mathrm{CH_4} = G^{\mathrm{CH_4}} \cdot p_\mathrm{CO_2} \cdot p_\mathrm{H_2}^2$ |
| 6 | N$_2$ | primary input | $p_\mathrm{N_2} = \mathrm{pin[N_2]}$ |
| 7 | NH$_3$ | requires N$_2$ + H$_2$ | $p_\mathrm{NH_3} = \sqrt{G^{\mathrm{NH_3}} \cdot p_\mathrm{N_2} \cdot p_\mathrm{H_2}^3}$ |
| 8 | O$_2$ | $f_{\mathrm{O}_2}$ buffer | $p_\mathrm{O_2} = 10^{\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T) + \Delta\mathrm{IW}}$ |
| 9 | S$_2$ | primary input | $p_\mathrm{S_2} = \mathrm{pin[S_2]}$ |
| 10 | SO$_2$ | requires S$_2$ + O$_2$ | $p_\mathrm{SO_2} = \sqrt{G^{\mathrm{SO_2}} \cdot p_\mathrm{S_2} \cdot p_\mathrm{O_2}^2}$ |
| 11 | H$_2$S | requires S$_2$ + H$_2$ | $p_\mathrm{H_2S} = \sqrt{G^{\mathrm{H_2S}} \cdot p_\mathrm{S_2} \cdot p_\mathrm{H_2}^2}$ |

After the walk, every `p_d[s]` is clipped to be non-negative. The function returns the dictionary that downstream pieces (atmospheric mass, dissolved mass) consume.

## Why this is the right level of detail

The choice to expand only six redox couples is a deliberate trade-off between completeness and well-posedness. Adding more species (e.g. HCN, HCl, CS$_2$) requires either (i) more elemental constraints (Cl, additional H atoms in HCN), or (ii) extra equilibrium reactions whose constants are calibrated outside the relevant $T$-range. In the present species set, every secondary species has a clean reduction back to one of the four primaries via a JANAF / IVTHANTHERMO fit that is validated for $1500 \le T \le 3000$ K, which matches the magma-ocean regime CALLIOPE is targeted at.

For application contexts that require Cl-bearing species, sub-ideal real-gas effects, or condensation, the [atmodeller](https://atmodeller.readthedocs.io/) JAX-based solver (Bower et al. 2025) is the supported alternative within the PROTEUS framework.

## See also

- [Oxygen fugacity](oxygen_fugacity.md): how the IW buffer enters the modified equilibrium constants
- [Mass balance & solver](mass_balance.md): how the four primary partial pressures are determined from the elemental conservation constraints
- [Solubility laws](solubility.md): how the dissolved-volatile masses close the system
