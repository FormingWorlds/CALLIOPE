# Oxygen fugacity

CALLIOPE parameterises the redox state of the magma ocean through the oxygen fugacity $f_{\mathrm{O}_2}$ at the surface, expressed in $\log_{10}$ units relative to the iron-wüstite (IW) mineral buffer:

$$
\Delta\mathrm{IW} \equiv \log_{10} f_{\mathrm{O}_2} - \log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T).
$$

The user supplies `fO2_shift_IW = ` $\Delta\mathrm{IW}$ as a scalar input; CALLIOPE computes the absolute $\log_{10} f_{\mathrm{O}_2}$ at $T = T_\mathrm{magma}$ by adding the buffer value to the shift. This page documents the buffer parameterisations, the physical meaning of $\Delta\mathrm{IW}$, and how $f_{\mathrm{O}_2}$ enters the chemistry.

## The IW mineral buffer

Iron in equilibrium with wüstite (FeO) in the presence of O$_2$ obeys

$$
\mathrm{Fe} + \tfrac{1}{2}\,\mathrm{O_2} \rightleftharpoons \mathrm{FeO},
$$

which fixes a single curve $\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T)$ in $T$-$f_{\mathrm{O}_2}$ space. CALLIOPE supports two parameterisations of this curve.

### O'Neill & Eggins (2002), `oneill` (default)

A thermochemically-constrained fit derived from low-temperature equilibrium data, expressed as Bower et al. (2022) Equation (7):

$$
\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T) = \frac{2\left[-244118 + 115.559\,T - 8.474\,T \ln T\right]}{\ln(10)\, R\, T},
$$

with $R = 8.31441$ J K$^{-1}$ mol$^{-1}$. Bower et al. (2022) adopted this as the "IW buffer to which $f_{\mathrm{O}_2}$ is referenced". This is the function `OxygenFugacity.oneill(T)` in `oxygen_fugacity.py`.

### Fischer et al. (2011), `fischer`

A simpler two-parameter fit calibrated against high-pressure ($\sim$25 GPa) experimental data:

$$
\log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T) = 6.94059 - \frac{28180.8}{T}.
$$

Implemented as `OxygenFugacity.fischer(T)`.

!!! note "Choice of buffer"
    The two parameterisations agree to $\le$0.2 dex over the full magma-ocean temperature range (1500-3500 K), and their numerical predictions at fixed $T$ differ by less than the typical uncertainty in $\Delta\mathrm{IW}$ inferred from petrological observations (Sossi et al. 2020 give $\Delta\mathrm{IW} = +3.5 \pm 0.5$ for Earth's modern surface mantle). The choice between them rarely matters for final partial pressures; CALLIOPE defaults to O'Neill & Eggins because Bower et al. (2022) used it and it was validated end-to-end against PROTEUS coupled runs.

## How $\Delta\mathrm{IW}$ enters the chemistry

The shift $\Delta\mathrm{IW}$ feeds the equilibrium chemistry through two channels:

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

For the H$_2$O-H$_2$ couple ($n_\mathrm{O_2} = +0.5$), reducing conditions ($f_{\mathrm{O}_2}$ smaller, $\Delta\mathrm{IW}$ more negative) drive $G_\mathrm{eq}$ larger, which in turn drives more H$_2$O to dissociate into H$_2$. This is the redox dependence visible in the redox-sweep tutorial, and it is the mechanism behind the H$_2$-dominated, long-lived magma-ocean atmospheres found in Nicholls et al. (2024) Figure 6 at $\Delta\mathrm{IW} \le -2$.

### 3. The S$_2$ Gaillard solubility

The Gaillard et al. (2022) sulfide solubility carries an explicit $\ln f_{\mathrm{O}_2}$ term, so $\Delta\mathrm{IW}$ enters the dissolved-S inventory directly. The implementation in `solubility.SolubilityS2.gaillard` calls back into `OxygenFugacity()` to compute the absolute $f_{\mathrm{O}_2}$.

### 4. The N$_2$ Dasgupta solubility

Similarly, the Dasgupta et al. (2022) N$_2$ solubility includes a $-1.6\,\Delta\mathrm{IW}$ term in its exponent, so reducing conditions sharply increase the dissolved-N inventory. This is one mechanism by which planet-scale N partitioning is tied to mantle redox; see Nicholls et al. (2026) Section 4 for an application to L 98-59 d, where the inferred SO$_2$/H$_2$ atmosphere implies $\Delta\mathrm{IW} \approx -1$.

## Reference values for $\Delta\mathrm{IW}$

| Reservoir | $\Delta\mathrm{IW}$ | Source |
|---|---|---|
| Mercury surface | $\sim -5$ | Cartier & Wood (2019) |
| Asteroidal material | $\sim -2$ | Doyle et al. (2019) |
| Mars mantle | $-3$ to $0$ | Wadhwa (2001) |
| Iron-wüstite buffer | $0$ | by definition |
| Earth modern surface | $+3.5$ | Sossi et al. (2020) |
| Earth modern surface (preferred) | $+3$ to $+5$ | Sossi et al. (2020) compilation |
| Modern Earth (deep mantle) | $\sim 0$ | Frost & McCammon (2008) |

CALLIOPE's PROTEUS-side default is `fO2_shift_IW = 4.0`, consistent with a near-surface terrestrial composition. Nicholls et al. (2024) explored $\Delta\mathrm{IW} \in \{-5, -3, -1, 0, +1, +3, +5\}$ on a 7-point grid and demonstrated that the resulting atmospheric composition spans the full range from H$_2$-dominated reduced atmospheres (TRAPPIST-1 c-like) to H$_2$O/CO$_2$-dominated oxidised atmospheres (Earth-like).

## Limitations

- **No $f_{\mathrm{O}_2}$ evolution**: $\Delta\mathrm{IW}$ is a constant input, not a state variable. In reality the mantle $f_{\mathrm{O}_2}$ should evolve with degree of crystallisation, fractional crystallisation depth, and atmospheric escape; CALLIOPE does not capture this and the user is responsible for choosing a representative value or sweeping over a grid.
- **No $f_{\mathrm{O}_2}$ depth profile**: the surface $f_{\mathrm{O}_2}$ alone enters the chemistry. Bower et al. (2022) §2.3 and Sossi et al. (2020) discuss why the *interface* fugacity (rather than the deep-mantle value) is the relevant choice; this assumption is consistent with CALLIOPE's ideal-gas, single-temperature treatment but breaks down if a Fe-FeO equilibrium curve in the deep mantle differs by more than a few dex.
- **No solid-FeO buffering**: when $\Phi_\mathrm{global} \to 0$, there is no melt to buffer $f_{\mathrm{O}_2}$ against the user-prescribed value. CALLIOPE keeps using the shift regardless of melt fraction, which is a reasonable bookkeeping choice but should not be over-interpreted physically.

## See also

- [Equilibrium chemistry](equilibrium_chemistry.md): the species-by-species speciation tree
- [Solubility laws](solubility.md): where $f_{\mathrm{O}_2}$ enters the S and N solubility paths
- [API reference for `calliope.oxygen_fugacity`](../Reference/api/calliope.oxygen_fugacity.md)
