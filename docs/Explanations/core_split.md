# Core and mantle split

How much of the H, C, N and S of a planet goes into its core when the core forms?

**The short answer.** The split divides the bulk mass of each element between the core (the metal) and the mantle (the silicate). No mass is lost, and the ratio of the two mass fractions equals the partition coefficient $D$ of the element. $D$ says how strongly the element prefers the metal. CALLIOPE takes $D$ from a law at one pressure, one temperature and one redox state. $D$ depends on the metal composition, so the code repeats the calculation until the metal stops changing. In the tutorial planet the core takes 0.871 to 0.992 of the S, depending on the [redox state](#the-redox-state). The split runs on its own: PROTEUS does not call it.

To run the split, follow the [core split tutorial](../Tutorials/core_split.md). Every argument, output flag and error is in the [partition file reference](../Reference/partition_file.md#named-laws-and-the-core-split).

## What the split does

When a planet forms its core, liquid iron sinks through the magma ocean and takes part of the light elements with it. `calliope.core_split.split_core_mantle` divides the bulk mass of H, C, N and S between the core (the metal) and the mantle (the silicate). It does this once, at one pressure, one temperature and one redox state. The [metal reservoir](metal.md) answers another question: how a melt and a metal share the elements at one moment of a run, together with the atmosphere.

## The mass balance

The mass balance of each element $i$ uses 4 quantities. The split gets the first 3 as inputs and computes $D_i$ from the law of the element:

- $M_i$, the bulk mass of the element, in kg;
- $M_\mathrm{metal}$, the mass of the metal that equilibrates with the silicate, in kg;
- $M_\mathrm{silicate}$, the mass of the silicate, in kg;
- $D_i$, the partition coefficient on the mass basis: the mass fraction of $i$ in the metal divided by its mass fraction in the silicate.

Each element is either in the metal or in the silicate. Let $w_i$ be its mass fraction in the silicate. The metal then holds the mass fraction $D_i w_i$, and the bulk is $M_i = M_\mathrm{silicate} w_i + M_\mathrm{metal} D_i w_i$. Solved for $w_i$:

$$
w_i = \frac{M_i}{M_\mathrm{silicate} + D_i\,M_\mathrm{metal}}.
$$

The mantle holds $M_\mathrm{silicate} w_i$ and the core holds $D_i M_\mathrm{metal} w_i$. With $x = M_\mathrm{metal}/M_\mathrm{silicate}$, the core takes the fraction $D_i x/(1 + D_i x)$ of the bulk.

**Worked example.** The Earth has $M_\mathrm{metal} = 1.94 \times 10^{24}$ kg and $M_\mathrm{silicate} = 4.03 \times 10^{24}$ kg, so $x = 0.481$. An element with a constant $D = 100$ puts $100 \times 0.481 / (1 + 100 \times 0.481) = 0.9796$ of its bulk into the core (Step 4 of the [core split tutorial](../Tutorials/core_split.md#step-4-your-own-s-law)). With the S law of the code, $D_\mathrm{S}$ depends on the conditions, so the core fraction does too; see [the redox state](#the-redox-state).

After each split, the code checks for each element that core plus mantle equals the bulk, to $10^{-12}$ of the bulk. It raises `RuntimeError` when this element ledger does not close.

## Why the split must iterate

$D_i$ is not a constant. It depends on what else is in the metal: the $\varepsilon_j$ terms of a law (one fitted coefficient for each other light element $j$ in the metal, see [The partition coefficient](metal.md#the-partition-coefficient)) change $D_i$ with the amount of the other light elements, and a law on the mole basis is converted to the mass basis with the moles per kg of the metal. The metal composition in turn depends on how much of each element the split puts into the metal. So the split repeats 3 steps until the metal stops changing: compute $D_i$ for the present metal, split the bulk, and compute the new metal from the split. One round of these 3 steps is one pass. This is a fixed point on the metal composition: a loop in which each pass uses the result of the pass before, until the result stops changing. The metal reservoir uses the same loop.

## The inputs

The split takes 5 kinds of input besides the masses.

### The equilibration point

The equilibration point $(P_\mathrm{eq}, T_\mathrm{eq})$ is the pressure and temperature at which the metal and the silicate equilibrate. The default is 60 GPa and 4000 K, the deep magma ocean column of Suer et al. (2023)[^cite-suer2023], Table 1. The S law is inside its pressure and temperature range there. The core-mantle boundary is not a good default. The Earth's core-mantle boundary is at about 136 GPa, above the data of every packaged law (at most 91 GPa), and larger planets have higher pressures there.

### The redox state

$\Delta\mathrm{IW}$ is the oxygen fugacity relative to the iron-wüstite buffer in log10 units. It has no default, and it is the one redox input of the split. The H, C and N laws read it through their $\Delta\mathrm{IW}$ term. The S law has no such term. Its redox dependence is the FeO of the melt: S enters the metal as FeS in exchange for the FeO of the melt, so $D_\mathrm{S}$ rises with the FeO. The split computes the FeO mole fraction of the melt, $X_\mathrm{FeO}$, from $\Delta\mathrm{IW}$ and the iron of the metal in each pass. It uses the definition of Suer et al. (2017, p. 90)[^cite-suer2017] for ideal mixing, $\Delta\mathrm{IW} = 2\log_{10}(X_\mathrm{FeO}/X_\mathrm{Fe})$:

$$
X_\mathrm{FeO} = x'_\mathrm{Fe}\,10^{\Delta\mathrm{IW}/2},
$$

with $x'_\mathrm{Fe}$ the Fe mole fraction over the metal without H and C. A more oxidised equilibrium therefore puts more S into the core. For the tutorial planet (the Earth's core and mantle masses, 2 wt% O and 3 wt% Si in the metal) the core takes 0.871, 0.974 and 0.992 of the S at $\Delta\mathrm{IW}$ = -3, -1.5 and -0.5.

$X_\mathrm{FeO}$ cannot be larger than 1. Above $\Delta\mathrm{IW} = -2\log_{10} x'_\mathrm{Fe}$ no melt fits, and a law with an FeO term raises `ValueError` that names $\Delta\mathrm{IW}$ and this limit. The first pass uses the metal without the split elements, which has the most Fe. The limit of the split is therefore the limit of the start metal `w_fixed`: +0.06 for 2 wt% O, and +0.11 for the tutorial metal with 2 wt% O and 3 wt% Si.

### The fixed metal elements

`w_fixed` holds the metal mass fractions of O and Si. The split does not change them (see [Limits](#limits)).

### The melt

`oxides` and `nbo_t` describe the melt. A law on the [mole basis](metal.md#terms-used-on-this-page) gives $D$ as a ratio of mole fractions. To turn it into a mass-basis $D$, the split multiplies by the moles of metal per kg and divides by the moles of melt per kg. The moles of melt per kg are one mole for each mole of cations of the oxides (the single-cation basis), plus the dissolved H and C. The melt composition `oxides` sets this number. `oxides` is also the melt that a law written as a Python function gets; it has no other role in the split. `nbo_t` is the NBO/T of the melt, the number of non-bridging O atoms per tetrahedral cation, a measure of how polymerised the melt is: it enters the NBO/T term of a law (the C law has one), and a function law gets it too. A function law that needs the melt FeO computes it from $\Delta\mathrm{IW}$ and the metal.

### The laws

The split takes one law per element, either by name from `calliope.metal.get_law` or as a loaded law entry. The defaults are `tagawa2021_eq2` (H), `fischer2020_eq1` (C), `suer2023_line` (N) and `suer2017_eq7` (S). The S law of the metal reservoir, `suer2023_line`, is a line in $\Delta\mathrm{IW}$ at one pressure and one temperature. The split uses the law of Suer et al. (2017)[^cite-suer2017] instead, because it has a pressure and a temperature term. To add your own law, see [Add a partition law to the core split](../How-to/core_split_law.md).

## What the output flags mean

The split evaluates each law at the given point. It does not move the point to the edge of the range of the law, because that would hide the temperature dependence that dominates $D_\mathrm{S}$, as Suer et al. (2017, p. 91)[^cite-suer2017] show. It flags the law instead, and a law outside its range logs one warning per process. `flags[e]` holds the flags of element `e`, and `flags['metal']` holds those of the metal:

| Flag | Meaning | Unit |
|---|---|---|
| `flags[e]['out_of_range']` | 1 when any of the 3 distances below is not 0, else 0 | none |
| `flags[e]['P_GPa']`, `['T_K']`, `['dIW']` | distance of the point outside each range of the law; 0 inside the range, NaN for an unknown range (`calliope.metal.range_distance`) | GPa, K, log10 units |
| `flags[e]['D_residual']` | relative difference between the metal mass fraction of `e` and $D_e$ times its silicate mass fraction, at the returned split | none |
| `flags['metal']['x_light']`, `['x_S']` | mole fraction of all light elements, and of S, over the whole metal | none |
| `flags['metal']['S_over_Fe']` | S to Fe mole ratio of the metal | none |
| `flags['metal']['beyond_FeS']` | 1 when the metal holds more S than Fe in moles, past FeS; with other light elements in the metal this happens below $x_\mathrm{S} = 0.5$. A flag only: the split does not clamp the metal | none |
| `flags['metal']['root_rule']` | 1 when the iron-rich root rule (next section) replaced the first root, else 0 | none |

## How the solver finds the split

The split starts from the metal of `w_fixed` and runs one pass after another until the metal stops changing. If the passes do not settle, it repeats them with smaller steps.

Some S-rich, iron-poor metals have 2 roots: 2 metal compositions that both close the mass balance, one with more Fe than the other. When the smaller steps end on a metal with more S than Fe, the split searches again with very small steps from the start metal and keeps the root with more Fe; `root_rule` is 1 when this changes the result.

The [Core split solver](../Reference/core_split_solver.md) page lists the steps, their tolerances and the errors they raise.

## How reliable it is

A sweep of 5500 random splits, from $10^{20}$ to $3 \times 10^{24}$ kg of metal, converges in 5074 cases and stops with a named error in the others. Where the metal keeps more than 1 % of Fe by mass, the returned metal meets each law to $2.6 \times 10^{-10}$. The [sweep of random inputs](../Validation/metal.md#core-and-mantle-split-over-random-inputs) gives the counts per input range and the known limits.

## Limits

Each limit below says what the code does when a value is outside it.

- **$\Delta\mathrm{IW}$ above $-2\log_{10} x'_\mathrm{Fe}$.** A law with an FeO term raises `ValueError` (see [the redox state](#the-redox-state)).
- **A law outside the range of its data.** The split flags it and logs one warning per process. It does not move the point.
- **A metal with no Fe left, or passes that do not settle.** The split raises `ValueError` when the iron runs out and `RuntimeError` otherwise. In the sweep of random inputs, 416 of 5500 splits raise `ValueError` and 10 raise `RuntimeError`.
- **An element ledger that does not close** to $10^{-12}$ of the bulk. The split raises `RuntimeError`.
- **O and Si are not split.** Their bulk is the O and Si of the rock, and CALLIOPE does not hold the rock. `split_core_mantle` refuses O and Si in its bulk and takes their metal mass fractions as fixed values (`w_fixed`). The code has no split of O and Si by their exchange $K_D$.
- **The melt basis of the S law is an assumption.** Suer et al. (2017)[^cite-suer2017] do not state the mole basis of the melt for $X_\mathrm{FeO}$ and $C_\mathrm{S}$. CALLIOPE takes $X_\mathrm{FeO}$ from $\Delta\mathrm{IW}$, uses $\log_{10} C_\mathrm{S} = -5.372$ as given and assumes the single-cation basis. The [Suer et al. (2017) S runs](../Validation/metal.md#suer-et-al-2017-s-runs) case shows how well this reproduces their Table 1.
- **One $C_\mathrm{S}$ for every melt.** The S law uses $\log_{10} C_\mathrm{S} = -5.372$ of a pyrolite melt, from Fig. 6 of Suer et al. (2017)[^cite-suer2017], for every melt composition.

**Where this is checked.** The [Core and mantle split](../Validation/metal.md#core-and-mantle-split) case checks the mass balance against its closed form, the [Suer et al. (2017) S runs](../Validation/metal.md#suer-et-al-2017-s-runs) case checks the S law, and the [sweep of random inputs](../Validation/metal.md#core-and-mantle-split-over-random-inputs) shows where the solver converges.

[^cite-suer2023]: T.-A. Suer, C. Jackson, D. S. Grewal, C. Dalou, T. Lichtenberg, *[The distribution of volatile elements during rocky planet formation](https://doi.org/10.3389/feart.2023.1159412)*, Frontiers in Earth Science, 11, 1159412, 2023. [SciX](https://scixplorer.org/abs/2023FrEaS..1159412S/abstract).
[^cite-suer2017]: T.-A. Suer, J. Siebert, L. Remusat, N. Menguy, G. Fiquet, *[A sulfur-poor terrestrial core inferred from metal-silicate partitioning experiments](https://doi.org/10.1016/j.epsl.2017.04.016)*, Earth and Planetary Science Letters, 469, 84-97, 2017. [SciX](https://scixplorer.org/abs/2017E%26PSL.469...84S/abstract).
