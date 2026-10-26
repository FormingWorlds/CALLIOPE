# Metal reservoir

**The question.** How much of the volatile elements H, C, N and S goes into the liquid iron that forms a planet's core, and how does that change the atmosphere?

**The short answer.** When liquid metal is in equilibrium with the silicate melt, each element splits between the two in a ratio, the partition coefficient $D$, that depends on pressure, temperature, oxygen fugacity and composition. CALLIOPE adds the metal as a third reservoir next to the gas and the melt, so an element that goes into the metal is no longer available to the atmosphere. The values of $D$ come from laboratory fits in a partition file.

The metal reservoir is off by default. It is on when `ddict['metal_mass']` is above zero; `ddict` is the input dictionary of a solve (see [Configuration](../How-to/configuration.md)). The code names the metal reservoir `core` in its keys (`<e>_kg_core`, `w_<e>_core`); this page says *metal* for the reservoir, and *core* only for a planet's core or in a key name. The [metal reservoir tutorial](../Tutorials/metal_reservoir.md) runs it step by step, and the [partition file reference](../Reference/partition_file.md) lists every input key, every output key and the file format.

## Terms used on this page

- **Metal reservoir.** Liquid iron, with light elements dissolved in it, that is in chemical equilibrium with the silicate melt. It stands for the metal that sinks to form the core.
- **`metal_mass`** (kg), $M_\mathrm{metal}$. The mass of metal that is in equilibrium with the present melt. It is an input; CALLIOPE does not compute it.
- **Melt mass** (kg), $M_\mathrm{melt} = M_\mathrm{mantle}\,\Phi_\mathrm{global}$. The molten part of the silicate mantle: $M_\mathrm{mantle}$ (`M_mantle`, kg) is the mass of the silicate mantle, and $\Phi_\mathrm{global}$ (`Phi_global`, no unit) the fraction of it that is molten. A metal reservoir needs melt: with `Phi_global` = 0 the solve raises `ValueError`.
- **Mass fraction** $w_i$ (kg per kg). The mass of element $i$ in a reservoir divided by the mass of that reservoir. $w_i^\mathrm{melt}$ is the value in the melt, $w_i^\mathrm{metal}$ the value in the metal.
- **Partition coefficient on the mass basis** $D_i = w_i^\mathrm{metal} / w_i^\mathrm{melt}$ (no unit). $D_i = 10$ means that 1 kg of metal holds 10 times as much of element $i$ as 1 kg of melt.
- **Mole basis.** Some laboratory studies fit the ratio of mole fractions, $x_i^\mathrm{metal} / x_i^\mathrm{melt}$, instead of mass fractions. CALLIOPE converts such a $D$ to the mass basis before it uses it (see below). In the melt, CALLIOPE counts moles on the single-cation basis: one unit per cation (SiO$_2$, MgO, AlO$_{1.5}$, ...).
- **$x'$ basis.** The mole fraction of an element over the metal without H and C, $x'_j = N_j / \sum_{k \ne \mathrm{H,C}} N_k$, with $N$ the number of moles. [Tagawa et al. (2021)](https://doi.org/10.1038/s41467-021-22035-0) and [Tsutsumi et al. (2025)](https://doi.org/10.1038/s41467-025-65729-5) use this basis; Tsutsumi et al. (2025, Eq. 10, p. 8) give the reason: H and C atoms fill gaps between the Fe atoms of the metal.
- **$\Delta\mathrm{IW}$ of the metal-melt equilibrium**, `dIW_eq` (log10 units). The oxygen fugacity at which the metal and the melt equilibrate, relative to the iron-wüstite (IW) buffer. It is a separate quantity from `fO2_shift_IW`, the oxygen fugacity of the gas and the melt at the surface.
- **$K_D$, the exchange coefficient.** Some studies fit the equilibrium of a reaction with FeO, for example $\mathrm{FeO} = \mathrm{Fe} + \mathrm{O}$, instead of $D$ itself. $K_D$ is the ratio of the mole fractions of the products to those of the reactants. For O and Si the file stores $K_D$. For H, C, N and S it stores a law for $D$; when a study fits $K_D$, the file stores the converted law.
- **The H law as an example of a converted law.** Its source fits $\log_{10} K_D = \log_{10} D + \tfrac{1}{4}\Delta\mathrm{IW}$. The file stores the law for $\log_{10} D = \log_{10} K_D - \tfrac{1}{4}\Delta\mathrm{IW}$, so its $\Delta\mathrm{IW}$ coefficient is $-0.25$.
- **KLB-1.** A peridotite composition (Takahashi 1986, 7 major oxides in wt%), the default melt composition `silicate_oxides`.
- **Authoritative-O solve.** The solve that takes the O of the gas and the melt as a budget and finds the oxygen fugacity, `equilibrium_atmosphere_authoritative_O` (see [Authoritative-oxygen mode](authoritative_oxygen.md)).
- **Buffered solve.** The solve that takes the oxygen fugacity `fO2_shift_IW` as an input, `equilibrium_atmosphere`.
- **Metal share.** The part of the total mass of an element (gas, melt and metal together, the bulk inventory) that is in the metal: `<e>_kg_core` divided by `<e>_kg_total + <e>_kg_core`.
- **NBO/T** (`silicate_nbo_t`, no unit). The number of non-bridging O atoms per tetrahedral cation of the melt: a measure of how polymerised the melt is. Some laws have a term in it.
- **Fixed point.** A loop in which each pass uses the result of the pass before, until the result stops changing. CALLIOPE uses it because $D$ depends on the metal composition, and the metal composition depends on $D$.
- **`tol_w_core`** (default $10^{-6}$). The stop of that loop: the passes end when each metal mass fraction changes by at most `tol_w_core` of its value between 2 passes.

## The element balance with a metal term

Without metal, CALLIOPE finds the partial pressures at which the gas and the melt together hold the budget of each element (see [Mass balance & solver](mass_balance.md)). With metal, the metal holds a third part:

$$
m_i^\mathrm{atm} + m_i^\mathrm{melt} + D_i\,M_\mathrm{metal}\,w_i^\mathrm{melt} = m_i^\mathrm{target},
\qquad i \in \{\mathrm{H}, \mathrm{C}, \mathrm{N}, \mathrm{S}\}.
$$

Here $m_i^\mathrm{atm}$ and $m_i^\mathrm{melt}$ are the masses of element $i$ in the gas and in the melt, and $m_i^\mathrm{target}$ is its budget, all in kg. The third term is the mass of the element in the metal: the mass fraction of the element in the metal is $D_i$ times its mass fraction in the melt, and there are $M_\mathrm{metal}$ kg of metal.

In the code, `calliope.metal.metal_masses` computes this term. `solve._core_masses` calls it, and both residual functions, `solve.func` (buffered solve) and `solve.func_authoritative_O` (authoritative-O solve), add it to the residual of H, C, N and S. The O residual of the authoritative-O solve gets no metal term (see [O and Si](#o-and-si-exchange-with-the-silicate-feo)).

**What this does to the element.** Take an element that is only in the melt and in the metal, with no gas. The melt holds $w^\mathrm{melt} M_\mathrm{melt}$ of it and the metal $D\,w^\mathrm{melt} M_\mathrm{metal}$. With the mass ratio $r = M_\mathrm{metal}/M_\mathrm{melt}$, the metal then holds the fraction $f = D r/(1 + D r)$ of the element. For Earth's core and a fully molten Earth mantle, $r$ is about 0.48 (1.94 over 4.03). Then $D = 1$ already puts a third of the element into the metal, and $D = 100$ puts 98 % there.

**Outputs.** `<e>_kg_core` is the mass of the element named in `<e>` in the metal, and `w_<e>_core` its mass fraction there. `<e>_kg_total` stays the gas plus the melt. The bulk inventory of an element is `<e>_kg_total + <e>_kg_core`.

## The partition coefficient

Each element has one law in the partition file, of the form

$$
\log_{10} D = a + \frac{b}{T} + c\,\frac{P}{T} + d \cdot \mathrm{dIW\_factor} \cdot \Delta\mathrm{IW} + e\,(\mathrm{NBO/T}) + \sum_j \varepsilon_j \log_{10}(1 - X_j),
$$

with $T$ = `T_eq` in K, $P$ = `P_eq` in GPa and $\Delta\mathrm{IW}$ = `dIW_eq`. $a$, $b$, $c$, $d$ and $e$ are the fitted coefficients of the law, and `dIW_factor` is a fixed multiplier of the $\Delta\mathrm{IW}$ term (1 in the packaged file). $X_j$ is the $x'$ mole fraction of element $j$ in the metal, and $\varepsilon_j$ its fitted coefficient. The last term is the composition term: other light elements in the metal change how much of element $i$ it takes. `calliope.metal.log10_D` evaluates this law.

The field `basis` of each law says what the law gives: a mass-basis $D$, a mole-basis $D$, or a mole-basis $D$ on the $x'$ basis. `calliope.metal.mass_D` converts a mole-basis $D$ to the mass basis with the moles per kg of metal and of melt: $D_w = D\,N_\mathrm{metal}/N_\mathrm{melt}$. The melt moles come from `silicate_oxides` (the KLB-1 composition by default), plus the H and C dissolved in the melt. With a mole basis, $D_w$ therefore depends on the metal composition, which is one more reason for the fixed point.

## The oxygen fugacity of the metal-melt equilibrium

`dIW_eq_mode` selects how CALLIOPE sets $\Delta\mathrm{IW}$ for the laws:

- `'parameter'` (default): you give `dIW_eq` directly.
- `'from_FeO'`: CALLIOPE computes it from the FeO of the melt and the Fe of the metal, $\Delta\mathrm{IW} = 2\log_{10}(x_\mathrm{FeO}/x'_\mathrm{Fe})$, with ideal activities. With the KLB-1 melt (8.10 wt% FeO) and pure Fe metal this gives $-2.48$ (step 3 of the [metal reservoir tutorial](../Tutorials/metal_reservoir.md#step-3-the-2-modes-of-deltamathrmiw) shows it). More O or Si in the metal lowers $x'_\mathrm{Fe}$ and raises $\Delta\mathrm{IW}$.

`dIW_eq` and `fO2_shift_IW` describe 2 different equilibria. `fO2_shift_IW` sets the absolute oxygen fugacity of the gas and the melt at the surface, through the IW buffer of the solve (`fO2_model`, see [Oxygen fugacity](oxygen_fugacity.md)). `dIW_eq` enters only the partition laws and the O and Si exchange, and CALLIOPE never converts it to an absolute oxygen fugacity. The buffer reaches the metal only through the composition that the solve returns. In `'from_FeO'` mode, a different buffer can give a different metal, and with it a different `dIW_eq`. With a fixed metal composition, `dIW_eq` is the same on every buffer. The [partition file reference](../Reference/partition_file.md#inputs-in-ddict) states which laws change with the buffer in each mode.

## The fixed point on the metal composition

$D$ depends on the metal (through the composition terms, the mole basis and, in `'from_FeO'` mode, through $\Delta\mathrm{IW}$), and the metal depends on $D$. `solve._with_metal` solves this in passes:

1. Compute $\Delta\mathrm{IW}$ and $D$ of H, C, N and S from the metal and the melt of the previous pass (pure Fe at the start, plus the start O and Si).
2. Solve the gas, melt and metal balance at these fixed $D$.
3. Take the new metal H, C, N and S from the metal masses of that solve.
4. Stop when every metal mass fraction changes by at most `tol_w_core` of its value (plus $10^{-15}$, an absolute floor for fractions near 0); else go to step 1.

For O and Si in `'partition'` mode (see the next section), a pass can take less than the step that the exchange asks for. The stop test for these 2 elements uses the full step, not the part that was taken.

At most 30 passes run; then the solve raises `RuntimeError`. `tol_w_core` bounds the last step, not the error of the result. After the loop, a closure check confirms that gas, melt and metal hold the budget of each of H, C, N and S, within the tolerance of the solve.

## O and Si: exchange with the silicate FeO

O and Si enter the metal by reactions with the silicate:

$$
\mathrm{FeO} = \mathrm{Fe} + \mathrm{O}, \qquad \mathrm{SiO_2} + 2\,\mathrm{Fe} = \mathrm{Si} + 2\,\mathrm{FeO}.
$$

`O_mode_core` and `Si_mode_core` select how CALLIOPE sets them:

- `'fixed'` (default): the metal holds the mass fractions `w_O_core` and `w_Si_core` (default 0).
- `'partition'`: `calliope.metal.exchange_metal` computes the mole fractions from the exchange $K_D$ of the file and the $\Delta\mathrm{IW}$ of each pass, with ideal activities. Each pass takes the whole step, with 2 exceptions. Let $s$ be the secant slope of the exchange value against the metal value over the last 2 passes. When $s$ is negative (the exchange value falls as the metal value rises), the pass takes the part $1/(1 - s)$ of the step: with $\Delta\mathrm{IW}$ from FeO, more metal Si raises $\Delta\mathrm{IW}$ and lowers the Si that the exchange gives, so full steps would swing back and forth about the solution. And while the O and Si moves of a pass together take at least half of the iron left in the metal, they are halved.

The O of the metal comes from the FeO of the rock, not from the volatile O of the gas and the melt. So `O_kg_core` never enters the O balance: in the authoritative-O solve, the gas and melt O close on the O budget alone.

CALLIOPE does not hold the FeO of the silicate as a reservoir. It reports what the metal took instead. `n_FeO_exchange` (mol) is $n_\mathrm{O} - 2\,n_\mathrm{Si}$ of the metal: the first reaction takes 1 FeO per O atom, the second gives back 2 FeO per Si atom. A negative value means that the metal gives FeO back to the melt. `n_FeO_silicate` (mol) is the FeO of the melt. When the exchange is more than 10 % of the melt FeO, `FeO_exchange_out_of_range` is 1: the melt FeO and the `dIW_eq` of this solve do not follow the exchange, so the result is not self-consistent. The caller applies the exchange to its own FeO budget.

## Worked example

This example puts an Earth-like volatile budget into a fully molten Earth mantle and lets the whole Earth core equilibrate with it once. The inputs:

- Planet and melt as in the [First run tutorial](../Tutorials/firstrun.md): `M_mantle` $4.03 \times 10^{24}$ kg, `Phi_global` 1, `T_magma` 2500 K, `fO2_shift_IW` +0.5, the H of 1 Earth ocean ($1.55 \times 10^{20}$ kg), C/H 0.1 by mass, N 2 ppmw and S 200 ppmw (parts per million of the mantle mass).
- Metal: `metal_mass` $1.94 \times 10^{24}$ kg (`calliope.constants.M_core_earth`), `P_eq` 50 GPa, `T_eq` 3500 K, `dIW_eq` $-1$ in `'parameter'` mode, no O and Si.

**Why this metal mass.** `metal_mass` is the metal that equilibrates with the present melt. Here it is the whole core, in equilibrium with the whole mantle at one pressure and temperature. This is the single-stage model of Tsutsumi et al. (2025, p. 4). In their words, it assumes "that the entire core and mantle reached chemical equilibrium under a certain P, T, and fO2 condition". Their list of published single-stage conditions includes 50 GPa and 3500 K, the values used here. It gives the effect of a metal of Earth's core mass that meets the whole molten mantle at once, at these conditions. It is not a model of how a core grows: in a planet, metal arrives in portions over time, and each portion meets the melt of its own time.

**Why these conditions.** 50 GPa, 3500 K and $\Delta\mathrm{IW} = -1$ lie inside the data ranges of the placeholder H law ([Tagawa et al. 2021](https://doi.org/10.1038/s41467-021-22035-0): 30-60 GPa, 3080-4560 K, $\Delta\mathrm{IW}$ $-1.33$ to $-0.88$) and of the C law ([Fischer et al. 2020](https://doi.org/10.1073/pnas.1919930117): 1-59 GPa, 1618-5200 K, $\Delta\mathrm{IW}$ $-3.9$ to $-0.53$). The surface value `fO2_shift_IW` = +0.5 and `dIW_eq` = $-1$ differ on purpose: they describe the gas at the surface and the metal at depth.

The result is below (with `np.random.seed(0)`, which sets the random start guesses of the solver). The flag is `<e>_D_out_of_range`: it is 1 when the law of the element runs outside the pressure, temperature or $\Delta\mathrm{IW}$ range of the data it was fitted to; [Limits](#limits) says what the solve does then.

| Element | $D$ (mass basis) | Metal share | Out of range flag |
|---|---|---|---|
| H | 22.14 | 91.4 % | 0 |
| C | 0.0828 | 0.55 % | 0 |
| N | 39 | 41.7 % | 1 |
| S | 173 | 98.8 % | 1 |

The surface pressure falls from 8.7 bar without metal to 7.67 bar with metal. The atmosphere here is mostly C species, and almost no C goes into the metal at these conditions, so the pressure changes little. The H in the gas falls by a factor of about 140, from $2.26 \times 10^{17}$ kg to $1.58 \times 10^{15}$ kg.

**How to read the flags.** N and S have flag 1. Each of their placeholder laws is a straight line in $\Delta\mathrm{IW}$ through 2 values of one column of a table ([Suer et al. 2023](https://doi.org/10.3389/feart.2023.1159412), Table 1). For N this column is at 20 GPa and 3000 K, for S at 60 GPa and 4000 K. Each law is in range only at its one pressure and temperature, so any other `P_eq` or `T_eq` is flagged. No `P_eq` and `T_eq` put N and S in range together. The values of $D_\mathrm{N}$ = 39 and $D_\mathrm{S}$ = 173 here are the table values at $\Delta\mathrm{IW} = -1$, carried to 50 GPa and 3500 K without a pressure or temperature term. Use such a result as an estimate only.

**O and Si in `'partition'` mode.** The same case with `O_mode_core` and `Si_mode_core` set to `'partition'` gives 13.7 wt% O and 0.27 wt% Si in the metal. Both flags are 1: the placeholder O and Si laws have no known temperature and $\Delta\mathrm{IW}$ range, so every use is flagged.

The net exchange (O minus 2 Si) takes $1.62 \times 10^{25}$ mol of FeO, 3.56 times the $4.54 \times 10^{24}$ mol of FeO in the melt, so `FeO_exchange_out_of_range` is 1 as well: there is not enough FeO in the melt for this exchange.

The O also enters the C law through its composition term, and $D_\mathrm{C}$ rises from 0.0828 to 430; 97 % of the C then goes into the metal, and the surface pressure falls to 1.1 bar. This case shows why the flags matter: one flagged input (the metal O) changes a law that is in range (C) by a factor of about 5000.

Steps 1, 2 and 4 of the [metal reservoir tutorial](../Tutorials/metal_reservoir.md) run this case and print these numbers.

## Limits

- **Placeholder laws.** The packaged file has `placeholder = true`, and each solve with it logs one warning per process. The N and S laws have no pressure or temperature term; the O and Si laws have no known temperature and $\Delta\mathrm{IW}$ range. Replace a law with a fitted law from a source you have read: see [Add a new partition law](../Tutorials/partition_law.md).
- **Out-of-range flags.** `<e>_D_out_of_range` is 1 when the law of element `e` runs outside the pressure, temperature or $\Delta\mathrm{IW}$ range of its source data, or when a range is unknown. The solve still returns a result and logs a warning once per process; it does not clamp the law.
- **Ideal activities for O and Si.** `'partition'` mode leaves out the activity coefficients that published fits of these exchanges include. A solve in this mode logs a warning once per process. Treat its O and Si as a first estimate.
- **No FeO reservoir.** The exchange of O and Si does not change the melt FeO, the melt composition or `dIW_eq` of the solve. Above 10 % of the melt FeO, `FeO_exchange_out_of_range` is 1.
- **Cases that do not converge.** With $\Delta\mathrm{IW}$ from FeO and the input of the worked example below:
    - O and Si both in `'partition'` mode converge at 5500 K and raise `ValueError` at 6000 to 7500 K: the O and Si mole fractions of the exchange fill the whole metal (no Fe left), and their sum did not fall since the pass before.
    - O alone in `'partition'` mode converges up to 8480 K. From 8490 to 9700 K the passes end at the 30-pass cap with `RuntimeError`, and from 9720 K the solve raises the same `ValueError`. For O alone the exchange asks for $x_\mathrm{O} = K_D(\mathrm{O})\,x_\mathrm{FeO}/(1 - x_\mathrm{O})$, which has a root only while $K_D(\mathrm{O})\,x_\mathrm{FeO} \le 1/4$ (here up to 9470 K). The [partition file reference](../Reference/partition_file.md#inputs-in-ddict) gives the exact stop rules.
- **`metal_mass` is held fixed.** Fe that moves into or out of the metal through the O and Si exchange does not change `metal_mass`.
- **PROTEUS does not pass `metal_mass`.** A PROTEUS run has no metal reservoir. The reservoir is available in standalone CALLIOPE runs (see [Coupling to PROTEUS](proteus_coupling.md)).

**Where this is checked.** The [metal reservoir validation page](../Validation/metal.md) checks the metal term against mass balance, the full solve against a solve with a larger melt, the O bookkeeping of the authoritative-O solve, and the laws against published laboratory values.
