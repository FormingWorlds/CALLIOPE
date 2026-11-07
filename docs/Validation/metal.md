# Validation: metal reservoir

## What this page checks

This page checks how CALLIOPE splits the volatile elements H, C, N and S between gas, silicate melt and metal.
The metal reservoir stands for the iron that forms a planet's core.
Validation here means that we run the code on a case where another source gives the right answer.
The test then checks that the code gives that answer within a stated tolerance.
The right answer comes from mass balance alone, from a second CALLIOPE solve that must agree, or from laboratory measurements and tables in published papers.
If this part of the code is wrong, the metal takes the wrong amount of each element.
The atmosphere of the model then gets the wrong amount too.
Each case below says why the check exists, what it tests, how we know the right answer and why it matters.
The box "Test details" at the end of each case gives the reference, the tolerance, the test id and the command.

**How the metal enters the calculation.**
A metal reservoir is on when `metal_mass` is above zero.
The partition coefficient $D_i$ of element $i$ is the ratio of its mass fraction in the metal to its mass fraction in the melt.
The element balance of H, C, N and S then gains a core term $D_i M_\mathrm{metal} w_i^\mathrm{melt}$.
Here $M_\mathrm{metal}$ is the metal mass and $w_i^\mathrm{melt}$ the mass fraction of the element in the melt.
The values of $D_i$ come from the partition file (`src/calliope/metal.py`); the packaged laws are placeholders.
The metal O comes from the FeO of the silicate and stays out of the volatile O budget.
`<e>_kg_total` stays gas plus melt and excludes the core; the bulk inventory of an element is `<e>_kg_total + <e>_kg_core`.
Inputs, outputs and the placeholder laws are on the [metal reservoir page](../Reference/partition_file.md).

## Gas-free core fraction

**Why this check exists.** When metal takes an element from the melt, the melt has less of it left. A code that forgets this loss puts too much of the element into the core. For a large $D$ it can put more than the whole budget there.

**What it tests.** A melt and a metal with no gas. The metal to melt mass ratio is $x = M_\mathrm{metal}/M_\mathrm{melt}$, here 0.325/0.675 as for the Earth's core and mantle. The core must then hold the fraction $f = Dx/(1 + Dx)$ of the element.

**How we know the right answer.** This result follows from mass balance alone. The element is either in the melt or in the metal, and $D$ sets their ratio.

**Why it matters.** The core fraction sets how much H, C, N and S stays in the silicate and can reach the atmosphere.

![Core fraction of an element against D](../assets/validation/metal_gas_free_core_fraction.svg#only-light){ width="80%" }
![Core fraction of an element against D](../assets/validation/metal_gas_free_core_fraction_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The x-axis is $D$ on a log scale and the y-axis is the fraction of the element in the core. The solid line is $Dx/(1 + Dx)$, the dashed line is $Dx$ without the loss from the melt, and the points are the values of the code. A pass means that the points lie on the solid line.

??? info "Test details"
    - Reference: the mass balance of melt and metal (analytical).
    - Tolerance: relative 1e-12; $f = Dx$, the law without the melt back-reaction, lies more than 10 % above it from $D = 0.3$ on.
    - Test id: `tests/test_metal.py::test_gas_free_core_fraction_curve`, `tests/test_metal.py::TestMetalReservoir` (unit tier)
    - Command: `pytest tests/test_metal.py -k "gas_free"`

## Equal D is a larger melt

**Why this check exists.** The metal term enters the full solve of gas, melt and metal. A wrong factor in that term, or an element left out of it, changes the partial pressures of the gas.

**What it tests.** H, C, N and S all get the same $D$, the metal holds no O or Si, and $\Delta\mathrm{IW}$ is a fixed parameter. $\Delta\mathrm{IW}$ is the oxygen fugacity of the metal-melt equilibrium in log10 units relative to the iron-wüstite (IW) buffer. In this case the metal acts exactly like extra melt of mass $D M_\mathrm{metal}$. So the solve with metal must give the same result as the solve without metal and a melt that is larger by $D M_\mathrm{metal}$. The test uses $D$ = 1 and 30. The gas holds the rest of each budget, so the core mass is $(M - M_\mathrm{gas})\,Dx/(1 + Dx)$. Here $M$ is the budget of the element and $M_\mathrm{gas}$ its mass in the gas.

**How we know the right answer.** The two solves have the same element balance, term by term, so they must give the same answer.

**Why it matters.** This case checks the metal term inside the full solve, with the gas, the solubility laws and the chemistry all active.

![Partial pressures with metal against a melt larger by D M_metal, buffered at D = 1 and 30 and with an O budget at D = 30](../assets/validation/metal_equal_d_larger_melt.svg#only-light){ width="80%" }
![Partial pressures with metal against a melt larger by D M_metal, buffered at D = 1 and 30 and with an O budget at D = 30](../assets/validation/metal_equal_d_larger_melt_dark.svg#only-dark){ width="80%" }

**How to read the plot.** Each point is the partial pressure of one gas species. The x-axis is the solve with the larger melt and the y-axis is the solve with metal, both on log scales. The 3 marker types are $D$ = 1, $D$ = 30, and $D$ = 30 with an O budget. The last one is the case [Equal D is a larger melt with an O budget](#equal-d-is-a-larger-melt-with-an-o-budget). A pass means that every point lies on the dashed 1:1 line.

??? info "Test details"
    - Reference: the identity of the element balance (analytical).
    - Tolerance: partial pressures and melt plus core masses relative 1e-6, with both solves from the same random state and the root accepted at relative 1e-10 plus 1 kg; core over gas-free core equals $1 - M_\mathrm{gas}/M$ to 2e-5.
    - Test id: `tests/test_metal_solve.py::test_metal_with_equal_d_is_a_larger_melt`, `tests/test_metal_solve.py::test_equal_d_identity_plot` (smoke tier)
    - Command: `pytest tests/test_metal_solve.py -k "larger_melt or identity_plot"`

## O budget with metal O

**Why this check exists.** The metal O comes from the FeO of the rock, not from the volatile O of the gas and the melt. A code that counts the metal O as volatile O gets the oxygen fugacity of the planet wrong.

**What it tests.** CALLIOPE can solve in 2 ways. The buffered solve takes the oxygen fugacity as an input, written as an offset from the IW buffer. The authoritative-O solve takes the total O of the gas and the melt as an input and finds the oxygen fugacity. First, a buffered solve at IW+2 with 1.9e24 kg of metal gives the O of its gas and melt. Its metal has $D$ = 30 for H, C, N and S and a metal-melt $\Delta\mathrm{IW}$ of -2. Then an authoritative-O solve gets that O as its budget and 0.1 wt% O in the metal. It must return IW+2, and the O of its gas and melt must equal the budget.

**How we know the right answer.** The second solve must recover the input of the first one (a round trip). The metal O (1.9e21 kg) is 22 times the gas and melt O (8.76e19 kg). If the code counted it in the budget, the result would be far from IW+2.

**Why it matters.** PROTEUS can use the authoritative-O solve to follow the oxygen of a planet over time (`planet.fO2_source = "from_O_budget"`). The metal must not change that bookkeeping.

??? info "Test details"
    - Reference: round trip between the buffered and the authoritative-O solve (self-consistency).
    - Tolerance: derived $\Delta\mathrm{IW}$ +2 to 1e-6; partial pressures relative 1e-6; gas and melt O equal the budget, and gas, melt and core of each of H, C, N and S equal theirs, to relative 1e-10 plus 10 kg. The metal O (1.9e21 kg) is 22 times the gas and melt O (8.76e19 kg), so charging it to the budget fails the test; partition-mode O (3.6 wt% in the metal, more than 100 times the budget) solves as well.
    - Test id: `tests/test_metal_solve.py::test_authoritative_o_budget_excludes_the_metal_o`, `tests/test_metal_solve.py::test_partition_o_far_above_the_o_budget_solves` (smoke tier)
    - Command: `pytest tests/test_metal_solve.py -k "excludes_the_metal_o or partition_o_far"`

## Equal D is a larger melt with an O budget

**Why this check exists.** The authoritative-O solve has its own path through the code. If that path leaves out the metal term or gives it a wrong factor, its result differs from the solve with the larger melt.

**What it tests.** The same comparison as in [Equal D is a larger melt](#equal-d-is-a-larger-melt), now with the authoritative-O solve, $D$ = 30 and no O in the metal. The larger melt also holds the O of its dissolved species. So its O budget is larger by $D M_\mathrm{metal}$ times the O mass fraction of the melt in the solve with metal.

**How we know the right answer.** As before, the two solves have the same element balance. So they must give the same partial pressures and the same oxygen fugacity.

**Why it matters.** With this case, both solve modes that PROTEUS can use are checked with a metal reservoir.

??? info "Test details"
    - Reference: the identity of the element balance (analytical).
    - Tolerance: partial pressures and melt plus core masses relative 1e-6, derived $\Delta\mathrm{IW}$ to 1e-6, roots accepted at relative 1e-10.
    - Test id: `tests/test_metal_solve.py::test_authoritative_o_with_equal_d_is_a_larger_melt` (smoke tier)
    - Command: `pytest tests/test_metal_solve.py -k authoritative_o_with_equal_d`

## Tagawa et al. (2021) H runs

**Why this check exists.** The placeholder law for H is a fit to laboratory data, and the paper gives more than one fit. Using the wrong one, for example their Eq. 10 instead of Eq. 2, gives a wrong $D$ for H with no other sign.

**What it tests.** The H law is Eq. 2 of [Tagawa et al. (2021)](https://doi.org/10.1038/s41467-021-22035-0). It gives $K_D$, the exchange coefficient of H between metal and silicate, from pressure and temperature. The test compares it with the 6 runs of their SI Table 1. Each run is an experiment in a laser-heated diamond anvil cell, with a measured $\log_{10} K_D$ and its 1 sigma error.

**How we know the right answer.** Tagawa et al. (2021) measured these 6 values in the laboratory.

**Why it matters.** The H law sets how much hydrogen the metal takes. Hydrogen in the metal is H that the atmosphere and the water inventory do not get.

![Measured against fitted log10 K_D of the Tagawa runs](../assets/validation/metal_tagawa2021_h_runs.svg#only-light){ width="80%" }
![Measured against fitted log10 K_D of the Tagawa runs](../assets/validation/metal_tagawa2021_h_runs_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The x-axis is $\log_{10} K_D$ from Eq. 2 and the y-axis is the measured value, with its 1 sigma error bar. A pass means that every error bar crosses the dashed 1:1 line.

??? info "Test details"
    - Reference: Tagawa et al. (2021), SI Table 1 (published table).
    - Tolerance: each run inside its 1 sigma (largest 0.61 sigma); their Eq. 10 misses every run by more than 2 sigma.
    - Test id: `tests/test_metal.py::test_tagawa_runs_lie_within_one_sigma_of_eq2` (unit tier)
    - Command: `pytest tests/test_metal.py -k tagawa`

## Tsutsumi et al. (2025) single-stage cores

**Why this check exists.** Papers give $D$ in different forms, for example per mole and over a metal without H and C. CALLIOPE must convert each form to its own mass basis. A wrong conversion changes the H and C content of the core.

**What it tests.** [Tsutsumi et al. (2025)](https://doi.org/10.1038/s41467-025-65729-5) give $D$ per mole over the metal without H and C. CALLIOPE converts it with the KLB-1 melt composition of Takahashi (1986). The test uses their 3 single-stage cases: 40 GPa and 3750 K, 50 GPa and 3500 K, 54 GPa and 3350 K. In each case the bulk silicate Earth holds 710 ppm H$_2$O and 120 ppm C. The code must give the core H and C of their SI Table 1.

**How we know the right answer.** Tsutsumi et al. (2025) list the core H and C of each case in their SI Table 1.

**Why it matters.** This case checks the conversion of $D$ with a full set of published inputs and outputs. The same conversion applies to each law in the partition file that gives $D$ per mole, such as the H law.

![Core H and C of the Tsutsumi single-stage cases](../assets/validation/metal_tsutsumi2025_core_h_c.svg#only-light){ width="80%" }
![Core H and C of the Tsutsumi single-stage cases](../assets/validation/metal_tsutsumi2025_core_h_c_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The x-axis has the 3 cases and the y-axis is the mass fraction of H and C in the core in wt%. The open markers are the table values and the filled markers are the values of the code. A pass means that each filled marker lies inside its open marker.

??? info "Test details"
    - Reference: Tsutsumi et al. (2025), SI Table 1 (published table).
    - Tolerance: core H and C equal to the table at its rounding, 0.005 wt%; D times the melt mass fraction without the basis conversion misses the 40 GPa case.
    - Test id: `tests/test_metal.py::TestBasis`, `tests/test_metal.py::test_tsutsumi_core_values_plot` (unit tier)
    - Command: `pytest tests/test_metal.py -k "tsutsumi or Basis"`

## Published law values

**Why this check exists.** Each placeholder law is a formula with coefficients from a paper. A wrong sign or a wrong coefficient gives a law that looks reasonable but is not the published one.

**What it tests.** Each law is evaluated at points where its source prints the result:

- Tagawa et al. (2021) Eq. 2 at 3 pressure and temperature points.
- Fischer et al. (2020) Eq. 1, $D_\mathrm{C}$ = 0.550 for their single-stage case (54 GPa, 4100 K, NBO/T 2.6, IW-2.2, $X_\mathrm{O}$ 0.06, $X_\mathrm{S}$ 0.03; p. 8745). NBO/T is the number of non-bridging oxygens per tetrahedral cation, a measure of how polymerised the melt is. $X_\mathrm{O}$ and $X_\mathrm{S}$ are the mole fractions of O and S in the metal.
- The lines through the values of Suer et al. (2023) Table 1: $D_\mathrm{N}$ 1 and 39 and $D_\mathrm{S}$ 10 and 173 at dIW -4 and -1.
- The ideal exchange of O and Si between metal and melt with the $K_D$ of Tagawa et al. (2021) SI Table 3.

**How we know the right answer.** The papers print these values.

**Why it matters.** These laws set $D$ for every element in the metal reservoir, so every other case on this page depends on them.

??? info "Test details"
    - Reference: the published equations and tables (published values).
    - Tolerance: Eq. 2 at its 4 printed decimals; $D_\mathrm{C}$ 5e-4; the Suer lines relative 1e-5; the O exchange relative 1e-12 and the Si exchange relative 1e-3 (the KLB-1 single-cation units by hand).
    - Test id: `tests/test_metal.py::TestLaw::test_tagawa_eq2_values`, `tests/test_metal.py::TestLaw::test_fischer_eq1_single_stage_value`, `tests/test_metal.py::TestLaw::test_suer_lines`, `tests/test_metal.py::TestExchange::test_ideal_exchange_mole_fractions` (unit tier)
    - Command: `pytest tests/test_metal.py -k "eq2_values or eq1_single or suer_lines or ideal_exchange"`

## Suer et al. (2017) S runs

**Why this check exists.** The core and mantle split takes S from the law of [Suer et al. (2017)](https://doi.org/10.1016/j.epsl.2017.04.016). A wrong term in that law, for example a natural log where the paper means log10, moves $D_\mathrm{S}$ by more than a factor of 10.

**What it tests.** The law is evaluated at the 9 diamond anvil cell runs of their Table 1 (46 to 91 GPa, 3100 to 4100 K) and compared with the measured $D_\mathrm{S}$. The law is their Eq. 7 with log10 in every term and the coefficients of their Table 3:

$$
\log_{10} D_\mathrm{S} = -3.50 + \frac{3000}{T} + 33\,\frac{P}{T} + \log_{10} X_\mathrm{FeO} - \log_{10} C_\mathrm{S} + 14 \log_{10}(1 - X_\mathrm{O}),
$$

with $T$ in K, $P$ in GPa, $X_\mathrm{FeO}$ the mole fraction of FeO in the silicate melt, $X_\mathrm{O}$ the mole fraction of O in the metal and $C_\mathrm{S}$ the sulfide capacity of the melt. $D_\mathrm{S}$ is a ratio of mole fractions (their Eq. 1). Each run uses its own $X_\mathrm{FeO}$ and $X_\mathrm{O}$ from Table 1, and $\log_{10} C_\mathrm{S} = -5.372$, the value of their Fig. 6 for a pyrolite mantle. Table 1 does not list the $C_\mathrm{S}$ of each run.

For run 7, Table 1 (p. 91) prints $X_\mathrm{FeO}$ = 0.02, which we read as a misprint of 0.20. The same row gives IW -1.06, and with $\mathrm{IW} = 2\log_{10}(X_\mathrm{FeO}/X_\mathrm{Fe})$ and $X_\mathrm{Fe}$ = 0.68 that needs $X_\mathrm{FeO}$ = 0.2007. The FeO of 25.62 wt% in their Table 2 (p. 92) gives 0.204 on the single-cation basis. The test uses 0.20 for run 7. For the other 8 runs the IW column and $X_\mathrm{FeO}$ agree to 0.005.

**How we know the right answer.** The measured $D_\mathrm{S}$ of each run is in Table 1. The tolerance is 0.47 dex, the 1-sigma uncertainty that Table 3 gives for the coefficient $a$.

**Why it matters.** The S law sets how much of the S goes into the core and how much stays in the mantle, where it can later reach the atmosphere.

![Measured log10 D_S of the 9 runs of Suer et al. (2017) against Eq. 7](../assets/validation/metal_suer2017_s_runs.svg#only-light){ width="80%" }
![Measured log10 D_S of the 9 runs of Suer et al. (2017) against Eq. 7](../assets/validation/metal_suer2017_s_runs_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The x-axis is $\log_{10} D_\mathrm{S}$ from Eq. 7 and the y-axis the measured value. The dashed line is 1:1 and the band is $\pm$0.47 dex around it. Filled points are runs within the band; the open point is run 7. A pass means that 8 runs lie in the band and run 7 lies outside it.

The residual (model minus measurement) is +0.30, -0.12, +0.42, +0.42, +0.03, -0.03, +0.50, +0.15 and +0.38 dex for runs 1 to 9, with a mean of +0.23 dex. Run 7 is a stated miss at +0.495 dex, 0.025 dex outside the band; with the IW value 0.2007 it would be +0.497 dex. With the printed $X_\mathrm{FeO}$ of 0.02 it would miss by 0.035 dex on the other side (-0.505 dex). A possible cause of the mean offset is the $C_\mathrm{S}$ value: it is for a pyrolite melt, while the runs use MORB and KLB-1 melts (Table 1). With $a = -3.30$, the value printed in their Eq. 8, every residual is 0.2 dex larger and 5 runs miss.

??? info "Test details"
    - Reference: Suer et al. (2017), Eq. 7, Table 1 (runs), Table 3 (coefficients), Fig. 6 ($\log_{10} C_\mathrm{S}$).
    - Tolerance: 0.47 dex per run for 8 runs, run 7 a stated miss pinned at +0.495 dex (5e-4); the law at the default point relative 1e-12; Eq. 8 read with natural logs lies more than 1 dex away.
    - Test id: `tests/test_metal.py::TestSuer2017::test_runs_of_table_1_against_eq7`, `tests/test_metal.py::TestSuer2017::test_eq7_value_at_the_default_point` (unit tier)
    - Command: `pytest tests/test_metal.py -k TestSuer2017`

## Core and mantle split

**Why this check exists.** `calliope.core_split.split_core_mantle` puts the bulk mass of each element into the core or the mantle. A wrong factor in the split, or a mantle mass that is not computed from the same melt fraction as the core mass, creates or destroys mass.

**What it tests.** With a constant mass-basis $D$ and $x = M_\mathrm{metal}/M_\mathrm{silicate}$ (here 1.94/4.03, the masses of Earth's core and mantle), the core must hold the fraction $Dx/(1 + Dx)$ of the bulk, for $D$ from 0.1 to 1000. A second test checks the element ledger, bulk = core + mantle with both parts at 0 or above, for random $D$ from $10^{-6}$ to $10^6$ and masses from $10^{15}$ to $10^{24}$ kg. The split itself checks the same ledger on every call and raises `RuntimeError` when it fails.

**How we know the right answer.** Mass balance: an element is either in the metal or in the silicate, and $D$ sets the ratio of its mass fractions.

**Why it matters.** The split sets the light-element mass of the core and the volatile mass left in the mantle, which can later reach the atmosphere.

![Core fraction of an element against D](../assets/validation/core_split_fraction.svg#only-light){ width="80%" }
![Core fraction of an element against D](../assets/validation/core_split_fraction_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The x-axis is $D$ on a log scale and the y-axis the core fraction of the bulk. The line is $Dx/(1 + Dx)$ and the points are the values of `split_core_mantle`. A pass means that the points lie on the line.

??? info "Test details"
    - Reference: the mass balance of metal and silicate (analytical).
    - Tolerance: relative 1e-12 for the core fraction and the ledger; $f = Dx$ without the back-reaction lies more than 10 % above it from $D = 1$ on.
    - Test id: `tests/test_core_split.py::test_constant_d_gives_the_closed_form_core_fraction`, `tests/test_core_split.py::test_ledger_closes_over_the_domain` (unit tier)
    - Command: `pytest tests/test_core_split.py -k "closed_form or ledger_closes"`
