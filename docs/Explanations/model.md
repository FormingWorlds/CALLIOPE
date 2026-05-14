# Model overview

CALLIOPE is a 0-D **equilibrium outgassing** solver for the magma-ocean atmosphere coupling. It treats the silicate mantle and the overlying gas-phase atmosphere as a single thermodynamic system in equilibrium at the surface, and asks: for a given total elemental inventory and a given magma-ocean state, what surface partial pressures and dissolved-volatile masses simultaneously satisfy (i) gas-phase chemical equilibrium, (ii) gas-melt solubility equilibrium, and (iii) elemental mass conservation?

This page summarises the model assumptions, the variables it solves for, and how it relates to the upstream papers ([Bower et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...631A.103B), [2022](https://ui.adsabs.harvard.edu/abs/2022PSJ.....3...93B); [Nicholls et al. 2024](https://ui.adsabs.harvard.edu/abs/2024JGRE..12908576N)). Each component has its own dedicated page.

## What is in the model

- **Eleven gas-phase species** (`calliope.constants.volatile_species`): H$_2$O, CO$_2$, O$_2$, H$_2$, CH$_4$, CO, N$_2$, S$_2$, SO$_2$, H$_2$S, NH$_3$.
- **Five elements**: H, C, N, S always solved. O is either a derived quantity set by the oxygen-fugacity buffer (buffered mode, `equilibrium_atmosphere`) or a fifth budgeted element with $\Delta\mathrm{IW}$ as the additional unknown (authoritative-O mode, `equilibrium_atmosphere_authoritative_O`). The two modes share all physics functions and differ only in their unknown set; see [Authoritative-oxygen mode](authoritative_oxygen.md).
- **Six equilibrium reactions** (gas-phase, surface temperature):

    | Reaction | Source |
    |---|---|
    | $\mathrm{H_2O} \rightleftharpoons \mathrm{H_2} + \tfrac{1}{2}\,\mathrm{O_2}$ | [JANAF](https://janaf.nist.gov/) (`janaf_H2`) and [Schaefer & Fegley 2017](https://ui.adsabs.harvard.edu/abs/2017ApJ...843..120S) (`schaefer_H`) |
    | $\mathrm{CO_2} \rightleftharpoons \mathrm{CO} + \tfrac{1}{2}\,\mathrm{O_2}$ | [JANAF](https://janaf.nist.gov/) (`janaf_CO`) and [Schaefer & Fegley 2017](https://ui.adsabs.harvard.edu/abs/2017ApJ...843..120S) (`schaefer_C`) |
    | $\mathrm{CO_2} + 2\,\mathrm{H_2} \rightleftharpoons \mathrm{CH_4} + \mathrm{O_2}$ | [Schaefer & Fegley 2017](https://ui.adsabs.harvard.edu/abs/2017ApJ...843..120S) (`schaefer_CH4`) |
    | $\tfrac{1}{2}\,\mathrm{S_2} + \mathrm{O_2} \rightleftharpoons \mathrm{SO_2}$ | JANAF, doubled form (`janaf_SO2`) |
    | $\tfrac{1}{2}\,\mathrm{S_2} + \mathrm{H_2} \rightleftharpoons \mathrm{H_2S}$ | JANAF, doubled form (`janaf_H2S`) |
    | $\tfrac{1}{2}\,\mathrm{N_2} + \tfrac{3}{2}\,\mathrm{H_2} \rightleftharpoons \mathrm{NH_3}$ | JANAF, doubled form (`janaf_NH3`) |

- **One oxygen-fugacity buffer**: [O'Neill & Eggins (2002)](https://ui.adsabs.harvard.edu/abs/2002ChGeo.186..151O) iron-wüstite (default), or [Fischer et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011E%26PSL.304..496F) IW. The shift $\Delta\mathrm{IW}$ sets $\log_{10} f_{\mathrm{O}_2}$ relative to the buffer; under the buffered mode it is a user-prescribed input, under the authoritative-O mode it is a solver unknown.
- **One solubility law per species** with multiple alternative compositions (peridotite, basalt, lunar glass, anorthite-diopside) selectable via constructor argument.

## What is *not* in the model

- **No interior structure**: gravity, radius, mantle mass come in as scalar inputs. Use [Zalmoxis](https://proteus-framework.org/Zalmoxis) for these.
- **No interior thermal evolution**: $T_\mathrm{magma}$ and $\Phi_\mathrm{global}$ come in as scalars. Use [Aragog](https://proteus-framework.org/aragog) or SPIDER.
- **No radiative transfer**: surface partial pressures come out, optical depths and surface temperature come from [AGNI](https://www.h-nicholls.space/AGNI/) or [JANUS](https://proteus-framework.org/JANUS/).
- **No atmospheric escape**: per-iteration mass loss is computed by the PROTEUS escape module ([ZEPHYRUS](https://proteus-framework.org/ZEPHYRUS/)).
- **No solid-phase partitioning**: dissolved-mass fields are written into `_kg_solid` slots that always read `0.0`; CALLIOPE only resolves melt and gas reservoirs. The PROTEUS atmosphere modules handle solid-phase trapping if any.
- **No real-gas EOS**: all species are treated as ideal gases, so partial pressure $\equiv$ fugacity. For non-ideal real-gas effects use [atmodeller](https://atmodeller.readthedocs.io/) ([Bower et al. 2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...59B)).
- **No condensation**: every species is in the gas phase. Condensation chemistry happens in AGNI / JANUS.

## Mathematical statement

CALLIOPE assembles one mass-conservation equation per solved element. Each equation has the structure

$$
\underbrace{m_e^{\mathrm{atm}}(p_{\mathrm{H_2O}}, p_{\mathrm{CO_2}}, p_{\mathrm{N_2}}, p_{\mathrm{S_2}})}_{\text{Bower 2019 Eq. 2 summed over all species}} + \underbrace{m_e^{\mathrm{melt}}(p_{\mathrm{H_2O}}, p_{\mathrm{CO_2}}, p_{\mathrm{N_2}}, p_{\mathrm{S_2}})}_{\text{Henry's law summed over all species}} = m_e^{\mathrm{target}}.
$$

Under the buffered mode the equation set spans $e \in \{\mathrm{H}, \mathrm{C}, \mathrm{N}, \mathrm{S}\}$, the four primary partial pressures are the unknowns, and the $4\times 4$ system is solved with `scipy.optimize.fsolve` (Powell hybrid). Under the authoritative-O mode the equation set spans $e \in \{\mathrm{H}, \mathrm{C}, \mathrm{N}, \mathrm{S}, \mathrm{O}\}$, the unknown vector is extended with $\Delta\mathrm{IW}$, and the $5\times 5$ system is solved with the same outer loop ([Authoritative-oxygen mode](authoritative_oxygen.md)).

The seven secondary partial pressures are *not* independent in either mode: they are algebraic functions of the primaries via the six equilibrium constants, evaluated at $T = T_\mathrm{magma}$ and $\log_{10} f_{\mathrm{O}_2} = \log_{10} f_{\mathrm{O}_2}^\mathrm{IW}(T) + \Delta\mathrm{IW}$.

The four pieces of physics decompose cleanly:

| Component | Page | Implementation |
|---|---|---|
| Speciation tree (primary $\to$ secondary) | [Equilibrium chemistry](equilibrium_chemistry.md) | `chemistry.ModifiedKeq`, `solve.get_partial_pressures` |
| Atmospheric column mass | [Mass balance & solver](mass_balance.md) | `solve.atmosphere_mass` |
| Dissolved mass via Henry / power-law / multi-arg solubility | [Solubility laws](solubility.md) | `solubility.SolubilityH2O`, ..., `solve.dissolved_mass` |
| O'Neill IW buffer | [Oxygen fugacity](oxygen_fugacity.md) | `oxygen_fugacity.OxygenFugacity` |

## Lineage

- **[Bower et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...631A.103B)** introduced the H$_2$O + CO$_2$ mass-balance + Henry's-law treatment that CALLIOPE inherits, including the molar-mass correction $\mu_v / \bar\mu$ in the column-mass relation that earlier studies ([Elkins-Tanton 2008](https://ui.adsabs.harvard.edu/abs/2008E%26PSL.271..181E); [Lebrun et al. 2013](https://ui.adsabs.harvard.edu/abs/2013JGRE..118.1155L); [Salvador et al. 2017](https://ui.adsabs.harvard.edu/abs/2017JGRE..122.1458S); [Nikolaou et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...875...11N)) had omitted.
- **[Bower et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022PSJ.....3...93B)** added the H$_2$, CO, CH$_4$ extensions and the explicit [Schaefer & Fegley (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...843..120S) IVTHANTHERMO / [Chase (1998)](https://janaf.nist.gov/) JANAF equilibrium constants for the H$_2$O–H$_2$, CO$_2$–CO, and CO$_2$+H$_2$–CH$_4$ couples; also adopted the [O'Neill & Eggins (2002)](https://ui.adsabs.harvard.edu/abs/2002ChGeo.186..151O) IW buffer (their Eq. 7) as the parameterisation of mantle redox state.
- **[Nicholls et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024JGRE..12908576N)** introduced N$_2$ via the [Libourel et al. (2003)](https://ui.adsabs.harvard.edu/abs/2003GeCoA..67.4123L) and [Dasgupta et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022GeCoA.336..291D) solubility laws, which is the species set in `calliope.solve.equilibrium_atmosphere` today.
- **[Nicholls et al. (2026)](https://ui.adsabs.harvard.edu/abs/2026NatAs.tmp...61N)** demonstrated the sulfur extension (S$_2$, SO$_2$, H$_2$S) on L 98-59 d, validating the equilibrium constants and the [Gaillard et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022E%26PSL.57717255G) S$_2$ solubility law against in-situ photochemical inferences.

!!! note "Why four primaries"
    CALLIOPE's prognostic *species* are the four primary partial pressures, not the eleven species partial pressures: the gas-phase chemistry collapses the eleven species into four independent mass-balance constraints. N has only one solved degree of freedom even though it appears in both N$_2$ and NH$_3$. O is either not solved (buffered mode) or carried as an additional scalar unknown $\Delta\mathrm{IW}$ alongside the four pressures (authoritative-O mode); in neither case is a new primary partial pressure introduced. Adding a new oxygen-bearing species (e.g. NO) would not require a new constraint, only a new entry in `get_partial_pressures()` and the corresponding contribution to atmospheric and dissolved mass.

## Validity range

CALLIOPE is calibrated for surface temperatures of roughly $1000 \le T_\mathrm{magma} \le 4000$ K and surface pressures of roughly $0.1 \le p_\mathrm{surf} \le 5000$ bar. The lower end of the pressure range is set by numerical stability of the speciation walk; the upper end is the loose envelope above which one or more solubility laws extrapolate. Individual solubility laws have tighter calibration windows than the envelope (Dixon CO$_2$: $\le$815 bar; Sossi H$_2$O: a few kbar; Ardia CH$_4$: 0.7-3 GPa total pressure), see the per-law table in [Solubility laws](solubility.md). Outside the envelope above:

- Below $T \sim 1000$ K the JANAF fits used for the equilibrium constants extrapolate beyond their validation range. The PROTEUS wrapper enforces a configurable `T_floor` (default 700 K), which clips temperatures below `T_floor` to this value, since thermochemical equilibrium does not necessarily hold at cooler temperatures.
- Above $T \sim 4000$ K the mantle-atmosphere partitioning approximation breaks down; switch to atmodeller ([Bower et al. 2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...59B)).
- At surface pressures above ~5 kbar, the H$_2$O solubility laws ([Sossi et al. 2023](https://ui.adsabs.harvard.edu/abs/2023E%26PSL.60117894S), [Newcombe et al. 2017](https://ui.adsabs.harvard.edu/abs/2017GeCoA.200..330N)) extrapolate beyond their experimental calibration window; results are still self-consistent but should be checked against atmodeller for robustness.
- Solid-phase partitioning is ignored; CALLIOPE strictly handles melt + gas. Use it only when $\Phi_\mathrm{global} > 0$, or accept that all dissolved masses will be zero.
