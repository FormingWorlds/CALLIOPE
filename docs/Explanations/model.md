# Model overview

CALLIOPE is a 0-D **equilibrium outgassing** solver for the magma-ocean atmosphere coupling. It treats the silicate mantle and the overlying gas-phase atmosphere as a single thermodynamic system in equilibrium at the surface, and asks: for a given total elemental inventory and a given magma-ocean state, what surface partial pressures and dissolved-volatile masses simultaneously satisfy (i) gas-phase chemical equilibrium, (ii) gas-melt solubility equilibrium, and (iii) elemental mass conservation?

This page summarises the model assumptions, the variables it solves for, and how it relates to the upstream papers (Bower et al. 2019[^cite-bower2019], 2022[^cite-bower2022]; Nicholls et al. 2024[^cite-nicholls2024]). Each component has its own dedicated page.

## What is in the model

- **Eleven gas-phase species** (`calliope.constants.volatile_species`): H$_2$O, CO$_2$, O$_2$, H$_2$, CH$_4$, CO, N$_2$, S$_2$, SO$_2$, H$_2$S, NH$_3$.
- **Five elements**: H, C, N, S always solved. O is either a derived quantity set by the oxygen-fugacity buffer (buffered mode, `equilibrium_atmosphere`) or a fifth budgeted element with $\Delta\mathrm{IW}$ as the additional unknown (authoritative-O mode, `equilibrium_atmosphere_authoritative_O`). The two modes share all physics functions and differ only in their unknown set; see [Authoritative-oxygen mode](authoritative_oxygen.md).
- **Six equilibrium reactions** (gas-phase, surface temperature):

    | Reaction | Source |
    |---|---|
    | $\mathrm{H_2O} \rightleftharpoons \mathrm{H_2} + \tfrac{1}{2}\,\mathrm{O_2}$ | [JANAF](https://janaf.nist.gov/) (`janaf_H2`) and Schaefer & Fegley 2017[^cite-schaeferfegley2017] (`schaefer_H`) |
    | $\mathrm{CO_2} \rightleftharpoons \mathrm{CO} + \tfrac{1}{2}\,\mathrm{O_2}$ | [JANAF](https://janaf.nist.gov/) (`janaf_CO`) and Schaefer & Fegley 2017[^cite-schaeferfegley2017] (`schaefer_C`) |
    | $\mathrm{CO_2} + 2\,\mathrm{H_2} \rightleftharpoons \mathrm{CH_4} + \mathrm{O_2}$ | Schaefer & Fegley 2017[^cite-schaeferfegley2017] (`schaefer_CH4`) |
    | $\tfrac{1}{2}\,\mathrm{S_2} + \mathrm{O_2} \rightleftharpoons \mathrm{SO_2}$ | JANAF, doubled form (`janaf_SO2`) |
    | $\tfrac{1}{2}\,\mathrm{S_2} + \mathrm{H_2} \rightleftharpoons \mathrm{H_2S}$ | JANAF, doubled form (`janaf_H2S`) |
    | $\tfrac{1}{2}\,\mathrm{N_2} + \tfrac{3}{2}\,\mathrm{H_2} \rightleftharpoons \mathrm{NH_3}$ | JANAF, doubled form (`janaf_NH3`) |

- **One oxygen-fugacity buffer**: Fischer et al. (2011)[^cite-fischer2011] iron-wüstite (default; close to atmodeller's Hirschmann composite across the magma-ocean range), or the legacy O'Neill & Eggins (2002)[^cite-oneilleggins2002] IW. The shift $\Delta\mathrm{IW}$ sets $\log_{10} f_{\mathrm{O}_2}$ relative to the buffer; under the buffered mode it is a user-prescribed input, under the authoritative-O mode it is a solver unknown.
- **One solubility law per species** with multiple alternative compositions (peridotite, basalt, lunar glass, anorthite-diopside) selectable via constructor argument.

## What is *not* in the model

- **No interior structure**: gravity, radius, mantle mass come in as scalar inputs. Use [Zalmoxis](https://proteus-framework.org/Zalmoxis) for these.
- **No interior thermal evolution**: $T_\mathrm{magma}$ and $\Phi_\mathrm{global}$ come in as scalars. Use [Aragog](https://proteus-framework.org/aragog) or SPIDER.
- **No radiative transfer**: surface partial pressures come out, optical depths and surface temperature come from [AGNI](https://www.h-nicholls.space/AGNI/) or [JANUS](https://proteus-framework.org/JANUS/).
- **No atmospheric escape**: per-iteration mass loss is computed by the PROTEUS escape module ([ZEPHYRUS](https://proteus-framework.org/ZEPHYRUS/)).
- **No solid-phase partitioning**: dissolved-mass fields are written into `_kg_solid` slots that always read `0.0`; CALLIOPE only resolves melt and gas reservoirs. The PROTEUS atmosphere modules handle solid-phase trapping if any.
- **No real-gas EOS**: all species are treated as ideal gases, so partial pressure $\equiv$ fugacity. For non-ideal real-gas effects use [atmodeller](https://atmodeller.readthedocs.io/) (Bower et al. 2025[^cite-bower2025]).
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
| Fischer IW buffer (default), O'Neill IW buffer (legacy) | [Oxygen fugacity](oxygen_fugacity.md) | `oxygen_fugacity.OxygenFugacity` |

## Lineage

- **Bower et al. (2019)[^cite-bower2019]** introduced the H$_2$O + CO$_2$ mass-balance + Henry's-law treatment that CALLIOPE inherits, including the molar-mass correction $\mu_v / \bar\mu$ in the column-mass relation that earlier studies (Elkins-Tanton 2008[^cite-elkinstanton2008]; Lebrun et al. 2013[^cite-lebrun2013]; Salvador et al. 2017[^cite-salvador2017]; Nikolaou et al. 2019[^cite-nikolaou2019]) had omitted.
- **Bower et al. (2022)[^cite-bower2022]** added the H$_2$, CO, CH$_4$ extensions and the explicit Schaefer & Fegley (2017)[^cite-schaeferfegley2017] IVTHANTHERMO / Chase (1998)[^cite-chase1998] JANAF equilibrium constants for the H$_2$O/H$_2$, CO$_2$/CO, and CO$_2$+H$_2$/CH$_4$ couples; also adopted the O'Neill & Eggins (2002)[^cite-oneilleggins2002] IW buffer (their Eq. 7) as the parameterisation of mantle redox state.
- **Nicholls et al. (2024)[^cite-nicholls2024]** introduced N$_2$ via the Libourel et al. (2003)[^cite-libourel2003] and Dasgupta et al. (2022)[^cite-dasgupta2022] solubility laws, which is the species set in `calliope.solve.equilibrium_atmosphere` today.
- **Nicholls et al. (2026)[^cite-nicholls2026]** demonstrated the sulfur extension (S$_2$, SO$_2$, H$_2$S) on L 98-59 d, validating the equilibrium constants and the Gaillard et al. (2022)[^cite-gaillard2022] S$_2$ solubility law against in-situ photochemical inferences.

!!! note "Why four primaries"
    CALLIOPE's prognostic *species* are the four primary partial pressures, not the eleven species partial pressures: the gas-phase chemistry collapses the eleven species into four independent mass-balance constraints. N has only one solved degree of freedom even though it appears in both N$_2$ and NH$_3$. O is either not solved (buffered mode) or carried as an additional scalar unknown $\Delta\mathrm{IW}$ alongside the four pressures (authoritative-O mode); in neither case is a new primary partial pressure introduced. Adding a new oxygen-bearing species (e.g. NO) would not require a new constraint, only a new entry in `get_partial_pressures()` and the corresponding contribution to atmospheric and dissolved mass.

## Validity range

CALLIOPE is calibrated for surface temperatures of roughly $1000 \le T_\mathrm{magma} \le 4000$ K and surface pressures of roughly $0.1 \le p_\mathrm{surf} \le 5000$ bar. The lower end of the pressure range is set by numerical stability of the speciation walk; the upper end is the loose envelope above which one or more solubility laws extrapolate. Individual solubility laws have tighter calibration windows than the envelope (Dixon CO$_2$: $\le$815 bar; Sossi H$_2$O peridotite: 1 atm; Hamilton H$_2$O basalt: 1-6 kbar; Ardia CH$_4$: 0.7-3 GPa total pressure), see the per-law table in [Solubility laws](solubility.md). Outside the envelope above:

- Below $T \sim 1000$ K the JANAF fits used for the equilibrium constants extrapolate beyond their validation range. The PROTEUS wrapper enforces a configurable `T_floor` (default 700 K), which clips temperatures below `T_floor` to this value, since thermochemical equilibrium does not necessarily hold at cooler temperatures.
- Above $T \sim 4000$ K the mantle-atmosphere partitioning approximation breaks down; switch to atmodeller (Bower et al. 2025[^cite-bower2025]).
- The Sossi (2023) peridotite and Newcombe (2017) anorthite-diopside / lunar-glass H$_2$O laws are calibrated at 1 atm. The Dixon (1995) basalt fit is calibrated to 717 bar $p_\mathrm{H_2O}$. For surface pressures above $\sim$1 kbar use the Hamilton (1964) basalt fit (1-6 kbar calibration range) or atmodeller for higher-pressure non-ideal behaviour. Applying the Sossi or Newcombe fits at kbar pressures extrapolates the partial-pressure input by 3 orders of magnitude beyond calibration; the resulting dissolved-mass error scales as $p^{0.5}$ for the power-law fits, so the error is a factor of $\sim$30 at 1 kbar.
- Solid-phase partitioning is ignored; CALLIOPE strictly handles melt + gas. Use it only when $\Phi_\mathrm{global} > 0$, or accept that all dissolved masses will be zero.

[^cite-bower2019]: D. J. Bower, D. Kitzmann, A. S. Wolf, P. Sanan, C. Dorn, A. V. Oza, *[Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations](https://doi.org/10.1051/0004-6361/201935710)*, Astronomy & Astrophysics, 631, A103, 2019. [SciX](https://scixplorer.org/abs/2019A%26A...631A.103B/abstract).
[^cite-bower2022]: D. J. Bower, K. Hakim, P. A. Sossi, P. Sanan, *[Retention of water in terrestrial magma oceans and carbon-rich early atmospheres](https://doi.org/10.3847/PSJ/ac5fb1)*, The Planetary Science Journal, 3(4), 93, 2022. [SciX](https://scixplorer.org/abs/2022PSJ.....3...93B/abstract).
[^cite-bower2025]: D. J. Bower, M. A. Thompson, K. Hakim, M. Tian, P. A. Sossi, *Diversity of low-mass planet atmospheres in the C-H-O-N-S-Cl system with interior dissolution, nonideality, and condensation: application to TRAPPIST-1e and sub-Neptunes*, The Astrophysical Journal, 995, 59, 2025. [SciX](https://scixplorer.org/abs/2025ApJ...995...59B/abstract).
[^cite-chase1998]: M. W. Chase, *[NIST-JANAF Thermochemical Tables, 4th edition](https://janaf.nist.gov/)*, Journal of Physical and Chemical Reference Data Monograph 9, 1998.
[^cite-dasgupta2022]: R. Dasgupta, E. Falksen, A. Pal, C. Sun, *[The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions](https://doi.org/10.1016/j.gca.2022.09.012)*, Geochimica et Cosmochimica Acta, 336, 291–307, 2022. [SciX](https://scixplorer.org/abs/2022GeCoA.336..291D/abstract).
[^cite-elkinstanton2008]: L. T. Elkins-Tanton, *[Linked magma ocean solidification and atmospheric growth for Earth and Mars](https://doi.org/10.1016/j.epsl.2008.03.062)*, Earth and Planetary Science Letters, 271, 181–191, 2008. [SciX](https://scixplorer.org/abs/2008E%26PSL.271..181E/abstract).
[^cite-fischer2011]: R. A. Fischer, A. J. Campbell, G. A. Shofner, O. T. Lord, P. Dera, V. B. Prakapenka, *[Equation of state and phase diagram of FeO](https://doi.org/10.1016/j.epsl.2011.02.025)*, Earth and Planetary Science Letters, 304, 496–502, 2011. [SciX](https://scixplorer.org/abs/2011E%26PSL.304..496F/abstract).
[^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
[^cite-lebrun2013]: T. Lebrun, H. Massol, E. Chassefière, A. Davaille, E. Marcq, P. Sarda, F. Leblanc, G. Brandeis, *[Thermal evolution of an early magma ocean in interaction with the atmosphere](https://doi.org/10.1002/jgre.20068)*, Journal of Geophysical Research: Planets, 118, 1155–1176, 2013. [SciX](https://scixplorer.org/abs/2013JGRE..118.1155L/abstract).
[^cite-libourel2003]: G. Libourel, B. Marty, F. Humbert, *[Nitrogen solubility in basaltic melt. Part I. Effect of oxygen fugacity](https://doi.org/10.1016/S0016-7037(03)00259-X)*, Geochimica et Cosmochimica Acta, 67(21), 4123–4135, 2003. [SciX](https://scixplorer.org/abs/2003GeCoA..67.4123L/abstract).
[^cite-nicholls2024]: H. Nicholls, T. Lichtenberg, D. J. Bower, R. Pierrehumbert, *[Magma ocean evolution at arbitrary redox state](https://doi.org/10.1029/2024JE008576)*, Journal of Geophysical Research: Planets, 129, e2024JE008576, 2024. [SciX](https://scixplorer.org/abs/2024JGRE..12908576N/abstract).
[^cite-nicholls2026]: H. Nicholls, T. Lichtenberg, R. D. Chatterjee, C. M. Guimond, E. Postolec, R. T. Pierrehumbert, *[Volatile-rich evolution of molten super-Earth L 98-59 d](https://doi.org/10.1038/s41550-026-02815-8)*, Nature Astronomy, 2026. [SciX](https://scixplorer.org/abs/2026NatAs.tmp...61N/abstract).
[^cite-nikolaou2019]: A. Nikolaou, N. Katyal, N. Tosi, M. Godolt, J. L. Grenfell, H. Rauer, *[What factors affect the duration and outgassing of the terrestrial magma ocean?](https://doi.org/10.3847/1538-4357/ab08ed)*, The Astrophysical Journal, 875, 11, 2019. [SciX](https://scixplorer.org/abs/2019ApJ...875...11N/abstract).
[^cite-oneilleggins2002]: H. St. C. O'Neill, S. M. Eggins, *[The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts](https://doi.org/10.1016/S0009-2541(01)00414-4)*, Chemical Geology, 186, 151–181, 2002. [SciX](https://scixplorer.org/abs/2002ChGeo.186..151O/abstract).
[^cite-salvador2017]: A. Salvador, H. Massol, A. Davaille, E. Marcq, P. Sarda, E. Chassefière, *The relative influence of H$_2$O and CO$_2$ on the primitive surface conditions and evolution of rocky planets*, Journal of Geophysical Research: Planets, 122, 1458–1486, 2017. [SciX](https://scixplorer.org/abs/2017JGRE..122.1458S/abstract).
[^cite-schaeferfegley2017]: L. Schaefer, B. Fegley, *[Redox states of initial atmospheres outgassed on rocky planets and planetesimals](https://doi.org/10.3847/1538-4357/aa784f)*, The Astrophysical Journal, 843(2), 120, 2017. [SciX](https://scixplorer.org/abs/2017ApJ...843..120S/abstract).
