# CALLIOPE

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/docs.yaml?branch=main&label=Docs)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/docs.yaml)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/CALLIOPE?label=coverage&logo=codecov&color=brightgreen)](https://app.codecov.io/gh/FormingWorlds/CALLIOPE)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/tests.yaml?branch=main&label=Unit%20Tests&color=brightgreen)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/nightly.yml?branch=main&label=Integration%20Tests&color=brightgreen)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/nightly.yml)

**CALLIOPE** is the equilibrium outgassing solver of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It computes the partitioning of volatile elements between a partially molten silicate mantle and an overlying gas-phase atmosphere, assuming both reservoirs are in thermochemical equilibrium at the planetary surface.

Given an elemental inventory and a magma ocean state ($T_\mathrm{magma}$, $\Phi_\mathrm{global}$), CALLIOPE solves a nonlinear mass-balance system for the surface partial pressures of the four primary species (H$_2$O, CO$_2$, N$_2$, S$_2$) and propagates the speciation to the seven secondary species (H$_2$, CH$_4$, CO, NH$_3$, SO$_2$, H$_2$S, O$_2$). The solver runs in either of two modes that share the same physics functions and differ only in their unknown set:

- **Buffered mode** ([`equilibrium_atmosphere`](Explanations/mass_balance.md)) takes an oxygen fugacity $f_{\mathrm{O}_2}$ (specified as a $\log_{10}$ shift from the iron-wüstite buffer) as input and solves a four-equation system for the H, C, N, S budget. Oxygen mass is derived.
- **Authoritative-oxygen mode** ([`equilibrium_atmosphere_authoritative_O`](Explanations/authoritative_oxygen.md)) takes a five-element budget including O and solves a five-equation system for the four pressures plus $\Delta\mathrm{IW}$. Oxygen fugacity is derived.

Named after the [Greek muse of eloquence and epic poetry](https://en.wikipedia.org/wiki/Calliope). Pronounced *kal-IGH-uh-pee*.

!!! tip "New to CALLIOPE?"
    See the **[Getting Started guide](getting_started.md)** for installation, first run, and basic usage.

## Features

- **Eleven volatile species**: H$_2$O, CO$_2$, N$_2$, S$_2$ as primary unknowns; H$_2$, CH$_4$, CO, NH$_3$, SO$_2$, H$_2$S, O$_2$ derived from gas-phase equilibrium
- **Five elemental conservation channels**: H, C, N, S always solved; O either derived from the $f_{\mathrm{O}_2}$ buffer (buffered mode) or supplied as a fifth budget (authoritative-O mode)
- **Configurable redox state**: Fischer et al. (2011)[^cite-fischer2011] iron-wüstite (IW) buffer with arbitrary $\Delta\mathrm{IW}$ shift (default; chosen to be close to atmodeller's Hirschmann composite across the magma-ocean range), or the legacy O'Neill & Eggins (2002)[^cite-oneilleggins2002] IW
- **Calibrated equilibrium constants**: JANAF[^cite-chase1998] and Schaefer & Fegley (2017)[^cite-schaeferfegley2017] fits for the H$_2$O-H$_2$, CO$_2$-CO, CO$_2$+H$_2$-CH$_4$, S$_2$-SO$_2$, S$_2$+H$_2$-H$_2$S, and N$_2$+H$_2$-NH$_3$ couples
- **Multiple solubility laws per species**: peridotite (default H$_2$O, Sossi et al. 2023[^cite-sossi2023]), basalt (Dixon et al. 1995[^cite-dixon1995], Wilson & Head 1981[^cite-wilsonhead1981], Hamilton et al. 1964[^cite-hamilton1964]), anorthite-diopside (Newcombe et al. 2017[^cite-newcombe2017]), lunar glass (Newcombe et al. 2017[^cite-newcombe2017]); CO$_2$ (Dixon et al. 1995[^cite-dixon1995]); CO (Armstrong et al. 2015[^cite-armstrong2015]); CH$_4$ (Ardia et al. 2013[^cite-ardia2013]); N$_2$ (Libourel et al. 2003[^cite-libourel2003] or Dasgupta et al. 2022[^cite-dasgupta2022]); S$_2$ (Gaillard et al. 2022[^cite-gaillard2022])
- **Robust hybrid solver**: alternating `scipy.optimize.fsolve` (Powell hybrid) and `trust-constr` minimisation, with Monte-Carlo restart on failure
- **PROTEUS-coupled or standalone**: the same equilibrium kernel powers both the in-loop call from PROTEUS and one-off scripts

!!! info "PROTEUS framework"
    When used within PROTEUS, CALLIOPE is called at every coupling timestep to update the surface partial pressures, atmosphere mass, and dissolved volatile masses, after the structure module (Zalmoxis or SPIDER) has updated the planetary radius and gravity, and before the interior energetics module (Aragog or SPIDER) advances the entropy. The PROTEUS-side documentation is at [proteus-framework.org/PROTEUS](https://proteus-framework.org/PROTEUS); the recipe for using CALLIOPE inside that pipeline is on the [Coupling to PROTEUS (how-to)](How-to/proteus_coupling.md) page; the per-iteration control flow and the `[outgas.calliope]` schema mapping is on the [Coupling to PROTEUS (theory)](Explanations/proteus_coupling.md) page.

## Quick links

<div class="grid cards" markdown>

-   :material-download: **Install**

    [Go to installation guide](How-to/installation.md)

-   :material-tune: **Configure**

    [Go to configuration](How-to/configuration.md)

-   :material-rocket-launch: **Use CALLIOPE**

    [Go to usage](How-to/usage.md)

-   :material-book-open-variant: **Understand the model**

    [Go to model overview](Explanations/model.md)

-   :material-code-braces: **Browse the API**

    [Go to API reference](Reference/api/index.md)

-   :material-github: **Contribute / browse code**

    [Go to source code](https://github.com/FormingWorlds/CALLIOPE)

-   :material-bug: **Raise an issue**

    [Go to issues](https://github.com/FormingWorlds/CALLIOPE/issues)

-   :material-email: **Get in touch**

    [Go to contact](Community/contact.md)

</div>

## Citation

If you use CALLIOPE in published work, please cite the original equilibrium-chemistry framework, the modern multi-species redox treatment, and the magma-ocean evolution study that introduced the present extended species list:

- Bower, D.J., Kitzmann, D., Wolf, A.S., Sanan, P., Dorn, C., & Oza, A.V. (2019). *Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations*. **Astronomy & Astrophysics**, 631, A103. \[[SciX](https://scixplorer.org/abs/2019A%26A...631A.103B/abstract) | [DOI](https://doi.org/10.1051/0004-6361/201935710)\]
- Bower, D.J., Hakim, K., Sossi, P.A., & Sanan, P. (2022). *Retention of water in terrestrial magma oceans and carbon-rich early atmospheres*. **The Planetary Science Journal**, 3, 93. \[[SciX](https://scixplorer.org/abs/2022PSJ.....3...93B/abstract) | [DOI](https://doi.org/10.3847/PSJ/ac5fb1)\]
- Nicholls, H., Lichtenberg, T., Bower, D.J., & Pierrehumbert, R. (2024). *Magma ocean evolution at arbitrary redox state*. **Journal of Geophysical Research: Planets**, 129, e2024JE008576. \[[SciX](https://scixplorer.org/abs/2024JGRE..12908576N/abstract) | [DOI](https://doi.org/10.1029/2024JE008576) | [arXiv](https://arxiv.org/abs/2411.19137)\]
- Nicholls, H., Lichtenberg, T., Chatterjee, R.D., Guimond, C.M., Postolec, E., & Pierrehumbert, R.T. (2026). *Volatile-rich evolution of molten super-Earth L 98-59 d*. **Nature Astronomy**. \[[SciX](https://scixplorer.org/abs/2026NatAs.tmp...61N/abstract) | [DOI](https://doi.org/10.1038/s41550-026-02815-8)\]

See the [Publications](Reference/publications.md) page for the full reference list, including the underlying solubility-law and equilibrium-constant sources.

## Code availability

- [GitHub repository (Forming Worlds)](https://github.com/FormingWorlds/CALLIOPE)
- [PyPI package `fwl-calliope`](https://pypi.org/project/fwl-calliope/)

If you plan to contribute to CALLIOPE, please read our [Code of Conduct](Community/CODE_OF_CONDUCT.md) and [contributing guidelines](Community/CONTRIBUTING.md).
If you are running into problems, please do not hesitate to raise an [Issue](https://github.com/FormingWorlds/CALLIOPE/issues).

## License

[Apache License 2.0](https://opensource.org/licenses/Apache-2.0). See [the included license](https://github.com/FormingWorlds/CALLIOPE/blob/main/LICENSE.txt).

[^cite-ardia2013]: P. Ardia, M. M. Hirschmann, A. C. Withers, B. D. Stanley, *[Solubility of CH$_4$ in a synthetic basaltic melt, with applications to atmosphere-magma ocean-core partitioning of volatiles and to the evolution of the Martian atmosphere](https://doi.org/10.1016/j.gca.2013.03.028)*, Geochimica et Cosmochimica Acta, 114, 52–71, 2013. [SciX](https://scixplorer.org/abs/2013GeCoA.114...52A/abstract).
[^cite-armstrong2015]: L. S. Armstrong, M. M. Hirschmann, B. D. Stanley, E. G. Falksen, S. D. Jacobsen, *[Speciation and solubility of reduced C-O-H-N volatiles in mafic melt: implications for volcanism, atmospheric evolution, and deep volatile cycles in the terrestrial planets](https://doi.org/10.1016/j.gca.2015.07.007)*, Geochimica et Cosmochimica Acta, 171, 283–302, 2015. [SciX](https://scixplorer.org/abs/2015GeCoA.171..283A/abstract).
[^cite-chase1998]: M. W. Chase, *[NIST-JANAF Thermochemical Tables, 4th edition](https://janaf.nist.gov/)*, Journal of Physical and Chemical Reference Data Monograph 9, 1998.
[^cite-dasgupta2022]: R. Dasgupta, E. Falksen, A. Pal, C. Sun, *[The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions](https://doi.org/10.1016/j.gca.2022.09.012)*, Geochimica et Cosmochimica Acta, 336, 291–307, 2022. [SciX](https://scixplorer.org/abs/2022GeCoA.336..291D/abstract).
[^cite-dixon1995]: J. E. Dixon, E. M. Stolper, J. R. Holloway, *[An experimental study of water and carbon dioxide solubilities in mid-ocean ridge basaltic liquids. Part I: Calibration and solubility models](https://doi.org/10.1093/oxfordjournals.petrology.a037267)*, Journal of Petrology, 36(6), 1607–1631, 1995. [SciX](https://scixplorer.org/abs/1995JPet...36.1607D/abstract).
[^cite-fischer2011]: R. A. Fischer, A. J. Campbell, G. A. Shofner, O. T. Lord, P. Dera, V. B. Prakapenka, *[Equation of state and phase diagram of FeO](https://doi.org/10.1016/j.epsl.2011.02.025)*, Earth and Planetary Science Letters, 304, 496–502, 2011. [SciX](https://scixplorer.org/abs/2011E%26PSL.304..496F/abstract).
[^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
[^cite-hamilton1964]: D. L. Hamilton, C. W. Burnham, E. F. Osborn, *[The solubility of water and effects of oxygen fugacity and water content on crystallization in mafic magmas](https://doi.org/10.1093/petrology/5.1.21)*, Journal of Petrology, 5(1), 21–39, 1964.
[^cite-libourel2003]: G. Libourel, B. Marty, F. Humbert, *[Nitrogen solubility in basaltic melt. Part I. Effect of oxygen fugacity](https://doi.org/10.1016/S0016-7037(03)00259-X)*, Geochimica et Cosmochimica Acta, 67(21), 4123–4135, 2003. [SciX](https://scixplorer.org/abs/2003GeCoA..67.4123L/abstract).
[^cite-newcombe2017]: M. E. Newcombe, A. Brett, J. R. Beckett, M. B. Baker, S. Newman, Y. Guan, J. M. Eiler, E. M. Stolper, *[Solubility of water in lunar basalt at low pH$_2$O](https://doi.org/10.1016/j.gca.2016.12.026)*, Geochimica et Cosmochimica Acta, 200, 330–352, 2017. [SciX](https://scixplorer.org/abs/2017GeCoA.200..330N/abstract).
[^cite-oneilleggins2002]: H. St. C. O'Neill, S. M. Eggins, *[The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts](https://doi.org/10.1016/S0009-2541(01)00414-4)*, Chemical Geology, 186, 151–181, 2002. [SciX](https://scixplorer.org/abs/2002ChGeo.186..151O/abstract).
[^cite-schaeferfegley2017]: L. Schaefer, B. Fegley, *[Redox states of initial atmospheres outgassed on rocky planets and planetesimals](https://doi.org/10.3847/1538-4357/aa784f)*, The Astrophysical Journal, 843(2), 120, 2017. [SciX](https://scixplorer.org/abs/2017ApJ...843..120S/abstract).
[^cite-sossi2023]: P. A. Sossi, P. M. E. Tollan, J. Badro, D. J. Bower, *[Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets](https://doi.org/10.1016/j.epsl.2022.117894)*, Earth and Planetary Science Letters, 601, 117894, 2023. [SciX](https://scixplorer.org/abs/2023E%26PSL.60117894S/abstract).
[^cite-wilsonhead1981]: L. Wilson, J. W. Head, *[Ascent and eruption of basaltic magma on the Earth and Moon](https://doi.org/10.1029/JB086iB04p02971)*, Journal of Geophysical Research, 86(B4), 2971–3001, 1981. [SciX](https://scixplorer.org/abs/1981JGR....86.2971W/abstract).
