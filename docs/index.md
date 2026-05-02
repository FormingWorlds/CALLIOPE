![Tests](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/fwl-calliope/badge/?version=latest)](https://fwl-calliope.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# CALLIOPE

**CALLIOPE** is the equilibrium outgassing solver of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It computes the partitioning of volatile elements between a partially molten silicate mantle and an overlying gas-phase atmosphere, assuming both reservoirs are in thermochemical equilibrium at the planetary surface.

Given an elemental inventory (H, C, N, S), a magma ocean temperature $T_\mathrm{magma}$, a melt fraction $\Phi_\mathrm{global}$, and an oxygen fugacity $f_{\mathrm{O}_2}$ (specified as a $\log_{10}$ shift from the iron-wüstite buffer), CALLIOPE solves a four-equation mass-balance system for the surface partial pressures of the four primary species (H$_2$O, CO$_2$, N$_2$, S$_2$) and propagates the speciation to the seven secondary species (H$_2$, CH$_4$, CO, NH$_3$, SO$_2$, H$_2$S, O$_2$).

Named after the [Greek muse of eloquence and epic poetry](https://en.wikipedia.org/wiki/Calliope). Pronounced *kal-IGH-uh-pee*.

!!! tip "New to CALLIOPE?"
    See the **[Getting Started guide](getting_started.md)** for installation, first run, and basic usage.

## Features

- **Eleven volatile species**: H$_2$O, CO$_2$, N$_2$, S$_2$ as primary unknowns; H$_2$, CH$_4$, CO, NH$_3$, SO$_2$, H$_2$S, O$_2$ derived from gas-phase equilibrium
- **Five elemental conservation channels**: H, C, N, S as solved constraints; O fixed by the $f_{\mathrm{O}_2}$ buffer
- **Configurable redox state**: O'Neill & Eggins (2002) iron-wüstite (IW) buffer with arbitrary $\Delta\mathrm{IW}$ shift, or Fischer et al. (2011) IW
- **Calibrated equilibrium constants**: JANAF and Schaefer & Fegley (2017) fits for the H$_2$O–H$_2$, CO$_2$–CO, CO$_2$+H$_2$–CH$_4$, S$_2$–SO$_2$, S$_2$+H$_2$–H$_2$S, and N$_2$+H$_2$–NH$_3$ couples
- **Multiple solubility laws per species**: peridotite (default H$_2$O, Sossi et al. 2022), basalt (Dixon et al. 1995, Wilson & Head 1981, Hamilton 1964), anorthite-diopside (Newcombe et al. 2017), lunar glass (Newcombe et al. 2017); CO$_2$ (Dixon et al. 1995); CO (Armstrong et al. 2015); CH$_4$ (Ardia et al. 2013); N$_2$ (Libourel et al. 2003 or Dasgupta et al. 2022); S$_2$ (Gaillard et al. 2022)
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

- Bower, D.J., Kitzmann, D., Wolf, A.S., Sanan, P., Dorn, C., & Oza, A.V. (2019). *Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations*. **Astronomy & Astrophysics**, 631, A103. \[[ADS](https://ui.adsabs.harvard.edu/abs/2019A%26A...631A.103B) | [DOI](https://doi.org/10.1051/0004-6361/201935710)\]
- Bower, D.J., Hakim, K., Sossi, P.A., & Sanan, P. (2022). *Retention of water in terrestrial magma oceans and carbon-rich early atmospheres*. **The Planetary Science Journal**, 3, 93. \[[ADS](https://ui.adsabs.harvard.edu/abs/2022PSJ.....3...93B) | [DOI](https://doi.org/10.3847/PSJ/ac5fb1)\]
- Nicholls, H., Lichtenberg, T., Bower, D.J., & Pierrehumbert, R. (2024). *Magma ocean evolution at arbitrary redox state*. **Journal of Geophysical Research: Planets**, 129, e2024JE008576. \[[ADS](https://ui.adsabs.harvard.edu/abs/2024JGRE..12908576N) | [DOI](https://doi.org/10.1029/2024JE008576) | [arXiv](https://arxiv.org/abs/2411.19137)\]

See the [Publications](Reference/publications.md) page for the full reference list, including the underlying solubility-law and equilibrium-constant sources.

## Code availability

- [GitHub repository (Forming Worlds)](https://github.com/FormingWorlds/CALLIOPE)
- [PyPI package `fwl-calliope`](https://pypi.org/project/fwl-calliope/)

If you plan to contribute to CALLIOPE, please read our [Code of Conduct](Community/CODE_OF_CONDUCT.md) and [contributing guidelines](Community/CONTRIBUTING.md).
If you are running into problems, please do not hesitate to raise an [Issue](https://github.com/FormingWorlds/CALLIOPE/issues).

## License

[Apache License 2.0](https://opensource.org/licenses/Apache-2.0). See [the included license](https://github.com/FormingWorlds/CALLIOPE/blob/main/LICENSE.txt).
