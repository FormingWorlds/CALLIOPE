# CALLIOPE

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/docs.yaml?branch=main&label=Docs)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/docs.yaml)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/CALLIOPE?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/CALLIOPE)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/CALLIOPE/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/nightly.yml)

**CALLIOPE** is the equilibrium outgassing solver of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It computes the partitioning of volatile elements (H, C, N, O, S) between a partially molten silicate mantle and an overlying gas-phase atmosphere, assuming both reservoirs are in thermochemical equilibrium at the planetary surface.

Given an elemental inventory, a magma ocean temperature, a melt fraction, and an oxygen fugacity (specified as a log<sub>10</sub> shift from the iron-wüstite buffer, defaulting to the Fischer et al. 2011 parameterisation), CALLIOPE returns the surface partial pressures of eleven volatile species, the dissolved volatile masses, and the atmospheric mass.

Two solver modes are available: `equilibrium_atmosphere` takes fO<sub>2</sub> as a control variable and derives O from the buffered chemistry, while `equilibrium_atmosphere_authoritative_O` takes total O mass as input and inverts to recover fO<sub>2</sub>. The buffered mode remains the default for standalone use; the authoritative-O mode is the chemistry side of whole-planet oxygen accounting on the PROTEUS side.

Named after the [Greek muse of eloquence and epic poetry](https://en.wikipedia.org/wiki/Calliope). Pronounced *kal-IGH-uh-pee*.

## Supported species

H<sub>2</sub>O, CO<sub>2</sub>, N<sub>2</sub>, S<sub>2</sub> (primary unknowns); H<sub>2</sub>, CH<sub>4</sub>, CO, NH<sub>3</sub>, SO<sub>2</sub>, H<sub>2</sub>S, O<sub>2</sub> (derived from gas-phase equilibrium).

## Documentation

Full documentation is at **[proteus-framework.org/CALLIOPE](https://proteus-framework.org/CALLIOPE)**, including:

- [Getting started](https://proteus-framework.org/CALLIOPE/getting_started.html): installation and a quick path to running.
- [Tutorials](https://proteus-framework.org/CALLIOPE/Tutorials/firstrun.html): first run, Earth and Mars fiducials, two-mode round-trip, coupled-loop driver, speciation phase diagram.
- [How-to guides](https://proteus-framework.org/CALLIOPE/How-to/installation.html): install, configure, run, couple to PROTEUS, use the authoritative-oxygen mode, test, release.
- [Explanations](https://proteus-framework.org/CALLIOPE/Explanations/model.html): model overview, equilibrium chemistry, solubility laws, oxygen fugacity, mass balance, authoritative-oxygen mode, [CALLIOPE-vs-atmodeller cross-backend comparison](https://proteus-framework.org/CALLIOPE/Explanations/cross_backend_comparison.html).
- [API reference](https://proteus-framework.org/CALLIOPE/Reference/api/index.html) and [validation anchors](https://proteus-framework.org/CALLIOPE/Validation/oxygen_fugacity.html): every public function with NumPy-style docstrings, plus the per-source reference-pinned test inventory.

## Installation

```console
pip install fwl-calliope
```

Or, for development:

```console
git clone https://github.com/FormingWorlds/CALLIOPE.git
cd CALLIOPE
pip install -e .[develop,docs]
```

The `docs` extra pulls in [Zensical](https://zensical.org/) so you can build this documentation locally with `zensical serve`.

## Quick start

```python
from calliope.solve import equilibrium_atmosphere, get_target_from_params
from calliope.constants import volatile_species

ddict = {
    'M_mantle': 4.03e24,                      # kg
    'gravity': 9.81, 'radius': 6.371e6,
    'T_magma': 2500.0, 'Phi_global': 1.0,
    'fO2_shift_IW': 0.5,                      # log10 shift relative to IW
    'hydrogen_earth_oceans': 1.0, 'CH_ratio': 0.1,
    'nitrogen_ppmw': 2.0, 'sulfur_ppmw': 200.0,
}
for sp in volatile_species:
    ddict[f'{sp}_included'] = 1
    ddict[f'{sp}_initial_bar'] = 0.0

target = get_target_from_params(ddict)
result = equilibrium_atmosphere(target, ddict)
print(f"P_surf = {result['P_surf']:.1f} bar, "
      f"H2O = {result['H2O_bar']:.1f} bar")
```

See the [first-run tutorial](https://proteus-framework.org/CALLIOPE/Tutorials/firstrun.html) for the full walkthrough.

## Citation

If you use CALLIOPE in published work, please cite the four methods papers below. The full reference list (chemistry constants, solubility laws, oxygen-fugacity buffers, applications) is on the [Publications page](https://proteus-framework.org/CALLIOPE/Reference/publications.html).

- Bower, D.J., Kitzmann, D., Wolf, A.S., Sanan, P., Dorn, C., & Oza, A.V. (2019). *Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations.* **A&A** 631, A103. [\[SciX\]](https://scixplorer.org/abs/2019A%26A...631A.103B/abstract) [\[DOI\]](https://doi.org/10.1051/0004-6361/201935710)
- Bower, D.J., Hakim, K., Sossi, P.A., & Sanan, P. (2022). *Retention of water in terrestrial magma oceans and carbon-rich early atmospheres.* **PSJ** 3, 93. [\[SciX\]](https://scixplorer.org/abs/2022PSJ.....3...93B/abstract) [\[DOI\]](https://doi.org/10.3847/PSJ/ac5fb1)
- Nicholls, H., Lichtenberg, T., Bower, D.J., & Pierrehumbert, R. (2024). *Magma ocean evolution at arbitrary redox state.* **JGR Planets** 129, e2024JE008576. [\[SciX\]](https://scixplorer.org/abs/2024JGRE..12908576N/abstract) [\[DOI\]](https://doi.org/10.1029/2024JE008576) [\[arXiv\]](https://arxiv.org/abs/2411.19137)
- Nicholls, H., Lichtenberg, T., Chatterjee, R.D., Guimond, C.M., Postolec, E., & Pierrehumbert, R.T. (2026). *Volatile-rich evolution of molten super-Earth L 98-59 d.* **Nature Astronomy**. [\[SciX\]](https://scixplorer.org/abs/2026NatAs.tmp...61N/abstract) [\[DOI\]](https://doi.org/10.1038/s41550-026-02815-8)

## License

[Apache License 2.0](LICENSE.txt). CALLIOPE is part of the [Forming Worlds Lab](https://formingworlds.space/) PROTEUS framework.
