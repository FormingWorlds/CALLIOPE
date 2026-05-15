# Getting started

!!! note "Usage within the PROTEUS framework"
    CALLIOPE is most commonly installed and used as part of the [PROTEUS framework](https://proteus-framework.org/PROTEUS). For coupled atmosphere-interior runs, the [PROTEUS Getting Started guide](https://proteus-framework.org/PROTEUS) is the right entry point; this site documents the standalone CALLIOPE API and the interface it exposes to PROTEUS.

## Quick path

1. **Install CALLIOPE**
   Editable install with the `docs` extra so you can also rebuild this site. <br>
   See [Installation guide](How-to/installation.md).

2. **Configure your problem**
   Build the input dictionary that CALLIOPE expects: planet mass, gravity, radius, $T_\mathrm{magma}$, $\Phi_\mathrm{global}$, $\Delta\mathrm{IW}$, plus the included-volatile flags. <br>
   See [Configuration](How-to/configuration.md).

3. **Run the equilibrium solver**
   Call `equilibrium_atmosphere(target, ddict)` from the `calliope.solve` module on your prepared input. <br>
   See [Usage](How-to/usage.md), or follow the end-to-end [First run tutorial](Tutorials/firstrun.md).

---

## What do you want to do?

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
