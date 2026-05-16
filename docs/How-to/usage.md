# Usage

This page shows the two canonical calling patterns for CALLIOPE: solving from an elemental inventory, and solving from a prescribed initial-pressure atmosphere. Both use the same `equilibrium_atmosphere()` entry point.

For an end-to-end runnable example with explanations, see the [First run tutorial](../Tutorials/firstrun.md). For the schema of each input field, see [Configuration](configuration.md).

## Pattern 1 - solve from elemental inventory

This is the workflow PROTEUS uses. You specify the total mass of H, C, N, S in the planet, plus the magma ocean state, and CALLIOPE returns the equilibrium speciation.

```python
from calliope.solve import equilibrium_atmosphere, get_target_from_params
from calliope.constants import volatile_species

# 1. Build ddict (planet + magma + inclusion flags)
ddict = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
    'Phi_global': 1.0,
    'T_magma': 2500.0,
    'fO2_shift_IW': 0.0,
    # elemental-budget inputs for get_target_from_params
    'hydrogen_earth_oceans': 1.0,    # 1 Earth ocean of H
    'CH_ratio': 0.1,
    'nitrogen_ppmw': 2.0,            # primitive Earth mantle estimate
    'sulfur_ppmw': 200.0,
}
for sp in volatile_species:
    ddict[f'{sp}_included'] = 1
    ddict[f'{sp}_initial_bar'] = 0.0

# 2. Build target from elemental budget
target = get_target_from_params(ddict)

# 3. Solve
result = equilibrium_atmosphere(target, ddict, print_result=True)

# 4. Read out
print(f"Surface pressure  : {result['P_surf']:.2f} bar")
print(f"Atmosphere mass   : {result['M_atm']:.2e} kg")
print(f"H2O bar           : {result['H2O_bar']:.2f}")
print(f"H2O VMR           : {result['H2O_vmr']:.3f}")
print(f"CO2 bar           : {result['CO2_bar']:.2f}")
print(f"H2  VMR           : {result['H2_vmr']:.3f}")
```

The `result` dictionary contains, for every species `s` in `volatile_species`:

- `s_bar`: surface partial pressure, bar
- `s_vmr`: volume mixing ratio (== mole fraction)
- `s_kg_atm`: atmospheric column mass, kg
- `s_kg_liquid`: dissolved mass, kg
- `s_kg_solid`: 0.0 (solid partitioning is handled by the caller, not CALLIOPE)
- `s_kg_total`: sum of atmosphere + liquid

and per-element totals `H_kg_atm`, `C_kg_atm`, ..., `S_res` (mass-balance residual in kg).

## Pattern 2 - solve from initial partial pressures

Use this when you know the initial composition of the atmosphere and want to back out the consistent total inventory plus the equilibrium update.

```python
from calliope.solve import equilibrium_atmosphere, get_target_from_pressures
from calliope.constants import volatile_species

ddict = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
    'Phi_global': 1.0,
    'T_magma': 2500.0,
    'fO2_shift_IW': 0.0,
}
for sp in volatile_species:
    ddict[f'{sp}_included'] = 1
    ddict[f'{sp}_initial_bar'] = 0.0

# Specify the initial atmosphere by its primary-species partial pressures
ddict['H2O_initial_bar'] = 220.0    # bar
ddict['CO2_initial_bar'] = 100.0
ddict['N2_initial_bar']  = 1.0
ddict['S2_initial_bar']  = 0.01

target = get_target_from_pressures(ddict)
result = equilibrium_atmosphere(target, ddict, print_result=True)
```

`get_target_from_pressures()` computes the implied total elemental masses by summing atmospheric column mass (Bower et al. (2019)[^cite-bower2019] Eq. 2) and dissolved mass (Henry's law) across every included species at the prescribed initial pressures.

## Warm-starting from a previous solve

For repeated calls (e.g. inside a time-evolution loop), pass the previous solution as the initial guess to dramatically reduce the number of `fsolve` iterations:

```python
p_guess = {
    'H2O': result['H2O_bar'],
    'CO2': result['CO2_bar'],
    'N2':  result['N2_bar'],
    'S2':  result['S2_bar'],
}
result_new = equilibrium_atmosphere(target_new, ddict_new, p_guess=p_guess, print_result=False)
```

This is exactly what the PROTEUS wrapper does. With a good warm start, convergence typically takes 1 to 3 `fsolve` iterations; without one, CALLIOPE may need 20+ Monte-Carlo restarts before it lands on a plausible basin.

## Choosing solubility laws

To override the default solubility law for a species, instantiate the law explicitly and pass it through your own wrapper rather than calling `equilibrium_atmosphere()` directly. The defaults are:

| Species | Default class | Default composition | Source |
|---|---|---|---|
| H$_2$O | `SolubilityH2O` | `peridotite` | Sossi et al. (2023)[^cite-sossi2023] |
| CO$_2$ | `SolubilityCO2` | `basalt_dixon` | Dixon et al. (1995)[^cite-dixon1995] |
| CO | `SolubilityCO` | `mafic_armstrong` | Armstrong et al. (2015)[^cite-armstrong2015] |
| CH$_4$ | `SolubilityCH4` | `basalt_ardia` | Ardia et al. (2013)[^cite-ardia2013] |
| N$_2$ | `SolubilityN2` | `dasgupta` (in `dissolved_mass`) | Dasgupta et al. (2022)[^cite-dasgupta2022] |
| S$_2$ | `SolubilityS2` | `gaillard` | Gaillard et al. (2022)[^cite-gaillard2022] |

Alternative compositions (e.g. `SolubilityH2O('basalt_dixon')`, `SolubilityH2O('lunar_glass')`) are documented in the [API reference](../Reference/api/calliope.solubility.md) and discussed in [Solubility laws](../Explanations/solubility.md).

## Next step

For the science behind these laws and the equilibrium constants, head to [Equilibrium chemistry](../Explanations/equilibrium_chemistry.md) and [Solubility laws](../Explanations/solubility.md). For the PROTEUS-side TOML recipe, head to [Coupling to PROTEUS](proteus_coupling.md). If you have an oxygen budget to enforce rather than a buffer offset to apply, switch to the [authoritative-O recipe](authoritative_oxygen.md).

[^cite-ardia2013]: P. Ardia, M. M. Hirschmann, A. C. Withers, B. D. Stanley, *[Solubility of CH$_4$ in a synthetic basaltic melt, with applications to atmosphere-magma ocean-core partitioning of volatiles and to the evolution of the Martian atmosphere](https://doi.org/10.1016/j.gca.2013.03.028)*, Geochimica et Cosmochimica Acta, 114, 52–71, 2013. [SciX](https://scixplorer.org/abs/2013GeCoA.114...52A/abstract).
[^cite-armstrong2015]: L. S. Armstrong, M. M. Hirschmann, B. D. Stanley, E. G. Falksen, S. D. Jacobsen, *[Speciation and solubility of reduced C-O-H-N volatiles in mafic melt: implications for volcanism, atmospheric evolution, and deep volatile cycles in the terrestrial planets](https://doi.org/10.1016/j.gca.2015.07.007)*, Geochimica et Cosmochimica Acta, 171, 283–302, 2015. [SciX](https://scixplorer.org/abs/2015GeCoA.171..283A/abstract).
[^cite-bower2019]: D. J. Bower, D. Kitzmann, A. S. Wolf, P. Sanan, C. Dorn, A. V. Oza, *[Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations](https://doi.org/10.1051/0004-6361/201935710)*, Astronomy & Astrophysics, 631, A103, 2019. [SciX](https://scixplorer.org/abs/2019A%26A...631A.103B/abstract).
[^cite-dasgupta2022]: R. Dasgupta, E. Falksen, A. Pal, C. Sun, *[The fate of nitrogen during parent body partial melting and accretion of the inner Solar System bodies at reducing conditions](https://doi.org/10.1016/j.gca.2022.09.012)*, Geochimica et Cosmochimica Acta, 336, 291–307, 2022. [SciX](https://scixplorer.org/abs/2022GeCoA.336..291D/abstract).
[^cite-dixon1995]: J. E. Dixon, E. M. Stolper, J. R. Holloway, *[An experimental study of water and carbon dioxide solubilities in mid-ocean ridge basaltic liquids. Part I: Calibration and solubility models](https://doi.org/10.1093/oxfordjournals.petrology.a037267)*, Journal of Petrology, 36(6), 1607–1631, 1995. [SciX](https://scixplorer.org/abs/1995JPet...36.1607D/abstract).
[^cite-gaillard2022]: F. Gaillard, F. Bernadou, M. Roskosz, M. A. Bouhifd, Y. Marrocchi, G. Iacono-Marziano, M. Moreira, B. Scaillet, G. Rogerie, *[Redox controls during magma ocean degassing](https://doi.org/10.1016/j.epsl.2021.117255)*, Earth and Planetary Science Letters, 577, 117255, 2022. [SciX](https://scixplorer.org/abs/2022E%26PSL.57717255G/abstract).
[^cite-sossi2023]: P. A. Sossi, P. M. E. Tollan, J. Badro, D. J. Bower, *[Solubility of water in peridotite liquids and the prevalence of steam atmospheres on rocky planets](https://doi.org/10.1016/j.epsl.2022.117894)*, Earth and Planetary Science Letters, 601, 117894, 2023. [SciX](https://scixplorer.org/abs/2023E%26PSL.60117894S/abstract).
