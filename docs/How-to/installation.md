# Installation

CALLIOPE is a pure-Python package with three runtime dependencies (`numpy`, `scipy`, `matplotlib`, plus `cmcrameri` for plotting). It supports CPython 3.11, 3.12, and 3.13.

## From PyPI (recommended for end users)

```console
pip install fwl-calliope
```

This is what the PROTEUS install pulls in by default; you do not need to install CALLIOPE separately if you have already installed `fwl-proteus`.

## From source (recommended for contributors)

Clone the repository and install in editable mode so that your changes are picked up without reinstalling:

```console
git clone https://github.com/FormingWorlds/CALLIOPE.git
cd CALLIOPE
pip install -e .
```

To also install the optional development tooling (pytest, ruff, pre-commit, bump-my-version, coverage):

```console
pip install -e .[develop]
```

To install the documentation tooling on top of that:

```console
pip install -e .[develop,docs]
```

The `docs` extra pulls in `mkdocs`, `mkdocs-material`, `mkdocstrings[python]`, and `markdown-include`, which are everything you need to build this site locally.

## Verifying the install

```python
import calliope
print(calliope.__version__)
```

Should print the current CalVer version (e.g. `25.05.04`).

A minimal smoke-test that exercises the equilibrium solver:

```python
from calliope.solve import equilibrium_atmosphere
from calliope.constants import volatile_species

ddict = {
    'M_mantle': 4.03e24,    # kg, Earth mantle mass
    'gravity': 9.81,        # m s-2
    'radius': 6.371e6,      # m, planetary surface radius
    'Phi_global': 1.0,      # fully molten
    'T_magma': 2500.0,      # K
    'fO2_shift_IW': 0.0,    # at IW buffer
}
for sp in volatile_species:
    ddict[f'{sp}_included'] = 1
    ddict[f'{sp}_initial_bar'] = 0.0

target = {'H': 1e20, 'C': 1e17, 'N': 1e17, 'S': 1e16}    # kg
result = equilibrium_atmosphere(target, ddict, print_result=False)
print(f"Surface pressure: {result['P_surf']:.2f} bar")
```

A successful run should print a non-zero surface pressure (typically a few hundred to a few thousand bar for these inputs). Any exception means the install or the runtime environment is broken; check that `numpy` and `scipy` import cleanly and that you are on Python 3.11+.

## Next step

Once the import works, head to [Configuration](configuration.md) to learn what each input field means, then to [Usage](usage.md) for the end-to-end calling pattern.
