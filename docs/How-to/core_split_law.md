# Add a partition law to the core split

The [core and mantle split](../Explanations/metal.md#core-and-mantle-split-at-an-equilibration-point) takes one partition law per element, by name. This guide shows how to add your own law for H, C, N or S and use it in `calliope.core_split.split_core_mantle`. It has 2 complete examples that you can copy and run:

- **Example A** adds a law as a **partition file entry**. Use this when the law has the standard form of CALLIOPE (below).
- **Example B** adds a law as a **Python function**. Use this for any other form.

You should already know what a partition coefficient is, and you should have done the [core split tutorial](../Tutorials/core_split.md). Run the code blocks in order in one Python session. The whole guide runs in a few seconds.

## Before you start: what you need from the paper

Write down these 6 things before you write any code. Each one maps onto a field or an argument below.

1. **The element**: H, C, N or S. O and Si are not split (see [Limits](../Explanations/metal.md#limits)).
2. **What the law gives**: $D$, the ratio of the concentration in the metal to the concentration in the silicate melt, or an exchange coefficient $K_D$ that contains $D$.
3. **The basis of $D$**: a ratio of mass fractions (`mass`), a ratio of mole fractions over the whole metal (`mole`), or a ratio of mole fractions over the metal without H and C (`mole_xprime_noHC`).
4. **The units**: CALLIOPE takes the pressure in GPa and the temperature in K.
5. **The range of the data**: the lowest and highest pressure, temperature and $\Delta\mathrm{IW}$ of the experiments that the law was fitted to. Take them from the table of runs of the paper.
6. **The source**: the paper, the equation, the table and the page, and the DOI.

## The standard form of a law

A partition file entry of the form `log10_linear` gives

$$
\log_{10} D = a + \frac{b}{T} + c\,\frac{P}{T} + d \cdot \mathrm{dIW\_factor} \cdot \Delta\mathrm{IW} + e \cdot \mathrm{NBO/T} + \sum_j \varepsilon_j \log_{10}(1 - X_j),
$$

with $T$ in K, $P$ in GPa, NBO/T the ratio of non-bridging O to tetrahedral cations in the melt (how polymerised the melt is), and $X_j$ the mole fraction of the light element $j$ over the metal without H and C (the $x'$ basis). An entry can add $f \log_{10} X_\mathrm{FeO} - \log_{10} C_\mathrm{S}$, with $X_\mathrm{FeO}$ the FeO mole fraction of the melt and $C_\mathrm{S}$ the sulfide capacity of the melt. If your law fits this form, use Example A. If it does not, use Example B.

## Example A: a law as a partition file entry

The law is the H law of [Tsutsumi et al. (2025)](https://doi.org/10.1038/s41467-025-65729-5), *Nature Communications* 16, 10038. The [partition law tutorial](../Tutorials/partition_law.md#step-1-the-published-law) explains how its coefficients map onto the file fields.

### Step A1: read the law from the paper

Their Eq. 3 (p. 3) is

$$
\log_{10} D_\mathrm{H} + \frac{1}{4}\Delta\mathrm{IW} = a + \frac{b}{T} + c\,\frac{P}{T} + d \cdot \log_{10}\left(1 - x'_\mathrm{C}\right),
$$

with $a = 2.31$, $b = -1542$, $c = -38.9$ and $d = 6.69$ (p. 3). $D_\mathrm{H}$ is a ratio of mole fractions over the metal without H and C (their Eq. 1 and Eq. 10). The ranges are those of the runs of the fit: runs 1 to 5 of their Table 1 (p. 2) and the runs of Tagawa et al. (2021), together 30 to 60 GPa, 3080 to 4760 K and $\Delta\mathrm{IW}$ -1.78 to -0.88.

Move $\frac{1}{4}\Delta\mathrm{IW}$ to the right side. Then $d = -0.25$ in the file, and the paper's $d$ goes into `eps`, the term for the C of the metal.

### Step A2: write the file

A partition file can hold one law. Each field is explained in its comment:

```python
from calliope.core_split import split_core_mantle
from calliope.metal import (
    get_law, law_names, load_partition_table, log10_D, register_law
)

law_file = '''schema_version = 1

[meta]
name = "my-laws"           # any name for your file
placeholder = false        # true makes each solve with this file log a warning

[H]                        # the element of the law
form = "log10_linear"      # the standard form above
quantity = "K_D"           # what the paper fits; the law still gives D
basis = "mole_xprime_noHC" # D over the metal without H and C (their Eq. 1 and 10)
silicate_basis = "single_cation"  # how the melt is counted for a mole basis
a = 2.31                   # constant term
b = -1542.0                # multiplies 1/T, T in K
c = -38.9                  # multiplies P/T, P in GPa
d = -0.25                  # multiplies dIW_factor * dIW
e = 0.0                    # multiplies NBO/T; this law has none
dIW_factor = 1.0
eps = { C = 6.69 }         # eps_j of the term eps_j log10(1 - X_j)
covariance = []            # the paper gives no covariance of the coefficients
P_range_GPa = [30.0, 60.0] # range of the data of the fit
T_range_K = [3080.0, 4760.0]
dIW_range = [-1.78, -0.88]
P_unit = "GPa"
source = "Tsutsumi et al. 2025, Eq. 3 (p. 3); ranges: Table 1 and Tagawa et al. 2021"
doi = "10.1038/s41467-025-65729-5"
'''
with open('my_laws.toml', 'w') as f:
    f.write(law_file)
```

The [partition file reference](../Reference/partition_file.md#file-schema) lists every field and its rule.

### Step A3: load the file and register the law

`load_partition_table` reads the file and checks every field. `register_law` stores the law under a name that the split can use. A name can be registered once per Python session; to run this step again in the same session, pass `replace=True`:

```python
table = load_partition_table('my_laws.toml')
register_law('H', 'tsutsumi2025_eq3', table['H'])
print(law_names('H'))
```

```console
['tagawa2021_eq2', 'tsutsumi2025_eq3']
```

### Step A4: check the law against the paper

A wrong number in the file is not a format error: it loads without an error and gives a wrong $D$. So compute Eq. 3 of the paper by hand at one point and compare it with the law. The point is 50 GPa, 3500 K, $\Delta\mathrm{IW} = -1$ and a metal without C:

```python
by_hand = 2.31 - 1542.0 / 3500.0 - 38.9 * 50.0 / 3500.0 - 0.25 * (-1.0)
law = get_law('H', 'tsutsumi2025_eq3')
by_law = log10_D(law, 50.0, 3500.0, -1.0, x_metal={'C': 0.0})
print(f'by hand {by_hand:.4f}, by the law {by_law:.4f}')
```

```console
by hand 1.5637, by the law 1.5637
```

### Step A5: use the law in the split

Pass the name in `laws`. `laws` also takes a loaded entry such as `table['H']` directly, without `register_law`. The other elements keep their default laws:

```python
bulk = {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}  # kg
M_metal, M_silicate = 1.94e24, 4.03e24  # kg, the Earth's core and mantle
metal = {'O': 0.02, 'Si': 0.03}  # metal mass fractions of O and Si
for name in ('tagawa2021_eq2', 'tsutsumi2025_eq3'):
    core, mantle, flags = split_core_mantle(
        bulk, M_metal, M_silicate, dIW=-1.5, w_fixed=metal, laws={'H': name}
    )
    share, flag = core['H'] / bulk['H'], flags['H']['out_of_range']
    print(f'{name}: H core share {share:.3f}, out of range {flag}')
```

```console
tagawa2021_eq2: H core share 0.949, out of range 1
tsutsumi2025_eq3: H core share 0.961, out of range 0
```

The split also logs one warning for each law outside its range, as in the [tutorial](../Tutorials/core_split.md#step-1-split-at-the-default-point). The new law puts a little more H into the core. It is inside its range at this point, because its data reach to $\Delta\mathrm{IW} = -1.78$; the default H law is 0.17 log10 units outside its $\Delta\mathrm{IW}$ range.

## Example B: a law as a Python function

A function law can have any form. The example is the S law of [Suer et al. (2017)](https://doi.org/10.1016/j.epsl.2017.04.016), Eq. 7. CALLIOPE already holds this law as the entry `suer2017_eq7`, so you can check your function against it. Example B runs on its own: its first block has all the imports it needs.

### Step B1: the arguments of a function law

CALLIOPE calls the function as `law(p_gpa, t_k, d_iw, x_metal, x_silicate, nbo_t)`:

| Argument | What it holds | Unit |
|---|---|---|
| `p_gpa` | pressure of the equilibration point | GPa |
| `t_k` | temperature of the equilibration point | K |
| `d_iw` | oxygen fugacity relative to the iron-wüstite buffer, $\Delta\mathrm{IW}$ | log10 units |
| `x_metal` | dict of the mole fractions of the light elements in the metal, over the metal without H and C (the $x'$ basis). Fe is not in the dict: Fe is the rest, $x'_\mathrm{Fe} = 1 - \sum x'_j$ over the elements other than H and C | none |
| `x_silicate` | dict of the oxide mole fractions of the melt on the single-cation basis, from the `oxides` of the split | none |
| `nbo_t` | NBO/T of the melt | none |

The function returns $\log_{10} D$ on the basis that you declare in Step B3. Where it cannot evaluate a point, it must raise `ValueError`.

### Step B2: write the function

Eq. 7 of Suer et al. (2017), with log10 in every term and the coefficients of their Table 3, is

$$
\log_{10} D_\mathrm{S} = -3.50 + \frac{3000}{T} + 33\,\frac{P}{T} + \log_{10} X_\mathrm{FeO} - \log_{10} C_\mathrm{S} + 14 \log_{10}(1 - X_\mathrm{O}),
$$

with $\log_{10} C_\mathrm{S} = -5.372$ (their Fig. 6, a pyrolite melt). The melt FeO comes from $\Delta\mathrm{IW}$ and the Fe of the metal, $X_\mathrm{FeO} = x'_\mathrm{Fe}\,10^{\Delta\mathrm{IW}/2}$ (see [the redox state](../Explanations/metal.md#the-inputs)):

```python
import math

from calliope.core_split import split_core_mantle
from calliope.metal import get_law, mass_D, register_law


def suer2017_s(p_gpa, t_k, d_iw, x_metal, x_silicate, nbo_t):
    """log10 D_S of Suer et al. (2017), Eq. 7, on the mole basis."""
    # Fe is not in x_metal: it is the rest of the metal without H and C
    x_fe = 1.0 - sum(x for el, x in x_metal.items() if el not in ('H', 'C'))
    x_feo = x_fe * 10.0 ** (d_iw / 2.0)  # melt FeO from dIW, ideal mixing
    if not 0.0 < x_feo <= 1.0:
        raise ValueError(f'X_FeO {x_feo:g} outside (0, 1] at dIW {d_iw:g}')
    x_o = x_metal.get('O', 0.0)
    return (
        -3.50 + 3000.0 / t_k + 33.0 * p_gpa / t_k
        + math.log10(x_feo) + 5.372 + 14.0 * math.log10(1.0 - x_o)
    )
```

### Step B3: register the function

A function has no file fields, so you declare them as keywords. `basis`, `P_range_GPa`, `source` and `doi` are required, and `silicate_basis` is required for a mole basis. If you leave out `T_range_K` or `dIW_range`, that range is unknown, and every use of the law is flagged out of range:

```python
register_law(
    'S',
    'my_suer2017',
    suer2017_s,
    basis='mole',  # D is a ratio of mole fractions over the whole metal (their Eq. 1)
    silicate_basis='single_cation',  # needed for a mole basis
    P_range_GPa=(1.0, 91.0),  # range of the data of the fit, GPa
    T_range_K=(1873.0, 4100.0),  # K
    dIW_range=(-1.7, -1.0),  # log10 units
    source='Suer et al. 2017, EPSL 469, Eq. 7, Table 3, Fig. 6',
    doi='10.1016/j.epsl.2017.04.016',
)
```

### Step B4: check the function

Compare the function with the packaged entry of the same law at one metal composition. `mass_D` gives $D$ on the mass basis, as the split uses it:

```python
w_metal = {'O': 0.02, 'Si': 0.03, 'S': 0.02}  # metal mass fractions; the rest is Fe
mine = mass_D(get_law('S', 'my_suer2017'), 60.0, 4000.0, -1.5, w_metal, {})
packaged = mass_D(get_law('S', 'suer2017_eq7'), 60.0, 4000.0, -1.5, w_metal, {})
print(f'mass-basis D_S: function {mine:.4f}, packaged {packaged:.4f}')
```

```console
mass-basis D_S: function 78.2195, packaged 78.2195
```

For a law that CALLIOPE does not hold, check the function against a value of the paper instead, as in Step A4.

### Step B5: use the function in the split

```python
bulk = {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}  # kg
M_metal, M_silicate = 1.94e24, 4.03e24  # kg, the Earth's core and mantle
metal = {'O': 0.02, 'Si': 0.03}  # metal mass fractions of O and Si
for name in ('suer2017_eq7', 'my_suer2017'):
    core, mantle, flags = split_core_mantle(
        bulk, M_metal, M_silicate, dIW=-1.5, w_fixed=metal, laws={'S': name}
    )
    print(f"{name}: S core share {core['S'] / bulk['S']:.4f}")
```

```console
suer2017_eq7: S core share 0.9744
my_suer2017: S core share 0.9744
```

## Errors you can meet

| What you did | Error | What to do |
|---|---|---|
| Registered an H entry as a law of S | `PartitionTableError: 'S': needs a law entry of 'S'` | Register the entry under its own element. |
| Used a name that is already held | `PartitionTableError: 'H': law name 'tagawa2021_eq2' is empty or taken` | Pick a new name, or pass `replace=True` to store over a law that you registered. A packaged law cannot be replaced. |
| Asked for a name that is not registered | `PartitionTableError: laws['H'] = 'tsutsumi': no law 'tsutsumi' for 'H'; known: ['tagawa2021_eq2']` | Register the law first; the message lists the names held. |
| Gave the split a law entry of another element | `PartitionTableError: laws['S'] = {...}: 'S': needs a law entry of 'S'` | Give each element a law of that element. |
| Gave a law for an element that is not in the bulk | `ValueError: laws names elements outside bulk_kg: ['H']` | Give laws only for the elements of `bulk_kg`. |
| Left out a required declaration of a function law | `PartitionTableError: 'S': 'source' must be a non-empty string` | Give `basis`, `P_range_GPa`, `source`, `doi` and, for a mole basis, `silicate_basis`. |
| Registered a function `nanlaw` that returns NaN | `ValueError: split of S: nanlaw law of S returned nan at P 60 GPa, T 4000 K, dIW -1.5` | Return a finite number. Where the law cannot be evaluated, raise `ValueError` with a clear message. |

A function must raise `ValueError` for a point that it cannot evaluate, as `suer2017_s` does above. What the split then does depends on the pass:

- in the first pass and in the relaxed passes, the split ends with `ValueError`, which names the element;
- in a later plain pass, the split goes on with relaxed passes;
- in the damped passes of the iron-rich root rule, the split keeps the first root.

Any other exception from the function ends the split at once.

## Where to go next

- [Partition file reference](../Reference/partition_file.md#named-laws-and-the-core-split): every argument of `register_law` and `split_core_mantle`.
- [Add a new partition law](../Tutorials/partition_law.md): a new law for the metal reservoir of the solve, and what to change in the code for a law in a new form.
- [Core and mantle split](../Explanations/metal.md#core-and-mantle-split-at-an-equilibration-point): the physics of the split.
