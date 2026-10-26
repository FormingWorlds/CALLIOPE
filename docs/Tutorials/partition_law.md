# Add a new partition law

The [metal reservoir](../Explanations/metal.md) reads its partition laws from a TOML file. CALLIOPE comes with a file of placeholder laws. This tutorial puts a published law into your own copy of that file, checks it, and uses it in a solve (Part 1). Part 2 explains the covariance, range and basis fields. Part 3 is for developers: it says what to change in the code to add a law in a new form or for a new element.

By the end of it you will:

- have mapped each coefficient of a published equation onto the fields of the file;
- have written your own partition file and loaded it with `load_partition_table`;
- know the error that a typical mistake in the file gives;
- have compared $D$ and the metal share of the new law with the placeholder law;
- know what to fill in the covariance, range and basis fields.

You should already have done the [metal reservoir tutorial](metal_reservoir.md). The file format is on the [partition file reference](../Reference/partition_file.md#file-schema). The whole tutorial runs in a few seconds.

## Part 1: replace the law of one element

### Step 1: the published law

The example law is the H law of [Tsutsumi et al. (2025)](https://doi.org/10.1038/s41467-025-65729-5), *Nature Communications* 16, 10038. They fit their own diamond anvil cell runs together with those of Tagawa et al. (2021), the source of the placeholder H law. Their law adds a term for the C in the metal, which lowers $D_\mathrm{H}$.

Their Eq. 3 (PDF p. 3) is

$$
\log_{10} K_D = \log_{10} D_\mathrm{H}^\mathrm{metal/silicate} + \frac{1}{4}\Delta\mathrm{IW} = a + \frac{b}{T} + c\,\frac{P}{T} + d \cdot \log_{10}\left(1 - x'^\mathrm{metal}_\mathrm{C}\right),
$$

with $a = 2.31(20)$, $b = -1542(500)$, $c = -38.9(78)$ and $d = 6.69(30)$ (p. 3; the numbers in brackets are the 1-sigma errors of the last digits). $T$ is in K and $P$ in GPa. The other definitions:

- $D_\mathrm{H} = x'^\mathrm{metal}_\mathrm{H} / x^\mathrm{silicate}_\mathrm{HO_{0.5}}$, a ratio of mole fractions (Eq. 1, p. 2). The note under their Table 1 (p. 2) says: "D$_\mathrm{H}$ and D$_\mathrm{C}$ are based on molar concentrations".
- $x'_i = N_i / \sum_{k \ne \mathrm{H,C}} N_k$, the mole fraction over the metal without H and C (Eq. 10, p. 8).

### Step 2: map the coefficients onto the file

The file needs $\log_{10} D$ alone on the left. Move the term $\frac{1}{4}\Delta\mathrm{IW}$ to the right:

$$
\log_{10} D_\mathrm{H} = 2.31 - \frac{1542}{T} - 38.9\,\frac{P}{T} - 0.25\,\Delta\mathrm{IW} + 6.69\,\log_{10}(1 - x'_\mathrm{C}).
$$

Compare it term by term with the law of the file, $\log_{10} D = a + b/T + c\,P/T + d \cdot \mathrm{dIW\_factor} \cdot \Delta\mathrm{IW} + e\,(\mathrm{NBO/T}) + \sum_j \varepsilon_j \log_{10}(1 - X_j)$:

| Paper (Eq. 3) | File field | Value | Why |
|---|---|---|---|
| paper $a$ | `a` | 2.31 | same term |
| paper $b$ | `b` | -1542.0 | same term, $T$ in K |
| paper $c$ | `c` | -38.9 | same term, $P$ in GPa |
| $\frac{1}{4}\Delta\mathrm{IW}$ on the left | `d`, `dIW_factor` | -0.25, 1.0 | moved to the right, so the sign changes |
| no NBO/T term | `e` | 0.0 | |
| paper $d \cdot \log_{10}(1 - x'_\mathrm{C})$ | `eps` | `{ C = 6.69 }` | CALLIOPE takes $X_j$ on the same $x'$ basis (Eq. 10) |
| fits $K_D$ | `quantity` | `"K_D"` | what the source fits; the law still gives $D$, because the $\Delta\mathrm{IW}$ term holds the conversion |
| $D$ over $x'$ in the metal | `basis` | `"mole_xprime_noHC"` | Eq. 1 and Eq. 10 |
| H in the melt as HO$_{0.5}$ | `silicate_basis` | `"single_cation"` | see below |

The letters of the paper are not the field names of the file: the paper's $d$ multiplies $\log_{10}(1 - x'_\mathrm{C})$ and goes into the file field `eps`, while the file field `d` multiplies $\Delta\mathrm{IW}$.

`silicate_basis`: Tsutsumi et al. count H in the melt as HO$_{0.5}$ units. CALLIOPE counts the melt on the single-cation basis, one unit per cation, and H as HO$_{0.5}$. This is the basis with which CALLIOPE reproduces the single-stage core H and C of Tsutsumi et al. (2025, SI Table 1) on the [validation page](../Validation/metal.md#tsutsumi-et-al-2025-single-stage-cores).

The ranges come from the runs of the fit:

- Tsutsumi et al. (2025) Table 1 (p. 2), the 5 runs with a $D_\mathrm{H}$ (runs 1 to 5): 42-56 GPa, 3630-4760 K, $\Delta\mathrm{IW}$ $-1.78$ to $-1.31$.
- The runs of Tagawa et al. (2021), the other data of the fit: 30-60 GPa, 3080-4560 K, $\Delta\mathrm{IW}$ $-1.33$ to $-0.88$ (the ranges of the placeholder H law).
- Together: 30-60 GPa, 3080-4760 K, $\Delta\mathrm{IW}$ $-1.78$ to $-0.88$.

Tsutsumi et al. give the 1-sigma error of each coefficient, but not the covariance between them. The covariance cannot be built from the errors alone, so `covariance` stays empty (see [Part 2](#part-2-covariance-ranges-and-basis)).

### Step 3: write your own file

Copy the packaged file and replace its `[H]` table. The other laws stay the placeholders, so `placeholder` stays `true`:

```python
import logging
from importlib import resources

import numpy as np

from calliope.constants import M_core_earth, volatile_species
from calliope.metal import PartitionTableError, load_partition_table
from calliope.solve import equilibrium_atmosphere, get_target_from_params

logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

packaged = resources.files('calliope.data').joinpath('partition_coefficients_placeholder.toml')
text = packaged.read_text()

h_law = '''[H]
form = "log10_linear"
quantity = "K_D"
basis = "mole_xprime_noHC"
silicate_basis = "single_cation"
a = 2.31
b = -1542.0
c = -38.9
d = -0.25
e = 0.0
dIW_factor = 1.0
eps = { C = 6.69 }
covariance = []
P_range_GPa = [30.0, 60.0]
T_range_K = [3080.0, 4760.0]
dIW_range = [-1.78, -0.88]
P_unit = "GPa"
source = "Tsutsumi et al. 2025, Nat. Commun. 16, 10038, Eq. 3 (p. 3); ranges of runs 1-5 of their Table 1 and of the Tagawa et al. 2021 runs"
doi = "10.1038/s41467-025-65729-5"

'''
start, end = text.index('[H]\n'), text.index('[C]\n')
mine = text[:start] + h_law + text[end:]
packaged_name = load_partition_table()['meta']['name']
mine = mine.replace(f'name = "{packaged_name}"', 'name = "tsutsumi2025-H"')
with open('my_partition_laws.toml', 'w') as f:
    f.write(mine)
```

The cut from `[H]` to `[C]` works because the packaged file has the H table first and the C table second. If you edit the file by hand, the result is the same.

### Step 4: load and check the file

`load_partition_table` reads the file and checks every field:

```python
table = load_partition_table('my_partition_laws.toml')
print(table['meta'])
h = table['H']
print({k: h[k] for k in ('a', 'b', 'c', 'd', 'eps', 'basis', 'P_range_GPa', 'T_range_K', 'dIW_range')})
```

```text
{'name': 'tsutsumi2025-H', 'placeholder': True}
{'a': 2.31, 'b': -1542.0, 'c': -38.9, 'd': -0.25, 'eps': {'C': 6.69}, 'basis': 'mole_xprime_noHC', 'P_range_GPa': (30.0, 60.0), 'T_range_K': (3080.0, 4760.0), 'dIW_range': (-1.78, -0.88)}
```

### Step 5: 2 typical mistakes

A field that breaks the file format raises `PartitionTableError`, a `ValueError` that names the element and the field. 2 mistakes that are easy to make:

- an `eps` key that is not an element name, for example `xC` for the $x'_\mathrm{C}$ of the paper;
- a basis name that the file does not know, for example `"molar"`.

```python
for wrong, right in (('eps = { xC = 6.69 }', 'eps = { C = 6.69 }'), ('basis = "molar"', 'basis = "mole_xprime_noHC"')):
    with open('broken.toml', 'w') as f:
        f.write(mine.replace(right, wrong, 1))
    try:
        load_partition_table('broken.toml')
    except PartitionTableError as err:
        print(f'{wrong}: {err}')
```

```text
eps = { xC = 6.69 }: 'H': 'eps' keys must be elements of ('H', 'C', 'N', 'S', 'O', 'Si')
basis = "molar": 'H': 'basis' must be one of ('mass', 'mole', 'mole_xprime_noHC'), got 'molar'
```

A mistake in a number is not a format error: a wrong sign of `d` loads without error and gives a wrong $D$. Check each coefficient against the paper. Where the paper gives $D$ at stated conditions, evaluate the law there with `calliope.metal.log10_D` and compare.

### Step 6: use the law in a solve

Pass the table as `ddict['partition_table']`, and solve the case of the metal reservoir tutorial with the placeholder law and with the new law:

```python
ddict = {
    'M_mantle': 4.03e24, 'gravity': 9.81, 'radius': 6.371e6,
    'T_magma': 2500.0, 'Phi_global': 1.0, 'fO2_shift_IW': 0.5,
    'hydrogen_earth_oceans': 1.0, 'CH_ratio': 0.1, 'nitrogen_ppmw': 2.0, 'sulfur_ppmw': 200.0,
}
for sp in volatile_species:
    ddict[f'{sp}_included'] = 1
    ddict[f'{sp}_initial_bar'] = 0.0
metal = {'metal_mass': M_core_earth, 'P_eq': 50.0, 'T_eq': 3500.0, 'dIW_eq': -1.0}
target = get_target_from_params(ddict)

np.random.seed(0)
for name, tab in (('placeholder', load_partition_table()), ('Tsutsumi 2025', table)):
    res = equilibrium_atmosphere(
        target, {**ddict, **metal, 'partition_table': tab}, print_result=False
    )
    share = res['H_kg_core'] / (res['H_kg_core'] + res['H_kg_total'])
    print(
        f"{name:13s}  D_H {res['H_D_metal']:.4g}  H metal share {share:.4f}  "
        f"H in gas {res['H_kg_atm']:.3e} kg  flag {res['H_D_out_of_range']}"
    )
```

```text
placeholder    D_H 22.14  H metal share 0.9142  H in gas 1.582e+15 kg  flag 0
Tsutsumi 2025  D_H 33.3  H metal share 0.9413  H in gas 7.413e+14 kg  flag 0
```

The new law gives $D_\mathrm{H}$ = 33.3, 1.5 times the placeholder value, and halves the H that is left in the gas. At these conditions the metal holds almost no C (the C law gives $D_\mathrm{C}$ = 0.08), so the C term of the new law is close to 0, and the difference comes from the new $a$, $b$ and $c$. In a metal with more C, the C term lowers $D_\mathrm{H}$. Both laws are in range here (flag 0), but the file is still marked `placeholder`, so the placeholder warning stays: the N, S, O and Si laws of the file have not changed.

## Part 2: covariance, ranges and basis

**Covariance.** `covariance` is the covariance matrix of the coefficients, in the order a, b, c, d, e and then the `eps` terms by element name, so its size is 5 plus the number of `eps` terms. It must be symmetric and positive semidefinite. Fill it only when the source gives it, as Fischer et al. (2020) do for the placeholder C law (their SI Appendix Table S3). Otherwise leave it as `[]`. CALLIOPE stores the covariance with the law and does not use it in the solve.

**Ranges.** `P_range_GPa`, `T_range_K` and `dIW_range` are the `[low, high]` ranges of the data that the law was fitted to, not the range where you want to use the law. Take them from the table of runs of the source. For a fit to several data sets, take the range over all of them. When the source does not give a temperature or $\Delta\mathrm{IW}$ range, write `[]`: the law is then flagged out of range at every use, and the result says so with `<e>_D_out_of_range` = 1. The pressure range cannot be empty.

**Basis.** `basis` must match the definition of $D$ in the source:

- `"mass"`: $D$ is a ratio of mass fractions (wt% or ppm by mass).
- `"mole"`: $D$ is a ratio of mole fractions, with the mole fraction over the whole metal.
- `"mole_xprime_noHC"`: a ratio of mole fractions, with the mole fraction over the metal without H and C.

For the 2 mole bases, `silicate_basis` must be `"single_cation"`. If the source counts the melt on another basis (for example per O atom), convert its law before you put it in the file.

**What CALLIOPE checks and what it does not.** `load_partition_table` checks the format: the field names, the types, the choices, `low <= high`, the covariance shape. It cannot check that a number is the number of the paper, or that the basis is the basis of the paper. Write the source, the equation and the page into the `source` field, so that the next reader can check them.

## Part 3: a law in a new form, or for a new element

This part is for developers. A law in the `log10_linear` form for one of H, C, N and S needs no code change: Part 1 is all. A law in a new functional form, or a law for an element that CALLIOPE does not partition yet, needs code. The checklist:

1. **`src/calliope/metal.py`**:
    - `FORMS`: add the name of the new form.
    - `KEYS` and `_entry`: add and check each new field of the form.
    - `_check_read_values`: state which forms and fields the evaluator of each element reads, so that a file with a value the code ignores is refused.
    - The evaluator: a function like `log10_D` for the new form, and the conversion to the mass basis in `mass_D` (or a new function like `exchange_metal` for an exchange reaction).
    - `BASES`, `QUANTITIES`, `SILICATE_BASES` and `COEFFICIENTS`: add a new basis, fitted quantity, melt basis or numeric coefficient; `_entry` refuses any value outside these tuples.
    - `ELEMENTS` (and `EXCHANGE_ELEMENTS` for an exchange with FeO): add a new element.
2. **`src/calliope/solve.py`**, for a new volatile element:
    - `func`, `func_authoritative_O` and `_with_metal`: a budget, a residual and a metal term for the element.
    - The element lists of `_with_metal`: the loops over H, C, N and S, and over O and Si.
    - The loops of `equilibrium_atmosphere` and `equilibrium_atmosphere_authoritative_O` that write `<e>_kg_core`: `_with_metal` reads those keys from each inner solve.
3. **The packaged file** `src/calliope/data/partition_coefficients_placeholder.toml`: a law for each new element, so the default table stays complete.
4. **The docs**: the field table and the placeholder table of [Reference/partition_file.md](../Reference/partition_file.md), and the [explanation page](../Explanations/metal.md) if the physics changes.
5. **The tests** (`tests/test_metal.py` for `metal.py`; `tests/test_metal_solve.py` and `tests/test_metal_fields.py` for `solve.py`; rules in `tests/AGENTS.md`):
    - a `reference_pinned` test (a test that runs the real code and compares it with a published value): evaluate the law at the conditions of a value that the source gives (a run of its table, or a value in its text) and assert that value within the error of the source;
    - a format test for each new field: a wrong type, NaN, a missing value and a value out of its range raise `PartitionTableError`;
    - a mutation check, which shows that the test can fail: change one coefficient sign or the basis in the code, run the test, confirm that it fails, and undo the change;
    - a case on the [validation page](../Validation/metal.md) and a row in `docs/Validation/index.md`, which the test `tests/test_validation_index.py` compares with the list of validation tests.

## Where to go next

- [Partition file reference](../Reference/partition_file.md): every field of the file and every input and output key.
- [Metal reservoir](../Explanations/metal.md): the physics and the limits.
- [Metal reservoir validation](../Validation/metal.md): the published values that the placeholder laws are checked against.
