# Core and mantle split

This tutorial splits the bulk H, C, N and S of a planet between its core and its mantle at one equilibration point, the step at which the core forms. It uses `calliope.core_split.split_core_mantle`.

By the end of it you will:

- have split an Earth-like bulk at the default point, 60 GPa and 4000 K;
- know how to read the out-of-range flag and the distance outside each range;
- have seen what a point at a core-mantle boundary pressure does to the S split;
- have added your own S law by name and used it in the split;
- have checked the element ledger, bulk = core + mantle.

You should already have done the [metal reservoir tutorial](metal_reservoir.md). The physics is in the section [Core and mantle split at an equilibration point](../Explanations/metal.md#core-and-mantle-split-at-an-equilibration-point), and every input and output is on the [partition file reference](../Reference/partition_file.md#named-laws-and-the-core-split). Run the steps in order in one Python session. The whole tutorial runs in a few seconds.

## Step 1: split at the default point

The bulk masses are example values. The metal and silicate masses are those of the Earth's core and mantle. The metal holds 2 wt% O and 3 wt% Si, which the split does not change.

```python
from calliope.core_split import split_core_mantle
from calliope.metal import register_law

bulk = {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}  # kg
M_metal, M_silicate = 1.94e24, 4.03e24  # kg, Earth's core and mantle
metal = {'O': 0.02, 'Si': 0.03}  # metal mass fractions of O and Si

core, mantle, flags = split_core_mantle(bulk, M_metal, M_silicate, dIW=-1.5, w_fixed=metal)
for el in bulk:
    print(
        f'{el}: core {core[el]:.3e} kg, mantle {mantle[el]:.3e} kg, '
        f'core share {core[el] / bulk[el]:.3f}, out of range {flags[el]["out_of_range"]}'
    )
```

```console
H: core 1.423e+20 kg, mantle 7.705e+18 kg, core share 0.949, out of range 1
C: core 1.396e+18 kg, mantle 1.360e+19 kg, core share 0.093, out of range 1
N: core 7.285e+18 kg, mantle 7.146e+17 kg, core share 0.911, out of range 1
S: core 7.795e+20 kg, mantle 2.047e+19 kg, core share 0.974, out of range 0
```

Most of the H, N and S goes into the core, and most of the C stays in the mantle.

The H, C and N laws are flagged, and each logs one warning with its distance outside the range:

```console
H partition law tagawa2021_eq2 evaluated outside its range: P_eq 60 GPa, T_eq 4000 K, dIW -1.5 (distance {'P_GPa': 0.0, 'T_K': 0.0, 'dIW': 0.17})
C partition law fischer2020_eq1 evaluated outside its range: P_eq 60 GPa, T_eq 4000 K, dIW -1.5 (distance {'P_GPa': 1.0, 'T_K': 0.0, 'dIW': 0.0})
N partition law suer2023_line evaluated outside its range: P_eq 60 GPa, T_eq 4000 K, dIW -1.5 (distance {'P_GPa': 40.0, 'T_K': 1000.0, 'dIW': 0.0})
```

**How to read the flags.** The H law comes from runs at dIW -1.33 to -0.88, so dIW -1.5 is 0.17 log10 units outside. The C law has data up to 59 GPa, 1 GPa below the point. The N law is a line at 20 GPa and 3000 K only, so every other point is outside. The split still uses each law at the given point: a flag tells you to read the result with care, it does not change it.

## Step 2: the redox state

$\Delta\mathrm{IW}$ is the one redox input of the split. The S law reads the FeO of the melt, and the split takes it from $\Delta\mathrm{IW}$ and the iron of the metal, $X_\mathrm{FeO} = x'_\mathrm{Fe}\,10^{\Delta\mathrm{IW}/2}$, so a more oxidised equilibrium puts more S into the core:

```python
for diw in (-3.0, -0.5):
    c, m, f = split_core_mantle(bulk, M_metal, M_silicate, dIW=diw, w_fixed=metal)
    print(f'dIW {diw}: S core share {c["S"] / bulk["S"]:.3f}, metal {f["metal"]["beyond_FeS"]}')
```

```console
S partition law suer2017_eq7 evaluated outside its range: P_eq 60 GPa, T_eq 4000 K, dIW -3 (distance {'P_GPa': 0.0, 'T_K': 0.0, 'dIW': 1.3})
dIW -3.0: S core share 0.871, metal 0
dIW -0.5: S core share 0.992, metal 0
```

At dIW -1.5 (Step 1) the share is 0.974. Both points lie outside the dIW range of the S law (-1.7 to -1.0), so each is flagged; the warning comes once per law and process. `flags['metal']` holds the light-element and S mole fractions of the metal; `beyond_FeS` is 1 when the metal holds more S than FeS, which a large bulk S on a small metal mass can give.

## Step 3: a point at a core-mantle boundary pressure

```python
core2, mantle2, flags2 = split_core_mantle(
    bulk, M_metal, M_silicate, dIW=-1.5, P_eq=300.0, T_eq=5500.0, w_fixed=metal
)
print('S at 300 GPa, 5500 K: core share', f'{core2["S"] / bulk["S"]:.3f}', flags2['S'])
```

```console
S at 300 GPa, 5500 K: core share 0.998 {'out_of_range': 1, 'P_GPa': 209.0, 'T_K': 1400.0, 'dIW': 0.0}
```

The S law now runs 209 GPa and 1400 K beyond its data. From 60 GPa and 4000 K to this point the $cP/T$ term of the law rises by 1.31 and the $b/T$ term falls by 0.20, so $D_\mathrm{S}$ rises and the core takes 99.8 % of the S. No measurement backs this number. This is why the default point is 60 GPa and 4000 K and not the core-mantle boundary.

## Step 4: your own S law

Register a law under a new name. This one is a constant, $D_\mathrm{S} = 100$ on the mass basis, so its result has a closed form.

```python
register_law(
    'S',
    'my_constant',
    lambda p, t, diw, x_metal, x_sil, nbo: 2.0,
    basis='mass',
    P_range_GPa=(0.0, 400.0),
    T_range_K=(1000.0, 6000.0),
    dIW_range=(-5.0, 0.0),
    source='a constant for this tutorial',
    doi='10.0000/none',
)
core3, mantle3, _ = split_core_mantle(
    bulk, M_metal, M_silicate, dIW=-1.5, laws={'S': 'my_constant'}, w_fixed=metal
)
x = M_metal / M_silicate
print(
    f'S with D = 100: core share {core3["S"] / bulk["S"]:.4f}, '
    f'closed form {100 * x / (1 + 100 * x):.4f}'
)
```

```console
S with D = 100: core share 0.9796, closed form 0.9796
```

A function law gets the pressure in GPa, the temperature in K, $\Delta\mathrm{IW}$, the metal mole fractions, the melt oxide fractions and the NBO/T of the melt, and returns $\log_{10} D$. A second `register_law` under the same name needs `replace=True`. To add a law from a paper as a file entry instead, see [Add a new partition law](partition_law.md) and pass the loaded entry to `register_law`.

## Step 5: the element ledger

```python
print(max(abs(core[el] + mantle[el] - bulk[el]) / bulk[el] for el in bulk))
```

```console
2.1845333333333333e-16
```

The split checks this ledger itself on every call and raises `RuntimeError` above $10^{-12}$; here the largest error is 2e-16 of the bulk, the round-off of one addition.
