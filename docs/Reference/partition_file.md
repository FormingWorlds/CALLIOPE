# Metal reservoir and partition file

With `metal_mass` above zero, `equilibrium_atmosphere` and `equilibrium_atmosphere_authoritative_O` add a metal reservoir in equilibrium with the melt. The element balance of H, C, N and S gains a core term $D_i M_\mathrm{metal} w_i^\mathrm{melt}$, with $D_i$ the mass-basis metal over melt partition coefficient. O and Si in the metal are set by a mass fraction or by their exchange reaction. The partition laws come from a TOML file; CALLIOPE ships a placeholder file. The validation cases are on the [solve validation page](../Validation/solve.md#metal-reservoir).

## Inputs in `ddict`

| Key | Unit | Default | Meaning |
|---|---|---|---|
| `metal_mass` | kg | absent (no metal) | Mass of the metal that equilibrates with the present melt; 0 is the solve without metal; above 0 needs melt (`Phi_global` above 0) |
| `P_eq`, `T_eq` | GPa, K | required | Pressure and temperature of the metal-melt equilibrium |
| `dIW_eq_mode` | | `'parameter'` | `'parameter'` reads `dIW_eq`; `'from_FeO'` uses $2\log_{10}(x_\mathrm{FeO}/x'_\mathrm{Fe})$ with ideal activities, $x_\mathrm{FeO}$ on the single-cation basis of the melt and $x'_\mathrm{Fe}$ over the metal without H and C (needs FeO above 0) |
| `dIW_eq` | log10 units | required for `'parameter'` | Oxygen fugacity of the equilibrium relative to the iron-wustite buffer |
| `O_mode_core`, `Si_mode_core` | | `'fixed'` | `'fixed'` reads `w_O_core`, `w_Si_core`; `'partition'` uses the exchange $K_D$ of the file at the dIW of each pass |
| `w_O_core`, `w_Si_core` | mass fraction | 0 | Metal O and Si in `'fixed'` mode, at least 0, each and their sum below 1 |
| `silicate_oxides` | wt% | KLB-1 (7 major oxides, Takahashi 1986) | Melt composition for the mole basis of $D$; sums to 95-101 wt% and is used as given, not normalised; Si `'partition'` needs SiO2 above 0 |
| `silicate_nbo_t` | | 2.6 | NBO/T of the melt |
| `partition_table` | | packaged placeholder | A table from `calliope.metal.load_partition_table` |
| `tol_w_core` | | 1e-6 | Relative stop of the fixed point on each metal mass fraction, above 0 and below 1 |

Every key is checked before the first solve: a value of the wrong type, NaN, infinity or a value out of range raises `ValueError` naming the key, and with `metal_mass` above zero `P_eq`, `T_eq` and, in `'parameter'` mode, `dIW_eq` are required. The solve runs a fixed point: it evaluates $D$ on the metal and melt of the previous pass with the O and Si of the present pass, solves at fixed $D$, and stops when each metal mass fraction moves by at most `tol_w_core` of its value plus 1e-15 (at most 30 passes). A partitioned O or Si moves by $1/(1 - s)$ of its exchange step, with $s$ the secant slope of the exchange value against the metal value of the last 2 passes when $s < 0$, and 1 otherwise: with dIW from FeO, more metal Si raises dIW and lowers the exchange Si, and whole steps alternate about the solution; for these elements the stop reads the whole step. A law whose $D$ or $K_D$ is beyond the float range raises `ValueError` with its P, T and dIW. The metal takes the H, C, N and S of the melt in proportion to `metal_mass`, so the caller sets `metal_mass` to the metal that equilibrates with the present melt.

The metal O comes from the silicate FeO, not from the volatile O budget: in the authoritative-O solve gas and melt O close on the O budget, and `O_kg_core` is reported apart. CALLIOPE does not hold the silicate FeO reservoir; the caller applies the exchange from `n_FeO_exchange`. `metal_mass` is held fixed while Fe moves into or out of the metal.

## Outputs

| Key | Unit | Meaning |
|---|---|---|
| `<e>_kg_core` | kg | Element mass in the metal, for H, C, N, S, O and Si |
| `w_<e>_core` | mass fraction | Element mass fraction of the metal |
| `<e>_D_metal` | | Mass-basis $D$ of H, C, N and S in the last pass |
| `<e>_D_out_of_range` | | 1 when the law of the element runs outside its (P, T, dIW) range or on an unknown range; O and Si only in `'partition'` mode |
| `dIW_eq` | log10 units | dIW of the equilibrium in the last pass |
| `n_FeO_exchange`, `FeO_kg_exchange` | mol, kg | Silicate FeO taken by the metal, $n_\mathrm{O} - 2 n_\mathrm{Si}$ of the metal (FeO = Fe + O takes 1 FeO per O; SiO$_2$ + 2 Fe = Si + 2 FeO gives back 2 per Si); negative when the metal returns FeO |
| `n_FeO_silicate` | mol | FeO of the melt, from `M_mantle * Phi_global` and the FeO of `silicate_oxides` |
| `FeO_exchange_out_of_range` | | 1 when the exchange is more than 10 % of the melt FeO, in either sign; the melt FeO and `dIW_eq` of the solve do not follow the exchange |

With `metal_mass` equal to 0 the result carries only `<e>_kg_core` = 0 for H, C, N and S; the other keys above come with `metal_mass` above zero. `<e>_kg_total` stays gas plus melt and excludes the core; the bulk inventory of an element is `<e>_kg_total + <e>_kg_core`. A law outside its range, the `'partition'` mode and an exchange above 10 % of the melt FeO each log a warning once per process.

## File schema

The file has `schema_version = 1`, a `[meta]` table with a non-empty `name` and a boolean `placeholder`, and one table per element of H, C, N, S, O and Si. A law is

$$
\log_{10} D = a + \frac{b}{T} + c\,\frac{P}{T} + d \cdot \mathrm{dIW\_factor} \cdot \Delta\mathrm{IW} + e\,(\mathrm{NBO/T}) + \sum_j \varepsilon_j \log_{10}(1 - X_j)
$$

with $T$ in K, $P$ in GPa, NBO/T the ratio of non-bridging O to tetrahedral cations of the melt, and $X_j$ the mole fraction of element $j$ over the metal without H and C, $x'_j = N_j / \sum_{k \ne \mathrm{H,C}} N_k$ (Tagawa et al. 2021, Methods Eq. 4, p. 5). Fischer et al. (2020, p. 8744) fit $X_\mathrm{O}$ and $X_\mathrm{S}$ as mole fractions in metal without H and do not name the elements of the denominator; $x'$ is the choice of CALLIOPE, the metal basis of the H law. The O and Si exchange uses the same $x'$ in the metal and the single-cation basis in the melt, as the fO2 path of the Tagawa et al. (2021) multi-stage model (p. 6).

| Field | Values | Meaning |
|---|---|---|
| `form` | `log10_linear` for H, C, N, S; `exchange_KD` for O, Si | The law above, or the exchange $K_D$ of O and Si |
| `quantity` | `D` or `K_D` for `log10_linear`; `K_D` for `exchange_KD` | What the source fits; a `log10_linear` law always gives $D$ in `basis` (the dIW term of a $K_D$ fit carries the conversion) |
| `basis` | `mass`, `mole`, `mole_xprime_noHC`; `mole` for `exchange_KD` | Metal basis of $D$; `mole_xprime_noHC` is the mole fraction over the metal without H and C |
| `silicate_basis` | `single_cation` | Melt basis of a mole-basis law (not read for `mass`) |
| `a`, `b`, `c`, `d`, `e`, `dIW_factor` | finite numbers | Coefficients of the law; an `exchange_KD` law uses a, b and c only, with d and e 0 and `dIW_factor` 1 |
| `eps` | table of element to number | The $\varepsilon_j$ terms; empty for `exchange_KD` |
| `covariance` | empty, or a square matrix | Covariance of a, b, c, d, e and the eps terms in name order; symmetric and positive semidefinite; stored, not used by the solve |
| `P_range_GPa`, `T_range_K`, `dIW_range` | `[low, high]` | Range of the source data, P and T at 0 or above; an empty T or dIW range means unknown, and every use is flagged out of range |
| `P_unit` | `GPa` | Pressure unit |
| `source`, `doi` | non-empty strings | Origin of the law |

`load_partition_table` raises `PartitionTableError` (a `ValueError`) that names the element and field of a value that breaks the schema or that the evaluator of the element does not read, and the keys of an element table outside the fields above.

## Placeholder laws

The packaged file `src/calliope/data/partition_coefficients_placeholder.toml` has `placeholder = true`.

| Element | Law | Range of the source data | Note |
|---|---|---|---|
| H | [Tagawa et al. (2021)](https://doi.org/10.1038/s41467-021-22035-0) Eq. 2, $K_D$ on the metal without H and C | 30-60 GPa, 3080-4560 K, dIW -1.33 to -0.88 | |
| C | [Fischer et al. (2020)](https://doi.org/10.1073/pnas.1919930117) Eq. 1, mass basis, with O and S interaction terms | 1-59 GPa, 1618-5200 K, dIW -3.9 to -0.53 | |
| N | Line in dIW through $D_N$ = 1 at dIW -4 and 39 at dIW -1, the 20 GPa, 3000 K column of [Suer et al. (2023)](https://doi.org/10.3389/feart.2023.1159412) Table 1 | 20 GPa, 3000 K, dIW -4 to -1 | The 5 GPa, 2000 K column gives 5 and 265 at the same dIW; Table 1 has no N value at 60 GPa |
| S | Line in dIW through $D_S$ = 10 at dIW -4 and 173 at dIW -1, the 60 GPa, 4000 K column of Suer et al. (2023) Table 1 | 60 GPa, 4000 K, dIW -4 to -1 | The oxidized column spans dIW -1 to about +1.3; the line uses -1 |
| O | Exchange $K_D$ of Tagawa et al. (2021) SI Table 3 for the data of Siebert et al. (2012, 2013) | 35-74 GPa; T and dIW unknown | Always flagged out of range; SI Table 3 states no basis, and the Siebert et al. papers are not read for it |
| Si | As O | 35-74 GPa; T and dIW unknown | Always flagged out of range |

Suer et al. (2023) give an uncertainty of up to 20 % on their Table 1 values. In `'partition'` mode O and Si use ideal activities in metal and melt. Badro et al. (2015) fit their Si law together with activity coefficients from interaction parameters; leaving those terms out changes $D_\mathrm{Si}$ by about an order of magnitude at the end of their accretion paths, so the Si mode gives a first estimate, not a calibrated value.
