# Metal reservoir and partition file

With `metal_mass` above zero, `equilibrium_atmosphere` and `equilibrium_atmosphere_authoritative_O` add a metal reservoir in equilibrium with the melt. The element balance of H, C, N and S gains a core term $D_i M_\mathrm{metal} w_i^\mathrm{melt}$, with $D_i$ the mass-basis metal over melt partition coefficient. O and Si in the metal are set by a mass fraction or by their exchange reaction. The partition laws come from a TOML file; CALLIOPE ships a placeholder file. The validation cases are on the [metal reservoir validation page](../Validation/metal.md).

## Inputs in `ddict`

| Key | Unit | Default | Meaning |
|---|---|---|---|
| `metal_mass` | kg | absent (no metal) | Mass of the metal that equilibrates with the present melt; 0 is the solve without metal; above 0 needs melt (`Phi_global` above 0) |
| `P_eq`, `T_eq` | GPa, K | required | Pressure and temperature of the metal-melt equilibrium |
| `dIW_eq_mode` | | `'parameter'` | `'parameter'` reads `dIW_eq`; `'from_FeO'` uses $2\log_{10}(x_\mathrm{FeO}/x'_\mathrm{Fe})$ with ideal activities, $x_\mathrm{FeO}$ on the single-cation basis of the melt and $x'_\mathrm{Fe}$ over the metal without H and C (needs FeO above 0) |
| `dIW_eq` | log10 units | required for `'parameter'` | Oxygen fugacity of the equilibrium relative to the iron-wustite buffer |
| `O_mode_core`, `Si_mode_core` | | `'fixed'` | `'fixed'` reads `w_O_core`, `w_Si_core`; `'partition'` uses the exchange $K_D$ of the file at the dIW of each pass |
| `w_O_core`, `w_Si_core` | mass fraction | 0 | Metal O and Si in `'fixed'` mode, and the start metal of the passes in `'partition'` mode; at least 0, each and their sum below 1 |
| `silicate_oxides` | wt% | KLB-1 (7 major oxides, Takahashi 1986) | Melt composition for the mole basis of $D$; sums to 95-101 wt% and is used as given, not normalised; Si `'partition'` needs SiO2 above 0 |
| `silicate_nbo_t` | | 2.6 | NBO/T of the melt |
| `partition_table` | | packaged placeholder | A table from `calliope.metal.load_partition_table` |
| `tol_w_core` | | 1e-6 | Relative stop of the fixed point on each metal mass fraction, above 0 and below 1 |

Every key except `partition_table` is checked before the first solve: a value of the wrong type, NaN, infinity or a value out of range raises `ValueError` naming the key, and with `metal_mass` above zero `P_eq`, `T_eq` and, in `'parameter'` mode, `dIW_eq` are required. The solve runs a fixed point: it evaluates $D$ on the metal and melt of the previous pass with the O and Si of the present pass, solves at fixed $D$, and stops when each metal mass fraction moves by at most `tol_w_core` of its value plus 1e-15 (at most 30 passes). A partitioned O or Si moves by $1/(1 - s)$ of its exchange step, with $s$ the secant slope of the exchange value against the metal value of the last 2 passes when $s < 0$, and 1 otherwise: with dIW from FeO, more metal Si raises dIW and lowers the exchange Si, and whole steps alternate about the solution; with O and Si both partitioned each secant also holds the move of the other element. The moves of a pass are halved until they take at most half of the iron left in the metal. A first exchange from pure Fe metal may leave no iron (Si at 4800 K with dIW from FeO, x_Si 1.10); the pass then takes its step, and with O and Si at 5000 K the halving keeps the metal in the range of the laws. An exchange that leaves no iron and does not fall from the previous pass raises `ValueError` naming the partitioned elements and their mode fields: no fixed point was reached from the start metal (`w_O_core`, `w_Si_core`), and another start metal can reach one; with Si from FeO at 6500 K the sum stays above 1 for 2 passes while it falls to its fixed point. This refusal is a stop rule, not a proof that no metal with iron solves the exchange; with dIW from FeO and the input of the worked example on the [metal reservoir page](../Explanations/metal.md#worked-example) (`metal_mass` 1.94e24 kg, P_eq 50 GPa, KLB-1 melt), O and Si both partitioned converge at 5500 K and raise this `ValueError` at 6000, 6500, 7000 and 7500 K. O alone partitioned converges up to 8480 K, ends at the 30-pass cap with `RuntimeError` from 8490 to 9700 K, and raises the `ValueError` from 9720 K (at 9720, 9740, 9780 and 9800-10000 K; other solves between 9710 and 9790 K ran longer than 300 s). For O alone the exchange asks for $x_\mathrm{O} = K_D(\mathrm{O})\,x_\mathrm{FeO}/(1 - x_\mathrm{O})$, which has a root only while $K_D(\mathrm{O})\,x_\mathrm{FeO} \le 1/4$ (here up to 9470 K); the passes stop converging already at $K_D(\mathrm{O})\,x_\mathrm{FeO} = 0.228$. For these elements the stop reads the whole step. A law whose $D$ or $K_D$ is beyond the float range raises `ValueError` that names the law and gives its P, T and dIW; an O and Si exchange that is not finite raises its own `ValueError`. The dIW term of each pass uses the metal before the O and Si of that pass move, one pass behind, which is below `tol_w_core` at convergence. `tol_w_core` bounds the last step of the passes, not the error of the result; when the passes contract slowly the error can be larger than the last step. The metal takes the H, C, N and S of the melt in proportion to `metal_mass`, so the caller sets `metal_mass` to the metal that equilibrates with the present melt.

`dIW_eq` and the `'from_FeO'` dIW enter only the partition laws and the O and Si exchange, and neither is converted to an absolute fO2. The IW buffer of the solve (`fO2_model`) turns `fO2_shift_IW` into the absolute fO2 of the gas and melt chemistry. It reaches the metal through the metal and melt composition that the solve returns, so a solved `'from_FeO'` dIW changes with the buffer, and with it the $D$ of each law with a dIW term; the $D$ of a law with a composition term or a mole basis changes with the buffer in either mode. With a fixed composition the dIW is the same on every buffer.

The metal O comes from the silicate FeO, not from the volatile O budget: in the authoritative-O solve gas and melt O close on the O budget, and `O_kg_core` is reported apart. CALLIOPE does not hold the silicate FeO reservoir; the caller applies the exchange from `n_FeO_exchange`. `metal_mass` is held fixed while Fe moves into or out of the metal.

`M_mantle` and `metal_mass` come from the caller. Inside PROTEUS `M_mantle` is the mantle mass that its interior steps set ([code architecture](../Explanations/code_architecture.md#structurepy)). PROTEUS does not pass `metal_mass`, so its runs have no metal reservoir. A caller that takes `M_mantle` from a structure model takes the metal from the same model: `metal_mass` is the part of its core that equilibrates with the present melt.

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
| `f`, `log10_C_S` | finite numbers, optional | Add $f \log_{10} X_\mathrm{FeO} - \log_{10} C_\mathrm{S}$, with $X_\mathrm{FeO} = x'_\mathrm{Fe}\,10^{\Delta\mathrm{IW}/2}$ of the metal and dIW (`mass_D`; dIW above $-2\log_{10} x'_\mathrm{Fe}$ raises `ValueError`); without the keys the law has no such terms; refused for `exchange_KD` |
| `eps` | table of element to number | The $\varepsilon_j$ terms; empty for `exchange_KD` |
| `covariance` | empty, or a square matrix | Covariance of a, b, c, d, e and the eps terms in name order; symmetric and positive semidefinite; stored, not used by the solve |
| `P_range_GPa`, `T_range_K`, `dIW_range` | `[low, high]` | Range of the source data, P and T at 0 or above; an empty T or dIW range means unknown, and every use is flagged out of range |
| `P_unit` | `GPa` | Pressure unit |
| `source`, `doi` | non-empty strings | Origin of the law |

`load_partition_table` raises `PartitionTableError` (a `ValueError`) that names the element and field of a value that breaks the schema or that the evaluator of the element does not read, and the keys of an element table outside the fields above.

## Placeholder laws

The packaged file `src/calliope/data/partition_coefficients_placeholder.toml` has `placeholder = true`; a solve with a table marked placeholder logs one warning per process.

| Element | Law | Range of the source data | Note |
|---|---|---|---|
| H | [Tagawa et al. (2021)](https://doi.org/10.1038/s41467-021-22035-0) Eq. 2, $K_D$ on the metal without H and C | 30-60 GPa, 3080-4560 K, dIW -1.33 to -0.88 | |
| C | [Fischer et al. (2020)](https://doi.org/10.1073/pnas.1919930117) Eq. 1, mass basis, with O and S interaction terms | 1-59 GPa, 1618-5200 K, dIW -3.9 to -0.53 | |
| N | Line in dIW through $D_N$ = 1 at dIW -4 and 39 at dIW -1, the 20 GPa, 3000 K column of [Suer et al. (2023)](https://doi.org/10.3389/feart.2023.1159412) Table 1 | 20 GPa, 3000 K, dIW -4 to -1 | The 5 GPa, 2000 K column gives 5 and 265 at the same dIW; Table 1 has no N value at 60 GPa |
| S | Line in dIW through $D_S$ = 10 at dIW -4 and 173 at dIW -1, the 60 GPa, 4000 K column of Suer et al. (2023) Table 1 | 60 GPa, 4000 K, dIW -4 to -1 | The oxidized column spans dIW -1 to about +1.3; the line uses -1 |
| O | Exchange $K_D$ of Tagawa et al. (2021) SI Table 3 for the data of Siebert et al. (2012, 2013) | 35-74 GPa; T and dIW unknown | Always flagged out of range; SI Table 3 states no basis, and the Siebert et al. papers are not read for it |
| Si | As O | 35-74 GPa; T and dIW unknown | Always flagged out of range |

Suer et al. (2023) give an uncertainty of up to 20 % on their Table 1 values. In `'partition'` mode O and Si use ideal activities in metal and melt. Badro et al. (2015) fit their Si law together with activity coefficients from interaction parameters; leaving those terms out changes $D_\mathrm{Si}$ by about an order of magnitude at the end of their accretion paths, so the Si mode gives a first estimate, not a calibrated value.

## Named laws and the core split

`calliope.metal.get_law(element, name)` returns a law by name; `law_names(element)` lists the names held. The packaged names:

| Element | Name | Law | File |
|---|---|---|---|
| H | `tagawa2021_eq2` | the H law above | `partition_coefficients_placeholder.toml` |
| C | `fischer2020_eq1` | the C law above | `partition_coefficients_placeholder.toml` |
| N | `suer2023_line` | the N law above | `partition_coefficients_placeholder.toml` |
| S | `suer2023_line` | the S law above | `partition_coefficients_placeholder.toml` |
| S | `suer2017_eq7` | [Suer et al. (2017)](https://doi.org/10.1016/j.epsl.2017.04.016) Eq. 7 with log10 in every term: a -3.50, b 3000, c 33 (Table 3), f 1, `log10_C_S` -5.372 (Fig. 6, pyrolite), eps O 14; mole basis (their Eq. 1), single-cation melt basis assumed; 1-91 GPa, 1873-4100 K (p. 91), dIW -1.7 to -1.0 (p. 90) | `partition_laws_core.toml` |

`calliope.metal.register_law(element, name, law, ...)` adds a law of H, C, N or S under a new name. The guide [Add a partition law to the core split](../How-to/core_split_law.md) shows both kinds of law in full.

| Argument | Rule |
|---|---|
| `element` | H, C, N or S; O and Si raise `PartitionTableError` |
| `name` | a non-empty name not yet held for the element; `replace=True` stores the law over a registered law of the same name; a packaged law cannot be replaced |
| `law` | a `log10_linear` entry of a loaded partition file, or a function `f(p_gpa, t_k, d_iw, x_metal, x_silicate, nbo_t)` that returns $\log_{10} D$ |
| `basis`, `silicate_basis` | function only: the metal basis of $D$ (`mass`, `mole`, `mole_xprime_noHC`) and, for a mole basis, `single_cation` |
| `P_range_GPa`, `T_range_K`, `dIW_range` | function only: the ranges of the source data; an empty T or dIW range is unknown |
| `source`, `doi` | function only: non-empty strings |

The arguments of a function law:

| Argument | Holds | Unit |
|---|---|---|
| `p_gpa`, `t_k`, `d_iw` | pressure, temperature and $\Delta\mathrm{IW}$ of the point | GPa, K, log10 units |
| `x_metal` | mole fractions of the light elements of the metal over the metal without H and C ($x'$ basis); Fe is not in the dict, it is the rest, $x'_\mathrm{Fe} = 1 - \sum x'_j$ over the elements other than H and C | none |
| `x_silicate` | oxide mole fractions of the melt on the single-cation basis | none |
| `nbo_t` | NBO/T of the melt | none |

Rules of the registry:

- The registry keeps a checked copy of an entry, and `get_law` returns a copy of the entry. A function is kept and returned as it is, so its state is shared.
- A taken name, an unknown element or basis, O or Si, an entry field that breaks the schema (as in a partition file), or a declaration that breaks its rule raises `PartitionTableError`.
- A function that returns a value that is not a finite real number (a bool, a complex or an int past the float range included) raises `ValueError` with the law name and the point.
- A function must raise `ValueError` for a point that it cannot evaluate. The core split handles that error; inside the damped passes of the iron-rich root rule it keeps the first root. Any other exception from the function ends the split, also inside the damped passes.
- `calliope.metal.check_law(entry, element)` gives the same check for an entry passed to the split.

`calliope.core_split.split_core_mantle(bulk_kg, M_metal, M_silicate, *, dIW, P_eq=60.0, T_eq=4000.0, laws=None, w_fixed=None, oxides=KLB1, nbo_t=2.6)` returns `(core_kg, mantle_kg, flags)`.

| Input | Unit | Rule |
|---|---|---|
| `bulk_kg` | kg | bulk mass per element of H, C, N, S; finite and >= 0; O and Si raise `ValueError` |
| `M_metal`, `M_silicate` | kg | finite; `M_metal` >= 0, `M_silicate` > 0; with `M_metal` 0 no law is evaluated, the core is 0 and the mantle the bulk |
| `dIW` | log10 units | finite, no default; with a law that has an FeO term (the default S law) at most $-2\log_{10} x'_\mathrm{Fe}$ of the start metal (`w_fixed`; +0.06 for 2 wt% O, +0.11 for 2 wt% O and 3 wt% Si), else `ValueError`; below about -600 `ValueError` (X_FeO underflows) |
| `P_eq`, `T_eq` | GPa, K | finite and > 0 |
| `laws` | | law per element of `bulk_kg`, a name or a D-law entry of that element; elements not given take `tagawa2021_eq2`, `fischer2020_eq1`, `suer2023_line`, `suer2017_eq7`; `laws` itself must be a dict, and a key outside `bulk_kg` raises `ValueError`; a value that is not a law name or an entry that passes `check_law` raises `PartitionTableError` |
| `w_fixed` | kg per kg | metal mass fractions of elements outside `bulk_kg` (O, Si), each from 0 to below 1, sum below 1 |
| `oxides`, `nbo_t` | wt%, none | melt composition (silicate units of the mole bases) and NBO/T; the FeO of a law with `f` comes from `dIW`, $X_\mathrm{FeO} = x'_\mathrm{Fe}\,10^{\mathrm{dIW}/2}$ |

| Output | Unit | Meaning |
|---|---|---|
| `core_kg[e]`, `mantle_kg[e]` | kg | core and mantle mass of element `e`; their sum is the bulk to 1e-12 relative (plus 1e-323 kg for the round-off of a subnormal bulk) |
| `flags[e]['out_of_range']` | | 1 when the law of `e` runs outside a range or on an unknown range |
| `flags['metal']['x_light']`, `['x_S']` | | mole fractions of all light elements and of S over the whole metal |
| `flags['metal']['S_over_Fe']` | | S to Fe mole ratio of the metal |
| `flags['metal']['beyond_FeS']` | | 1 when the metal holds more S than Fe, past FeS (flag only, no clamp) |
| `flags['metal']['root_rule']` | | 1 when the damped passes of a relaxed split past FeS replaced the first root with one that holds more Fe, else 0 |
| `flags[e]['P_GPa']`, `['T_K']`, `['dIW']` | GPa, K, log10 units | distance outside each range, 0 inside, NaN for an unknown range |
| `flags[e]['D_residual']` | | relative difference between the core mass fraction of `e` in the metal and $D$ at the returned split times its silicate mass fraction, 0 for no core mass |

Errors of the split ([how the solver finds the split](../Explanations/metal.md#how-the-solver-finds-the-split) gives the steps):

- A bool, a string, NaN or inf as an input raises `ValueError` with the name of the input.
- Relaxed passes that leave Fe below $10^{-8}$ by mass start again as common passes.
- The split raises `ValueError` (the iron runs out) when the common passes leave Fe below $10^{-8}$, or when they do not settle after the first passes ran out of iron.
- A relaxed return with more S than Fe in moles runs up to 20000 damped passes from the start metal (20000 more evaluations of the laws when they do not settle). The split takes their root when it holds more Fe and is another metal.
- A metal composition that does not settle in 50 plain and 2 x 1000 relaxed passes, or a broken ledger, raises `RuntimeError`.
