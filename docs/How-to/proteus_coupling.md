# Coupling to PROTEUS

This page is the **practical recipe** for using CALLIOPE inside a PROTEUS coupled run: which TOML blocks matter, which knobs are documented where, and which combinations are known to be safe.

For the underlying control flow (per-iteration sequence diagram, what the wrapper passes to `equilibrium_atmosphere()`, how warm-starts are constructed, what gets written back into `hf_row`), see the [theory page](../Explanations/proteus_coupling.md).

## Minimal `[outgas]` block

```toml
[outgas]
module       = "calliope"     # or "atmodeller" or "dummy"
fO2_shift_IW = 0.0            # log10 shift relative to IW buffer
mass_thresh  = 1e16           # kg; doubles as the absolute mass-balance
                              # tolerance and the per-element zero threshold
T_floor      = 700.0          # K, magma temperatures below this are clipped
solver_rtol  = 1e-4           # relative tolerance on mass balance
solver_atol  = 1e-6           # fsolve step tolerance (mapped to xtol, not atol)
h2_binodal   = false          # apply Rogers+2025 H2-MgSiO3 binodal override

  [outgas.calliope]
  include_H2O = true
  include_CO2 = true
  include_N2  = true
  include_S2  = true
  include_SO2 = true
  include_H2S = true
  include_NH3 = true
  include_H2  = true
  include_CH4 = true
  include_CO  = true
  solubility  = true          # if false, force Phi_global = 0 (no melt dissolution)
```

The four primary species (`H2O`, `CO2`, `N2`, `S2`) **must** be included; the wrapper will refuse to run with any of them set to `false`. The seven secondary species can be disabled individually if you want a constrained sub-system (e.g. `include_NH3 = false` to suppress ammonia chemistry).

## Choosing the elemental-budget mode

CALLIOPE needs total H, C, N, S inventories in kilograms. PROTEUS gives you two ways to specify them via `[planet]`:

```toml
[planet]
volatile_mode       = "elements"        # or "gas_prs"
volatile_reservoir  = "mantle"          # or "mantle+core"; sets the reservoir for ppmw
```

### `volatile_mode = "elements"` (most common)

```toml
[planet.elements]
use_metallicity = false      # if true, scales C/N/S to H from solar abundances

H_mode   = "oceans"          # "oceans" | "ppmw" | "kg"
H_budget = 1.0               # interpreted per H_mode

C_mode   = "C/H"             # "C/H" | "ppmw" | "kg"
C_budget = 0.1

N_mode   = "ppmw"
N_budget = 2.0

S_mode   = "ppmw"
S_budget = 200.0
```

The wrapper translates these into the four budget fields CALLIOPE expects:

- `H_mode = "oceans"` → `hydrogen_earth_oceans = H_budget`
- `H_mode = "ppmw"`   → `H_kg = H_budget * 1e-6 * M_reservoir`
- `H_mode = "kg"`     → `H_kg = H_budget`

C, N, S follow the same pattern; `C/H` is interpreted as a mass ratio relative to H.

### `volatile_mode = "gas_prs"` (initial-pressure mode)

```toml
[planet]
volatile_mode = "gas_prs"

[planet.gas_prs]
H2O = 220.0     # bar
CO2 = 100.0
N2  = 1.0
S2  = 0.01
# secondary species default to 0.0 unless explicitly set
```

When `volatile_mode = "gas_prs"`, the wrapper calls `get_target_from_pressures(ddict)` to back out the elemental inventory implied by the prescribed initial atmosphere; on every subsequent iteration the same elemental inventory is preserved.

## Selecting the fO2 dispatch

The PROTEUS schema field `[planet].fO2_source` selects which CALLIOPE entry point the wrapper calls:

```toml
[planet]
fO2_source = "user_constant"   # default; calls equilibrium_atmosphere
# fO2_source = "from_O_budget"  # calls equilibrium_atmosphere_authoritative_O
```

Under `user_constant` (the default) CALLIOPE buffers the redox state to the configured `outgas.fO2_shift_IW` and solves the four-equation H/C/N/S mass balance. Under `from_O_budget` the wrapper passes the running whole-planet O total (maintained by the PROTEUS element-budget bookkeeping) as a fifth elemental target, uses `outgas.fO2_shift_IW` only as an initial-guess hint, and solves a five-equation system that returns the derived $\Delta\mathrm{IW}$ in `hf_row['fO2_shift_IW_derived']`. See [Coupling to PROTEUS (theory)](../Explanations/proteus_coupling.md#step-5-call-the-solver) for the per-iteration control flow and [Authoritative-oxygen mode](../Explanations/authoritative_oxygen.md) for the augmented mass balance.

## Redox state

`fO2_shift_IW` is in $\log_{10}$ units relative to the [O'Neill & Eggins (2002)](https://ui.adsabs.harvard.edu/abs/2002ChGeo.186..151O) IW buffer. Under `fO2_source = "user_constant"` this value is the buffer offset; under `fO2_source = "from_O_budget"` it is only the initial-guess seed. Common reference values:

| $\Delta\mathrm{IW}$ | Description |
|---|---|
| $-5$ | Highly reduced (Mercury-like; sulphur-derived estimate IW-5.4, [Cartier & Wood 2019](https://ui.adsabs.harvard.edu/abs/2019Eleme..15...39C)) |
| $-3$ | Reduced (e.g. enstatite-chondrite-like; Mercury Fe-based estimate IW-2.8 to IW-4.5, [Cartier & Wood 2019](https://ui.adsabs.harvard.edu/abs/2019Eleme..15...39C)) |
| $-1$ | Moderately reduced; near the Mars-mantle source range ([Wadhwa 2001](https://ui.adsabs.harvard.edu/abs/2001Sci...291.1527W) places the shergottite-source mantle at $\approx$ IW) |
| $0$  | At iron-wüstite buffer (core formation equilibrium at depth) |
| $+3.5$ | [Sossi et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020SciA....6.1387S) preferred Earth's mantle $f_{\mathrm{O}_2}$ |
| $+4$ | CALLIOPE PROTEUS-side default; near-modern Earth upper mantle (within FMQ$\,\pm\,2$ per [Frost & McCammon 2008](https://ui.adsabs.harvard.edu/abs/2008AREPS..36..389F)) |

[Sossi et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020SciA....6.1387S) place Earth's modern upper mantle at $\Delta\mathrm{IW} \approx +3.5$; [Frost & McCammon (2008)](https://ui.adsabs.harvard.edu/abs/2008AREPS..36..389F) report a broader FMQ$\,\pm\,2$ range across mantle settings (approximately IW+1.5 to IW+5.5). CALLIOPE defaults sit at $\Delta\mathrm{IW} = 4.0$, consistent with a modern terrestrial composition.

## Solver tolerances

The PROTEUS wrapper hard-codes a few solver knobs that override the CALLIOPE library defaults:

| Wrapper override | Value | Reason |
|---|---|---|
| `nguess` | $10^3$ | The PROTEUS warm-start almost always lands on the right basin in the first attempt; $10^3$ Monte-Carlo restarts is plenty of margin without spending wall-time on degenerate cases. |
| `nsolve` | $3 \times 10^3$ | Per-restart `fsolve` iteration cap. |
| `opt_solver` | `False` | Skip the alternating-solver fallback; if `fsolve` cannot converge from a warm start, it almost always means an upstream bug rather than a basin-of-attraction problem. |
| `print_result` | `False` | Suppress the per-call INFO log; the PROTEUS main loop already prints the partial pressures. |

The TOML field names map to `equilibrium_atmosphere` keyword arguments as follows:

| TOML field | Solver argument | Role |
|---|---|---|
| `solver_rtol` | `rtol` | Relative tolerance on the elemental mass-balance residual. |
| `solver_atol` | `xtol` | Step tolerance on the inner `fsolve` Powell-hybrid iteration. The name is historical; it does *not* set the absolute mass-balance tolerance. |
| `mass_thresh` | `atol` | Absolute tolerance on mass balance, in kg. Also the per-element threshold below which a species' inventory is treated as zero. |

## Warm-start behaviour

For `Time > 1` yr, the wrapper builds `p_guess` from the previous-iteration `H2O_bar`, `CO2_bar`, `N2_bar`, `S2_bar` rows. For `Time = 0` (first call) the wrapper passes `p_guess = None` and lets CALLIOPE's Monte-Carlo initial-guess generator find the basin. If an element's target inventory drops below `mass_thresh`, the corresponding `p_guess` is forced to zero so that the solver does not waste iterations chasing a vanishing species.

## Common pitfalls

- **`Missing required volatile`** at startup: one of H2O / CO2 / N2 / S2 is set to `include = false`. Fix: re-enable it (it can still be set to a tiny initial pressure if you want to suppress its mass).
- **`Could not find solution for volatile abundances`** mid-run: CALLIOPE exhausted its `nguess` Monte-Carlo restarts. Almost always caused by an upstream NaN or unphysical state in `hf_row` (e.g. `T_magma < T_floor` after AGNI failed, `M_mantle = 0` after a structure-solve failure). Inspect `proteus_00.log` for the iteration immediately before the CALLIOPE failure; do not bump `nguess` blindly.
- **Atmosphere mass collapses to zero** without escape accounting it: the `check_desiccation` gate in `proteus.outgas.wrapper` will refuse to mark the planet desiccated if the loss is unexplained by cumulative escape. This is by design: a real desiccation event must be supported by tracked atmospheric escape, so the gate firing means an upstream module silently zeroed the atmosphere and the failure is upstream of CALLIOPE.
- **`Hydrogen inventory must be > 0`**: `H_budget = 0` or `H_mode` mistyped. CALLIOPE refuses to solve a zero-H system because the speciation is degenerate.
- **Solubility off but `Phi_global > 0`**: not an error, but the wrapper will silently force `Phi_global = 0` when `[outgas.calliope].solubility = false`, so dissolved masses will all be zero regardless of the live melt fraction.

## Next step

For what the wrapper actually does on each iteration (sequence diagram, mapping table, hf_row keys), read [Coupling to PROTEUS (theory)](../Explanations/proteus_coupling.md).
