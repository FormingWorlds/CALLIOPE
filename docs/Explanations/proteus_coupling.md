# Coupling to PROTEUS

This page is the **theory** of how CALLIOPE plugs into a PROTEUS coupled run: the per-iteration sequence, the dictionary the wrapper builds, the warm-start convention, the elemental-budget translation, and the helpfile keys CALLIOPE is responsible for. For the practical TOML recipe, see the [how-to page](../How-to/proteus_coupling.md).

## Where CALLIOPE sits in the iteration

PROTEUS runs an iterative loop in which structure, outgassing, atmosphere radiative transfer, interior thermal evolution, and atmospheric escape each take turns updating their portion of the planet's state, encoded in the `hf_row` dictionary. CALLIOPE occupies the **outgassing** slot. A typical iteration:

```
   ┌──────────────────────────────────────────────────────────────────────────┐
   │  while not done:                                                         │
   │                                                                          │
   │      # 1. Static structure (Zalmoxis or SPIDER)                          │
   │      run_structure(config, hf_row)            # writes R_int, M_mantle,  │
   │                                               # gravity, R_cmb, P_cmb    │
   │                                                                          │
   │      # 2. Outgassing (this is CALLIOPE)                                  │
   │      run_outgassing(dirs, config, hf_row)     # writes <s>_bar,          │
   │                                               # <s>_kg_atm/liquid/total, │
   │                                               # P_surf, M_atm,           │
   │                                               # atm_kg_per_mol           │
   │                                                                          │
   │      # 3. Atmosphere radiative transfer (AGNI or JANUS)                  │
   │      run_atmosphere(config, hf_row)           # writes T_surf, F_atm,    │
   │                                               # outgoing fluxes          │
   │                                                                          │
   │      # 4. Interior thermal evolution (Aragog or SPIDER)                  │
   │      run_interior(config, hf_row)             # writes T_magma, S_mantle,│
   │                                               # Phi_global, T_core,      │
   │                                               # F_int                    │
   │                                                                          │
   │      # 5. Atmospheric escape (ZEPHYRUS)                                  │
   │      run_escape(dirs, config, hf_row)         # writes esc_kg_*,         │
   │                                               # adjusts <s>_kg_total     │
   │                                                                          │
   │      # 6. Time advance                                                   │
   │      hf_row['Time'] += dt                                                │
   └──────────────────────────────────────────────────────────────────────────┘
```

CALLIOPE is the cheapest of the five (a few tens of milliseconds per call vs. seconds-to-minutes for atmosphere and interior), so it runs every iteration without gating.

## What the wrapper does

`proteus.outgas.calliope.calc_surface_pressures(dirs, config, hf_row)` is the single entry point. Its job is to translate from PROTEUS's TOML config + `hf_row` state vector into the two dicts CALLIOPE expects, call the equilibrium solver, and write the result back.

### Step 1 - build `solvevol_inp`

`construct_options(dirs, config, hf_row)` assembles the `ddict` argument:

| `solvevol_inp[…]` | Source |
|---|---|
| `M_mantle` | `hf_row['M_mantle']` (set by structure module) |
| `gravity` | `hf_row['gravity']` |
| `radius` | `hf_row['R_int']` |
| `Phi_global` | `hf_row['Phi_global']` if `config.outgas.calliope.solubility` else `0.0` |
| `T_magma` | `hf_row['T_magma']`, clipped to `config.outgas.T_floor` |
| `fO2_shift_IW` | `config.outgas.fO2_shift_IW` |
| `<s>_initial_bar` | `config.planet.gas_prs.get_pressure(s)` (only used when `volatile_mode == 'gas_prs'`) |
| `<s>_included` | `config.outgas.calliope.is_included(s)` (per-species TOML flag) |

For `volatile_mode == 'elements'`, the wrapper additionally translates the elemental-budget config:

| Field | Computation | TOML knob |
|---|---|---|
| `hydrogen_earth_oceans` | $H_\text{kg} / (N_\text{ocean} \cdot \mu_\mathrm{H_2})$ | `[planet.elements] H_mode = "oceans"\|"ppmw"\|"kg"`, `H_budget = …` |
| `CH_ratio` | $C_\text{kg} / H_\text{kg}$ | `C_mode = "C/H"\|"ppmw"\|"kg"`, `C_budget = …` |
| `nitrogen_ppmw` | $10^6 \cdot N_\text{kg} / M_\text{mantle}$ | `N_mode = "ppmw"\|"kg"`, `N_budget = …` |
| `sulfur_ppmw` | $10^6 \cdot S_\text{kg} / M_\text{mantle}$ | `S_mode = "ppmw"\|"kg"`, `S_budget = …` |

If `[planet.elements].use_metallicity = true`, C, N, S are scaled to H using the solar abundances `C_solar`, `N_solar`, `S_solar` from `proteus.utils.constants`, with `metallicity` as the multiplier.

The reservoir mass for `ppmw` conversions is `M_mantle` by default, switchable to `M_int` (= `M_mantle + M_core`) via `config.planet.volatile_reservoir = "mantle+core"`.

### Step 2 - build `target`

`calc_target_masses` (called from `wrapper.calc_target_elemental_inventories`) decides between two paths based on `volatile_mode`:

- `elements`: `target = get_target_from_params(solvevol_inp)` (translates `hydrogen_earth_oceans`, `CH_ratio`, `nitrogen_ppmw`, `sulfur_ppmw` into kg).
- `gas_prs`: `target = get_target_from_pressures(solvevol_inp)` (sums atmospheric + dissolved mass at the prescribed initial pressures).

The result is stored in `hf_row['<element>_kg_total']` for downstream code to read; CALLIOPE itself reads it back for the conservation constraint.

### Step 3 - build `p_guess`

`construct_guess(hf_row, target, mass_thresh)`:

```
if hf_row['Time'] < 1 yr:
    return None                # let CALLIOPE Monte-Carlo guess
else:
    p_guess = {
        'H2O': hf_row['H2O_bar'],   # previous-iteration partial pressures
        'CO2': hf_row['CO2_bar'],
        'N2':  hf_row['N2_bar'],
        'S2':  hf_row['S2_bar'],
    }
    # zero out species whose elemental inventory dropped below mass_thresh
    for s in ('H2O', 'CO2', 'N2', 'S2'):
        if any(target[e] < mass_thresh
               for e in 'HCNS' if e in s and e != 'O'):
            p_guess[s] = 0.0
    return p_guess
```

The previous-iteration partial pressures are an excellent warm start because the planet evolves slowly compared to the equilibrium-chemistry timescale: typically the new partial pressures differ from the old by less than 1%, well within `fsolve`'s basin of attraction.

### Step 4 - flag included volatiles

`flag_included_volatiles(p_guess, config)` reconciles the static `[outgas.calliope].include_<s>` config with the runtime state: a species that is configured-on but whose previous-iteration pressure is exactly zero is treated as off for this iteration. This stops the solver from chasing a species that has already been fully removed by escape or a previous depletion event.

### Step 5 - call the solver

The wrapper dispatches between the two CALLIOPE entry points based on `config.planet.fO2_source`:

```python
if config.planet.fO2_source == 'user_constant':
    solvevol_result = equilibrium_atmosphere(
        target,                            # H, C, N, S
        opts,
        xtol=config.outgas.solver_atol,    # fsolve step tolerance (despite the TOML name)
        rtol=config.outgas.solver_rtol,    # relative mass-balance tolerance
        atol=config.outgas.mass_thresh,    # absolute mass-balance tolerance
        nguess=config.outgas.calliope.nguess,
        nsolve=config.outgas.calliope.nsolve,
        p_guess=p_guess,
        print_result=False,
        opt_solver=False,
    )
elif config.planet.fO2_source == 'from_O_budget':
    target['O'] = hf_row['O_kg_total']     # add the fifth element budget
    solvevol_result = equilibrium_atmosphere_authoritative_O(
        target,                            # H, C, N, S, O
        opts,
        fO2_hint=config.outgas.fO2_shift_IW,
        xtol=config.outgas.solver_atol,
        rtol=config.outgas.solver_rtol,
        atol=config.outgas.mass_thresh,
        nguess=config.outgas.calliope.nguess,
        nsolve=config.outgas.calliope.nsolve,
        p_guess=p_guess,
        print_result=False,
        opt_solver=False,
    )
```

`config.outgas.calliope.nguess` and `nsolve` default to $10^3$ and $3\times 10^3$ respectively in the PROTEUS schema; override them in `[outgas.calliope]` if a particular run needs more attempts.

Under `user_constant` (the default) CALLIOPE uses the configured `outgas.fO2_shift_IW` as the buffer offset and solves for the four H/C/N/S pressures. Under `from_O_budget` the wrapper passes `hf_row['O_kg_total']` (the whole-planet oxygen total maintained by the PROTEUS element-budget bookkeeping) as the fifth target, uses `outgas.fO2_shift_IW` only as an initial-guess hint, and solves for the four pressures plus $\Delta\mathrm{IW}$. The two dispatches return result dicts with the same key schema; the authoritative-O dict additionally carries `fO2_shift_derived` and `O_res`, which the wrapper writes to `hf_row['fO2_shift_IW_derived']` and `hf_row['O_res']`. After the writeback the wrapper restores `hf_row['O_kg_total']` to the user-supplied target (the value carried by the PROTEUS element-budget bookkeeping), overriding the solver-derived total. The two values agree to within `O_res` by construction, so restoring the authoritative input prevents accumulated solver-tolerance drift from leaking into the running budget across iterations.

If either call raises `RuntimeError` (Monte-Carlo restarts exhausted), the wrapper writes status code 27 to the run's status file and re-raises; the PROTEUS main loop then handles cleanup.

### Step 6 - write back to `hf_row`

For each key in `expected_keys()` (defined in `proteus.outgas.common`), the wrapper copies the matching field from the CALLIOPE result dict into `hf_row`. The full set of CALLIOPE-owned `hf_row` keys:

| Key pattern | Meaning |
|---|---|
| `<s>_bar` | Surface partial pressure of species `s`, bar |
| `<s>_vmr` | Volume mixing ratio (== mole fraction) of `s` |
| `<s>_kg_atm` | Atmospheric column mass of `s`, kg |
| `<s>_kg_liquid` | Mass of `s` dissolved in melt, kg |
| `<s>_kg_solid` | Always 0.0 (CALLIOPE does not partition into solid) |
| `<s>_kg_total` | `_kg_atm + _kg_liquid + _kg_solid` |
| `<s>_mol_atm/liquid/solid/total` | Same in mol |
| `P_surf` | Total surface pressure, bar |
| `M_atm` | Total atmospheric mass, kg |
| `atm_kg_per_mol` | Mean molecular mass, kg/mol |
| `<e>/<e'>_atm` | Atmospheric elemental mass ratios |
| `<e>_res` | Mass-balance residual for element `e`, kg |

After the writeback, `wrapper.run_outgassing` recomputes `M_atm` from the per-species `_kg_atm` sum (defensive: catches any internal inconsistency in the result dict).

## Where the binodal lives

If `config.outgas.h2_binodal = true`, `wrapper.run_outgassing` calls `apply_binodal_h2` *after* CALLIOPE returns. This applies the [Rogers et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.544.3496R) H$_2$-MgSiO$_3$ miscibility correction as a post-processing step on top of the CALLIOPE equilibrium. CALLIOPE itself does not know about miscibility; it produces the ideal-mixing baseline that the binodal then perturbs.

When `config.interior_struct.zalmoxis.global_miscibility = true` (Zalmoxis radial binodal), the bulk binodal step is skipped because Zalmoxis has already done a per-radial-shell partition during the structure update. See the [Zalmoxis binodal page](https://proteus-framework.org/Zalmoxis/Explanations/binodal.html) for that mechanism.

## Crystallised vs. desiccated states

`wrapper.run_crystallized` and `wrapper.run_desiccated` substitute for `run_outgassing` when the planet has either solidified ($\Phi_\mathrm{global} < \phi_\text{crit}$) or lost all volatiles. Neither calls CALLIOPE: the former preserves the existing reservoirs, the latter zeros them. The desiccation gate is inactive while any element is still above `mass_thresh`, and is additionally guarded by an escape-balance check that refuses to claim desiccation if the implied mass loss exceeds `1.5 * esc_kg_cumulative + 1000 kg`. This guard prevents an upstream module failure from silently zeroing the atmosphere via CALLIOPE: a real desiccation event must be backed by accounted-for atmospheric escape.

## See also

- [Coupling to PROTEUS (how-to)](../How-to/proteus_coupling.md) for the TOML recipe and pitfalls.
- [Mass balance & solver](mass_balance.md) for what `equilibrium_atmosphere` does inside.
- [Authoritative-oxygen mode](authoritative_oxygen.md) for the augmented mass balance that `from_O_budget` dispatches to.
- The PROTEUS-side wrapper code: [`src/proteus/outgas/calliope.py`](https://github.com/FormingWorlds/PROTEUS/blob/main/src/proteus/outgas/calliope.py), [`src/proteus/outgas/wrapper.py`](https://github.com/FormingWorlds/PROTEUS/blob/main/src/proteus/outgas/wrapper.py), [`src/proteus/config/_outgas.py`](https://github.com/FormingWorlds/PROTEUS/blob/main/src/proteus/config/_outgas.py).
