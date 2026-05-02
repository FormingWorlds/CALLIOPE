# Configuration

CALLIOPE is a library, not an application: there is no TOML or YAML configuration file at the CALLIOPE level. Instead, the public function `calliope.solve.equilibrium_atmosphere(target, ddict, ...)` takes two dictionaries that the caller assembles. This page documents what those dictionaries must contain.

For the PROTEUS-side TOML schema (`[outgas]` and `[outgas.calliope]`), which is the form most users will encounter, see [Coupling to PROTEUS (how-to)](proteus_coupling.md).

## The `ddict` dictionary

`ddict` collects everything CALLIOPE needs to know about the planet, the magma ocean state, and which volatiles to include. Every key listed below is **required** unless flagged otherwise.

### Planet properties

| Key | Type | Units | Meaning |
|---|---|---|---|
| `M_mantle` | float | kg | Mass of the (molten + solid) silicate mantle. Used as the reservoir for ppmw conversions. |
| `gravity` | float | m s$^{-2}$ | Surface gravity. Used to convert column mass to partial pressure. |
| `radius` | float | m | Planetary surface radius. Used to convert column mass to total atmospheric mass. |

### Magma ocean state

| Key | Type | Units | Meaning |
|---|---|---|---|
| `T_magma` | float | K | Magma ocean temperature at the surface. Used in equilibrium constants and solubility laws. |
| `Phi_global` | float | dimensionless | Global melt fraction of the mantle, $0 \le \Phi \le 1$. Set to $0$ to disable solubility (no melt to dissolve into). |
| `fO2_shift_IW` | float | $\log_{10}$ units | Oxygen fugacity expressed as $\Delta\mathrm{IW} = \log_{10} f_{\mathrm{O}_2} - \log_{10} f_{\mathrm{O}_2}^{\mathrm{IW}}(T)$. Negative is reducing; positive is oxidising. |

### Per-species inclusion flags

For every species in `calliope.constants.volatile_species` (i.e. H$_2$O, CO$_2$, O$_2$, H$_2$, CH$_4$, CO, N$_2$, S$_2$, SO$_2$, H$_2$S, NH$_3$):

| Key | Type | Meaning |
|---|---|---|
| `<species>_included` | int (0 or 1) | If `1`, the species participates in the speciation. If `0`, its partial pressure is forced to zero. |
| `<species>_initial_bar` | float, bar | Optional. Used by `get_target_from_pressures()` when the elemental inventory is specified as initial partial pressures rather than masses. |

The four primary species (H$_2$O, CO$_2$, N$_2$, S$_2$) **must** be included; CALLIOPE will raise if any of them is flagged `0`. The seven secondary species (H$_2$, CH$_4$, CO, NH$_3$, SO$_2$, H$_2$S, O$_2$) are optional, but excluding them removes their contribution to gas-phase equilibrium and to dissolved-mass mass balance.

### Optional: elemental-budget inputs

If the calling code uses `calliope.solve.get_target_from_params(ddict)` to build the `target` dictionary (the typical workflow), the following four extra keys are also required:

| Key | Type | Units | Meaning |
|---|---|---|---|
| `hydrogen_earth_oceans` | float | Earth-ocean equivalents | Total H inventory, expressed as multiples of one present-day Earth ocean ($1.39 \times 10^{21}\,\mathrm{kg}$ of H$_2$O, equivalently $1.55 \times 10^{20}\,\mathrm{kg}$ of H). |
| `CH_ratio` | float | mass ratio | Total C / total H mass ratio (Earth value $\approx 0.1$). |
| `nitrogen_ppmw` | float | ppmw | Total N inventory in parts per million by mass, relative to `M_mantle`. |
| `sulfur_ppmw` | float | ppmw | Total S inventory in ppmw relative to `M_mantle`. |

When CALLIOPE is called from PROTEUS, the wrapper in `proteus.outgas.calliope` sets these four fields automatically based on the `[planet.elements]` block of the PROTEUS config.

## The `target` dictionary

`target` carries the elemental conservation constraints, in kilograms:

```python
target = {
    'H': 1.0e20,    # total H mass [kg]
    'C': 1.0e17,    # total C mass [kg]
    'N': 1.0e17,    # total N mass [kg]
    'S': 1.0e16,    # total S mass [kg]
}
```

CALLIOPE solves a four-equation system that requires the sum of dissolved + atmospheric mass of each element to equal the target value. Oxygen is **not** a free constraint: it is set by $f_{\mathrm{O}_2}$ and the speciation of the H, C, N, S equilibria.

You can build `target` two ways:

- `get_target_from_params(ddict)`: derives it from `hydrogen_earth_oceans`, `CH_ratio`, `nitrogen_ppmw`, `sulfur_ppmw`. Used when you specify the planet's bulk composition as an inventory.
- `get_target_from_pressures(ddict)`: derives it from `<species>_initial_bar` by computing the implied dissolved + atmospheric masses at $t=0$. Used when you want to specify the planet by its initial atmospheric composition.

The PROTEUS-side `volatile_mode` config knob switches between these two paths.

## Solver tuning knobs

Optional keyword arguments to `equilibrium_atmosphere()`:

| Argument | Default | Meaning |
|---|---|---|
| `rtol` | `1e-5` | Relative mass-balance tolerance. Solution is accepted if max residual $< \max(\text{target}) \cdot r_\text{tol} + a_\text{tol} + 10\,\mathrm{kg}$. |
| `atol` | `1e10` | Absolute mass-balance tolerance, kg. The default is loose; PROTEUS uses `mass_thresh` ($10^{16}$ kg) instead. |
| `xtol` | `1e-8` | Relative tolerance on partial-pressure updates passed to `fsolve`. |
| `nsolve` | `1500` | Maximum iterations per `fsolve` / `trust-constr` invocation. |
| `nguess` | `7500` | Maximum number of Monte-Carlo restarts before giving up. |
| `p_guess` | `None` | Optional dict `{'H2O': p, 'CO2': p, 'N2': p, 'S2': p}` used as the first initial guess. PROTEUS warm-starts from the previous-iteration partial pressures via this argument. |
| `opt_solver` | `True` | If `True`, alternate between `fsolve` and `trust-constr` on failure. If `False` (PROTEUS default), use `fsolve` only. |
| `hide_warnings` | `True` | Suppress the `RuntimeWarning`/`UserWarning` floods from `scipy` when the solver evaluates physically nonsensical guesses. |
| `print_result` | `True` | If `True`, log the final partial pressures at INFO level. PROTEUS sets this to `False`. |

## Next step

Now you know what to put in. See [Usage](usage.md) for the calling patterns, or jump straight to the [First run tutorial](../Tutorials/firstrun.md) for an end-to-end walkthrough.
