# Authoritative-oxygen mode

CALLIOPE offers two equilibrium-chemistry entry points that share the same physics but differ in which quantity is the input and which is the unknown:

- `equilibrium_atmosphere` (the [buffered mode](mass_balance.md)) takes the oxygen fugacity $f_{\mathrm{O}_2}$ as an input via `ddict['fO2_shift_IW']` and solves for the four primary partial pressures $(p_\mathrm{H_2O}, p_\mathrm{CO_2}, p_\mathrm{N_2}, p_\mathrm{S_2})$ against the H, C, N, S elemental budgets. Oxygen mass is *implicit*: it is whatever the chemistry requires at the prescribed buffer, and it can change freely between calls.
- `equilibrium_atmosphere_authoritative_O` takes the total oxygen mass as a fifth budget alongside H, C, N, S, and solves for the four pressures *plus* $\Delta\mathrm{IW}$ as a fifth unknown.

This page documents the augmented mass-balance system, the additional solver controls, and the contract the two modes share.

## When the two modes agree

The buffered mode and the authoritative-O mode are dual formulations of the same underlying chemistry. For any pair of $\Delta\mathrm{IW}$ and target H, C, N, S budgets, the buffered mode produces an equilibrium atmosphere whose total O mass can be read off the result dict as `O_kg_total`. Feeding $(H, C, N, S, \texttt{O\_kg\_total})$ to the authoritative-O mode reproduces the same partial pressures and recovers the same $\Delta\mathrm{IW}$ to within the solver tolerance. The regression test `tests/test_authoritative_O.py::TestRoundTrip::test_round_trip_recovers_fO2_within_tolerance` enforces this round-trip and is the operational definition of "same physics, different unknowns".

The two modes differ only in their degrees-of-freedom accounting:

| Mode | Inputs | Unknowns | Equations |
|---|---|---|---|
| `equilibrium_atmosphere` | $(H, C, N, S)$ targets, $\Delta\mathrm{IW}$ | $(p_\mathrm{H_2O}, p_\mathrm{CO_2}, p_\mathrm{N_2}, p_\mathrm{S_2})$ | 4 mass-balance |
| `equilibrium_atmosphere_authoritative_O` | $(H, C, N, S, O)$ targets | $(p_\mathrm{H_2O}, p_\mathrm{CO_2}, p_\mathrm{N_2}, p_\mathrm{S_2}, \Delta\mathrm{IW})$ | 5 mass-balance |

## The augmented mass-balance system

Adding an O budget closes the system with one extra constraint. For each element $e \in \{\mathrm{H}, \mathrm{C}, \mathrm{N}, \mathrm{S}, \mathrm{O}\}$:

$$
r_e(\mathbf{x}) = m_e^\mathrm{atm}(\mathbf{x}) + m_e^\mathrm{melt}(\mathbf{x}) - m_e^\mathrm{target} = 0,
$$

where the unknown is the five-vector $\mathbf{x} = (p_\mathrm{H_2O}, p_\mathrm{CO_2}, p_\mathrm{N_2}, p_\mathrm{S_2}, \Delta\mathrm{IW})$. The per-species column masses $m_v^\mathrm{atm}$ and dissolved masses $m_v^\mathrm{melt}$ are computed from the [same physics functions](mass_balance.md#atmospheric-column-mass) as the buffered mode; $\Delta\mathrm{IW}$ is read from $\mathbf{x}[4]$, and the residual vector has five entries (one per element in H, C, N, S, O).

The private versions `_atmosphere_mass` and `_dissolved_mass` accept `fO2_shift` as a positional argument, which lets both solver modes consume identical physics through one set of functions. The five-residual function is `solve.func_authoritative_O`:

```python
def func_authoritative_O(x_arr, ddict, mass_target_d):
    pin_dict = {'H2O': x_arr[0], 'CO2': x_arr[1],
                'N2':  x_arr[2], 'S2':  x_arr[3]}
    fO2_shift = x_arr[4]
    mass_atm_d = _atmosphere_mass(pin_dict, fO2_shift, ddict)
    mass_int_d = _dissolved_mass(pin_dict, fO2_shift, ddict)
    return [mass_atm_d[v] + mass_int_d[v] - mass_target_d[v]
            for v in ('H', 'C', 'N', 'S', 'O')]
```

## How $\Delta\mathrm{IW}$ enters as an unknown

$\Delta\mathrm{IW}$ enters the chemistry through the four channels listed on the [oxygen fugacity page](oxygen_fugacity.md#how-deltamathrmiw-enters-the-chemistry): the free $\mathrm{O_2}$ partial pressure, the modified equilibrium constants $G_\mathrm{eq}$, the Gaillard S$_2$ solubility, and the Dasgupta N$_2$ solubility. In the buffered mode $\Delta\mathrm{IW}$ is an input held fixed during the solve; in the authoritative-O mode it is a variable the solver adjusts in concert with the four pressures to satisfy the five-residual system.

The closure mechanism is intuitive: increasing $\Delta\mathrm{IW}$ (more oxidising) drives more H$_2$O at fixed H, more CO$_2$ at fixed C, and more SO$_2$ and dissolved S at fixed S. The Dasgupta N$_2$ solubility carries an explicit $-1.6\,\Delta\mathrm{IW}$ term in its exponent, so dissolved N drops by roughly an order of magnitude per dex of oxidation, while N$_2$ in the atmosphere takes up the slack. The aggregated atmospheric + dissolved O mass therefore monotonically increases with $\Delta\mathrm{IW}$ over the physically relevant range, which gives the 5-residual system a unique root for any feasible O target. The `tests/test_authoritative_O_monotonicity.py::TestMonotonicity::test_O_kg_total_strictly_increasing_with_fO2` property test pins this monotonicity over $\Delta\mathrm{IW} \in [-4, +6]$ at two $T_\mathrm{magma}$ values, with a sibling check at $\Phi_\mathrm{global} = 0.3$ over $\Delta\mathrm{IW} \in [-2, +4]$.

## Solver: per-element residual gate

`equilibrium_atmosphere_authoritative_O` uses the same hybrid `fsolve` + `trust-constr` outer loop as the buffered mode (see [Mass balance & solver](mass_balance.md#solver-hybrid-powell-trust-region-with-monte-carlo-restart)), with three differences in the acceptance test and the initial-guess generation.

### Per-element tolerance

The buffered mode accepts a solution if the scalar bound

$$
\max_e |r_e| < r_\mathrm{tol} \cdot \max_e m_e^\mathrm{target} + a_\mathrm{tol} + 10\,\mathrm{kg}
$$

is satisfied. The buffered mode never sees O in its residual vector, so the four H/C/N/S targets dominate $\max_e m_e^\mathrm{target}$. In the authoritative-O mode O is now in the residual vector, and at planet scale $m_\mathrm{O}^\mathrm{target} \sim 10^{22}\,\mathrm{kg}$ can be five orders of magnitude larger than $m_\mathrm{N}^\mathrm{target} \sim 10^{17}\,\mathrm{kg}$. A tolerance scaled to the largest target would silently admit N residuals of order $10^{17}\,\mathrm{kg}$, which is the whole N budget.

The authoritative-O mode replaces the scalar bound with a per-element gate:

$$
|r_e| < \max\left(r_\mathrm{tol} \cdot m_e^\mathrm{target}, \tfrac{a_\mathrm{tol}}{5}, \mathrm{TRUNC\_MASS}\right) \quad \forall e \in \{\mathrm{H}, \mathrm{C}, \mathrm{N}, \mathrm{S}, \mathrm{O}\}.
$$

The $a_\mathrm{tol}$ budget is split evenly across the five elements so the per-element floor stays in the kilogram range, and the `TRUNC_MASS` constant ($10\,\mathrm{kg}$, from `solve.py`) handles benign sub-kilogram noise that should not gate convergence.

### The fO2 hint and the restart redraw

The five-dimensional problem has a larger initial-guess space than the four-dimensional buffered problem, so the user supplies `fO2_hint` (default $+4.0$) as a starting value for $\Delta\mathrm{IW}$. On a cold start (`p_guess=None`), the four pressures are drawn log-uniformly over $[10^{-12}, 10^5]$ bar (same as the buffered mode) and the fifth unknown is initialised to `fO2_hint`. On every Monte-Carlo restart after the initial attempt, the four pressures are redrawn and $\Delta\mathrm{IW}$ is redrawn from $\mathrm{Uniform}(-6, +8)$. This window covers the physically relevant mantle redox states (reducing through Mercury-like at $-5$ to highly oxidised at $+5$) so a restart can pick a fundamentally different basin from the initial hint.

### Bounds

The trust-region solver uses bounds $[0, 10^7]$ bar on each pressure (same as the buffered mode) and $[-12, +12]$ on $\Delta\mathrm{IW}$. The pressure bounds are physical; the fO2 bounds are wider than the $[-6, +8]$ redraw window so the solver has slack to escape a poor cold start without immediately hitting the constraint, but tight enough to reject runaway trajectories into thermodynamically meaningless regions.

### Reproducibility

A `random_seed` argument (default `None`) seeds a `np.random.default_rng` for the restart draws. Passing an integer makes the solver deterministic across calls at fixed inputs, which is required for regression tests and for diffing two runs. The default `None` preserves the non-deterministic global-state behaviour shared with the buffered mode.

## Buffer convention

Like the buffered mode, the authoritative-O mode references $f_{\mathrm{O}_2}$ to the iron-wüstite buffer of [O'Neill & Eggins (2002)](https://ui.adsabs.harvard.edu/abs/2002ChGeo.186..151O) by default, with the [Fischer et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011E%26PSL.304..496F) alternative selectable through `OxygenFugacity()` instantiation. The returned `fO2_shift_derived` is the $\Delta\mathrm{IW}$ relative to whichever buffer was chosen.

!!! note "Cross-backend buffer divergence"
    PROTEUS supports a second outgassing backend, [atmodeller](https://atmodeller.readthedocs.io/), whose authoritative-O implementation uses the Hirschmann combined IW buffer. The Hirschmann and O'Neill & Eggins parameterisations differ by ${\sim}0.95$ dex at $T = 3000$ K. PROTEUS records both backends' derived offsets under the helpfile column `fO2_shift_IW_derived`, and the discrepancy is documented in the column's schema comment. The two backends agree on the underlying physics (same chemistry of FeO-O$_2$ equilibrium); they disagree on the numerical parameterisation of the buffer curve. Choose one backend per run and stay with it for any cross-time-step comparison.

## When to use this mode

The authoritative-O mode is the right tool when atmospheric + dissolved O is a budgeted quantity rather than a derived one. Three typical use cases:

1. **Whole-planet O accounting.** When PROTEUS treats O as one of the five tracked elements (so escape debits it, the mantle initial inventory budgets it, and the structure solver weights it into M_planet), the chemistry must invert against the running O total rather than buffer to a user-prescribed $\Delta\mathrm{IW}$. The PROTEUS `planet.fO2_source = "from_O_budget"` config switch flips outgassing into this mode at runtime.
2. **Inferring $\Delta\mathrm{IW}$ from observed composition.** Given a measured or assumed elemental inventory including O, the authoritative-O mode returns the $\Delta\mathrm{IW}$ consistent with that inventory at the chosen $T_\mathrm{magma}$ and $\Phi_\mathrm{global}$. This complements the buffered-mode forward calculation when the modelling question is "what redox state does this composition imply?".
3. **Sensitivity studies that vary the O reservoir.** When sweeping over a mantle FeO inventory grid, the authoritative-O mode lets you specify the O budget directly rather than having to translate each grid point through a buffer offset first.

For the routine "set $\Delta\mathrm{IW}$, get a self-consistent atmosphere" workflow, the buffered mode is faster and has fewer unknowns. The authoritative-O mode adds one residual equation and one unknown; convergence is empirically a few Monte-Carlo restarts with a good `fO2_hint` and PROTEUS warm-starting, compared to typically a single attempt for the buffered mode under the same conditions.

## Limitations

- **Solubility-law calibration.** The Dasgupta N$_2$ and Gaillard S$_2$ solubility laws have explicit $\ln f_{\mathrm{O}_2}$ terms. Outside their calibration footprints the solver still produces a solution, but the implied O partitioning between dissolved S/N and atmospheric SO$_2$/N$_2$ extrapolates the underlying experimental data. The [solubility-laws validity envelope](solubility.md#validity-envelope) gives the temperature, pressure, and fO$_2$ ranges over which each law was fit. The authoritative-O entry point itself does not flag extrapolation; the caller is expected to consult the envelope.
- **Solid-FeO buffering not modelled.** The chemistry assumes a single surface $\Delta\mathrm{IW}$ controls O speciation everywhere in the atmosphere. As $\Phi_\mathrm{global} \to 0$ the melt cannot buffer O against the gas-phase composition, but the entry point still produces a solution: at zero melt fraction the dissolved-mass channel is identically zero and the O constraint reduces to a pure atmospheric balance.
- **Monotonicity not guaranteed pathologically.** For target O budgets so small that O is dominated by dissolved species at strongly reducing $\Delta\mathrm{IW}$, the residual surface can develop a second root in the $\Delta\mathrm{IW}$ direction. The Monte-Carlo restart with fO2 redraw is designed to escape such basins, but for adversarial inputs (sub-trace O at very low $T_\mathrm{magma}$) it may need many restarts or fail. Failure raises `RuntimeError` with the final attempt's pressures and $\Delta\mathrm{IW}$ for diagnosis; the right response is usually to inspect whether the upstream O budget is physical, not to bump `nguess`.
- **No O kinetics.** Like the buffered mode, the authoritative-O mode is a snapshot equilibrium computation. Time-resolved O evolution is the responsibility of the calling code (escape, ingassing, diffusion).

## See also

- [How-to: authoritative oxygen mode](../How-to/authoritative_oxygen.md): the recipe for calling the entry point.
- [Mass balance & solver](mass_balance.md): the four-residual buffered mode.
- [Oxygen fugacity](oxygen_fugacity.md): the IW-buffer parameterisations and the four channels through which $\Delta\mathrm{IW}$ enters the chemistry.
- [Coupling to PROTEUS (theory)](proteus_coupling.md): how the PROTEUS wrapper selects between the two modes.
- [API reference for `calliope.solve`](../Reference/api/calliope.solve.md).
