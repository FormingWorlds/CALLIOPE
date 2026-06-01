# Solve with an oxygen budget (authoritative-O mode)

This page shows how to call `equilibrium_atmosphere_authoritative_O` directly from Python, the inverse of the buffered-mode workflow on the [Usage page](usage.md). Use this entry point when you have a total oxygen mass to enforce as a budget rather than an oxygen fugacity to apply as a buffer.

For the science behind the mode and the augmented mass-balance system, see [Authoritative-oxygen mode](../Explanations/authoritative_oxygen.md).

## Recipe 1 - solve from a 5-element inventory

The minimum call signature: supply a `target_d` with five element keys (H, C, N, S, O in kilograms), the same `ddict` you would build for the buffered mode, and a sensible `fO2_hint`.

```python
from calliope.solve import equilibrium_atmosphere_authoritative_O
from calliope.constants import volatile_species

ddict = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
    'Phi_global': 1.0,
    'T_magma': 2500.0,
    # NOTE: fO2_shift_IW is ignored by this entry point; the solver
    # determines it. Provide it only if you also intend to call the
    # buffered entry point against the same ddict.
}
for sp in volatile_species:
    ddict[f'{sp}_included'] = 1
    ddict[f'{sp}_initial_bar'] = 0.0

target_d = {
    'H': 1.55e20,    # ~1 Earth ocean of H
    'C': 1.55e19,    # C/H ~ 0.1
    'N': 8.06e18,    # 2 ppmw vs mantle
    'S': 8.06e20,    # 200 ppmw vs mantle
    'O': 1.4e21,     # the budgeted O
}

result = equilibrium_atmosphere_authoritative_O(
    target_d, ddict, fO2_hint=4.0, print_result=True,
)

print(f"Derived fO2_shift : {result['fO2_shift_derived']:+.3f} dex")
print(f"O residual        : {result['O_res']:.2e} kg")
print(f"Surface pressure  : {result['P_surf']:.2f} bar")
print(f"H2O bar           : {result['H2O_bar']:.2f}")
```

The returned `result` dict carries every field `equilibrium_atmosphere` returns (per-species `_bar`, `_vmr`, `_kg_atm`, `_kg_liquid`, `_kg_solid`, `_kg_total`; per-element `_kg_atm`, `_res`; total `P_surf`, `M_atm`, `atm_kg_per_mol`) plus two new keys: `fO2_shift_derived` for the converged $\Delta\mathrm{IW}$ and `O_res` for the fifth mass-balance residual.

## Recipe 2 - round-trip from buffered mode

The cleanest sanity check on a new dataset is to confirm that the buffered and authoritative-O modes are dual: solve once in buffered mode, take the resulting O total as the new O target, and re-solve in authoritative-O mode. The recovered $\Delta\mathrm{IW}$ should match the original to within solver tolerance.

```python
from calliope.solve import (
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
    get_target_from_params,
)

target_HCNS = get_target_from_params(ddict)

# Step 1: buffered mode at fO2_shift_IW = +4
ddict['fO2_shift_IW'] = 4.0
result_buffered = equilibrium_atmosphere(target_HCNS, ddict, print_result=False)

# Step 2: take the implied O total and re-solve in authoritative-O mode
target_HCNSO = dict(target_HCNS, O=result_buffered['O_kg_total'])
result_auth = equilibrium_atmosphere_authoritative_O(
    target_HCNSO, ddict, fO2_hint=4.0, print_result=False,
)

assert abs(result_auth['fO2_shift_derived'] - 4.0) < 0.01
```

This pattern is the basis of the `test_authoritative_O_round_trip` regression test. Failing it at non-extreme $\Delta\mathrm{IW}$ typically indicates a bug in one of the solubility laws or the speciation tree.

## Recipe 3 - warm-start the solver via `p_guess`

When you have a converged solution from a previous call (the typical pattern inside an iterative loop), pass it in via `p_guess` to skip the cold-start Monte-Carlo draw:

```python
p_guess = {
    'H2O': result['H2O_bar'],
    'CO2': result['CO2_bar'],
    'N2':  result['N2_bar'],
    'S2':  result['S2_bar'],
    'fO2_shift_IW': result['fO2_shift_derived'],   # optional, defaults to fO2_hint
}

result_new = equilibrium_atmosphere_authoritative_O(
    target_d_new, ddict_new, p_guess=p_guess, print_result=False,
)
```

If you omit the `'fO2_shift_IW'` key in `p_guess`, the solver falls back to `fO2_hint` for the fifth unknown. Supplying it from the previous call's `fO2_shift_derived` is the right choice when the inputs change slowly between calls (small time step, small budget perturbation): convergence then needs ${\sim}1$ to ${\sim}3$ `fsolve` iterations.

## Recipe 4 - reproducibility for regression tests

The Monte-Carlo restart draws nondeterministic guesses by default. Pass `random_seed=<int>` to make the solver reproducible across runs to solver tolerance:

```python
import math

result_a = equilibrium_atmosphere_authoritative_O(
    target_d, ddict, fO2_hint=4.0, random_seed=42, print_result=False,
)
result_b = equilibrium_atmosphere_authoritative_O(
    target_d, ddict, fO2_hint=4.0, random_seed=42, print_result=False,
)
assert math.isclose(
    result_a['fO2_shift_derived'], result_b['fO2_shift_derived'], abs_tol=1e-8,
)
```

The seed fixes the Monte-Carlo guess sequence; the inner `fsolve` and `trust-constr` calls are deterministic in their own right, so two seeded runs converge to the same root within `xtol`. For regression tests that compare against checked-in golden values, this is the right level of reproducibility. For production runs, leave `random_seed=None` (the default) so restart draws use the global `np.random` state.

## Pitfalls

- **`KeyError: target_d is missing required element keys: ['O']`** - you forgot to include `'O'` in `target_d`. The four-element dict that `get_target_from_params` returns is the buffered-mode input; in authoritative-O mode you have to add the O budget yourself.
- **`ValueError: fO2_hint=...` outside `[-12, +12]`** - the hint must be a finite real number in the solver-bounds window. Common values lie in $[-4, +6]$; values outside $[-12, +12]$ are rejected up front because they almost always indicate a unit confusion (e.g. passing absolute $\log_{10} f_{\mathrm{O}_2}$ instead of the IW-buffer offset).
- **`RuntimeError: Could not find solution ... (max attempts: N)`** - the Monte-Carlo restart loop exhausted. The integer `N` is whatever `nguess` was set to ($7500$ if the entry point is called directly with defaults; $1000$ if called through the PROTEUS wrapper). The failure message includes the final-attempt partial pressures and $\Delta\mathrm{IW}$. The diagnostic question is whether the target O budget is physically reachable at the supplied $(H, C, N, S, T_\mathrm{magma}, \Phi_\mathrm{global})$. Bumping `nguess` rarely helps; investigate whether the upstream inventory is consistent first.
- **Extrapolated solubility laws** - the entry point does not flag extrapolation. The Dasgupta N$_2$ and Gaillard S$_2$ laws carry $\ln f_{\mathrm{O}_2}$ terms with finite calibration footprints; outside those footprints the dissolved-N and dissolved-S branches of the O balance are uncertain. Consult the [validity envelope](../Explanations/solubility.md#validity-envelope) before interpreting results at unusual $T_\mathrm{magma}$ or $\Delta\mathrm{IW}$.
- **`p_guess` missing one of `'H2O'`, `'CO2'`, `'N2'`, `'S2'`** - the four primary keys are required. A missing `'fO2_shift_IW'` is fine and just falls back to `fO2_hint`. Non-finite values in any key raise `ValueError`.

## Next step

For the role of this entry point inside a PROTEUS-coupled simulation, see [Coupling to PROTEUS](proteus_coupling.md). For the mathematical setup, see [Authoritative-oxygen mode](../Explanations/authoritative_oxygen.md).
