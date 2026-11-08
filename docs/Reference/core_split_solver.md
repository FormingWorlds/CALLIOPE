# Core split solver

This page lists the steps with which `calliope.core_split.split_core_mantle` finds the metal composition of a split. The [Core and mantle split](../Explanations/core_split.md) page explains why the split must iterate and what the steps do in plain words; the [sweep of random inputs](../Validation/metal.md#core-and-mantle-split-over-random-inputs) shows how often each outcome occurs.

## Terms

- **Pass.** One update of the metal composition: compute $D$ for the present metal, split the bulk, and compute the new metal from the split.
- **$w$ and $g$.** In a pass, $w$ is the present metal composition and $g$ is the new composition that the pass computes from $w$. A root is a metal with $g = w$.

## Steps

1. **Plain passes.** Each pass sets the metal to $g$. The split stops when no metal mass fraction and no melt mass fraction of H and C changes by more than $10^{-12}$ of its value plus $10^{-15}$. The melt H and C are part of the test because they enter the moles of melt per kg of the mole basis.
2. **When the plain passes fail.** Where the metal takes so much S that $D_\mathrm{S}$ falls fast with the metal S ($X_\mathrm{FeO}$ falls with the Fe of the metal), the plain passes swing about the root. After 50 plain passes that do not settle, or after a pass that leaves no Fe, the split starts again with relaxed passes.
3. **Relaxed passes.** Each element moves $1/(1 + |s|)$ of its step, with $s$ the slope of $g$ against $w$ from its last 2 passes. A pass adds at most half of the Fe that is left.
4. **Common passes.** If 1000 relaxed passes do not settle, or a pass leaves Fe below $10^{-8}$ by mass, the relaxed passes start again with the smallest of these factors for every element. This keeps coupled elements in step.
5. **A root that floating-point numbers cannot resolve.** On a map so steep that the floating-point numbers of the metal cannot resolve the root, the relaxed passes accept an element when all of these hold: its $g$ falls with $w$ ($s < 0$); its move at its own factor is within 2 ulp (units in the last place of the number); its step is below $10^{-9}$ of $g$; and $g - w$ changes sign within 4 ulp on each side of the evaluated metal. The split is then returned at that metal for these bracketed elements, while an element that met the test of step 1 keeps $g$.
6. **The iron-rich root rule.** A relaxed split can settle or return on a root with little Fe while a root with more Fe exists. An S-C-H-N split with a silicate of $1.8 \times 10^{-5}$ of the metal mass has roots at Fe $5.0 \times 10^{-2}$ and $6.2 \times 10^{-7}$ by mass, both past FeS; round-off decides which one the relaxed passes reach. So when the metal of a relaxed return holds more S than Fe in moles, the split runs damped passes from the start metal: every element moves 0.02 of its step (the melt H and C undamped), for at most 20000 passes. The split returns their root when it holds more Fe and is another metal: an element is apart by more than $10^{-9}$ of its value plus $10^{-15}$, the grid tolerance of the relaxed passes. It then sets `root_rule` to 1. This rule prefers the root with more Fe; it does not prove that the split has one root.
7. **Errors.** When the common passes leave Fe below $10^{-8}$ by mass, the iron runs out and the split raises `ValueError`. When the common passes do not settle, the split raises that `ValueError` if the first relaxed passes ran out of iron, and `RuntimeError` otherwise.
8. **Residual.** Each returned split reports in `flags[e]['D_residual']` how well the returned metal meets the law of each element (see [What the output flags mean](../Explanations/core_split.md#what-the-output-flags-mean)).

## Errors

The errors of the split, with the inputs that cause them, are listed in the [partition file reference](partition_file.md#named-laws-and-the-core-split).
