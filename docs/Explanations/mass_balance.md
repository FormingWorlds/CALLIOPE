# Mass balance & solver

CALLIOPE's prognostic equations are four nonlinear elemental mass-conservation constraints, one per solved element (H, C, N, S). This page documents the residual function, the solver strategy, the mass-from-pressure relations, and the convergence criterion of this "buffered" mode, where the oxygen fugacity is an input and oxygen mass is derived. CALLIOPE also offers an [authoritative-oxygen mode](authoritative_oxygen.md) where the system is closed by adding O as a fifth budget and treating $\Delta\mathrm{IW}$ as an unknown; the two modes share all the physics functions and differ only in their unknown set.

## The conservation system

For each element $e \in \{\mathrm{H}, \mathrm{C}, \mathrm{N}, \mathrm{S}\}$ the conservation equation reads

$$
m_e^\mathrm{atm}(\mathbf{p}) + m_e^\mathrm{melt}(\mathbf{p}) - m_e^\mathrm{target} = 0,
$$

where $\mathbf{p} = (p_\mathrm{H_2O}, p_\mathrm{CO_2}, p_\mathrm{N_2}, p_\mathrm{S_2})$ is the four-vector of primary partial pressures in bar. The seven secondary partial pressures (and $p_\mathrm{O_2}$) are not independent variables: they are algebraic functions of $\mathbf{p}$ via the [equilibrium chemistry](equilibrium_chemistry.md) speciation tree. Oxygen is *not* conserved; it is set by the $f_{\mathrm{O}_2}$ buffer and is allowed to leak in or out of the system as the speciation requires.

The four-vector residual function is `solve.func`:

```python
def func(pin_arr, ddict, mass_target_d):
    pin_dict = {'H2O': pin_arr[0], 'CO2': pin_arr[1],
                'N2':  pin_arr[2], 'S2':  pin_arr[3]}
    mass_atm_d = atmosphere_mass(pin_dict, ddict)
    mass_int_d = dissolved_mass(pin_dict, ddict)
    return [mass_atm_d[v] + mass_int_d[v] - mass_target_d[v]
            for v in ('H', 'C', 'N', 'S')]
```

## Atmospheric column mass

The relation between a species' surface partial pressure and its column mass follows directly from hydrostatic equilibrium under the assumption of a well-mixed atmosphere. Bower et al. (2019) [^cite-bower2019] Equation (2) writes it as

$$
m_v^\mathrm{atm} = 4\pi R_p^2 \cdot \frac{\mu_v}{\bar\mu} \cdot \frac{p_v}{g},
$$

where $R_p$ is the planetary surface radius, $\mu_v$ is the species molar mass, $\bar\mu$ is the atmospheric mean molar mass weighted by partial pressure (`atmosphere_mean_molar_mass()` in `solve.py`), and $g$ is the surface gravity. CALLIOPE stores partial pressures in bar, so the implementation carries the conversion factor $1.0 \times 10^5$ Pa/bar:

```python
mass_atm_d[key] = value * 1.0e5 / ddict['gravity']
mass_atm_d[key] *= 4.0 * np.pi * ddict['radius'] ** 2.0
mass_atm_d[key] *= molar_mass[key] / mu_atm
```

The `mu_v / mu_atm` ratio is the part Bower et al. (2019) [^cite-bower2019] §4.1.1 emphasises was missing from the pre-2019 mass-balance formulations of Elkins-Tanton (2008) [^cite-elkinstanton2008], Lebrun et al. (2013) [^cite-lebrun2013], Salvador et al. (2017) [^cite-salvador2017], and Nikolaou et al. (2019) [^cite-nikolaou2019]. Without it, multi-species atmospheres receive an unphysical bias in the inferred reservoir partitioning.

After computing per-species column masses, `atmosphere_mass()` aggregates them into per-element atomic masses by stoichiometric atom-counting:

$$
m_\mathrm{H}^\mathrm{atm} = 2\,\frac{m_\mathrm{H_2O}^\mathrm{atm}}{\mu_\mathrm{H_2O}} + 2\,\frac{m_\mathrm{H_2}^\mathrm{atm}}{\mu_\mathrm{H_2}} + 4\,\frac{m_\mathrm{CH_4}^\mathrm{atm}}{\mu_\mathrm{CH_4}} + 2\,\frac{m_\mathrm{H_2S}^\mathrm{atm}}{\mu_\mathrm{H_2S}} + 3\,\frac{m_\mathrm{NH_3}^\mathrm{atm}}{\mu_\mathrm{NH_3}}
$$

times $\mu_\mathrm{H}$, and analogously for C, N, O, and S. The factor-3 in NH$_3$ and factor-4 in CH$_4$ are exactly the cases the `tests/test_stoichiometry.py::TestAtmosphericStoichiometry` tests pin down (e.g. `test_NH3_contributes_1_N_not_3` verifies that NH$_3$ contributes 1 N atom not 3).

## Dissolved mass

For each solubility-supported species, `dissolved_mass()` evaluates the per-species ppmw concentration via the chosen [solubility law](solubility.md), then converts to absolute mass via the prefactor

$$
m_i^\mathrm{melt} = 10^{-6}\, M_\mathrm{mantle} \cdot \Phi_\mathrm{global} \cdot X_i^\mathrm{melt}\,[\text{ppmw}],
$$

where $M_\mathrm{mantle}$ is the (molten + solid) silicate mantle mass and $\Phi_\mathrm{global}$ is the global melt fraction. Setting $\Phi_\mathrm{global} = 0$ disables solubility entirely; setting $\Phi_\mathrm{global} = 1$ (fully molten) gives the maximum dissolved-mass contribution.

Like for atmospheric mass, the per-species dissolved masses are aggregated into per-element atomic masses. Note the asymmetry with the atmospheric path: CALLIOPE only includes a subset of species in the dissolved-mass tally (H$_2$O, CO$_2$, CO, CH$_4$, N$_2$, S$_2$); the remaining species (H$_2$, NH$_3$, SO$_2$, H$_2$S, O$_2$) are assumed to have negligible solubility, consistent with Bower et al. (2022) [^cite-bower2022] §2.2.3.

## Solver: hybrid Powell + trust-region with Monte-Carlo restart

`equilibrium_atmosphere()` wraps the residual evaluation in a robust outer loop:

1. **Initial guess**: either the user-supplied `p_guess` (if non-`None`), or a Monte-Carlo draw via `get_initial_pressures()` that samples each primary partial pressure log-uniformly in $[10^{-12}, 10^5]$ bar.
2. **Inner solve**: alternates between
    - `scipy.optimize.fsolve` (default), which uses the MINPACK Powell hybrid algorithm and is the fastest path when the initial guess is in the right basin;
    - `scipy.optimize.minimize` with `method='trust-constr'` and bounds $[0, 10^7]$ bar, which is more robust to bad guesses but slower;
    - the alternation is gated by `opt_solver=True`. PROTEUS sets `opt_solver=False` and stays on `fsolve`.
3. **Acceptance test**: even if the solver flags `success`, CALLIOPE re-evaluates the residual and accepts only if $\max_e |r_e| < r_\text{tol} \cdot \max_e m_e^\mathrm{target} + a_\text{tol} + 10$ kg. This catches the "false convergence" failure mode where `fsolve` lands on a stationary point of the residual norm rather than a true zero.
4. **Restart on rejection**: if either the inner solver or the acceptance test fails, draw a new Monte-Carlo guess and try again, up to `nguess` times.
5. **Hard failure**: if `nguess` restarts all fail, raise `RuntimeError`. The PROTEUS wrapper catches this and writes a status code 27 ("outgassing failure") into the run's status file.

## Why the Monte-Carlo restart is necessary

The residual function $\mathbf{r}(\mathbf{p})$ has multiple physically valid roots when the elemental inventory is small enough that some species can be driven to numerical zero, and multiple physically *invalid* roots that come from the implicit positivity constraints not being explicitly enforced inside `fsolve`. The Powell hybrid algorithm has no notion of physical bounds and can happily walk off into negative-pressure territory if the initial guess is too far from the basin.

The Monte-Carlo restart cures both pathologies: by sweeping log-uniformly over 17 orders of magnitude, the solver eventually lands in the basin of the *physically* correct root from a starting point close enough that `fsolve` converges before any negative excursion occurs. With a good warm start (PROTEUS pattern), the first attempt succeeds in 99%+ of cases; without one, $\sim 10$-50 restarts are typical.

!!! note "On the iteration caps"
    The 1500-iteration `nsolve` cap and the 7500-restart `nguess` cap in the library defaults are both deliberately loose: they bound the wall-time at $\sim$10 seconds per call (a conservative ceiling) without cutting off pathological cases that genuinely need more attempts. The PROTEUS wrapper tightens both ($n_\text{solve} = 3000$, $n_\text{guess} = 1000$) because the warm-start strategy makes the long tails irrelevant in normal operation. If you see the wrapper hit `nguess = 1000`, the upstream physics is broken, not the solver.

## Convergence diagnostics

The `result` dictionary returned by `equilibrium_atmosphere()` includes `H_res`, `C_res`, `N_res`, `S_res` (residuals in kg). Their absolute values bound how well the elemental conservation was satisfied; their relative values $|r_e| / m_e^\mathrm{target}$ should be $\lesssim 10^{-5}$ for an `rtol = 1e-5` solve.

`solve.equilibrium_atmosphere` also logs the chosen restart count `count` at DEBUG level. If `count > 100` consistently across iterations, the warm-start is failing or the basin is genuinely degenerate; switch on `print_result=True` and re-run the failing case in isolation to inspect the convergence trajectory.

## See also

- [Authoritative-oxygen mode](authoritative_oxygen.md) for the dual five-residual formulation where O is an input budget and $\Delta\mathrm{IW}$ is the additional unknown.
- [Equilibrium chemistry](equilibrium_chemistry.md) for the speciation tree that maps $\mathbf{p}$ to all eleven partial pressures.
- [Solubility laws](solubility.md) for the form of $X_i^\mathrm{melt}(p_i)$ in `dissolved_mass()`.
- [Coupling to PROTEUS (theory)](proteus_coupling.md) for how the wrapper builds `target` and `ddict` from `hf_row`.
- [API reference for `calliope.solve`](../Reference/api/calliope.solve.md).

 [^cite-bower2019]: D. J. Bower, D. Kitzmann, A. S. Wolf, P. Sanan, C. Dorn, A. V. Oza, *[Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations](https://doi.org/10.1051/0004-6361/201935710)*, Astronomy & Astrophysics, 631, A103, 2019. [SciX](https://scixplorer.org/abs/2019A%26A...631A.103B/abstract).
 [^cite-bower2022]: D. J. Bower, K. Hakim, P. A. Sossi, P. Sanan, *[Retention of water in terrestrial magma oceans and carbon-rich early atmospheres](https://doi.org/10.3847/PSJ/ac5fb1)*, The Planetary Science Journal, 3(4), 93, 2022. [SciX](https://scixplorer.org/abs/2022PSJ.....3...93B/abstract).
 [^cite-elkinstanton2008]: L. T. Elkins-Tanton, *[Linked magma ocean solidification and atmospheric growth for Earth and Mars](https://doi.org/10.1016/j.epsl.2008.03.062)*, Earth and Planetary Science Letters, 271, 181–191, 2008. [SciX](https://scixplorer.org/abs/2008E%26PSL.271..181E/abstract).
 [^cite-lebrun2013]: T. Lebrun, H. Massol, E. Chassefière, A. Davaille, E. Marcq, P. Sarda, F. Leblanc, G. Brandeis, *[Thermal evolution of an early magma ocean in interaction with the atmosphere](https://doi.org/10.1002/jgre.20068)*, Journal of Geophysical Research: Planets, 118, 1155–1176, 2013. [SciX](https://scixplorer.org/abs/2013JGRE..118.1155L/abstract).
 [^cite-salvador2017]: A. Salvador, H. Massol, A. Davaille, E. Marcq, P. Sarda, E. Chassefière, *The relative influence of H$_2$O and CO$_2$ on the primitive surface conditions and evolution of rocky planets*, Journal of Geophysical Research: Planets, 122, 1458–1486, 2017. [SciX](https://scixplorer.org/abs/2017JGRE..122.1458S/abstract).
 [^cite-nikolaou2019]: A. Nikolaou, N. Katyal, N. Tosi, M. Godolt, J. L. Grenfell, H. Rauer, *[What factors affect the duration and outgassing of the terrestrial magma ocean?](https://doi.org/10.3847/1538-4357/ab08ed)*, The Astrophysical Journal, 875, 11, 2019. [SciX](https://scixplorer.org/abs/2019ApJ...875...11N/abstract).
