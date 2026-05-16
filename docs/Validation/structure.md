# Validation: `src/calliope/structure.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.structure` against a published source.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_structure.py::test_calculate_mantle_mass_recovers_wang_2018_earth_core_fraction` | Wang, Lineweaver & Ireland (2018)[^cite-wang2018], arxiv:1708.08718: Earth core mass fraction 32.5 +/- 0.3 wt% | [arxiv:1708.08718](https://arxiv.org/abs/1708.08718) | Pins the Earth-like mantle mass against the published core mass fraction and verifies the result lies in the 1e24 to 1e25 kg envelope expected for an Earth-mass planet. |

## Re-derivation note

`structure.calculate_mantle_mass` computes mantle mass by subtracting the
core mass from the total planetary mass. The core mass is the volume of a
sphere of radius `core_frac * R_planet` multiplied by a core density derived
from Earth: `core_rho = 3 * earth_fm * M_earth / (4 pi * (earth_fr * R_earth)^3)`
with `earth_fm = 0.325` and `earth_fr = 0.55` from Wang et al. (2018)[^cite-wang2018].

For Earth-like input (`R = R_earth`, `M = M_earth`, `core_frac = 0.55`), the
construction is degenerate: the core density times the Earth-radius core
volume reproduces the cited 0.325 mass fraction, so the mantle equals
`(1 - 0.325) * M_earth = 0.675 * M_earth ~ 4.03e24 kg`.

Scale: a regression that swaps the `r^3` core-volume factor for `r^2` would
land at ~3e18 kg, six orders of magnitude below the correct value. The
test's scale guard `1e24 < mantle < 1e25` brackets the correct order and
fails on any factor-of-10 unit slip.

## Anchor type

Published benchmark + analytical limit. The Wang+2018 cite is the
published-benchmark anchor; the conservation closure `M_mantle + M_core = M_planet`
is asserted separately in `test_calculate_mantle_mass_closure_holds_for_earth_like`
as the analytical-limit second-line check.

## Cross-references

- `src/calliope/structure.py` line 50: cites arxiv:1708.08718 for Earth core mass fraction 0.325.
- PROTEUS `src/proteus/orbit/satellite.py` uses the same Earth core fraction in the satellite angular-momentum decomposition; both codes are pinned against the same source.

## References

[^cite-wang2018]: H. S. Wang, C. H. Lineweaver, T. R. Ireland, *[The elemental abundances (with uncertainties) of the most Earth-like planet](https://doi.org/10.1016/j.icarus.2017.08.024)*, Icarus, 299, 460-474, 2018. [SciX](https://scixplorer.org/abs/2018Icar..299..460W/abstract).
