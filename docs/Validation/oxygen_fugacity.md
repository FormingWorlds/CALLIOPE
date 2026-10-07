# Validation: `src/calliope/oxygen_fugacity.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.oxygen_fugacity` against published sources.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_fischer_value_at_2000K_matches_published_fit` | Fischer et al. (2011) [^cite-fischer2011], EPSL 304, 496, Eq. 2; cross-checked against O'Neill & Eggins (2002) [^cite-oneilleggins2002] | [doi:10.1016/j.epsl.2011.02.025](https://doi.org/10.1016/j.epsl.2011.02.025) | Pins the Fischer-vs-O'Neill cross-calibration offset (0.258 dex at T = 2000 K) as the independent anchor, with a secondary regression check on the coded Fischer fit `6.94059 - 28.1808e3 / T` and a wrong-buffer discrimination guard against O'Neill & Eggins (2002) at the same T. |
| `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_hirschmann_value_at_1473K_matches_published_figure` | Hirschmann (2021) [^cite-hirschmann2021], GCA 313, 74, Table 1 and Fig. 5 | [doi:10.1016/j.gca.2021.08.039](https://doi.org/10.1016/j.gca.2021.08.039) | Pins the 1 bar value at T = 1473 K (-11.957) against the IW value read from Fig. 5 (~-11.94, 0.05 dex tolerance), with wrong-buffer guards against Fischer (2011) and O'Neill & Eggins (2002). |
| `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_hirschmann_low_pressure_slope_matches_volume_change` | Analytical limit: d log10 fO2 / dP = 2 dV / (R T ln 10) | Hirschmann (2021) Appendix A for V(FeO) | The 0-1 GPa slope of the fit matches the room-temperature volume change of Fe + 1/2 O2 = FeO within 15% at 1200, 2000 and 3000 K. The fit runs 7-10% higher than this estimate. |

## Re-derivation note

`OxygenFugacity('fischer')` implements Fischer et al. (2011) Eq. 2 for
the iron-wüstite (IW) buffer:

```
log10(fO2)_IW(T) = 6.94059 - 28.1808e3 / T
```

At T = 2000 K this evaluates to:

```
log10(fO2) = 6.94059 - 14.0904 = -7.14981
```

Cross-check: O'Neill & Eggins (2002) IW at the same T:

```
log10(fO2) = 2 * (-244118 + 115.559*T - 8.474*T*ln(T)) / (ln(10) * 8.31441 * T)
           ≈ -7.4078 at T = 2000 K
```

The 0.26 dex offset between the two buffers at 2000 K is the discrimination
guard's anchor: a regression that silently dispatches to the wrong buffer
would land 0.26 dex away from the expected value.

### Hirschmann (2021)

`OxygenFugacity('hirschmann')` evaluates Table 1 of Hirschmann (2021) at
P = 1 bar (1e-4 GPa). At T = 1473 K:

```
log10(fO2) = -11.957   (Fischer: -12.191, O'Neill & Eggins: -11.699)
```

The minus signs in Table 1 do not survive text extraction from the
published PDF, so the signs were taken from the rendered table. They are
checked independently by the continuity of the fcc/bcc and hcp branches
at the Eq. 18 boundary: the jump stays below 0.019 dex from 1000 to
3000 K, against a published maximum fit mismatch of 0.028 dex.

## Default buffer

The CALLIOPE default IW buffer is Fischer (2011). Tests that pin an IW
value MUST carry a discrimination guard so a change of default (or an
accidental config-side override) does not silently move the test's
reference point.

## Anchor type

Cross-implementation cross-check plus published benchmark. The
independent anchor is the cross-calibration offset between the Fischer
(2011) and O'Neill & Eggins (2002) IW fits, which are coded from separate
published formulae; the offset is independent of either single fit, so a
coefficient error in either moves it and fails the test. The coded
Fischer value is pinned as a secondary regression check, and the O'Neill
value also serves as the wrong-buffer discrimination guard.

## Cross-references

- `src/calliope/oxygen_fugacity.py` lines 26-34: implementations of both
  buffers, with a comment at line 32 pinning the `8.31441` constant to
  the literal O'Neill & Eggins (2002) value (do not replace with
  `constants.R_gas`).
- `docs/Explanations/oxygen_fugacity.md`: user-facing concept page with
  the empirical anchor (Earth ~ IW+3.5, Mars ~ IW, Mercury IW-3 to IW-5)
  for both buffers.
- `docs/Explanations/cross_backend_comparison.md`: empirical comparison
  of CALLIOPE Fischer vs atmodeller Hirschmann combined, showing the
  0.16 dex residual at Earth fiducial after the buffer-default flip.

## References

 [^cite-fischer2011]: R. A. Fischer, A. J. Campbell, G. A. Shofner, O. T. Lord, P. Dera, V. B. Prakapenka, *[Equation of state and phase diagram of FeO](https://doi.org/10.1016/j.epsl.2011.02.025)*, Earth and Planetary Science Letters, 304, 496-502, 2011. [SciX](https://scixplorer.org/abs/2011E%26PSL.304..496F/abstract).
 [^cite-hirschmann2021]: M. M. Hirschmann, *[Iron-wüstite revisited: A revised calibration accounting for variable stoichiometry and the effects of pressure](https://doi.org/10.1016/j.gca.2021.08.039)*, Geochimica et Cosmochimica Acta, 313, 74-84, 2021.
 [^cite-oneilleggins2002]: H. St. C. O'Neill, S. M. Eggins, *[The effect of melt composition on trace element partitioning: an experimental investigation of the activity coefficients of FeO, NiO, CoO, MoO$_2$ and MoO$_3$ in silicate melts](https://doi.org/10.1016/S0009-2541(01)00414-4)*, Chemical Geology, 186, 151-181, 2002. [SciX](https://scixplorer.org/abs/2002ChGeo.186..151O/abstract).
