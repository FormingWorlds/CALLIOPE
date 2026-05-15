# Validation: `src/calliope/oxygen_fugacity.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the
behaviour of `calliope.oxygen_fugacity` against published sources.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_oxygen_fugacity.py::test_oxygen_fugacity_fischer_value_at_2000K_matches_published_fit` | Fischer et al. (2011), EPSL 304, 496, Eq. 2 | [doi:10.1016/j.epsl.2011.02.011](https://doi.org/10.1016/j.epsl.2011.02.011) | Pins Fischer IW value at T = 2000 K against the closed-form fit `6.94059 - 28.1808e3 / T`; includes a wrong-buffer discrimination guard against O'Neill & Eggins (2002) at the same T. |

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

## Default buffer history

The CALLIOPE default flipped from `'oneill'` to `'fischer'` in the
2026-05 sweep that introduced the `from_O_budget` Path C entry point.
Tests that pin an IW value MUST carry a discrimination guard so a future
default flip (or accidental config-side override) does not silently
change the test's reference point. This is the buffer-flip trap from
`.github/.claude/rules/calliope-tests.md` Section 16.

## Anchor type

Published benchmark + cross-buffer discrimination. The Fischer 2011 cite
is the published-benchmark anchor; the O'Neill 2002 value at the same T
is the discrimination guard required by Section 2 rule 4 of
`calliope-tests.md`.

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
