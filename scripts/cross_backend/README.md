# Cross-backend comparison harness

This directory contains the reusable harness that produces the figures
embedded in `docs/Explanations/cross_backend_comparison.md`.

The harness is investigation tooling. It is not shipped with the
installed `fwl-calliope` distribution and not part of the test suite.

## What's here

| File | Purpose |
|---|---|
| `inventories.py` | Earth BSE H / C / N / S inventory (Krijt et al. 2023 PPVII Tables 1+2); volatile-O reference derived self-consistently at Delta-IW = +3.5 (Sossi 2020) |
| `buffers.py` | Analytical IW buffer formulae: O'Neill & Eggins 2002, Fischer 2011, Hirschmann composite |
| `runners.py` | Backend-uniform call wrappers: `run_calliope`, `run_atmodeller` |
| `verification.py` | Pre-flight per-backend round-trip checks (callable as `python -m scripts.cross_backend.verification`) |
| `plot_style.py` | Shared matplotlib styling and output paths |
| `fig1_buffers.py` ... `fig5_earth_anchor.py` | One script per figure |
| `run_all.sh` | One-shot regenerator |
| `data/` | Raw CSV output from each figure script (created on first run) |

## How to re-run

```bash
# from the repo root, with the proteus conda env active
bash scripts/cross_backend/run_all.sh
```

To regenerate a single figure:

```bash
python3 -m scripts.cross_backend.fig3_grid
```

To run the verification harness without producing any figure:

```bash
python3 -m scripts.cross_backend.verification
```

## Wall time

Approximate timings on a 2024 M-series Mac (proteus conda env):

| Figure | Wall time |
|---|---|
| Fig 1 (analytical) | < 1 s |
| Fig 2 (round-trip, 4 T x 4 dIW x 2 backends) | 5-10 min |
| Fig 3 (grid, 4 T x 5 O-factor x 2 backends) | 5-10 min |
| Fig 4 (attribution, 3 calls) | 1 min |
| Fig 5 (Earth anchor, 2 calls) | 30 s |

The dominant cost is the first atmodeller call per Python process
(JAX compile, ~60 s); subsequent calls are ~15 s warm.

## Re-using on a different fiducial

The pattern is:

```python
from scripts.cross_backend.inventories import Inventory, EARTH_BSE_KRIJT23, scale_O
from scripts.cross_backend.runners import run_calliope, run_atmodeller

# Build a different inventory
mars_inv = Inventory(name='Mars (placeholder)',
                     H=..., C=..., N=..., O=..., S=...,
                     citation='Wanke & Dreibus 1988')

# Run both backends
cal = run_calliope(mars_inv, T_magma=1800.0, fO2_hint=-1.0)
atm = run_atmodeller(mars_inv, T_magma=1800.0)
print(cal.fO2_shift_derived, atm.fO2_shift_derived)
```

For sensitivity sweeps, the `scale_O` helper varies only the O budget
while keeping H/C/N/S fixed (used by `fig3_grid.py`).

## Provenance

Every figure script saves its raw inputs and outputs to
`data/fig<N>_*.csv`. The committed CSVs are the provenance for the
PDFs and PNGs in `docs/assets/figures/cross_backend/`. Re-running a
figure script overwrites the CSV; re-running at a different
calliope / atmodeller commit will change the numbers.
