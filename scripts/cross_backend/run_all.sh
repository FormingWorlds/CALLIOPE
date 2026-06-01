#!/usr/bin/env bash
# Regenerate every figure on the cross-backend comparison page.
#
# Run from the CALLIOPE repository root with the `proteus` conda env
# active (or any env with calliope + atmodeller installed).
#
# Each figure script is independent and writes its own PDF + PNG into
# docs/assets/figures/cross_backend/ and its raw CSV data into
# scripts/cross_backend/data/. Total wall time is dominated by
# atmodeller solver calls (~15 s each warm; ~60 s cold for the first
# JAX compile).
#
# Approximate total runtime on a 2024 M-series Mac: ~20 minutes.
set -euo pipefail

cd "$(dirname "$0")/../.."  # repository root

echo "=== Fig 1 (buffers) ==="
python3 -m scripts.cross_backend.fig1_buffers

echo "=== Fig 2 (round-trip) ==="
python3 -m scripts.cross_backend.fig2_roundtrip

echo "=== Fig 3 (grid) ==="
python3 -m scripts.cross_backend.fig3_grid

echo "=== Fig 4 (attribution) ==="
python3 -m scripts.cross_backend.fig4_attribution

echo "=== Fig 5 (Earth anchor) ==="
python3 -m scripts.cross_backend.fig5_earth_anchor

echo
echo "All figures written to docs/assets/figures/cross_backend/"
echo "Raw data written to scripts/cross_backend/data/"
