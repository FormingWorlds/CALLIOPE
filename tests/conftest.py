"""Shared test configuration for CALLIOPE."""

from __future__ import annotations

from contextlib import contextmanager

import numpy as np
import pytest

GLOBAL_RNG_SEED = 0


@contextmanager
def seeded_global_rng(seed: int = GLOBAL_RNG_SEED):
    """Seed the global NumPy RNG for the block and restore the caller's state after it.

    Parameters
    ----------
    seed : int
        Seed passed to ``np.random.seed``.
    """
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(state)


@pytest.fixture(autouse=True)
def _seed_global_rng():
    """Run every test from the same global NumPy RNG state, so its solver restarts repeat."""
    with seeded_global_rng():
        yield
