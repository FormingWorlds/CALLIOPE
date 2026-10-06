"""Shared test configuration for CALLIOPE."""

from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture(autouse=True)
def _seed_global_rng():
    """Start every test from global NumPy RNG seed 0, so its solver restarts repeat."""
    np.random.seed(0)
