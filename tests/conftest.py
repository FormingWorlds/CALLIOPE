"""Shared test configuration for CALLIOPE."""

from __future__ import annotations

import math

import numpy as np
import pytest


def _positive_budget(mark) -> bool:
    """True when a timeout mark carries exactly one finite positive numeric budget."""
    if mark is None:
        return False
    budgets = [*mark.args[:1], *(v for k, v in mark.kwargs.items() if k == 'timeout')]
    b = budgets[0] if len(budgets) == 1 else None
    return isinstance(b, (int, float)) and not isinstance(b, bool) and 0 < b < math.inf


def pytest_collection_modifyitems(items):
    """Stop the run when a test, selected or not, has no positive time limit.

    Reads the closest ``timeout`` mark, which is the mark pytest-timeout applies.
    """
    bad = [i.nodeid for i in items if not _positive_budget(i.get_closest_marker('timeout'))]
    if bad:
        raise pytest.UsageError('tests without a positive timeout mark: ' + ' '.join(bad))


@pytest.fixture(autouse=True)
def _seed_global_rng():
    """Start every test from global NumPy RNG seed 0, so its solver restarts repeat."""
    np.random.seed(0)
