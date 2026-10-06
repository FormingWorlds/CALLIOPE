"""Shared test configuration for CALLIOPE."""

from __future__ import annotations

import math

import numpy as np
import pytest


def _positive_budget(mark) -> bool:
    """True when a timeout mark carries exactly one finite positive budget."""
    if mark is None:
        return False
    budgets = list(mark.args[:1]) + (
        [mark.kwargs['timeout']] if 'timeout' in mark.kwargs else []
    )
    return len(budgets) == 1 and type(budgets[0]) in (int, float) and 0 < budgets[0] < math.inf


def pytest_collection_modifyitems(config, items):
    """Stop the run when a test has no positive time limit.

    Reads ``item.get_closest_marker('timeout')``, the mark pytest-timeout applies, and takes
    its budget from the first positional argument or the ``timeout`` keyword.
    """
    bad = [i.nodeid for i in items if not _positive_budget(i.get_closest_marker('timeout'))]
    if bad:
        raise pytest.UsageError('tests without a positive timeout mark: ' + ' '.join(bad))


@pytest.fixture(autouse=True)
def _seed_global_rng():
    """Start every test from global NumPy RNG seed 0, so its solver restarts repeat."""
    np.random.seed(0)
