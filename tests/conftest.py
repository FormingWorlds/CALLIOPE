"""Shared test configuration for CALLIOPE."""

from __future__ import annotations

import math
from pathlib import Path

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


def _janaf_log_kf(species):
    """Map T [K] to log10 Kf in the NIST-JANAF table copy ``tools/janaf_<species>.dat``."""
    tab = np.loadtxt(Path(__file__).resolve().parents[1] / 'tools' / f'janaf_{species}.dat')
    return dict(zip(tab[:, 0], tab[:, 7]))


@pytest.fixture(scope='session')
def janaf_log_k():
    """Return ``f(product) -> (T, log10 K)`` per mole of product from the JANAF rows.

    0.5 S2 + H2 = H2S: log Kf(H2S) - 0.5 log Kf(S2). Both tables share the sulfur
    reference state, which cancels, so the combination holds on both sides of 882 K.
    """

    def per_mole(product):
        h2s, s2 = _janaf_log_kf('H2S'), _janaf_log_kf('S2')
        rows = {'H2S': {T: h2s[T] - 0.5 * s2[T] for T in h2s.keys() & s2.keys()}}[product]
        T = np.array(sorted(rows))
        return T, np.array([rows[t] for t in T])

    return per_mole
