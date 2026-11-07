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


def pytest_sessionstart(session):
    """Build the matplotlib font cache before the first plotting test and its time limit."""
    from matplotlib import font_manager  # noqa: F401


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
    0.5 N2 + 1.5 H2 = NH3: log Kf(NH3), since N2 and H2 are reference states.
    """

    def per_mole(product):
        rows = _janaf_log_kf(product)
        if product == 'H2S':
            s2 = _janaf_log_kf('S2')
            rows = {T: y - 0.5 * s2[T] for T, y in rows.items() if T in s2}
        T = np.array(sorted(rows))
        return T, np.array([rows[t] for t in T])

    return per_mole


def _legend_clear(fig):
    """Assert that no legend box covers a line, a marker or a filled band; close the figure."""
    try:
        fig.canvas.draw()
        for ax in fig.axes:
            if ax.get_legend() is not None:
                _check_legend(ax, ax.get_legend().get_window_extent(fig.canvas.get_renderer()))
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


def _check_legend(ax, box):
    def covered(p):
        return (
            (p[:, 0] >= box.x0)
            & (p[:, 0] <= box.x1)
            & (p[:, 1] >= box.y0)
            & (p[:, 1] <= box.y1)
        ).any()

    for line in ax.get_lines():
        xy = line.get_transform().transform(line.get_xydata())
        for a, b in zip(xy[:-1], xy[1:]):
            pts = a + np.linspace(0.0, 1.0, int(np.hypot(*(b - a))) + 2)[:, None] * (b - a)
            assert not covered(pts), f'legend covers {line.get_label()!r}'
        if line.get_marker() not in ('None', None, '') and len(xy):
            assert not covered(xy), f'legend covers {line.get_label()!r}'
    grid = np.meshgrid(np.linspace(box.x0, box.x1, 20), np.linspace(box.y0, box.y1, 20))
    corners = np.reshape(grid, (2, -1)).T
    for coll in ax.collections:
        offsets = coll.get_offsets()
        if len(offsets) and coll.get_offset_transform() is not None:
            pts = coll.get_offset_transform().transform(offsets)
            assert not covered(pts), f'legend covers {coll.get_label()!r}'
        for path in coll.get_paths():
            inside = path.contains_points(corners, transform=coll.get_transform())
            assert not inside.any(), f'legend covers {coll.get_label()!r}'


@pytest.fixture
def legend_clear():
    """The check of :func:`_legend_clear` for a test that makes a plot."""
    return _legend_clear
