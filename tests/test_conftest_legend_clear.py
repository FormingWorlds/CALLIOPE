"""Tests for the legend_clear check of `tests/conftest.py`.

Contract: the check fails when a legend covers a line, a marker, a scatter point or a
filled band, passes when the legend is clear, and closes the figure. Commands:
docs/How-to/build_tests.md.
"""

from __future__ import annotations

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _figure(kind):
    fig, ax = plt.subplots()
    ax.set(xlim=(0, 1), ylim=(0, 1))
    if kind == 'line':
        ax.plot([0, 1], [0.5, 0.5], label='line')
    elif kind == 'band':
        ax.fill_between([0, 1], [0, 0], [1, 1], alpha=0.2, label='band')
    elif kind == 'scatter':
        ax.scatter([0.5], [0.5], label='point')
    else:
        ax.plot([0, 1], [0.02, 0.02], label='low line')
    ax.legend(loc='center')
    return fig


@pytest.mark.parametrize('kind', ['line', 'band', 'scatter'])
def test_legend_over_a_mark_fails(legend_clear, kind):
    """A legend in the centre over a line, a band or a scatter point fails the check."""
    fig = _figure(kind)
    with pytest.raises(AssertionError, match='legend covers'):
        legend_clear(fig)
    assert not plt.fignum_exists(fig.number)


def test_clear_legend_passes_and_closes_the_figure(legend_clear):
    """A legend in the centre above a line near the bottom passes and closes the figure."""
    fig = _figure('clear')
    legend_clear(fig)
    assert not plt.fignum_exists(fig.number)
