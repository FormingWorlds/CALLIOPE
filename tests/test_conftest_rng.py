"""The autouse fixture seeds the global NumPy RNG with 0 before every test."""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.parametrize('run', [1, 2])
def test_each_test_starts_from_the_seeded_state(run):
    """Every test, also the second of two, draws the seed-0 sequence before it seeds."""
    expected = np.random.RandomState(0).random_sample(5)
    np.testing.assert_array_equal(np.random.random_sample(5), expected)
    assert not np.array_equal(expected, np.random.RandomState(run).random_sample(5))
