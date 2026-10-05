"""The autouse fixture seeds the global NumPy RNG and restores the caller's state."""

from __future__ import annotations

import numpy as np
import pytest
from conftest import GLOBAL_RNG_SEED, seeded_global_rng

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_each_test_starts_from_the_seeded_state():
    """A test that does not seed draws the sequence of GLOBAL_RNG_SEED."""
    expected = np.random.RandomState(GLOBAL_RNG_SEED).random_sample(5)
    np.testing.assert_array_equal(np.random.random_sample(5), expected)
    assert not np.array_equal(
        expected, np.random.RandomState(GLOBAL_RNG_SEED + 1).random_sample(5)
    )


def test_two_seeded_blocks_draw_the_same_sequence_and_restore_the_state():
    """Two seeded blocks give the same draws, and the caller's sequence resumes after each."""
    np.random.seed(123)
    caller = np.random.RandomState(123)
    draws = []
    for _ in range(2):
        with seeded_global_rng(7):
            draws.append(np.random.random_sample(4))
        np.testing.assert_array_equal(np.random.random_sample(3), caller.random_sample(3))
    np.testing.assert_array_equal(draws[0], draws[1])
    np.testing.assert_array_equal(draws[0], np.random.RandomState(7).random_sample(4))
