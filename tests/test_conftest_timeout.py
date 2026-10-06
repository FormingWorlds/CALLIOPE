"""The conftest collection hook stops a run whose tests lack a positive time limit."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

_GOOD = 'pytestmark = pytest.mark.timeout(30)\ndef test_a(): pass'
_BAD = {
    'first_mark_wins': 'pytestmark = [pytest.mark.timeout(0), pytest.mark.timeout(30)]\ndef test_a(): pass',
    'zero_after_tier': 'pytestmark = [pytest.mark.unit, pytest.mark.timeout(0)]\ndef test_a(): pass',
    'both_forms': 'pytestmark = pytest.mark.timeout(30, timeout=0)\ndef test_a(): pass',
    'function_override': _GOOD.replace('def', '@pytest.mark.timeout(0)\ndef'),
    'class_mark': 'class TestA:\n    pytestmark = pytest.mark.timeout(-1)\n    def test_a(self): pass',
    'missing': 'def test_a(): pass',
}


def _collect(tmp_path, sources):
    """Collect the given test modules in a fresh pytest run with this conftest."""
    (tmp_path / 'conftest.py').write_text(Path(__file__).with_name('conftest.py').read_text())
    for name, body in sources.items():
        (tmp_path / f'test_{name}.py').write_text(f'import pytest\n{body}\n')
    cmd = [sys.executable, '-m', 'pytest', '--collect-only', '-q', '-p', 'no:cacheprovider']
    return subprocess.run([*cmd, str(tmp_path)], capture_output=True, text=True, cwd=tmp_path)


def test_hook_names_every_test_without_a_positive_budget(tmp_path):
    """Each bad module fails collection by name; the module with a 30 s mark is not named."""
    out = _collect(tmp_path, {**_BAD, 'good': _GOOD})
    assert out.returncode == pytest.ExitCode.USAGE_ERROR, out.stdout + out.stderr
    named = out.stderr.split('positive timeout mark:')[-1].split()
    expected = [f'test_{n}.py::test_a' for n in _BAD if n != 'class_mark']
    assert sorted(named) == sorted([*expected, 'test_class_mark.py::TestA::test_a'])
