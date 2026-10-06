"""The conftest collection hook stops a run whose tests lack a positive time limit."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

_T = '\ndef test_a(): pass'
_GOOD = ['timeout(30)', 'timeout(timeout=30)', "timeout(0.5, 'thread')"]
_BAD = {
    'first_mark_wins': 'pytestmark = [pytest.mark.timeout(0), pytest.mark.timeout(30)]' + _T,
    'zero_after_tier': 'pytestmark = [pytest.mark.unit, pytest.mark.timeout(0)]' + _T,
    'both_forms': 'pytestmark = pytest.mark.timeout(30, timeout=0)' + _T,
    'function_override': 'pytestmark = pytest.mark.timeout(30)\n@pytest.mark.timeout(0)' + _T,
    'class_mark': 'class TestA:\n    pytestmark = pytest.mark.timeout(-1)\n    def test_a(self): pass',
    'infinite': "pytestmark = pytest.mark.timeout(float('inf'))" + _T,
    'boolean': 'pytestmark = pytest.mark.timeout(True)' + _T,
    'string': "pytestmark = pytest.mark.timeout('30')" + _T,
    'missing': _T,
}


def test_hook_names_every_test_without_a_positive_budget(tmp_path):
    """Each bad module stops collection by name; the valid budget forms are not named."""
    (tmp_path / 'conftest.py').write_text(Path(__file__).with_name('conftest.py').read_text())
    good = {f'good{k}': f'pytestmark = pytest.mark.{m}{_T}' for k, m in enumerate(_GOOD)}
    for name, body in {**_BAD, **good}.items():
        (tmp_path / f'test_{name}.py').write_text(f'import pytest\n{body}\n')
    cmd = [sys.executable, '-m', 'pytest', '--collect-only', '-q', '-p', 'no:cacheprovider']
    out = subprocess.run([*cmd, str(tmp_path)], capture_output=True, text=True, cwd=tmp_path)
    assert out.returncode == pytest.ExitCode.USAGE_ERROR, out.stdout + out.stderr
    named = out.stderr.split('positive timeout mark:')[-1].split()
    assert sorted(n.split('::')[0] for n in named) == sorted(f'test_{n}.py' for n in _BAD)
