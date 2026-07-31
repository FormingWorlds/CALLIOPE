"""Tests for the hook exclusions declared in `.pre-commit-config.yaml`.

Exercises the one property that pre-commit itself does not check on an ordinary
run: an `exclude` regex that stops matching any file is not an error, it is a
silent no-op. The excluded file quietly rejoins the set the rewriting hooks act
on, and the next contributor sees the style check go red for a reason that looks
unrelated to the path change that caused it.

Scope of the check. Every `exclude` in the config, both the config-level one and
each hook's own, must match at least one file from `git ls-files`. That file set
is a superset of what any single hook inspects, because the check does not narrow
by a hook's `files`, `types`, `types_or` or `exclude_types` the way pre-commit
does when it selects work. The check is therefore one-sided: it never fires
merely because a hook would have filtered a match out by type, and it does not
catch an exclusion that is dead only within one hook's narrower file set.

The check also requires the config to declare at least one exclusion, so it
cannot pass by having nothing left to inspect. Removing the last exclusion from
the config means removing this file along with it.

The CALLIOPE and MORS repositories keep this file identical, since both vendor
the same font assets behind the same exclusions. A change to one belongs in the
other.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO_ROOT / '.pre-commit-config.yaml'

# pre-commit's own "match nothing" pattern, and the value it applies to any hook
# that declares no exclusion. Matching nothing is what it is for, so it is never
# a stale pattern.
NULL_EXCLUDE = '^$'


def _declared_excludes(config):
    """Yield `(label, pattern)` for the config-level and per-hook exclusions."""
    candidates = [('<config>', config.get('exclude'))]
    for repo in config.get('repos') or ():
        for hook in repo.get('hooks') or ():
            candidates.append((hook.get('id', '<unnamed hook>'), hook.get('exclude')))
    for label, pattern in candidates:
        if pattern and pattern != NULL_EXCLUDE:
            yield label, pattern


def _stale_excludes(config, files):
    """Return `(label, pattern)` for every exclusion matching none of `files`.

    pre-commit selects files with `re.search`, not `re.match`, so the same call
    is used here to keep the two in agreement.
    """
    stale = []
    for label, pattern in _declared_excludes(config):
        regex = re.compile(pattern)
        if not any(regex.search(name) for name in files):
            stale.append((label, pattern))
    return stale


def _tracked_files():
    """Return the repository's tracked paths, or None when git cannot be run."""
    if shutil.which('git') is None:
        return None
    try:
        completed = subprocess.run(
            ['git', 'ls-files', '-z'],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            # Pinned rather than left to the locale, so a non-ASCII path under a
            # C locale yields surrogates instead of raising mid-call.
            encoding='utf-8',
            errors='surrogateescape',
            timeout=30,
            check=True,
        )
    except (subprocess.SubprocessError, OSError):
        return None
    return [name for name in completed.stdout.split('\0') if name]


def test_every_exclude_matches_a_tracked_file():
    """Every exclusion in the pre-commit config still applies to a real file.

    An exclusion that matches nothing has stopped protecting whatever it was
    written for, most often because the file moved or was renamed, and
    pre-commit reports no error when that happens.
    """
    files = _tracked_files()
    if files is None:
        pytest.skip('git unavailable, cannot enumerate tracked files')

    config = yaml.safe_load(CONFIG_PATH.read_text(encoding='utf-8'))
    declared = list(_declared_excludes(config))
    # Without this the test would pass vacuously on a config that had lost its
    # exclusions altogether, which is the very state being guarded.
    assert declared, f'{CONFIG_PATH.name} declares no exclude'

    stale = _stale_excludes(config, files)
    assert not stale, 'exclude matches no tracked file: ' + ', '.join(
        f'{label} -> {pattern}' for label, pattern in stale
    )


def test_stale_exclude_is_reported():
    """A dead exclusion is reported while live ones beside it are not.

    Models the realistic break: the vendored font assets move to a new
    directory and one hook's regex is left pointing at the old path. The
    config-level exclusion and the match-nothing pattern are both in the
    fixture so that neither is mistaken for a stale one.
    """
    files = ['docs/stylesheets/fonts/OFL.txt', 'README.md']
    config = {
        'exclude': r'^docs/stylesheets/',
        'repos': [
            {
                'hooks': [
                    {
                        'id': 'trailing-whitespace',
                        'exclude': r'^docs/stylesheets/fonts/OFL\.txt$',
                    },
                    {
                        'id': 'end-of-file-fixer',
                        'exclude': r'^docs/assets/fonts/OFL\.txt$',
                    },
                    # Unanchored, so it matches only under `search`. Were the
                    # check to anchor at position 0 the way `match` does, this
                    # live exclusion would be reported as stale.
                    {'id': 'check-yaml', 'exclude': r'fonts/OFL\.txt$'},
                    {'id': 'ruff', 'exclude': NULL_EXCLUDE},
                ]
            }
        ],
    }

    stale = _stale_excludes(config, files)
    assert stale == [('end-of-file-fixer', r'^docs/assets/fonts/OFL\.txt$')]
    # The surviving exclusions must stay out of the report, so the check
    # discriminates on the regex rather than flagging every hook it sees.
    assert {label for label, _ in stale} == {'end-of-file-fixer'}

    # A config-level exclusion is checked too, and goes stale the same way.
    config['exclude'] = r'^docs/assets/'
    assert ('<config>', r'^docs/assets/') in _stale_excludes(config, files)


def test_tracked_files_stands_down_instead_of_raising(monkeypatch):
    """`_tracked_files` answers None wherever git cannot supply a file list.

    The check has to stand down when the tree cannot be enumerated, so both
    routes to that state return None and the caller skips rather than
    reporting every exclusion as stale. A path that is not valid UTF-8
    survives the round trip as surrogates, so an undecodable name in the tree
    leaves the rest of the check working instead of raising mid-call.
    """
    monkeypatch.setattr(shutil, 'which', lambda name: None)
    assert _tracked_files() is None

    monkeypatch.setattr(shutil, 'which', lambda name: '/usr/bin/git')

    def _fail(*args, **kwargs):
        raise subprocess.CalledProcessError(128, 'git')

    monkeypatch.setattr(subprocess, 'run', _fail)
    assert _tracked_files() is None

    undecodable = b'docs/f\xffle.txt'.decode('utf-8', 'surrogateescape')
    seen = {}

    def _listing(*args, **kwargs):
        seen.update(kwargs)
        return subprocess.CompletedProcess(args, 0, f'README.md\0{undecodable}\0', '')

    monkeypatch.setattr(subprocess, 'run', _listing)
    assert _tracked_files() == ['README.md', undecodable]
    # Decoding is pinned rather than taken from the locale. Left to the
    # locale, a C-locale runner raises on the first non-ASCII path and the
    # exception escapes the clause above, which only covers process failures.
    assert (seen.get('encoding'), seen.get('errors')) == ('utf-8', 'surrogateescape')
    # The trailing NUL must not become an empty path, which would match any
    # unanchored exclusion and hide a stale one.
    assert '' not in _tracked_files()
