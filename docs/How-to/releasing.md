# Releasing

CALLIOPE uses [setuptools-scm](https://setuptools-scm.readthedocs.io/) for automatic versioning. The version is derived from git tags, so there are no version strings to edit manually anywhere in the codebase.

## Version scheme

CALLIOPE follows [Calendar Versioning](https://calver.org/) (CalVer) with the format `YY.MM.DD`:

- **Release versions** (from tags): `26.05.10`, `26.06.01`
- **Development versions** (between tags): `26.05.10.post1.dev5+gabc1234` (5 commits after the `26.05.10` tag, at commit `abc1234`)

Dev builds carry the `.post1.devN` suffix, not `.devN` of a guessed next version. This is deliberate: the default `setuptools-scm` `guess-next-dev` scheme would label an untagged commit as `26.5.11.devN+...`, claiming a future date that would collide with the actual next daily release. The `version_scheme = "no-guess-dev"` setting in `pyproject.toml` keeps dev builds anchored on the most recent real tag.

## How to make a release

### 1. Ensure main is ready

All CI checks should pass on the `main` branch. Verify locally:

```console
pytest -m "(unit or smoke) and not skip"
ruff check src/ tests/ tools/
```

### 2. Tag the release

From the `main` branch (after merging any feature branches):

```console
git checkout main
git pull
git tag 26.05.10        # Use today's date in YY.MM.DD format
git push origin 26.05.10
```

!!! warning "Tag format"
    Tags are bare CalVer dates with no leading `v` (`26.05.10`, not `v26.05.10`). setuptools-scm consumes the tag verbatim, so prefixed tags would produce a non-PEP-440 version string and fail the build.

### 3. Create a GitHub Release

```console
gh release create 26.05.10 --title "26.05.10" --generate-notes
```

Or create the release manually at [github.com/FormingWorlds/CALLIOPE/releases/new](https://github.com/FormingWorlds/CALLIOPE/releases/new):

1. Choose the tag you just pushed (`26.05.10`)
2. Set the title to the version (`26.05.10`)
3. Click "Generate release notes" for an automatic changelog
4. Click "Publish release"

### 4. PyPI publication (automatic)

Publishing the GitHub Release triggers the `publish.yaml` workflow, which:

1. Checks out the repository with full git history (for tag resolution)
2. Builds the package with `python -m build` (setuptools-scm writes the version from the tag)
3. Publishes to PyPI via trusted publishing (no API tokens needed)

The package appears on [pypi.org/project/fwl-calliope](https://pypi.org/project/fwl-calliope/) within a few minutes.

### 5. Verify

```console
pip install --upgrade fwl-calliope
python -c "from calliope import __version__; print(__version__)"
# Should print: 26.05.10
```

## How versioning works

setuptools-scm reads the git history to determine the version:

| Scenario | Example version |
|----------|----------------|
| Exactly on a tag (`26.05.10`) | `26.05.10` |
| 5 commits after tag, dirty | `26.05.10.post1.dev5+gabc1234.d20260515` |
| No tags in history | `0.0.0` (fallback) |

The version is written to `src/calliope/_version.py` at install time (this file is gitignored). At runtime, `__init__.py` imports it:

```python
from calliope import __version__
```

## Multiple releases per day

If you need a second release on the same day, append a patch number:

```console
git tag 26.05.10.1
```

This produces version `26.05.10.1`.

## CITATION.cff

`CITATION.cff` carries a `version:` field that is not auto-updated by setuptools-scm. Bump it by hand when cutting a release if you want the citation metadata to track the package version.

## Troubleshooting

### Version shows `0.0.0`

setuptools-scm cannot find any tags. Check:

```console
git tag --list
git describe --tags
```

If no tags appear, you may need to fetch them:

```console
git fetch --tags
```

### CI build shows wrong version

The CI checkout must use `fetch-depth: 0` (full history). This is already configured in `publish.yaml`.
