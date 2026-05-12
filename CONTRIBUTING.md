# Contributing guidelines

## Development

### Building the documentation

The documentation is written in [Markdown](https://www.markdownguide.org/basic-syntax/) and built with [Zensical](https://zensical.org/), the modern static site generator from the Material for MkDocs team. Zensical reads the existing `mkdocs.yml` directly so the source format is unchanged from a stock MkDocs project.

To build and serve the documentation for yourself:

```console
pip install -e .[docs]
zensical serve
```

For a one-shot production build:

```console
zensical build --clean
```

The build artefacts land in `site/` (gitignored). You can find the documentation source in the [docs](https://github.com/FormingWorlds/CALLIOPE/tree/main/docs) directory. If you are adding new pages, update the listing in [`mkdocs.yml`](https://github.com/FormingWorlds/CALLIOPE/blob/main/mkdocs.yml) under the `nav` entry.

The documentation is hosted at [proteus-framework.org/CALLIOPE](https://proteus-framework.org/CALLIOPE).

### Running tests

CALLIOPE uses [pytest](https://docs.pytest.org/en/latest/) with a four-tier marker scheme (`unit`, `smoke`, `integration`, `slow`). Common selections:

```console
pytest -m unit                              # fast in-process tests
pytest -m smoke                             # minimal-config solver tests
pytest -m integration                       # full multi-species CHNS solves
pytest -m "(unit or smoke) and not skip"    # PR-gate selection
pytest -m "not skip"                        # full nightly selection
```

To check coverage:

```console
pytest --cov=src/calliope --cov-report=term -m "not skip"
pytest --cov=src/calliope --cov-report=html -m "not skip"   # htmlcov/
```

For details on the marker scheme, badge system, and coverage gate, see the [testing suite](https://proteus-framework.org/CALLIOPE/Explanations/testing.html) explanation. To add a new test, see the [build a new test](https://proteus-framework.org/CALLIOPE/How-to/build_tests.html) how-to.

### Making a release

The versioning scheme is [CalVer](https://calver.org/).

0. Update requirements files:

   ```console
   python tools/generate_requirements_txt.py
   pip-compile -o requirements_full.txt pyproject.toml
   ```

1. Tag the release on `main` (CalVer `YY.MM.DD`, bare tag with no leading `v`):

   ```console
   git checkout main
   git pull
   git tag 26.05.10
   git push origin 26.05.10
   ```

   `setuptools-scm` derives the package version from the tag at build time; no source files need editing.

2. Create a new [release](https://github.com/FormingWorlds/CALLIOPE/releases) on GitHub against that tag, e.g. `26.05.10`.

3. The [upload to PyPI](https://pypi.org/project/fwl-calliope) is triggered automatically when the release is published, handled by [this workflow](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/publish.yaml). See [How-to: releasing](docs/How-to/releasing.md) for the full procedure.
