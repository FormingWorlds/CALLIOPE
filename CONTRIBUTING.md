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

CALLIOPE uses [pytest](https://docs.pytest.org/en/latest/) to run the tests. You can run them with:

```console
pytest
```

To check coverage:

```console
coverage run -m pytest
coverage report   # text summary in the terminal
coverage html     # HTML report under htmlcov/
```

### Making a release

The versioning scheme is [CalVer](https://calver.org/).

0. Update requirements files:

   ```console
   python tools/generate_requirements_txt.py
   pip-compile -o requirements_full.txt pyproject.toml
   ```

1. Bump the version (`release` / `patch` as needed):

   ```console
   bump-my-version bump release
   # e.g. 25.05.04 → 26.05.02
   ```

2. Commit and push your changes.

3. Create a new [release](https://github.com/FormingWorlds/CALLIOPE/releases) on GitHub. Set the tag to the specified version, e.g. `26.05.02`.

4. The [upload to PyPI](https://pypi.org/project/fwl-calliope) is triggered when a release is published, handled by [this workflow](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/publish.yaml).
