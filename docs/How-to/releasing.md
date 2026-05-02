# Releasing

CALLIOPE uses [CalVer](https://calver.org/) versioning (`YY.0M.0D[.patch]`) and is published to [PyPI](https://pypi.org/project/fwl-calliope/) automatically when a GitHub release is tagged.

## Release checklist

1. **Regenerate requirements files**

    ```console
    python tools/generate_requirements_txt.py
    pip-compile -o requirements_full.txt pyproject.toml
    ```

2. **Bump the version**

    ```console
    bump-my-version bump release    # e.g. 25.05.04 → 26.05.02
    # or
    bump-my-version bump patch      # e.g. 26.05.02 → 26.05.02.1
    ```

    `bump-my-version` updates `src/calliope/__init__.py`, `pyproject.toml`, and `CITATION.cff` in one go.

3. **Commit and push** the version bump and the regenerated requirements files. Open a PR if the change touches anything substantive.

4. **Create a GitHub release**:
    - Go to [Releases](https://github.com/FormingWorlds/CALLIOPE/releases) and click *Draft a new release*.
    - Set the tag to the new version (e.g. `26.05.02`) and target `main`.
    - Auto-generate release notes from merged PRs since the previous release; keep the title plain (no leading `v`).
    - Publish.

5. **PyPI upload** is triggered automatically by [`publish.yaml`](https://github.com/FormingWorlds/CALLIOPE/actions/workflows/publish.yaml) when the release is published. Verify the new version appears on [PyPI](https://pypi.org/project/fwl-calliope/) within ~5 minutes.

6. **Update PROTEUS** if the bump introduces an API change: PROTEUS pins `fwl-calliope` in its `pyproject.toml`; bump the lower bound there and re-run the PROTEUS test suite.

## Hotfixes

For an urgent fix on top of an already-released version, bump the patch component (`bump-my-version bump patch`) so the new release inherits the previous date string. This keeps the release timeline aligned with the calendar without losing the date prefix.

## Yanking a release

If a release is broken, [yank it on PyPI](https://pypi.org/help/#yanked) immediately to prevent new installs from picking it up, then publish a follow-up patch release with the fix. Do not delete the GitHub release; mark it pre-release so users can still trace the timeline.
