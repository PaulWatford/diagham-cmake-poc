# Releasing DiagHam

The full checklist is `docs/how-to/develop/release.md`. In short:

1. `main` is green on every CI job.
2. Regenerate the generated files (program reference, coverage page, manifests, site navigation).
3. Set the version in `CMakeLists.txt`, `CITATION.cff` and `codemeta.json`, and date the release in `CHANGELOG.md`. `python3 scripts_cmake/check_version.py` must pass.
4. Write the release notes in `docs/history/release-notes/`.
5. Commit, tag `vYEAR.MONTH.PATCH`, push the tag, and publish the release.
6. For a stable release (a version 1.0, or a version described in a paper), archive it on Zenodo. On GitHub, the Zenodo integration does this when the release is published; on GitLab, upload the release archive to Zenodo by hand. `.zenodo.json` holds the record's metadata.
7. After the first archived release, add the Zenodo concept DOI to `CITATION.cff` (`identifiers:`) and to the README, and commit.

Each paper names the DiagHam version it describes and links to that version's Zenodo record.

The creators of every Zenodo record are the DiagHam authors, as listed in `CITATION.cff`; `scripts_cmake/check_version.py` checks that `.zenodo.json` agrees with it. Who maintains the Zenodo releases after the current maintainer is for the project lead to name.
