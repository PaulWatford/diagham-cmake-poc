# Bring new DiagHam revisions in from Subversion

Purpose: update `upstream` and `main` when the DiagHam authors commit to their Subversion repository, without ever editing `upstream` by hand.
Source: new (2026-09-28); `scripts_cmake/sync_upstream.sh`, `scripts_cmake/svn-authors.txt`, the weekly *upstream drift* job in `.github/workflows/ci.yml`. The procedure was exercised by re-importing r4491–r4493 onto an `upstream` reset to r4490 and checking the trees.

## How the two repositories relate

The canonical DiagHam is Subversion at `https://www.nick-ux.org/diagham/svn/DiagHam` (trunk, branches, tags). This repository carries it on the `upstream` branch as a `git svn` conversion — every trunk revision as one commit, author by name only (`scripts_cmake/svn-authors.txt`), the revision in the `git-svn-id:` trailer — and never edits that branch by hand. `main` merges `upstream` and adds the build, tests, docs and a few fixes. `scripts_cmake/upstream-revision.txt` records the revision `main` carries.

```
SVN trunk --(git svn fetch)--> conversion clone --(rebase onto upstream)--> upstream --(merge)--> main
```

## Once: trust the server's certificate

nick-ux.org uses a certificate no public CA signs. Subversion (which `git svn` uses) must be told once to trust it, on the machine that runs the sync:

```
svn info --non-interactive --trust-server-cert-failures=unknown-ca,cn-mismatch,expired,not-yet-valid,other https://www.nick-ux.org/diagham/svn/DiagHam/trunk
```

writes nothing; for `git svn` the trust must be stored: run `svn info https://www.nick-ux.org/diagham/svn/DiagHam/trunk` interactively once and answer `p` (permanently), or copy the entry under `~/.subversion/auth/svn.ssl.server/` from a machine where it exists.

## Once: the conversion clone

The script keeps a `git svn` clone next to this repository (`../DiagHam-svn-full`, or `DIAGHAM_SVN_CLONE`). If it is absent the script creates it — a full clone of 4,500 revisions takes about an hour; afterwards a fetch takes seconds.

## Each time

```
git checkout main && git status        # must be clean
scripts_cmake/sync_upstream.sh --dry-run   # what is new
scripts_cmake/sync_upstream.sh             # import, merge, regenerate, record, build, test
```

What the script does, and refuses to do:

1. Fetches new revisions into the clone with the names-only authors file.
2. Finds the commit of the recorded revision on both sides and checks their **trees are identical** (our `upstream` has author corrections, so SHAs differ; content must not).
3. Rebases the new conversion commits onto `upstream`, checks the resulting tree equals the clone's trunk, and fast-forwards `upstream`. Nothing is ever written to `upstream` except this.
4. Merges `upstream` into `main`. A conflict means upstream touched a file this repository patched: resolve it (if upstream fixed the same thing, drop our fix and its row in `docs/explanation/upstream-fixes.md`), commit, rerun the script.
5. Regenerates the per-directory `CMakeLists.txt` (`extract_autotools.py`), writes the new revision to `upstream-revision.txt`, commits.
6. Configures, builds and runs `ctest` (default preset). A failing golden after a sync is information: either a golden was wrong, or upstream changed behaviour. Read the failure; never adjust an expected value to make it pass.

Then, if programs or options changed: regenerate the program reference and the coverage page (the script prints the commands), commit, push.

## As a pull request (recommended)

```
scripts_cmake/sync_upstream.sh --pr
```

does the same on a branch `sync/r<rev>`, pushes it and opens the pull request through your `gh` login, so the commits, the push and the pull request are yours and CI runs before anything reaches `main`. Review the upstream commits in the PR, wait for the checks, merge. The weekly job never does this by itself: an automatic sync would put unreviewed physics changes on `main` under a bot's name.

## The weekly drift check

The scheduled CI job compares trunk's revision with `upstream-revision.txt` and, when trunk has moved on, opens (or updates) an issue labelled `upstream-drift` naming the revisions. Bringing them in is deliberate: this script with `--pr`, run by a person, reviewed as a pull request; never automatic.

## Branches and tags of the Subversion repository

Only trunk is carried on `upstream`. The Subversion branches (`testing_nr`, `testing_zp`, `regnault`) and the tag `first-release` exist in the conversion clone; they were not part of the migration and are not synchronised.
