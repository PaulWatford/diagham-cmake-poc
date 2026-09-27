#!/usr/bin/env bash
# sync_upstream.sh: bring new DiagHam Subversion revisions onto `upstream`, merge them
# into `main`, regenerate the CMake files, record the revision, build and test.
#
# The canonical DiagHam is Subversion at nick-ux.org. `upstream` is its faithful
# conversion (names-only authors, never hand-edited); `main` = upstream + this
# repository's build, tests, docs and fixes. New revisions therefore travel
#     SVN trunk --(git svn fetch)--> conversion clone --(rebase onto upstream)--> upstream --(merge)--> main
#
# Usage (from the repository root, on a clean `main`):
#   scripts_cmake/sync_upstream.sh [--dry-run] [--no-build] [--to REV]
#
#   --dry-run   fetch and report what would change; touch nothing in this repository
#   --no-build  skip the configure/build/ctest step after merging
#   --to REV    stop at SVN revision REV instead of trunk HEAD (for tests)
#
# Environment:
#   DIAGHAM_SVN_CLONE   the git-svn conversion clone (default: ../DiagHam-svn-full, created
#                       with `git svn clone --stdlayout -A scripts_cmake/svn-authors.txt` if absent
#                       -- the initial clone of 4,500 revisions takes about an hour)
#   DIAGHAM_SVN_URL     default https://www.nick-ux.org/diagham/svn/DiagHam
#
# The server's certificate is not signed by a public CA: Subversion must have been
# told once to trust it (see docs/how-to/develop/sync-upstream.md).
set -euo pipefail

DRY_RUN=0; BUILD=1; TO_REV=""
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY_RUN=1 ;;
        --no-build) BUILD=0 ;;
        --to) TO_REV="$2"; shift ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
    esac
    shift
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
SVN_URL="${DIAGHAM_SVN_URL:-https://www.nick-ux.org/diagham/svn/DiagHam}"
CLONE="${DIAGHAM_SVN_CLONE:-$REPO_ROOT/../DiagHam-svn-full}"
AUTHORS="$REPO_ROOT/scripts_cmake/svn-authors.txt"
RECORDED=$(tr -dc 0-9 < scripts_cmake/upstream-revision.txt)

say() { echo "== $*"; }

if [ -n "$(git status --porcelain)" ]; then
    echo "working tree not clean; commit or stash first" >&2; exit 1
fi
if [ "$(git rev-parse --abbrev-ref HEAD)" != main ]; then
    echo "run this on main" >&2; exit 1
fi

# --- 1. fetch ---------------------------------------------------------------------
if [ ! -d "$CLONE/.git" ]; then
    say "no conversion clone at $CLONE: cloning the Subversion repository (long)"
    git svn clone --stdlayout -A "$AUTHORS" "$SVN_URL" "$CLONE"
fi
say "fetching new revisions into $CLONE"
git -C "$CLONE" svn fetch -A "$AUTHORS"
TRUNK=$(git -C "$CLONE" rev-parse refs/remotes/svn/trunk)
rev_of() { git -C "$1" log -1 --format=%B "$2" | sed -n 's/^git-svn-id: .*@\([0-9]*\) .*/\1/p'; }
HEAD_REV=$(rev_of "$CLONE" "$TRUNK")
say "trunk is at r$HEAD_REV; this repository carries r$RECORDED"

# the conversion commit for the recorded revision, and the target commit
find_rev() {   # find_rev <repo> <ref> <rev> -> commit whose git-svn-id is @rev
    git -C "$1" log --format=%H --grep="^git-svn-id: .*@$3 " "$2" | tail -1
}
BASE=$(find_rev "$CLONE" "$TRUNK" "$RECORDED")
if [ -z "$BASE" ]; then
    echo "cannot find r$RECORDED in the conversion clone" >&2; exit 1
fi
TARGET="$TRUNK"; TARGET_REV="$HEAD_REV"
if [ -n "$TO_REV" ]; then
    TARGET=$(find_rev "$CLONE" "$TRUNK" "$TO_REV"); TARGET_REV="$TO_REV"
    [ -n "$TARGET" ] || { echo "cannot find r$TO_REV" >&2; exit 1; }
fi
NEW=$(git -C "$CLONE" rev-list --count "$BASE".."$TARGET")
if [ "$NEW" -eq 0 ]; then
    say "up to date: nothing after r$RECORDED"; exit 0
fi
say "$NEW new revision(s), r$((RECORDED + 1))..r$TARGET_REV:"
git -C "$CLONE" log --reverse --format="   %h %ad %an: %s" --date=short "$BASE".."$TARGET" | cut -c1-120
if [ "$DRY_RUN" -eq 1 ]; then
    say "dry run: stopping here"; exit 0
fi

# --- 2. upstream: rebase the new conversion commits onto our upstream ---------------------
# (our upstream may carry author corrections, so its SHAs differ from the clone's;
#  the trees are identical, which the rebase relies on)
OUR_BASE=$(find_rev . upstream "$RECORDED")
[ -n "$OUR_BASE" ] || { echo "cannot find r$RECORDED on upstream" >&2; exit 1; }
if [ "$(git rev-parse "$OUR_BASE^{tree}")" != "$(git -C "$CLONE" rev-parse "$BASE^{tree}")" ]; then
    echo "tree of r$RECORDED differs between upstream and the conversion clone; refusing" >&2; exit 1
fi
say "importing onto upstream"
git fetch -q "$CLONE" "$TARGET:refs/svn/import"
git branch -f svn-import refs/svn/import
git checkout -q svn-import
git rebase -q --onto upstream "$BASE" svn-import
NEW_TIP=$(git rev-parse HEAD)
if [ "$(git rev-parse "$NEW_TIP^{tree}")" != "$(git -C "$CLONE" rev-parse "$TARGET^{tree}")" ]; then
    echo "tree after import differs from the conversion clone; refusing" >&2; git checkout -q main; exit 1
fi
git branch -f upstream "$NEW_TIP"
git checkout -q main
git branch -D svn-import > /dev/null
git update-ref -d refs/svn/import
say "upstream is now r$TARGET_REV ($(git rev-parse --short upstream)); tree verified against the clone"

# --- 3. main: merge, regenerate, record ------------------------------------------------------
say "merging upstream into main"
if ! git merge --no-edit -m "Merge DiagHam SVN r$((RECORDED + 1))-r$TARGET_REV (upstream) into main" upstream; then
    cat >&2 <<EOF
merge conflicts: an upstream change touched a file this repository patched.
Resolve them (the fix may be obsolete if upstream fixed the same thing -- then drop
it and its row in docs/explanation/upstream-fixes.md), commit, and rerun this script
with the same arguments to finish the regeneration and record steps.
EOF
    exit 1
fi
say "regenerating the per-directory CMakeLists.txt"
python3 scripts_cmake/extract_autotools.py "$REPO_ROOT" > /dev/null
echo "$TARGET_REV" > scripts_cmake/upstream-revision.txt
if [ -n "$(git status --porcelain)" ]; then
    git add -A
    git commit -q -m "Generated CMakeLists.txt and recorded revision for DiagHam r$TARGET_REV"
fi
say "recorded r$TARGET_REV"

# --- 4. build and test ------------------------------------------------------------------------
if [ "$BUILD" -eq 1 ]; then
    say "configure, build, test (default preset)"
    cmake --preset default > /dev/null
    cmake --build --preset default --parallel "$(nproc)" > build/sync-build.log 2>&1 || { echo "build failed: see build/sync-build.log" >&2; exit 1; }
    ctest --preset default --parallel "$(nproc)" || { echo "tests failed: a golden may now be wrong, or upstream changed behaviour -- read the failures before touching any expected value" >&2; exit 1; }
    say "regenerate the program reference and the coverage page if programs or options changed:"
    echo "   python3 scripts_cmake/gen_program_reference.py build/default docs/reference/programs"
    echo "   python3 tests/coverage.py build/default --manifest tests/manifest.txt --write docs/reference/test-coverage.md"
fi
say "done: main now carries DiagHam r$TARGET_REV; review 'git log upstream -$NEW' and push when satisfied"
