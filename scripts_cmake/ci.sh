#!/usr/bin/env bash
# ci.sh: the whole CI pipeline, host-independent.
#
# Fetches a pinned upstream DiagHam revision, overlays this repository's
# CMake build onto it, then configures, builds and runs ctest with one
# preset. GitHub Actions (.github/workflows/ci.yml) and GitLab CI
# (.gitlab-ci.yml, for the Kent GitLab move) both just call this script,
# so moving hosts changes nothing about what is tested.
#
# Usage: scripts_cmake/ci.sh [PRESET]            (default preset: "default")
#
# Environment:
#   DIAGHAM_UPSTREAM_URL  upstream git repository (default: the guysoft mirror)
#   DIAGHAM_UPSTREAM_REF  commit to build (default: pinned below; bump it
#                         deliberately, in its own change, so an upstream
#                         change can never silently alter what CI tests)
#   DIAGHAM_CI_WORK       scratch directory (default: ./_ci)
#   CMAKE_GENERATOR       honoured by cmake (e.g. Ninja)
#   CTEST_PARALLEL_LEVEL  honoured by ctest
set -euo pipefail

PRESET="${1:-default}"
UPSTREAM_URL="${DIAGHAM_UPSTREAM_URL:-https://github.com/guysoft/DiagHam.git}"
UPSTREAM_REF="${DIAGHAM_UPSTREAM_REF:-ed78a30e2ec1ae013c6115fc6c9700854acaf0a5}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${DIAGHAM_CI_WORK:-$REPO_ROOT/_ci}"
TREE="$WORK/DiagHam"

echo "== upstream: $UPSTREAM_URL @ $UPSTREAM_REF"
if [ ! -d "$TREE/.git" ]; then
    mkdir -p "$TREE"
    git -C "$TREE" init -q
    git -C "$TREE" remote add origin "$UPSTREAM_URL"
fi
git -C "$TREE" fetch -q --depth 1 origin "$UPSTREAM_REF"
git -C "$TREE" checkout -q --force FETCH_HEAD
# Start from pristine sources every time: the patch series is applied at
# configure time and must never be stacked on a previous run's result.
git -C "$TREE" clean -q -fdx -e build/

echo "== overlay"
python3 "$REPO_ROOT/scripts_cmake/overlay.py" "$TREE"

cd "$TREE"
# A cached build directory from an earlier run is fine to reuse (it makes
# rebuilds incremental), but its patch sentinel went with `git clean`.
echo "== configure ($PRESET)"
cmake --preset "$PRESET"
echo "== build ($PRESET)"
cmake --build --preset "$PRESET" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-$(nproc)}"
echo "== test ($PRESET)"
ctest --preset "$PRESET" --parallel "${CTEST_PARALLEL_LEVEL:-$(nproc)}"
