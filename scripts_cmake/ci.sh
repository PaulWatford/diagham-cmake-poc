#!/usr/bin/env bash
# ci.sh: the whole CI pipeline, host-independent.
#
# Configures, builds and tests this repository in place with one preset.
# GitHub Actions (.github/workflows/ci.yml) and GitLab CI (.gitlab-ci.yml)
# both just call this script, so moving hosts changes nothing about what
# is tested. The DiagHam sources are part of the tree (branch `main` carries
# the r4493 sources merged from `upstream`); nothing is fetched.
#
# Usage: scripts_cmake/ci.sh [PRESET]            (default preset: "default")
#
# Environment honoured:
#   CMAKE_GENERATOR              e.g. Ninja
#   CMAKE_BUILD_PARALLEL_LEVEL   build jobs (default: nproc)
#   CTEST_PARALLEL_LEVEL         test jobs (default: nproc)
#   DIAGHAM_CI_EXTRA_CMAKE_ARGS  extra -D arguments for the configure step
set -euo pipefail

PRESET="${1:-default}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "== environment"
echo "commit:   $(git rev-parse --short HEAD 2>/dev/null || echo n/a)"
echo "preset:   $PRESET"
echo "cmake:    $(cmake --version | head -1)"
echo "compiler: ${CXX:-c++} -> $(${CXX:-c++} --version | head -1)"

echo "== configure ($PRESET)"
# shellcheck disable=SC2086
cmake --preset "$PRESET" ${DIAGHAM_CI_EXTRA_CMAKE_ARGS:-}

echo "== build ($PRESET)"
cmake --build --preset "$PRESET" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-$(nproc)}"

echo "== configuration as built"
conf="$(find "build/$PRESET" -type f -name TestDiagHamConf -perm -u+x | head -1 || true)"
if [ -n "$conf" ]; then "$conf" || true; fi

echo "== test ($PRESET)"
ctest --preset "$PRESET" --parallel "${CTEST_PARALLEL_LEVEL:-$(nproc)}"
