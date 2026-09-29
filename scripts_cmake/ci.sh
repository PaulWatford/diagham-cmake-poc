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

echo "== environment as built (build/$PRESET/environment.txt)"
scripts_cmake/environment.sh "build/$PRESET" | tee "build/$PRESET/environment.txt"

echo "== test ($PRESET)"
ctest --preset "$PRESET" --parallel "${CTEST_PARALLEL_LEVEL:-$(nproc)}"

if [ "$PRESET" = default ] && command -v python3 >/dev/null; then
    echo "== version consistent?"
    python3 scripts_cmake/check_version.py
    echo "== README program count matches the coverage page?"
    python3 - <<'PY'
import re, sys
from pathlib import Path
cov = Path("docs/reference/test-coverage.md").read_text().splitlines()
n_phys = n_union = n_any = 0
for line in cov:
    m = re.match(r"^\| (\S+) \| (\S+) \| (✓?) \| (✓?) \| (✓?) \| `", line)
    if not m:
        continue
    p, c, r = (bool(x) for x in m.groups()[2:])
    n_phys += p; n_union += (p or c); n_any += (p or c or r)
want = f"{n_union} of 603 programs have a physics or cross-check test today ({n_phys} with an independently known answer) and {n_any - n_union} more a regression test ({n_any} with any test beyond --help)."
if want not in Path("README.md").read_text():
    print("README.md program-count sentence is stale; it should read:", want, file=sys.stderr)
    sys.exit(1)
print("README program count is current")
PY
    echo "== coverage page up to date?"
    python3 tests/coverage.py "build/$PRESET" --write "build/$PRESET/test-coverage.md" > /dev/null
    if ! diff -q <(tail -n +5 "build/$PRESET/test-coverage.md") <(tail -n +5 docs/reference/test-coverage.md) > /dev/null; then
        echo "docs/reference/test-coverage.md is stale: regenerate with python3 tests/coverage.py build/default --write docs/reference/test-coverage.md" >&2
        exit 1
    fi
fi
