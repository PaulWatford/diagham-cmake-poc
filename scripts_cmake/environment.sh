#!/usr/bin/env bash
# environment.sh: record what a build was made with, so a result can be reproduced.
#
# Usage: scripts_cmake/environment.sh [BUILD_DIR]      (default build/default)
# Prints, in order: the repository commit and DiagHam revision, the build's
# CMake options, the toolchain and library versions, and the output of
# TestDiagHamConf (what the compiled code believes it was built with).
# CI writes this to build/<preset>/environment.txt and keeps it as an artefact.
set -u
BUILD_DIR="${1:-build/default}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "# DiagHam build environment ($(date -u +%Y-%m-%dT%H:%M:%SZ))"
echo "commit:            $(git rev-parse HEAD 2>/dev/null || echo unknown)"
echo "diagham revision:  r$(tr -dc 0-9 < scripts_cmake/upstream-revision.txt)"
echo "build directory:   $BUILD_DIR"
echo
echo "## CMake options"
if [ -f "$BUILD_DIR/CMakeCache.txt" ]; then
    grep -E "^(DIAGHAM_[A-Z_]+|CMAKE_BUILD_TYPE|CMAKE_CXX_COMPILER|CMAKE_CXX_FLAGS|BLAS_LIBRARIES|LAPACK_LIBRARIES|MPI_CXX_COMPILER|GSL_LIBRARY|GMP_LIBRARY|FFTW3_LIBRARY|BZIP2_LIBRARIES):[A-Z]+=" "$BUILD_DIR/CMakeCache.txt" | sort
else
    echo "(no CMakeCache.txt in $BUILD_DIR)"
fi
echo
echo "## Toolchain"
echo "os:                $(. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME" || uname -a)"
echo "kernel:            $(uname -r) $(uname -m)"
echo "cmake:             $(cmake --version 2>/dev/null | head -1)"
CXX_BIN=$(grep -E "^CMAKE_CXX_COMPILER:" "$BUILD_DIR/CMakeCache.txt" 2>/dev/null | cut -d= -f2)
echo "c++:               $("${CXX_BIN:-c++}" --version 2>/dev/null | head -1)"
echo "generator:         $(grep -E "^CMAKE_GENERATOR:" "$BUILD_DIR/CMakeCache.txt" 2>/dev/null | cut -d= -f2)"
echo "python3:           $(python3 --version 2>/dev/null) $(python3 -c 'import numpy; print("numpy", numpy.__version__)' 2>/dev/null)"
echo "mpirun:            $(mpirun --version 2>/dev/null | head -1 || echo none)"
echo
echo "## Libraries (package manager)"
if command -v dpkg-query >/dev/null; then
    dpkg-query -W -f='${Package} ${Version}\n' liblapack-dev liblapack3 libblas-dev libblas3 libopenblas-dev libgsl-dev libgmp-dev libfftw3-dev libbz2-dev libopenmpi-dev libscalapack-openmpi-dev 2>/dev/null
elif command -v rpm >/dev/null; then
    rpm -q lapack-devel blas-devel gsl-devel gmp-devel fftw-devel bzip2-devel openmpi-devel 2>/dev/null
else
    echo "(no dpkg or rpm; list your modules here)"
fi
command -v module >/dev/null 2>&1 && { echo; echo "## Loaded modules"; module list 2>&1; }
echo
echo "## TestDiagHamConf (what the compiled code reports)"
conf=$(find "$BUILD_DIR" -type f -name TestDiagHamConf -perm -u+x 2>/dev/null | head -1)
if [ -n "$conf" ]; then "$conf" 2>&1; else echo "(TestDiagHamConf not built)"; fi
