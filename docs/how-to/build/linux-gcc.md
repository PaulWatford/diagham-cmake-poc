# Build on Linux with GCC

Purpose: get a working DiagHam build on a Linux workstation with GCC.
Source: upstream wiki "Install" (Requirements, First installation, Checking your installation; as of 2026-09-24) and `README.md`. Changed: `svn checkout` + `bootstrap.sh` + `configure` replaced by `git clone` + CMake presets; the module flags (`--enable-fqhe`, `--enable-spin`, `--enable-fti`) are on by default here; verification now uses `ctest`.

## Requirements

- GCC (any C++11-capable version; CI uses 13, the migration was verified with 15)
- CMake ≥ 3.21 and a generator (`make` or `ninja`)
- pthreads (always used: DiagHam's "SMP" mode is pthreads)
- Python 3, with numpy if you want the Python cross-check test

Optional, each enabled by a preset or an option:

- BLAS/LAPACK — not mandatory, but they provide better performance and enable
  the `LapackDiagonalize` paths (`liblapack-dev libopenblas-dev` on
  Debian/Ubuntu)
- GSL, GMP (`libgsl-dev libgmp-dev`), FFTW3 (`libfftw3-dev`), bzip2
  (`libbz2-dev`)
- MPI and ScaLAPACK — only for a cluster, see [hpc-cluster.md](hpc-cluster.md)

## First build

```
git clone <this repository> DiagHam
cd DiagHam
cmake --preset default              # every module, SMP, no optional libraries
cmake --build --preset default -j   # a few minutes on a workstation
ctest --preset default              # physics goldens, --help on every program, install test
```

With system LAPACK/BLAS installed, use `--preset lapack` instead; with the
full set of serial libraries, `--preset full` (the equivalent of the old
`--enable-fqhe --enable-lapack --enable-gsl --enable-bz2 --enable-gmp
--enable-spin --enable-fti` line). All presets: [presets.md](presets.md).

Unlike the autotools build, every module (core, FQHE, FTI, Spin,
QuantumDots) is built by default. To build a subset, set the module options,
e.g. `cmake --preset default -DDIAGHAM_BUILD_SPIN=OFF`; FTI requires FQHE.

## Checking the installation

Everything went fine? Try a program's help:

```
build/default/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator --help
```

Every program answers `--help` (or `-h`) with its options; this is the same
text the [program reference](../../reference/programs/README.md) is
generated from. `ctest --preset default -L physics` runs the physics
goldens alone (seconds). To see which compile-time features the build
enabled, run `build/default/src/Programs/TestDiagHamConf`.

Keep the source, the build and your data separate: build trees live under
`build/<preset>`; run programs from a separate working directory, because
they write their output files into the current directory.

## Updating the code

`git pull`, then rebuild with the same preset; CMake re-runs itself when the
`CMakeLists.txt` files changed. If upstream added a program or library
(`Makefile.am` changed on the `upstream` branch), regenerate the
per-directory files first: `python3 scripts_cmake/extract_autotools.py .`.

## Installing

`cmake --install build/default --prefix $HOME/opt/diagham` installs the
programs to `bin/`, the static libraries to `lib/diagham/`, the headers to
`include/diagham/` and a CMake package to `lib/cmake/DiagHam/`; see
[../install-and-use-as-library.md](../install-and-use-as-library.md).
