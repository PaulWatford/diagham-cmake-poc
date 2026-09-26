# Build on an HPC cluster (MPI + ScaLAPACK)

Purpose: build DiagHam on a cluster with a module system for a distributed-memory (MPI) run.
Source: `HPC.md` (2026-09-24) and the upstream wiki "MPI" and "Scalapack" pages (as of 2026-09-24). Changed: the illustrative configure line has been replaced by the wiki's documented ones; the `hpc` preset is the equivalent. Running the programs under MPI is a separate guide: [../run-mpi.md](../run-mpi.md).

## 1. Load a toolchain

Module names differ per site; the pieces needed are a C++11 compiler,
CMake ≥ 3.21, an MPI, a BLAS/LAPACK, ScaLAPACK, and GSL, for example

```
module load gcc cmake openmpi openblas scalapack gsl
```

For an Intel stack (`intel-oneapi-compilers`, `-mpi`, `-mkl`) use the `mkl`
preset with `-DDIAGHAM_USE_MPI=ON` instead — see [intel-mkl.md](intel-mkl.md),
including the `-DMPICH_IGNORE_CXX_SEEK` flag that Intel MPI needs.

## 2. Configure, build, test

```
git clone <this repository> DiagHam && cd DiagHam
cmake --preset hpc                  # MPI, LAPACK, ScaLAPACK, GSL
cmake --build --preset hpc -j 16
ctest --preset hpc -j 16            # physics goldens + --help smoke test (serial)
```

`--enable-scalapack` requires MPI and LAPACK, exactly as in `configure.ac`;
the preset turns all three on. Library discovery follows the loaded modules
(`CMAKE_PREFIX_PATH`, `MPI_HOME`, `MKLROOT`, …). When a library is not
discoverable, pass the link line directly, as `--with-*-libs` did:

```
cmake --preset hpc \
  -DDIAGHAM_LAPACK_LIBS="-L$OPENBLAS_ROOT/lib -lopenblas" \
  -DDIAGHAM_SCALAPACK_LIBS="-L$SCALAPACK_ROOT/lib -lscalapack"
```

The wiki records one ScaLAPACK link problem with Intel MPI
(`undefined reference to zhemv_ / zher2_` from `libscalapack.a`) whose fix
was `--with-scalapack-libs="-lscalapack -llapack -lblas"`; the equivalent is
`-DDIAGHAM_SCALAPACK_LIBS="-lscalapack -llapack -lblas"`.

To use a specific MPI wrapper (`--with-mpi-cxx=mpicxx`), set
`-DMPI_CXX_COMPILER=mpicxx`; or simply `CXX=mpicxx cmake --preset hpc`.

## 3. Install for the group

```
cmake --install build/hpc --prefix /path/to/shared/diagham
```

`bin/` holds the programs, `lib/diagham/` the static libraries,
`include/diagham/` the headers and `lib/cmake/DiagHam/` a CMake package so
other codes can `find_package(DiagHam)`; see
[../install-and-use-as-library.md](../install-and-use-as-library.md).
Install only the programs with `--component Runtime`.

## What has been checked

The `hpc` preset builds and passes the full ctest suite in CI (OpenMPI 4.1
and ScaLAPACK 2.2 from Ubuntu 24.04, GCC 13). It has not been run on a real
multi-node cluster under this build, and no ctest launches a program under
`mpirun`: the MPI and ScaLAPACK code paths are compiled and linked but not
exercised by the test suite. The Kent site recipe is a draft awaiting the site's values:
[../../drafts/kent-cluster-build.md](../../drafts/kent-cluster-build.md).
