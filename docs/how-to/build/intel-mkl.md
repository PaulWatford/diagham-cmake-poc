# Build with the Intel compilers and MKL

Purpose: use Intel MKL for BLAS/LAPACK (and FFTW3), optionally with Intel MPI.
Source: upstream wiki "Install", section "Using Intel compiler and MKL", and the wiki "MPI" page's Intel-MPI notes (as of 2026-09-24); `CONFIGURE_FLAGS.md`. Changed: the configure line is mapped to the `mkl` preset. **Not yet verified on a machine with MKL** — read the caveat.

DiagHam is developed with GCC as the default compiler but works with the
Intel compilers. The wiki's recipe

```
CXX="icc" CC="icc" FC="ifort" ../configure --enable-fti --enable-fqhe --enable-mpi --enable-intelmkl --with-intelmkl-libdir=$MKL_HOME/lib/intel64
```

maps to

```
CC=icx CXX=icpx cmake --preset mkl -DDIAGHAM_USE_MPI=ON -DMPI_CXX_COMPILER=mpiicpx
cmake --build --preset mkl -j
ctest --preset mkl
```

- `--enable-intelmkl` ⇒ `DIAGHAM_USE_LAPACK=ON DIAGHAM_USE_MKL=ON`. As
  upstream, MKL *is* the LAPACK provider (`HAVE_LAPACK` is defined); you do
  not also point CMake at a separate LAPACK. `DIAGHAM_USE_FFTW=ON` (in the
  preset) takes FFTW3 from MKL's FFTW interface.
- `--with-intelmkl-libdir=DIR` ⇒ the `MKL_ROOT` cache variable or the
  `MKLROOT` environment variable (set by Intel's `setvars.sh` / the cluster
  module). `cmake/FindMKL.cmake` (from Gunnar Möller's BDMC_UFL project)
  looks there.
- `icc`/`icpc` are the classic compilers, discontinued by Intel in 2023;
  `icx`/`icpx` are the current ones and take the same role.

## Intel MPI

With Intel MPI the wiki warns of

```
mpicxx.h(95): error: #error directive: "SEEK_SET is #defined but must not be for the C++ binding of MPI. Include mpi.h before stdio.h"
```

and prescribes `CPPFLAGS="-DMPICH_IGNORE_CXX_SEEK"`. The CMake equivalent is
`-DCMAKE_CXX_FLAGS=-DMPICH_IGNORE_CXX_SEEK` on the configure line. The
wiki's full Intel example (MKL + Intel MPI + ScaLAPACK, ILP64 interface):

```
CPPFLAGS="-DMPICH_IGNORE_CXX_SEEK" ../configure --enable-fqhe --enable-fti --enable-spin --enable-intelmkl --enable-mpi --enable-debug CC="icc" CXX="icpc" --enable-lapack --with-intelmkl-libdir=.../mkl/lib/intel64 --with-blas-libs="-lmkl_blas95_ilp64" --with-lapack-libs="-lmkl_lapack95_ilp64" --with-mpi-cxx="mpiicpc" --enable-scalapack --with-scalapack-libs="-lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64"
```

## Caveat: LP64 versus ILP64 — an open decision

Upstream's `configure.ac` links the **ILP64** MKL interface
(`-lmkl_intel_ilp64`, `-DMKL_ILP64`, 64-bit integers in every BLAS/LAPACK
call). `cmake/FindMKL.cmake` supports **LP64** only (32-bit integers, the
interface that matches `int` arguments), and adds no `-DMKL_ILP64`. DiagHam
passes `int` dimensions to LAPACK through `FORTRAN_NAME`, so LP64 is the
consistent choice and upstream's ILP64 line is suspect — but nobody has yet
compiled either variant against a real MKL under this build. Until someone
does (the Kent cluster is the obvious place), treat the `mkl` preset as
untested and check `TestDiagHamConf` and the physics goldens on first use.
If the link line needs to be given by hand, use
`-DDIAGHAM_LAPACK_LIBS="<link line>"`, which replaces `FindMKL` entirely.
