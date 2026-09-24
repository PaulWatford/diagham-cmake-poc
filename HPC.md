# Building DiagHam on an HPC cluster

The autotools recipe most clusters carry around looks something like

```bash
./configure --enable-fqhe --enable-fti --enable-lapack --enable-gsl \
            --enable-mpi --with-mpi-cxx=mpicxx --enable-scalapack \
            --with-blas-libs="-lopenblas" --with-lapack-libs="-lopenblas"
make -j8
```

The CMake equivalent is the `hpc` preset plus the site's modules.

## 1. Load a toolchain

Module names differ per site; the pieces needed are a C++11 compiler,
CMake >= 3.21 (presets), an MPI, a LAPACK/BLAS, ScaLAPACK, and GSL:

```bash
module load gcc/13 cmake/3.28 openmpi/4.1 openblas/0.3 scalapack/2.2 gsl/2.7
```

(An Intel stack works the same way with `intel-oneapi-compilers`,
`intel-oneapi-mpi` and `intel-oneapi-mkl`, then the `mkl` preset plus
`-DDIAGHAM_USE_MPI=ON`.)

## 2. Get the sources and overlay the CMake build

```bash
git clone https://github.com/guysoft/DiagHam.git        # or the canonical repository
git clone https://github.com/PaulWatford/diagham-cmake-poc.git
python3 diagham-cmake-poc/scripts_cmake/overlay.py DiagHam
cd DiagHam
```

## 3. Configure, build, test

```bash
cmake --preset hpc
cmake --build --preset hpc -j 16
ctest --preset hpc -j 16          # physics goldens + --help smoke test
```

The `hpc` preset turns on MPI, LAPACK, ScaLAPACK and GSL. Library
discovery follows the loaded modules (`CMAKE_PREFIX_PATH`, `MPI_HOME`,
`MKLROOT`, ...). When a module is not discoverable, pass the link line
directly, as `--with-*-libs` did:

```bash
cmake --preset hpc \
  -DDIAGHAM_LAPACK_LIBS="-L$OPENBLAS_ROOT/lib -lopenblas" \
  -DDIAGHAM_SCALAPACK_LIBS="-L$SCALAPACK_ROOT/lib -lscalapack"
```

Every configure flag and its CMake replacement is listed in
`CONFIGURE_FLAGS.md`.

## 4. Install

```bash
cmake --install build/hpc --prefix $HOME/opt/diagham
```

puts the programs in `bin/`, the static libraries in `lib/diagham/`, the
headers in `include/diagham/`, and a CMake package in
`lib/cmake/DiagHam/`, so your own code can use DiagHam as a library:

```cmake
find_package(DiagHam REQUIRED)            # with CMAKE_PREFIX_PATH=$HOME/opt/diagham
target_link_libraries(mycode PRIVATE DiagHam::DiagHam)
```

Install only the programs with `--component Runtime`, or only the
libraries, headers and CMake package with `--component Development`.

## What has been checked

The `hpc` preset was configured, built and tested (OpenMPI 4.1 and
ScaLAPACK 2.2 from Ubuntu 24.04, GCC 13) on 24/09: all 61 libraries and
561 programs build, and all 572 ctest tests pass. Doing so found an
upstream uninitialised-member bug in `SpinChainHamiltonianWithTranslations`
that only crashed in the MPI build (patch 15). It has not been
run on a real multi-node cluster, and no ctest yet launches a program
under `mpirun`, so the MPI code paths are compiled and linked but not
exercised by the test suite.
