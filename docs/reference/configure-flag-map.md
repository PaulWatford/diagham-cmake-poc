# configure.ac → CMake flag map

Source: `CONFIGURE_FLAGS.md` (2026-09-24), moved to the reference section 2026-09-25 with two corrections (MKL interface note; --enable-debug mapping).

Every `AC_ARG_ENABLE` / `AC_ARG_WITH` in upstream `configure.ac` (34
distinct flags, pinned upstream revision `ed78a30`), and what to pass to
CMake instead. The "Tested" column says where the option was built and the
full ctest suite passed with it: in CI (a GitHub Actions job runs that
preset), or by hand on Ubuntu 24.04 / GCC 13 (24/09). Options with nothing
in the column are wired but untested, because their libraries have no
distribution package to test against (MPACK) or no hardware/licence here
(MKL).

## Packages (modules)

| configure | CMake | Default | Notes |
|---|---|---|---|
| `--enable-fqhe` | `DIAGHAM_BUILD_FQHE` | ON | |
| `--enable-fti` | `DIAGHAM_BUILD_FTI` | ON | needs FQHE (enforced at configure time) |
| `--enable-spin` | `DIAGHAM_BUILD_SPIN` | ON | new in this PoC |
| `--enable-quantumdots` | `DIAGHAM_BUILD_QUANTUMDOTS` | ON | new in this PoC |
| `--enable-anyons` | — | — | `Anyons/` has no sources in the upstream tree and its `AC_CONFIG_FILES` line is commented out, so the autotools flag does nothing either |

autotools builds none of the modules unless asked; the CMake default builds
all of them, which is what "first-class" means here. `cmake --preset core`
gives the core-only build.

## Optional libraries

| configure | CMake | Tested | Notes |
|---|---|---|---|
| `--enable-lapack` | `DIAGHAM_USE_LAPACK=ON` | CI (`full`) | `find_package(LAPACK)` |
| `--with-lapack-libs="..."`, `--with-blas-libs="..."` | `DIAGHAM_LAPACK_LIBS="..."` | by hand (`core` + `-llapack -lblas`) | replaces `find_package(LAPACK)`; give the whole LAPACK+BLAS link line |
| `--with-lapack-libdir=DIR`, `--with-blas-libdir=DIR` | `CMAKE_PREFIX_PATH=DIR/..` or `CMAKE_LIBRARY_PATH=DIR` | | standard CMake search paths |
| `--enable-lapack-only` | `DIAGHAM_USE_LAPACK_ONLY=ON` | | requires LAPACK |
| `--enable-intelmkl` | `DIAGHAM_USE_LAPACK=ON DIAGHAM_USE_MKL=ON` (preset `mkl`) | | `cmake/FindMKL.cmake`, from Gunnar Möller's BDMC_UFL. **Interface differs from upstream**: configure links the ILP64 MKL interface with `-DMKL_ILP64`; FindMKL is LP64-only. Untested against a real MKL — see `docs/how-to/build/intel-mkl.md` |
| `--with-intelmkl-libdir=DIR` | `MKLROOT` environment variable | | what `FindMKL.cmake` reads |
| `--enable-gsl` | `DIAGHAM_USE_GSL=ON` | CI (`full`) | `find_package(GSL)` |
| `--with-gsl-libs="..."` | `DIAGHAM_GSL_LIBS="..."` | | |
| `--enable-gmp` | `DIAGHAM_USE_GMP=ON` | CI (`full`) | links `gmpxx` and `gmp`, as configure does |
| `--with-gmp-libs="..."` | `DIAGHAM_GMP_LIBS="..."` | | |
| `--enable-mpack` | `DIAGHAM_USE_MPACK=ON` | | requires LAPACK; no system package exists to test against |
| `--with-mpack-libs="..."` | `DIAGHAM_MPACK_LIBS="..."` | | default `mlapack_gmp;mblas_gmp`, as configure |
| `--enable-fftw` | `DIAGHAM_USE_FFTW=ON` | CI (`full`) | FFTW3; taken from MKL when `DIAGHAM_USE_MKL=ON` |
| `--enable-bz2` | `DIAGHAM_USE_BZ2=ON` | CI (`full`) | `find_package(BZip2)` |
| `--enable-mpi` | `DIAGHAM_USE_MPI=ON` | CI (`hpc`, OpenMPI 4.1) | `find_package(MPI COMPONENTS CXX)`; no test runs under `mpirun` yet |
| `--with-mpi-cxx=mpicxx` | `CMAKE_CXX_COMPILER=mpicxx` | | or leave the compiler alone: `find_package(MPI)` adds the flags |
| `--with-mpi-libs`, `--with-mpi-incdir`, `--with-mpi-libdir` | `MPI_HOME=...` or `MPI_CXX_COMPILER=...` | | FindMPI's own hints |
| `--enable-scalapack` | `DIAGHAM_USE_SCALAPACK=ON` | CI (`hpc`, ScaLAPACK 2.2) | requires MPI and LAPACK (same rule as configure); tries `find_package(scalapack CONFIG)` then `find_library` |
| `--with-scalapack-libs="..."` | `DIAGHAM_SCALAPACK_LIBS="..."` | | |
| `--with-scalapack-libdir=DIR` | `CMAKE_LIBRARY_PATH=DIR` | | |
| `--with-fortran2c-libs="..."` | `DIAGHAM_FORTRAN2C_LIBS="..."` | | appended to every program's link line |

## Build modes

| configure | CMake | Notes |
|---|---|---|
| `--enable-debug` | `CMAKE_BUILD_TYPE=RelWithDebInfo` (`-O2 -g`) | configure keeps the default `-O2 -Wall` and appends `-g -fPIC`, i.e. an optimised build with symbols; `Debug` (preset `debug`) is `-g` without optimisation, which is *not* what `--enable-debug` did |
| `--enable-profile` | `DIAGHAM_USE_PROFILE=ON` | `-g -pg` on compile and link |
| `--enable-m64` | `CMAKE_CXX_FLAGS=-m64` | obsolete on every 64-bit toolchain; not given its own option |
| `--with-cmd-log=FILE` | `DIAGHAM_COMMAND_LOG=FILE` | defines `HAVE_GLOBAL_COMMAND_LOG` and `GLOBAL_COMMAND_LOG` |

## Macros that are not options

`src/config.h` defines `__SMP__` and `__DEBUG__` unconditionally, whatever
configure (or CMake) says. `DIAGHAM_USE_SMP=OFF` therefore only drops the
pthread link, and `DIAGHAM_USE_DEBUG=OFF` changes nothing. Fixing that
needs an upstream `config.h` change (gate both on `HAVE_*` like the other
features); see the "FLAGGED, NOT FIXED" note in `CMakeLists.txt`.
