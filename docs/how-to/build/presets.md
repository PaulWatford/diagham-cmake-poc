# CMake presets

Purpose: what each preset in `CMakePresets.json` turns on, and how presets replace the old configure lines.
Source: `CMakePresets.json` (2026-09-25); descriptions checked against the top-level `CMakeLists.txt`.

A preset fixes the build directory (`build/<preset>`), the install prefix
(`install/<preset>`) and a set of options. Use the same name for the three
steps:

```
cmake --preset <name>
cmake --build --preset <name> -j
ctest --preset <name>
```

| Preset | Modules | Optional libraries | Build type | Old configure equivalent |
|---|---|---|---|---|
| `default` | all | none | Release (`-O2`) | `--enable-fqhe --enable-fti --enable-spin --enable-quantumdots` |
| `core` | Base/src + src only | none | Release | no `--enable-*` module flag |
| `lapack` | all | LAPACK/BLAS | Release | `… --enable-lapack` |
| `full` | all | LAPACK, GSL, GMP, FFTW3, bzip2 | Release | `… --enable-lapack --enable-gsl --enable-gmp --enable-fftw --enable-bz2` (the wiki's LLVM line) |
| `hpc` | all | MPI, LAPACK, ScaLAPACK, GSL | Release | `… --enable-mpi --enable-lapack --enable-scalapack --enable-gsl` |
| `mkl` | all | Intel MKL (BLAS/LAPACK/FFTW3) | Release | `… --enable-intelmkl` (untested; see intel-mkl.md) |
| `debug` | all | none | Debug (`-g`, no optimisation) | closest: `--enable-debug` *plus* `-O0`; upstream's `--enable-debug` keeps `-O2` and adds `-g -fPIC`, which is `-DCMAKE_BUILD_TYPE=RelWithDebInfo` here |

Test presets exist with the same names, plus `physics` (only the physics
goldens, seconds). Any option can be added on the command line after the
preset, e.g. `cmake --preset default -DDIAGHAM_BUILD_SPIN=OFF` or
`-DDIAGHAM_USE_PROFILE=ON`. The complete option list is the
[configure-flag map](../../reference/configure-flag-map.md).

Build types: `Release` is `-O2` (upstream's default `CXXFLAGS=" -O2 -Wall "`,
no `NDEBUG`, so DiagHam's asserts stay live); `RelWithDebInfo` and
`RelWithChecks` are `-O2 -g`; `Debug` is CMake's `-g`.
