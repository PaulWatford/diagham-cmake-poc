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
| `core-programs` | Base, src, FQHE, FTI, Spin; **only the core programs** | LAPACK (OpenBLAS where installed), GSL, GMP | Release | none: the recommended build of the core programs |
| `core-programs-portable` | as `core-programs` | none | Release | none: the core programs where nothing can be installed |
| `debug` | all | none | Debug (`-g`, no optimisation) | closest: `--enable-debug` *plus* `-O0`; upstream's `--enable-debug` keeps `-O2` and adds `-g -fPIC`, which is `-DCMAKE_BUILD_TYPE=RelWithDebInfo` here |

Test presets exist with the same names, plus `physics` (only the physics
goldens, seconds). Any option can be added on the command line after the
preset, e.g. `cmake --preset default -DDIAGHAM_BUILD_SPIN=OFF` or
`-DDIAGHAM_USE_PROFILE=ON`. The complete option list is the
[configure-flag map](../../reference/configure-flag-map.md).

## Build only the programs you need

Two options choose which programs are built. The libraries of every enabled module are always built.

| Option | Values | Effect |
|---|---|---|
| `DIAGHAM_PROGRAM_SET` | `all` (default) | every program, as before |
| | `core` | the core programs listed in `cmake/core-programs.txt` |
| | `custom` | only the programs named in `DIAGHAM_PROGRAMS` |
| `DIAGHAM_PROGRAMS` | upstream binary names, separated by `;` or `,` | the whole set with `custom`; added to the core set with `core` |

Examples:

```
cmake --preset core-programs                       # the core programs, with LAPACK, GSL, GMP
cmake --preset core-programs -DDIAGHAM_PROGRAMS="FCIDiceLatticeModel"     # core plus one more
cmake --preset default -B build/mine -DDIAGHAM_PROGRAM_SET=custom -DDIAGHAM_PROGRAMS="HubbardSquareLatticeModel,FQHESphereFermionsTwoBodyGeneric"
```

A name that matches no program in the enabled modules stops the configure step and names it. The test suite
registers only the tests of the programs that are built.
- **`core` set:** checked against its own locked list, `tests/manifest-core-programs.txt`.
- **Any other selection:** every test it registers must be in the full list, `tests/manifest.txt`.

The core list is provisional: the 42 programs of the candidate list of 19 October 2026. The team settles the
final list by how developed each program's wiki page is. After editing `cmake/core-programs.txt`, regenerate the
locked list from a `core-programs-portable` build:

```
cmake --preset core-programs-portable && cmake --build --preset core-programs-portable -j
LC_ALL=C python3 tests/coverage.py build/core-programs-portable --manifest tests/manifest-core-programs.txt
```

**What it saves.** Measured on 25 October 2026 on a 32-core machine, fresh build directories:

| Build | Time | Programs | Tests |
|---|---|---|---|
| `core-programs-portable` | 50 s | 42 | 266 |
| `core-programs` | 48 s | 42 | 272 |
| custom, two programs | 45 s | 2 | 50 |
| `default`, every program | 60 s | 603 | 1013 |

The shared libraries take most of the build time, so the gain is mainly a smaller install and a shorter test
run.

Build types: `Release` is `-O2` (upstream's default `CXXFLAGS=" -O2 -Wall "`,
no `NDEBUG`, so DiagHam's asserts stay live); `RelWithDebInfo` and
`RelWithChecks` are `-O2 -g`; `Debug` is CMake's `-g`.
