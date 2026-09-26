# Build system reference

Purpose: every option, file and mechanism of the CMake build, in one place; the numbers other pages quote come from here.
Source: `MODULE_MAP.md` (2026-09-18, task-oriented map), the top-level `CMakeLists.txt`, `cmake/DiagHamHelpers.cmake`, `cmake/DiagHamInstall.cmake` and `scripts_cmake/extract_autotools.py` at the 2026-09-26 state of `main`. Changed: MODULE_MAP's entries for configure-time patching and the "79-check suite" removed (patches are commits; `verify_build.sh` is described as what it is); options listed from the file, not from memory.

## Numbers (DiagHam r4493, `default` preset)

61 static libraries · 603 programs (599 with an option parser) · 88 generated `CMakeLists.txt` · 34 upstream configure flags mapped · 613 ctest tests (13 physics goldens). Five sources are deliberately not built (see "Exclusions and extras").

## Files

| File | Role |
|---|---|
| `CMakeLists.txt` | top-level: options, external libraries, the config header, module `add_subdirectory` calls in upstream's `SUBDIRS` order (`Base src QuantumDots Spin FQHE FTI`), install and tests |
| `CMakePresets.json` | configure/build/test presets — see [../how-to/build/presets.md](../how-to/build/presets.md) |
| `cmake/DiagHamHelpers.cmake` | the functions every generated file calls (below) |
| `cmake/DiagHamInstall.cmake`, `cmake/DiagHamConfig.cmake.in` | install rules and the `find_package(DiagHam)` package |
| `cmake/config_ac.h.in` | replaces autoconf's `config_ac.h`: one `#cmakedefine HAVE_<X>` per `AC_DEFINE` in `configure.ac`; upstream's `src/config.h` includes it and derives the `__X__` macros itself |
| `cmake/CompilerDefaults.cmake` | build-type validation and per-compiler warning flags (adapted from Gunnar Möller's BDMC_UFL) |
| `cmake/KentDefaults.cmake` | site-defaults scaffold for Kent, `DIAGHAM_KENT_DEFAULTS=ON` (placeholders; see [../how-to/build/kent.md](../how-to/build/kent.md)) |
| `cmake/Find{MKL,NAG,GSL,GMP,MPFR,MPC,FFTW3,DLR}.cmake`, `gsl_external_cblas.patch` | Find modules from BDMC_UFL; MKL, GSL, GMP and FFTW3 are used, the rest staged |
| `cmake/verify_build.sh` | CMake-vs-autotools build parity (library existence, `nm` symbol counts, one program's output); needs an autotools build of the same tree; 79 checks pass in the default configuration; not part of ctest |
| `scripts_cmake/extract_autotools.py` | generates the per-directory `CMakeLists.txt` from upstream `Makefile.am` |
| `scripts_cmake/ci.sh` | the CI pipeline: pinned upstream, build, ctest |
| `scripts_cmake/gen_program_reference.py`, `attach_wiki_manuals.py` | generate `docs/reference/programs/` |
| `tests/` | the ctest suite — see `docs/reference/tests.md` |
| `<module>/src/**/CMakeLists.txt` | generated; do not edit |

## Options

Module and build options (`option(...)` in `CMakeLists.txt`; defaults shown):

| Option | Default | Effect |
|---|---|---|
| `DIAGHAM_BUILD_FQHE` | ON | FQHE module (sphere, disk, torus, cylinder, lattice) |
| `DIAGHAM_BUILD_FTI` | ON | FTI module (Hubbard, FCI, FTI); requires FQHE (fatal error otherwise) |
| `DIAGHAM_BUILD_SPIN` | ON | Spin module |
| `DIAGHAM_BUILD_QUANTUMDOTS` | ON | QuantumDots module |
| `DIAGHAM_BUILD_PROGRAMS` | ON | build the executables, not only the libraries |
| `DIAGHAM_BUILD_TESTS` | ON | register the ctest suite (needs programs) |
| `DIAGHAM_USE_SMP` | ON | link pthreads (DiagHam's "SMP" mode). The `__SMP__` macro itself is defined unconditionally by upstream `config.h` |
| `DIAGHAM_USE_DEBUG` | ON | no effect on the code: `__DEBUG__` is defined unconditionally by upstream `config.h`; kept for compatibility |
| `DIAGHAM_KENT_DEFAULTS` | OFF | include `cmake/KentDefaults.cmake` |

Optional libraries (each defines `HAVE_<X>` in the config header; the old flag in brackets):

| Option | Default | Notes |
|---|---|---|
| `DIAGHAM_USE_LAPACK` (`--enable-lapack`) | OFF | `find_package(LAPACK)`, or `DIAGHAM_USE_MKL`, or `DIAGHAM_LAPACK_LIBS` |
| `DIAGHAM_USE_MKL` (`--enable-intelmkl`) | OFF | requires `DIAGHAM_USE_LAPACK`; `cmake/FindMKL.cmake`; also provides FFTW3 — untested, see [../how-to/build/intel-mkl.md](../how-to/build/intel-mkl.md) |
| `DIAGHAM_USE_LAPACK_ONLY` (`--enable-lapack-only`) | OFF | requires LAPACK |
| `DIAGHAM_USE_MPI` (`--enable-mpi`) | OFF | `find_package(MPI COMPONENTS CXX)` |
| `DIAGHAM_USE_SCALAPACK` (`--enable-scalapack`) | OFF | requires MPI and LAPACK; config package, then `find_library`, then `DIAGHAM_SCALAPACK_LIBS` |
| `DIAGHAM_USE_GSL` (`--enable-gsl`) | OFF | |
| `DIAGHAM_USE_GMP` (`--enable-gmp`) | OFF | links `gmpxx` and `gmp` |
| `DIAGHAM_USE_MPACK` (`--enable-mpack`) | OFF | requires LAPACK; link line from `DIAGHAM_MPACK_LIBS` (default `mlapack_gmp;mblas_gmp`) |
| `DIAGHAM_USE_FFTW` (`--enable-fftw`) | OFF | from MKL when `DIAGHAM_USE_MKL`, else `find_package(FFTW3)` |
| `DIAGHAM_USE_BZ2` (`--enable-bz2`) | OFF | `find_package(BZip2)` |
| `DIAGHAM_USE_PROFILE` (`--enable-profile`) | OFF | `-g -pg` on compile and link |
| `DIAGHAM_COMMAND_LOG` (`--with-cmd-log=FILE`) | "" | path; defines `HAVE_GLOBAL_COMMAND_LOG` and `GLOBAL_COMMAND_LOG` |

Link-line overrides (the `--with-<pkg>-libs` equivalents; a `;`-list or a space-separated string): `DIAGHAM_LAPACK_LIBS`, `DIAGHAM_GSL_LIBS`, `DIAGHAM_GMP_LIBS`, `DIAGHAM_MPACK_LIBS`, `DIAGHAM_SCALAPACK_LIBS`, `DIAGHAM_FORTRAN2C_LIBS`. The complete flag-by-flag map is [configure-flag-map.md](configure-flag-map.md).

Every external library is attached to one interface target,
`diagham_external_deps`; programs link it, and the exported
`DiagHam::DiagHam` forwards it to downstream projects.

## Compile definitions and flags

- `HAVE_CONFIG_H` and `MACHINE_PRECISION=1e-14` are defined for every target;
  everything else comes from `config_ac.h` → upstream `src/config.h`. The
  build deliberately does **not** inject `__64_BITS__`, `__SMP__` or
  `__DEBUG__` (config.h derives or defines them).
- C++11 with GNU extensions (`-std=gnu++11`), as upstream's g++ default.
- Release is `-O2`, RelWithDebInfo and RelWithChecks `-O2 -g`, Debug `-g`;
  no `NDEBUG` (autotools' `CXXFLAGS=" -O2 -Wall "`). Warning flags per
  compiler from `CompilerDefaults.cmake`; `-Wno-deprecated -Wno-write-strings
  -Wno-narrowing -Wno-format` silence 2003-era idioms.
- In-source builds are refused.

## Helper functions (what a generated file calls)

| Function | Does |
|---|---|
| `diagham_add_library(<Name> SOURCES …)` | `add_library(<Name> STATIC …)`, output `lib<Name>.a` as autotools; registers the library so programs can link "everything"; installs it (`Development`) |
| `diagham_add_programs(PREFIX <p> PROGRAMS <name> <source.cc> …)` | one executable per pair, target `<p>_<name>`, binary named `<name>`; skips `DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS`; installs to `bin/` (`Runtime`) |
| `diagham_add_program(<target> <source>)` | the single-program primitive: links every registered library inside `--start-group`/`--end-group` (GNU ld; the libraries are mutually recursive) plus `diagham_external_deps` |
| `diagham_add_programs_in_directory()` | legacy glob form; no generated file uses it since the single-policy generator |

Target names are `<directory>_<Program>` (`FQHEOnSphere_FQHESphereJackGenerator`, `HubbardModels_HubbardSquareLatticeModel`); for a directory literally named `Programs` the module name is used (`Spin_Cobalt`, `QuantumDots_QuantumDot`), and `src/Programs` keeps `Programs_` (`Programs_TestDiagHamConf`).

## The generator

`python3 scripts_cmake/extract_autotools.py .` walks every `Makefile.am`
under `Base/src src FQHE/src FTI/src Spin/src QuantumDots/src` and writes
one `CMakeLists.txt` per directory:

- `noinst_LIBRARIES` + `lib<X>_a_SOURCES` → `diagham_add_library` (comment
  lines stripped, duplicates collapsed);
- `SUBDIRS` → `add_subdirectory` in upstream's order;
- `bin_PROGRAMS` + `<prog>_SOURCES` → `diagham_add_programs` — the one rule
  for every module; a directory is a programs directory iff it has
  `bin_PROGRAMS`, whatever its name (so `QuantumDots/src/Tools/Analysis`
  is included);
- a documented allow-list adds sources upstream never listed
  (`EXTRA_PROGRAMS`) and recovers three orphaned library sources needed by
  the Dice-lattice program (`ORPHAN_RECOVERY_ALLOWLIST`).

Re-run it after any upstream `Makefile.am` change and commit the result;
it is path-separator-safe and runs on Windows.

## Exclusions and extras

Not built, on purpose (the only sources in scope that are not):

| Source | Why |
|---|---|
| `FQHE/src/Programs/FQHEOnTorus/QHEFermionsTorusWithSpin.cc` | dead duplicate of `FQHETorusFermionsWithSpin` that never compiled (`../explanation/deferred-code.md`); `DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS` |
| `Spin/src/Programs/SUNSpinsOnLatticeCorrelations.cc` | not in upstream `bin_PROGRAMS`; includes a non-existent header |
| `QuantumDots/src/Tools/Analysis/PeriodicOscillatorForce.cc`, `QuantumDots/src/Programs/ExplicitPeriodic3DQuantumDots.cc`, `…/VisualPeriodic2D.cc` | not in upstream `bin_PROGRAMS` |

Built although upstream never lists them (`EXTRA_PROGRAMS`, each compiles and passes the `--help` smoke test): `QHEBosonsDiskDelta` (a 2006 duplicate of `FQHEDiskBosonsDelta`; retirement candidate), `FQHELatticeEntanglementSpectrum`, `FQHESphereFermionMonteCarloEnergy`, `FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition`, `FQHETorusWithSU2SpinSingleModeApproximation`, `FCIDiceLatticeModel`, `FCIWannierConstruction`.

Renamed in this repository (name collisions at install): `QHEBosonsDelta` → `FQHEDiskBosonsDelta`, `QHEFermionsOverlap` → `FQHEDiskLaughlinMonteCarloOverlap` (the FQHEOnDisk copies; the FQHEOnSphere programs keep their names).

## Install layout

`cmake --install build/<preset> --prefix <dir>`:

```
bin/                      every program            (component Runtime)
lib/diagham/lib<NAME>.a   every static library      (component Development)
include/diagham/...       every public header, one merged tree (no path is duplicated across modules)
lib/cmake/DiagHam/        DiagHamConfig.cmake, exported targets, the Find modules the config needs
```

Downstream: `find_package(DiagHam REQUIRED)` then
`target_link_libraries(myprog PRIVATE DiagHam::DiagHam)` — see
[../how-to/install-and-use-as-library.md](../how-to/install-and-use-as-library.md).
