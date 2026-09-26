# Tutorial: your first build, and how to know it is right

Purpose: take a newcomer from a clone to a verified DiagHam build in a few minutes, and show what "verified" means here.
Source: upstream wiki "Install" (sections "Checking your DiagHam installation" and "Testing the configuration", as of 2026-09-24) and `TESTING.md`. Changed: autotools steps replaced by the CMake preset; the `TestDiagHamConf` output shown is from the r4493 default build (the wiki's showed a build with LAPACK, bzip2, GMP and GSL on, and an older `config.h` that printed "Pathscale compiler found").

You need a C++ compiler, CMake ≥ 3.21 and Python 3 (see
[../how-to/build/](../how-to/build/README.md) for platform detail).

## 1. Build

```
git clone <this repository> DiagHam
cd DiagHam
cmake --preset default
cmake --build --preset default -j
```

`default` builds every module with no optional libraries; a few minutes on
a workstation, about a minute on 32 cores. The build tree is
`build/default`; the source tree is untouched (an in-source configure is
refused).

## 2. Ask a program for its help

Every DiagHam program answers `--help`. The wiki's example:

```
build/default/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator --help
```

prints the program's option groups (system, output, precalculation,
parallelisation, misc). The same text, for all 599 programs with an option
parser, is the [program reference](../reference/programs/README.md).

## 3. Check what the build enabled

`TestDiagHamConf` prints the compile-time configuration the binaries were
built with:

```
build/default/src/Programs/TestDiagHamConf
```

```text
__LITTLEENDIAN__ defined
__SSTREAM_STYLE__ defined
MACHINE_PRECISION 1e-14 defined
__SMP__ defined
__DEBUG__ defined
__64_BITS__ defined
x86-64 architecture detected
__64_BITS__ defined
__128_BIT_LONGLONG__ defined
LONGLONG int128_t defined
ULONGLONG uint128_t defined
USE_OUTPUT defined
USE_CLUSTER_ARCHITECTURE defined
USE_HILBERT_SPACE defined
```

With the `lapack` or `full` preset you will also see `__LAPACK__ defined`
(and `__GSL__`, `__GMP__`, `__BZ2LIB__` for `full`). If a feature you
expected is missing here, it is missing from the build, not from the
program.

## 4. Run the physics goldens

```
ctest --preset physics
```

runs the 13 physics tests in a few seconds. They are not "does it link"
checks: each runs a real program in a clean directory and compares the
spectrum it writes with an answer that is known independently —

- the 2×2 Hubbard model at U = 4 must give exactly −4√2 (to 4 units in the
  last place of a double; 0 observed);
- the Laughlin state must be the unique zero-energy state of the V₁
  pseudopotential on the sphere, and three-fold degenerate on the torus;
- Heisenberg rings must give their closed-form energies;
- the spinful torus Coulomb spectrum must agree with two independently
  written programs.

`TESTING.md` lists every test and why its answer is known. The full suite
(`ctest --preset default`) adds a `--help` run of every program and an
install-and-link test; expect 613 tests, all passing.

## 5. Run something

Programs write their output files into the current directory, so work in a
separate directory:

```
mkdir run && cd run
../build/default/FTI/src/Programs/HubbardModels/HubbardSquareLatticeModel \
    -p 4 -x 2 -y 2 --u-potential 4 --nn-t 1.0 --full-diag 1000
sort -k4 -g fermions_hubbard_square_x_2_y_2_n_4_ns_4_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat | head -2
```

```text
# kx ky sz E
0 0 0 -5.6568542494923806
```

That number is −4√2 rounded to double precision. What it means and how to
go further is the [Hubbard walkthrough](hubbard-walkthrough.md); the FQHE
equivalent is [your first FQHE run](first-fqhe-run.md).
