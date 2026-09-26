# Why CMake, and why a generated build

Purpose: the reasoning behind replacing the autotools build with a generated CMake build.
Source: `README.md` sections "Why CMake" and "Methodology" (proof-of-concept era, 2026), lightly expanded 2026-09-26.

## What the autotools build is

DiagHam's build is a 2003-era GNU autotools system: `bootstrap.sh`
(`aclocal`, `automake`, `autoconf`, `autoheader`, twice), a `configure.ac`
with 34 `--enable`/`--with` flags — LAPACK, MKL, MPI, ScaLAPACK, GSL, GMP,
MPACK, FFTW, bzip2, profiling, the module switches — and per-directory
`Makefile.am` files kept in step with the source tree by two custom Perl
scripts (`scripts/genmake.pl`, `scripts/genam.pl`). Each flag is some
cluster administrator's hard-won dependency chain. It works, but it is
brittle (out-of-tree includes are fragile, programs silently fall out of
`bin_PROGRAMS`, compilers newer than the code reject idioms autotools never
checked), slow to configure, and unfamiliar to anyone who learned to build
software after 2010.

## Why CMake

CMake is what current compilers, IDEs, package managers and HPC module
systems expect. It gives, for no extra work: out-of-source builds and
several configurations of one tree; presets that replace folklore configure
lines; `ctest` for a test suite; `install()` and an exported package so
DiagHam can be used as a library from another project; and a build that
Windows, macOS and Linux toolchains all read. The Kent group's own projects
(BDMC_UFL) already use it, and their Find modules and defaults were adopted
directly.

## Why generate the build instead of writing it

Hand-writing CMake for ~1,450 classes and 603 programs would be wrong the
day upstream added a file. Instead `scripts_cmake/extract_autotools.py`
reads upstream's own `Makefile.am` files — the library source lists, the
`SUBDIRS` order, the `bin_PROGRAMS` lists — and writes the per-directory
`CMakeLists.txt`. The CMake build therefore *tracks* DiagHam rather than
forking it: after an upstream change, one script run gives the updated
build. All complexity lives in one helpers file; each generated file is a
few lines anyone can read. The generated files are committed so the
repository builds without running the script.

The one deliberate departure from "exactly what autotools builds" is a
short, documented allow-list of sources upstream never listed
(`docs/reference/build-system.md`, "Exclusions and extras").

## What "correct" means here

Three levels, each checked: the build (every library and program upstream
builds, we build); parity with autotools (`cmake/verify_build.sh`, library
coverage and symbol counts against an autotools build of the same tree);
and physics (ctest goldens with exact or independently computed answers,
plus a from-scratch Python diagonaliser). Compiling is not correctness; the
Hubbard 2×2 ground state matching −4√2 to the last bit is.
