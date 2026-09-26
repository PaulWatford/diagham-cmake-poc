# Build on Linux with Clang / LLVM

Purpose: build DiagHam with clang instead of GCC.
Source: upstream wiki "Install", section "Using LLVM" (last upstream test: September 2024 on Arch Linux; as of 2026-09-24). Changed: configure line replaced by the equivalent preset; CI now builds the `core` preset with Clang 18 on every change.

DiagHam compiles without any problem with LLVM (clang, or the macOS gcc
front end to LLVM). On a Linux machine the equivalent of the wiki's

```
CC="clang" CXX="clang++" ../configure --enable-fqhe --enable-lapack --enable-gsl --enable-bz2 --enable-gmp --enable-spin --enable-fti
```

is

```
CC=clang CXX=clang++ cmake --preset full
cmake --build --preset full -j
ctest --preset full
```

(`full` = LAPACK, GSL, GMP, FFTW3, bzip2; all modules are on by default.
None of the optional libraries is mandatory — use `--preset default` for a
build with no external libraries.)

Notes:

- The build uses the GNU dialect (`-std=gnu++11`), as upstream's g++ does;
  clang accepts it.
- Static libraries are linked inside `--start-group`/`--end-group` on GNU ld
  and on clang with GNU ld. Apple's linker does not take these flags; see
  [macos.md](macos.md) for the Mac case.
- Expect the same warnings as with GCC (unused parameters, overloaded
  virtuals) — DiagHam is 2003-era C++; none of them is an error.
