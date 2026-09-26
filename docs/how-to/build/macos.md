# Build on an Apple-silicon Mac

Purpose: build DiagHam on an M1/M2/… Mac.
Source: upstream wiki "Install", sections "Compiling on M1/M2/... Mac", "Common issues" (fortran2c note) and the Mac note under "First installation" (as of 2026-09-24). Changed: configure lines mapped to CMake options. **Not yet verified under this build on a Mac.**

Two options, as on the wiki.

## 1. Homebrew GCC + OpenBLAS

Install `gcc`, `openblas`, `lapack` (and `gmp` if wanted) with Homebrew. The
wiki's configure line

```
./configure --enable-fqhe --enable-fti --enable-gmp --enable-lapack --with-blas-libs=-lopenblas --with-lapack-libs=-llapack 'CFLAGS=-I/opt/homebrew/include ...' 'LDFLAGS=-L/opt/homebrew/lib ...'
```

becomes

```
cmake --preset lapack -DDIAGHAM_USE_GMP=ON -DBLA_VENDOR=OpenBLAS \
      -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/openblas;/opt/homebrew/opt/lapack;/opt/homebrew"
cmake --build --preset lapack -j
ctest --preset lapack
```

If `find_package(LAPACK)` picks the wrong library, give the link line
directly, as `--with-*-libs` did:
`-DDIAGHAM_LAPACK_LIBS="-L/opt/homebrew/opt/openblas/lib -lopenblas"`.

The drawback the wiki notes still applies: OpenBLAS/LAPACK from Homebrew do
not benefit from Apple's optimisations.

## 2. Apple's Accelerate BLAS through R's libraries

The wiki's second recipe uses the BLAS/LAPACK shipped with R (which can be
switched to Apple's vecLib):

```
../configure --enable-fqhe --enable-fti --enable-gmp --enable-lapack --with-lapack-libs="-L/Library/Frameworks/R.framework/Resources/lib/ -lRlapack" --with-blas-libs="-L/Library/Frameworks/R.framework/Resources/lib/ -lRblas" ...
```

becomes `-DDIAGHAM_LAPACK_LIBS="-L/Library/Frameworks/R.framework/Resources/lib -lRlapack -lRblas"`
on the `lapack` preset.

## Notes

- Case: on older macOS the wiki reports LAPACK/BLAS found only as
  `-lLAPACK -lBLAS`; the same `DIAGHAM_LAPACK_LIBS` override covers it.
- Fortran runtime: if the link complains about `g2c`/`f2c`/`gfortran`, the
  code runs fine without them; the old `--with-fortran2c-libs=""` is
  `-DDIAGHAM_FORTRAN2C_LIBS=""` (the default).
- Linker: with Apple's own compiler (`AppleClang`) the static libraries are
  linked without `--start-group` (Apple's linker rescans archives). A
  Homebrew LLVM `clang` reports itself as `Clang` and would be given the GNU
  `--start-group` flags, which Apple's linker rejects — known issue; use
  Apple's compiler or GCC for now.
- Apple silicon is 64-bit ARM: upstream's `src/config.h` (r4493) has an
  `__aarch64__` branch for the 128-bit integer support; the 2020 snapshot did
  not, so an older tree will not build here.
