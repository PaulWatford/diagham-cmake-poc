# Open defect U42: RealSymmetricMatrix::Householder overruns its arrays for a 1x1 matrix

Purpose: reproducer for the DiagHam authors; no fix is carried in the tree.
Source: `docs/reference/known-defects.md` (U42), found 2026-09-30 while auditing the default-build crash of `HaahCodeEntropy`; r4493 sources (`src/Matrix/RealSymmetricMatrix.cc` is byte-identical to SVN), GCC 15.

## What happens

`RealSymmetricMatrix::Householder(RealTriDiagonalSymmetricMatrix& M, double err)` (line 1861 onwards) and the eigenvector overload (starts at line 1876) begin with

```
int ReducedNbrRow = this->NbrRow - 1;
double* TmpUpperDiagonalElements = new double[ReducedNbrRow];   // zero elements for n = 1
...
M.UpperDiagonalElement(ReducedNbrRow - 1) = ...                 // index -1
```

so for n = 1 the loop bodies write and read at index -1 of arrays of length 0. `Diagonalize` on a 1x1 matrix therefore corrupts the heap. `HermitianMatrix` handles n = 1 correctly, and a LAPACK build never reaches this code.

Reproducer (`t2.cc`, linked against the built libraries with `-fsanitize=address`): `RealSymmetricMatrix m(1, true); m.SetMatrixElement(0, 0, 0.5); RealDiagonalMatrix d(1, true); m.Diagonalize(d);` reports a heap-buffer-overflow in `RealSymmetricMatrix::Householder`; without the sanitizer the process aborts with `free(): invalid size`. The same happens with the eigenvector overload.

In the default build `HaahCodeEntropy -x 2 -y 2 -z 2 --nbra-sitex 1 --nbra-sitey 1 --nbra-sitez 1 --export-entspectrum` (a one-site region has 1x1 parity blocks) aborts this way. The LAPACK build gives the entropy 2 ln 2, equal to the GF(2) rank value.

## Fix direction (validated, not applied)

In both overloads, before `int ReducedNbrRow = this->NbrRow - 1;`:

```
if (this->NbrRow == 1)
  {
    M.DiagonalElement(0) = this->DiagonalElements[0];
    return M;
  }
```

On a copy of the sources under AddressSanitizer this gives correct eigenvalues (and eigenvectors) for n = 1, 2 and 5, with and without the eigenvector matrix.

## Consequence for the tests

`knownbug.haah_code.householder_1x1_crash.2x2x2_region_1x1x1` (default build, WILL_FAIL) documents the crash; with LAPACK the same case is the physics test `physics.entanglement.haah_code.2x2x2_region_1x1x1`. The other callers of `RealSymmetricMatrix::Diagonalize` that could receive a 1x1 matrix were only listed statically and are not tested.
