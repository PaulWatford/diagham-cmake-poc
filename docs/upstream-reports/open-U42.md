# Open finding U42: RealSymmetricMatrix::Householder overruns its arrays for a 1x1 matrix

Purpose: reproducer for the DiagHam authors; no fix is carried in the tree.
Source: `docs/reference/open-findings.md` (U42), found 2026-09-30 while auditing the default-build crash of `HaahCodeEntropy`; r4493 sources (`src/Matrix/RealSymmetricMatrix.cc` is byte-identical to SVN), GCC 15.

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

`openfinding.haah_code.householder_1x1_crash.2x2x2_region_1x1x1` (default build, WILL_FAIL) documents the crash; with LAPACK the same case is the physics test `physics.entanglement.haah_code.2x2x2_region_1x1x1`. The other callers were swept dynamically: the library was instrumented to log every 1x1 call (a scratch copy, not in the tree) and the whole suite (895 tests, default build) run, plus about 30 targeted runs (sphere and torus Coulomb with dimension-1 sectors, with and without eigenstates; toric-code, X-cube and Haah regions of one and more sites). Only `HaahCodeEntropy` reached it. The sphere, torus, disk and generic real main tasks guard the call with `GetHilbertSpaceDimension() > 1`. A static scan finds about 330 further `Householder`/`Diagonalize` sites on `RealSymmetricMatrix` (entanglement-entropy, MPS and DMRG programs and tools), which the guard-free library would let through for a 1x1 input; they are unverified, and the library guard above covers all of them at once.

A second caller, found 2026-10-26 while writing its first independent tests: `FQHESphereMPSEntanglementSpectrumParticlePartition` diagonalises each sector's overlap matrix with `RealSymmetricMatrix::Diagonalize` when LAPACK is off (`FQHESphereMPSEntanglementSpectrumParticlePartition.cc`, around line 1215), and some sectors are 1x1. Every bosonic Laughlin run tried (N = 5, nu = 1/2, with and without `--orbital-es`, `--normalize-cylinder` and `--trim-qsector`) crashes there in the default build: AddressSanitizer reports the heap-buffer-overflow in `Householder` (line 2025, the eigenvector overload); without it the run ends with a segmentation fault, `munmap_chunk(): invalid pointer` or `double free or corruption`. The same on a pristine autotools build. With LAPACK its orbital entanglement spectra match the Schmidt decomposition of the exact Laughlin state to 5e-15 (`tests/oracles/mps_oes.py`, four physics tests); without LAPACK `openfinding.fqhe.cylinder.mps_orbital_entanglement_spectrum.householder_1x1_crash.boson_laughlin_m2_n5_perimeter8` (WILL_FAIL) documents the crash.

A third caller, found the same day: `FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity` with a single input state. Upstream, this program calls `LapackDiagonalize` unconditionally and does not build without LAPACK; this repository's patch 03 (`patches/PATCHES.md`) falls back to `RealSymmetricMatrix::Diagonalize` there, and the integrated-charge and imbalance matrices are then 1x1. The run writes its density file, prints "Integrated charge in the upper layer / lower layer:" and aborts (`munmap_chunk(): invalid pointer`). Two or more states (`--degenerate-states`) work in every build; with LAPACK a single state works too. `openfinding.fqhe.sphere.quasiholes_time_reversal_density.householder_1x1_crash.pairing_2s5_sz1_lz1_single_state` (WILL_FAIL, builds without LAPACK) documents it; with LAPACK the same case is a physics test.
