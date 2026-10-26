# Open finding U51: SpinChainMultipleEntanglementSpectra crashes in a build without LAPACK

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U51), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The same crash on a pristine autotools build of r4493 (configured without LAPACK).

## What happens

In a build without LAPACK, the program segfaults on every input tried, including DiagHam's own eigenvector matrices:

```
GenericOpenSpinChain -p 8 --nbr-sz 1 --full-diag 5000 --all-eigenstates
SpinChainMultipleEntanglementSpectra --multiple-states spin_1_2_openchain_n_8_j_1.000000_sz_0.eigenvec.mat
```

(and the same with a `PeriodicSpinChainAKLT --all-eigenstates` matrix, the wiki's kind of input). An address-sanitizer build locates it: `Spin/src/Architecture/ArchitectureOperation/SpinChainMultipleEntanglementSpectrumOperation.cc:195`, which squares the values returned by `RealMatrix::SingularValueDecomposition()`. That function's body is entirely inside `#ifdef HAVE_LAPACK` (`src/Matrix/RealMatrix.cc`, from line 1415), so without LAPACK it returns a null pointer. With LAPACK the program is right: its spectra agree with an independent Schmidt decomposition to about 1e-16 (`tests/oracles/spin_tools.py`).

## Fix direction

Refuse to run without LAPACK, with a message (as for U40), or fall back on diagonalising rho_A = M M^T, which needs no SVD.

## Consequence for the tests

In builds without LAPACK, `openfinding.spin.multiple_entanglement_spectra.svd_null_without_lapack.spin_1_2_n8_sz0` (WILL_FAIL) fails while the crash is there; in builds with LAPACK three cases are physics tests.
