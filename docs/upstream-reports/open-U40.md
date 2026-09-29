# Open defect U40: SpinChainEntanglementEntropy aborts in a build without LAPACK

Purpose: reproducer for the DiagHam authors; no fix is carried here (`knownbug.spin_chain_entanglement.use_svd_ungated.*` in ctest hold it in the default build).
Source: `docs/reference/known-defects.md` (U40), found 2026-09-29 while writing the entanglement goldens; r4493 sources, `default` preset (no LAPACK), GCC 15.

## What happens

```
SpinChainEntanglementEntropy spin_1_2_n_2_sz_0.0.vec -o bell.ent
Option 'use-svd' was requested, but is not implemented!
```

exit status 255, for every input vector and with or without `-s`. `tests/data/entanglement/spin_1_2_n_2_sz_0.0.vec` is the Bell pair (1, -1)/sqrt 2 as int32 dimension plus doubles.

## Why

`Spin/src/Programs/SpinChainEntanglementEntropy.cc` declares `use-lapack`, `use-svd` and `diag-precision` inside `#ifdef __LAPACK__` (lines 85-90), but reads `Manager.GetBoolean("use-svd")` (line 115) unconditionally, and the non-LAPACK branch of the diagonalisation calls `Manager.GetDouble("diag-precision")` (line 688), an option that does not exist there. Reading an undeclared option aborts.

## What was expected

The internal (non-LAPACK) diagonaliser is the documented default ("convergence precision in non LAPACK mode"). In a LAPACK build the program is right: the Bell pair gives S = ln 2 (0.6931471805599453) with two reduced-density-matrix levels of 1/2, the product state gives S = 0 with one level equal to 1 (the four `physics.entanglement.spin_chain.*` tests, `DIAGHAM_USE_LAPACK=ON`).

## Fix direction

Guard line 115 with `#ifdef __LAPACK__` (`SVDFlag = false` otherwise) and read `diag-precision` only where it is declared, or declare `diag-precision` outside the guard.
