# Open finding U62: complex block Lanczos without --use-lapack never diagonalises its reduced matrix

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U62), found 2026-11-03 while turning Nicolas Regnault's `FTIGenericInteractionFromFileSingleBand` example into a test; r4493 sources, GCC 15.

## What happens

`FullReorthogonalizedComplexBlockLanczosAlgorithm` (`src/LanczosAlgorithm/FullReorthogonalizedComplexBlockLanczosAlgorithm.cc`) diagonalises its reduced (block tridiagonal) matrix only through LAPACK:

- in `Diagonalize()` (from line 529) the branch for `LapackFlag == false`, or a build without `__LAPACK__`, is commented out, so `DiagonalizedMatrix` is never filled;
- `GetGroundState()` (line 294) prints "error, Lapack is required" and goes on.

So `--block-lanczos` on a complex Hamiltonian, run without `--use-lapack` (in any build) or in a build without LAPACK, prints nan from the first iteration, iterates up to `--iter-max`, and ends in a segmentation fault (exit 139). The run gives no other sign that the option combination is unsupported.

Nicolas Regnault's 5 x 5 example (8 fermions, sector kx = 0, ky = 0, dimension 43,263; `tests/data/fti_single_band_regnault`), his options without `--use-lapack`:

| Build | Exit | Time | Lines with nan |
|---|---|---|---|
| default (no LAPACK) | segmentation fault | about 3 minutes | from the first step |
| LAPACK, without `--use-lapack` | 139 | 232 s | 1,494 |
| LAPACK, with `--use-lapack` | 0 | 4.9 s | 0 (E1 = 0.36694040849738, his value to all 14 printed digits) |

A small case shows the same in 0.1 s: four particles, sector (0,0), 506 states, `--full-diag 101 --block-lanczos --block-size 2 --iter-max 100`.

## Fix direction

Stop with an error when the block algorithm is selected without LAPACK (or without `--use-lapack`), or restore the commented-out tridiagonalisation branch.

## Consequence for the tests

`openfinding.fti.generic_interaction_single_band.block_lanczos_without_lapack.n4_k0_0` (WILL_FAIL, every build) compares the two lowest levels of that small case with plain Lanczos. The physics test with block Lanczos (`physics.fti.generic_interaction_single_band.regnault_5x5_n8_k0_0_block_lanczos`) runs with `--use-lapack`, in LAPACK builds only.
