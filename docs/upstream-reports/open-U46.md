# Open finding U46: the built-in full diagonalisation skips Householder steps below 1e-7, so saved eigenvectors can be wrong at about 1e-8

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U46), found 2026-10-26 while writing the entanglement tests; r4493 sources, GCC 15. The program's output files are byte-identical on a pristine autotools build of r4493.

## What happens

Without `--use-lapack`, `GenericRealMainTask` diagonalises a small Hamiltonian with `RealSymmetricMatrix::Householder(TmpTriDiag, 1e-7, Q)` followed by `RealTriDiagonalSymmetricMatrix::Diagonalize(Q)` (`src/MainTask/GenericRealMainTask.cc`, lines 651 and 663). Inside `Householder` (`src/Matrix/RealSymmetricMatrix.cc`, from line 1876), a reflection step is applied only if the squared norm of the column remainder exceeds `SquareErr = err * err`, and `SquareErr` is raised to `MACHINE_PRECISION` (1e-14, `src/config.h`) when it is smaller. A step whose remainder has a norm below 1e-7 is therefore skipped, whatever `err` the caller passes, and that remainder is lost from the matrix.

Eigenvalues change only at second order, so every spectrum agrees with an independent calculation. Eigenvectors change at first order. Example:

```
printf "Pseudopotentials = 1.0\n" > v0.dat
FQHESphereBosonsTwoBodyGeneric -p 4 -l 6 --interaction-file v0.dat --interaction-name x --nbr-lz 1 --eigenstate -n 1 --full-diag 5000
```

The Lz = 0 ground state is the bosonic Laughlin state, unique (gap 1.59) at energy 0. Checked against the same Hamiltonian built independently (`tests/oracles/sphere_ed.py`), the saved vector has a residual |Hv - Ev| = 1.2e-8. With `--use-lapack` the residual is 2.1e-15. With the floor in `Householder` lowered (a local experiment, not carried), it is 2.9e-15.

The error is visible downstream. In an L = 0 state of 4 bosons in 7 orbitals, the pole orbital is occupied with probability exactly 4/7, so the reduced density matrix of the first orbital has eigenvalues 3/7 and 4/7. `FQHESphereBosonEntanglementEntropy` on the saved vector gives 0.42857142857601 instead of 0.428571428571429, and its entropies are off by up to 3e-10. On the vector from `--use-lapack` it agrees to the 14 printed digits.

How often it happens depends on the matrix. Of five small sphere cases checked, two were affected: the Laughlin case above (1.2e-8), and two-body pseudopotentials V0 = 1, V2 = 0.5 with N = 4, 2S = 6 (1.6e-13). The other three were at 1e-14.

## Fix direction

In `Householder` (both overloads), compare the squared norm with `err * err` floored at the square of a machine-precision value, not at `MACHINE_PRECISION` itself, and pass a much smaller `err` from `GenericRealMainTask` (and from any other caller that keeps eigenvectors). Alternatively, use LAPACK for full diagonalisation whenever it is available.

## Consequence for the tests

`openfinding.full_diagonalisation.householder_cutoff.laughlin_eigenvector` (WILL_FAIL) computes the Laughlin state with the built-in full diagonalisation and compares its entanglement with the exact answer at the agreed tolerance. It fails today by 3e-10. When the finding is fixed it passes, and becomes a physics test. The physics tests of `FQHESphereBosonEntanglementEntropy` use states whose saved vectors are accurate (residuals of 1e-14).
