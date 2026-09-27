# Open defect U30: plain Lanczos (-n 1) returns an energy below the ground state on a degenerate spectrum

Purpose: reproducer for the DiagHam authors; no fix is carried here (`knownbug.hubbard.2x4.U0.plain_lanczos_below_ground_state` in ctest holds it).
Source: `docs/reference/known-defects.md` (U30), found 2026-09-27; r4493 sources, `default` preset, GCC 15.

## What happens

```
HubbardSquareLatticeModel -p 8 -x 2 -y 4 --u-potential 0 -n 1 --full-diag 100
```

writes, for the momentum sectors (1,2) and (0,0), the "ground energies" **−18.412247481626451** and **−17.867733398521771**.

## What was expected

At U = 0 this is free fermions on the 2×4 periodic square lattice: ε(k) = −2(cos kx + cos ky), 8 electrons fill the four lowest levels twice, E₀ = 2·(−4 −2 −2 + 0) = **−16**, and no eigenvalue of the Hamiltonian is below −16. Full diagonalisation (`--full-diag 5000`), Lanczos with `-n 4`, and `--force-reorthogonalize` all give −16.000000000000 in both sectors; these are `physics.hubbard.2x4.U0.tight_binding` in this repository's suite.

## What is known

The value only appears with the plain Lanczos path (`-n 1`, no reorthogonalisation), and only on this highly degenerate spectrum (U = 0). It looks like a Krylov breakdown — the invariant subspace of the start vector is reached in a few steps, β → 0, and the recurrence continues on noise — that the algorithm does not detect. A Ritz value below the true minimum is not possible in exact arithmetic, so a check on β (and/or on the Ritz values against a full diagonalisation of the small tridiagonal matrix) would catch it. For non-degenerate cases (the U = 4 golden) plain Lanczos agrees with full diagonalisation to 10⁻¹³.
