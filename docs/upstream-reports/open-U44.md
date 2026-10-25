# Open finding U44: FCIHofstadterModel --real-space crashes for spinless fermions

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U44), found 2026-10-25 while writing the core-program tests; r4493 sources, GCC 15. The same on a pristine autotools build of r4493.

## What happens

```
FCIHofstadterModel --real-space -p 2 -x 4 -y 2 -X 1 -Y 3 -q 1 --full-diag 5000
```

prints `case not implented yet` for the first momentum sector and then exits with a segmentation fault (exit 139). It does so with or without interaction (`--u-potential 0`) and with `--redundant-kx --redundant-ky`.

In `FTI/src/Programs/FCI/FCIHofstadterModel.cc`, the real-space branch with translations handles bosons and spinful fermions. For spinless fermions it reaches the `else` at line 896, which only prints the message: `Space` and `Hamiltonian` are never set, and the program goes on to use them.

## Fix direction

Either implement the spinless real-space case (the class `FermionOnLatticeRealSpaceAnd2DTranslation` is already included by the program), or stop with an error message and a non-zero exit before using `Space`.

## Consequence for the tests

`openfinding.fci.hofstadter.realspace_fermions_segfault` (WILL_FAIL) runs the command above.
