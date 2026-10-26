# Open finding U52: SpinChainComputeS2 --degenerate-states uses only the first state when the file name has no momentum

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U52), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15.

## What happens

With `--degenerate-states LIST`, the program should diagonalise S^2 in the span of the listed states and print one value per state, as its wiki page shows. It reads the list into `InputStateNames[0..NbrStates-1]` (`Spin/src/Programs/SpinChainComputeS2.cc`, around line 141), but only the branch for states with momentum quantum numbers (from line 506) uses them all. For an open chain, or any state whose file name carries no momentum, the program takes the other branch (around line 285), reads `InputStateNames[0]` only, and prints <S^2> of that first state. The other states are ignored without a message.

Example: three orthonormal mixtures of S = 0, 1 and 2 states of six spins 1/2 in the Sz = 0 sector (`tests/oracles/spin_tools.py`, case `6 1 0 3`): the program prints one value (0.676, the expectation in the first mixture) instead of 0, 2 and 6.

## Fix direction

Build and diagonalise the S^2 matrix of the listed states in the non-translation branch too, as the translation branch does, or reject `--degenerate-states` there with a message.

## Consequence for the tests

`openfinding.spin.compute_s2.degenerate_states_without_momentum_first_only` (WILL_FAIL) expects the three values and fails while only one is printed. The single-state mode is a physics test (three cases).
