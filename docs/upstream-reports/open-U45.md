# Open finding U45: FCIHofstadterModel --real-space accepts -Y other than 1 and gives a wrong spectrum

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U45), found 2026-10-25 while writing the core-program tests; r4493 sources, GCC 15. The program's output file is byte-identical on a pristine autotools build of r4493. The program's wiki page already warns about this; what is missing is a check in the program.

## What happens

The program's wiki page says: "Note that the code is only working if -Y is set to one, another choice might lead to erroneous results." The program does not enforce this.

**With `-Y 1` the real-space mode is right.** Two and three hard-core bosons on 12 x 2 sites, at flux 1/3 and 1/4, match an independent exact diagonalisation (`tests/oracles/hofstadter.py`) to about 1e-14: `physics.fci.hofstadter.realspace_hardcore_bosons.*`. The one-particle lowest band is right with `-Y 3` and `-Y 4` too (`physics.fci.hofstadter.one_particle_lowest_band.*`).

**With `-Y 3` the result is wrong.** For example:

```
FCIHofstadterModel --real-space --redundant-kx --redundant-ky -p 2 -x 4 -y 2 -X 1 -Y 3 -q 1 --boson --hardcore --full-diag 5000
```

gives 276 levels (the right count, C(24, 2)) with ground state -6.1290191820741. That is impossible for the program's own one-particle Hamiltonian, whose lowest energy is e_0 = -(1 + sqrt 3) = -2.7320508. The hard-core space is a subspace of the two-boson space, so its lowest energy cannot be below 2 e_0 = -5.4641016. The consistent model gives -5.3042053409. Three hard-core bosons with `-Y 4` disagree in the same way.

## Fix direction

Refuse `--real-space` with `-Y` other than 1, with a message, as the wiki describes. Alternatively, make the real-space hopping (`GetRealSpaceTightBindingHamiltonian`) agree with the Bloch Hamiltonian for any unit cell; the band mode already works for any `-Y`.

## Consequence for the tests

`openfinding.fci.hofstadter.realspace_unit_cell_y_not_1_accepted.n2_x4_y2` and `.n3_x3_y2` (WILL_FAIL) carry the spectra the program should give with `-Y 3` and `-Y 4` (`tests/data/hofstadter/hardcore_bosons_*_Y3_*` and `_Y4_*`). If the fix is to refuse the input, these tests become tests that the program stops with an error.
