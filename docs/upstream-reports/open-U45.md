# Open finding U45: FCIHofstadterModel --real-space hard-core bosons lie below their own lower bound

Purpose: reproducer for the DiagHam authors; no fix is carried here. This is a question for the authors: which Hamiltonian the real-space mode is meant to build.
Source: `docs/reference/open-findings.md` (U45), found 2026-10-25 while writing the core-program tests; r4493 sources, GCC 15. The program's output file is byte-identical on a pristine autotools build of r4493.

## What happens

The one-particle spectrum of the program agrees with an independent Hofstadter exact diagonalisation (`tests/oracles/hofstadter.py`, Landau gauge, the hop to +x carrying e^{2 pi i phi y}). Its lowest band is a physics test: `physics.fci.hofstadter.one_particle_lowest_band.*`, three lattices, to about 1e-14. For 4 x 2 cells of 1 x 3 sites with one flux quantum per cell (24 sites, phi = 1/3), its lowest one-particle energy is e_0 = -(1 + sqrt 3) = -2.7320508.

Two hard-core bosons in the real-space mode,

```
FCIHofstadterModel --real-space --redundant-kx --redundant-ky -p 2 -x 4 -y 2 -X 1 -Y 3 -q 1 --boson --hardcore --full-diag 5000
```

give 276 levels (the correct count, C(24, 2)) with ground state -6.1290191820741. This is impossible for the same one-particle Hamiltonian. The hard-core space is a subspace of the two-boson space, so its lowest energy cannot be below the two-boson minimum 2 e_0 = -5.4641016. The consistent model gives -5.3042053409.

The spectrum does not change with `--u-potential` or `--v-potential`, as expected for hard-core bosons. No uniform flux on this lattice (phi = k/24 in either Landau gauge) reproduces it. Three hard-core bosons on 3 x 2 cells of 1 x 4 sites (phi = 1/4) disagree in the same way.

## Where to look

The real-space path builds its hopping from `TightBindingModel->GetRealSpaceTightBindingHamiltonian()`, while the band projection uses the Bloch Hamiltonian. A difference between the two (a phase or gauge convention in the real-space embedding, or the flux assigned per plaquette) would explain it. This has not been checked in the source.

## Consequence for the tests

`openfinding.fci.hofstadter.realspace_hardcore_bosons_below_bound.n2_x4_y2` and `.n3_x3_y2` (WILL_FAIL) carry the spectra of the consistent model (`tests/data/hofstadter/hardcore_bosons_*.dat`). If the authors say the real-space mode intends a different model, these references change, not the program.
