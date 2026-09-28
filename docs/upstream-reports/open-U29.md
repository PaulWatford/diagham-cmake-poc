# Open defect U29: FQHEDiskFermionsTwoBodyGeneric hangs after printing "start"

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/known-defects.md` (U29), found 2026-09-27 while writing the disk goldens; seen on the r4493 sources, `default` and `full` presets, GCC 15 (Ubuntu) and GCC 13 (Ubuntu 24.04, CI).

## What happens

```
printf 'Pseudopotentials = 0 1\n' > v1.dat
FQHEDiskFermionsTwoBodyGeneric -p 4 --interaction-file v1.dat --interaction-name v1 \
    --minimum-momentum 18 --maximum-momentum 18 --full-diag 5000
```

prints `memory requested for Hilbert space`, `nbr interaction = 624`, a line of `=`, `start` — and never returns (killed after 60 s and after 10 min). The same with `--force-maxmomentum 12`, with Lz 20 or 22, and with `-n 1`.

## What was expected

The Lz = 18 sector of 4 fermions with the V₁ pseudopotential holds the Laughlin ν = 1/3 state as its unique zero mode (dimension 1); at Lz = 18 + ΔL there are p(ΔL) = 1, 2, 3, 5 zero modes for ΔL = 1…4. The bosonic twin does exactly this in milliseconds:

```
printf 'Pseudopotentials = 1\n' > v0.dat
FQHEDiskBosonsTwoBodyGeneric -p 4 --interaction-file v0.dat --interaction-name v0 \
    --minimum-momentum 12 --maximum-momentum 16 --full-diag 5000
```

gives 1, 1, 2, 3, 5 zero modes at Lz = 12…16 (`physics.fqhe.disk.*` in this repository's ctest suite).

## What is known

Nothing beyond the symptom: the hang is after the interaction coefficients are built, before any eigenvalue is written. The fermionic disk Hamiltonian path (`ParticleOnDiskGenericHamiltonian` with `FermionOnDisk`) is the suspect; the bosonic path with the same options works.

## Update 2026-09-29

At N = 2 the program runs at L_z = 3 and segfaults at L_z = 5 ("Hilbert space dimension = 3", then the crash); at N = 3 it hangs at L_z = 6 and 9. `FQHEThinAnnulusFermionsTwoBodyGeneric` (same `FermionOnDisk` space) hangs the same way. The bosonic twin runs, so the fermionic Hilbert-space class is the suspect: `FermionOnDisk` reuses the state generation and look-up tables of `FermionOnSphere` with a shifted L_z convention. The pseudopotential Hamiltonian it would use is itself wrong for N >= 3 (U34).
