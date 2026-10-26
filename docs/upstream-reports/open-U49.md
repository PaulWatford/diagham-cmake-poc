# Open finding U49: FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing crashes (SIGFPE) after its first sector that uses fast multiplication

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U49), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The same crash on a pristine autotools build of r4493.

## What happens

With the default `--memory` (fast multiplication on), the program computes and writes the first Lz sector large enough to use fast multiplication, then dies with a floating-point exception (signal 8, exit 136). Every later sector is lost. The wiki's own example crashes:

```
FQHESphereQuasiholeMatrixElements -l 6
FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 6 --charging-energy 1 --average-nbrparticles 4 --interaction-file onebodypotentials_2s_6.dat --interaction-name test --full-diag 5000 --directory nphi_6_sphere
```

(with the wiki's `onebodypotentials_2s_6.dat`; note that the matrix-element program now writes `nphi_6_sphere`, not `nphi_6` as the wiki says). Small sectors at the edge of the Lz range, which do not use fast multiplication, run to the end.

A build with GCC's `-fsanitize=integer-divide-by-zero` locates it: `FQHE/src/Hamiltonian/AbstractQHEOnSphereHamiltonian.cc:92`, in the destructor, `EffectiveHilbertSpaceDimension / this->FastMultiplicationLargeStep`. The base constructor sets `FastMultiplicationLargeStep` to 0 (line 69). `ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing` sets up fast multiplication its own way and never sets it, so destroying the Hamiltonian divides by zero.

## Workaround

`--memory 0` disables fast multiplication; the program then runs every sector. The spectra are right: with `--memory 0` the full two-layer spectrum (one-body potentials, pairing, charging energy) agrees with an independent calculation to 4e-14 (`tests/oracles/quasiholes_trs.py`).

## Fix direction

Set `FastMultiplicationLargeStep` (and the related fields the base destructor reads) when this Hamiltonian enables fast multiplication, or give it a destructor that frees its own fast-multiplication arrays without calling into the base class's.

## Consequence for the tests

`openfinding.fqhe.sphere.quasiholes_time_reversal_pairing.default_memory_crash` (WILL_FAIL) runs the default `--memory` and fails while the crash is there. The physics tests of the program use `--memory 0`.
