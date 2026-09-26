# Tutorial: your first FQHE run — the Laughlin state on the sphere

Purpose: generate the ν = 1/3 Laughlin state two independent ways and check they agree; on the way, meet pseudopotential files, binary vectors, the squeezed (Haldane) basis and the Jack generator.
Source: upstream wiki "Laughlin_state" (as of 2026-09-24). Changed: the program `QHEFermionsTwoBodyGeneric` is now `FQHESphereFermionsTwoBodyGeneric`; paths use the CMake build tree; every command and number was re-run on the r4493 build (`lapack` preset) on 2026-09-26.

We take 8 fermions at filling 1/3 on the sphere: 2S = 3(N−1) = 21 flux
quanta. Work in a scratch directory; the programs write into the current
directory.

## 1. Exact diagonalisation of the hollow-core interaction

The Laughlin state is the exact zero-energy ground state of the
interaction with only the V₁ pseudopotential. Create
`pseudopotential_laughlin3_2s_21.dat`:

```
Pseudopotentials = 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

(22 entries for 2S = 21; V₀ is irrelevant for fermions, V₁ = 1 is the
hollow-core interaction.) Then diagonalise in the Lz = 0 sector:

```
../build/lapack/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsTwoBodyGeneric \
    -p 8 -l 21 --interaction-name laughlin --interaction-file pseudopotential_laughlin3_2s_21.dat \
    --use-lapack --eigenstate -n 1 --nbr-lz 1
```

(`--use-lapack` exists only in a LAPACK-enabled build; drop it with the
`default` preset.) The Hilbert space has 8,512 states; the run takes a
fraction of a second. The spectrum file `fermions_laughlin_n_8_2s_21_lz.dat`
shows the zero-energy state:

```text
# Lz E
0 1.4210854715202004e-14
```

and `--eigenstate` wrote it as the binary vector
`fermions_laughlin_n_8_2s_21_lz_0.0.vec` (format:
[data formats](../reference/data-formats.md)).

## 2. The same state from the Jack generator

A more accurate way is the Jack polynomial: the Laughlin state is the Jack
with α = −2 and root configuration 100100…1. Describe the root
configuration in `laughlin_n_8_2s_21.dat`:

```
NbrParticles=8
LzMax=21
ReferenceState=1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1
```

and generate:

```
../build/lapack/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator \
    -a -2 -o fermions_haldane_laughlinjack_n_8_2s_21_lz_0.0.vec -n \
    --reference-file laughlin_n_8_2s_21.dat --fermion
```

`-a -2` is α, `--fermion` multiplies by the Vandermonde, `-n` normalises,
`-o` writes a binary vector. The output lives in the **squeezed (Haldane)
basis** — the 5,302 configurations reachable from the root by squeezing,
not the full 8,512-state Lz = 0 basis. Convert it:

```
../build/lapack/FQHE/src/Programs/FQHEOnSphere/FQHESphereConvertHaldaneBasis \
    fermions_haldane_laughlinjack_n_8_2s_21_lz_0.0.vec \
    --reference-file laughlin_n_8_2s_21.dat -o fermions_laughlinjack_n_8_2s_21_lz_0.0.vec
```

## 3. Compare the two

```
../build/lapack/src/Programs/GenericOverlap \
    fermions_laughlinjack_n_8_2s_21_lz_0.0.vec fermions_laughlin_n_8_2s_21_lz_0.0.vec
```

```text
File 0  fermions_laughlinjack_n_8_2s_21_lz_0.0.vec
File 1  fermions_laughlin_n_8_2s_21_lz_0.0.vec
Overlap |<0|1>|^2 = 0.999999999999994
```

Two independent constructions — a diagonalisation and an algebraic
recursion — agree to 6×10⁻¹⁵. That agreement is what the ctest golden
`physics.fqhe.sphere.laughlin_1_3.unique_zero_mode` checks in miniature
(N = 6: exactly one zero mode among 338 states).

## Where next

- The manual pages for the programs used:
  [FQHESphereFermionsTwoBodyGeneric](../reference/programs/manuals/FQHESphereFermionsTwoBodyGeneric.md),
  [FQHESphereJackGenerator](../reference/programs/manuals/FQHESphereJackGenerator.md),
  [FQHESphereConvertHaldaneBasis](../reference/programs/manuals/FQHESphereConvertHaldaneBasis.md).
- Other geometries (torus, cylinder, disk) and the FCI programs: the
  [program reference](../reference/programs/README.md).
- Running larger systems in parallel: [../how-to/run-mpi.md](../how-to/run-mpi.md).
