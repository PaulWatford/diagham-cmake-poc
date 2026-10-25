# Cross-benchmarking with ALPS: which programs, and how the conventions map

Purpose: a first note for comparing DiagHam with ALPS on the same models. It proposes the spin-chain and Hubbard programs to keep for that comparison and maps each model's parameters between the two codes, so that a result from one can be checked against the other without guessing signs or factors of 2.
Source: written 2026-10-25 for the team's decision of 20 October 2026 (cross-benchmarking with ALPS where useful; whether spin-chain codes are a formal test priority, or a subset is kept mainly for ALPS, is still open). DiagHam conventions: the tests named in the tables, which fix each program's convention against an independent calculation. ALPS conventions: the ALPS 2.0 paper (Bauer et al., J. Stat. Mech. P05001 (2011), arXiv:1101.2646, appendix D.2) and the ALPS model library. No ALPS run has been made yet; the items marked "to confirm" are checked on the first ALPS installation used.

## Which programs

ALPS can only check a DiagHam program whose model ALPS can also write down. These six programs qualify. Each already has an independent physics test, so their conventions below are measured, not assumed.

| DiagHam program | Model | ALPS model and lattice | DiagHam convention checked by |
|---|---|---|---|
| `GenericOpenSpinChain` | XXZ chain with a field along z, any spin, open (or `--use-periodic`) | `spin` on `open chain lattice` | `physics.spin.ed.open_heisenberg_L8`, `open_xxz_L8_field`, `open_spin1_heisenberg_L6` |
| `GenericPeriodicSpinChain` | XXZ ring, any spin | `spin` on `chain lattice` | `physics.spin.heisenberg_ring_L4` to `L10`, `ring_xxz_L8`, `ring_xxz_L10`, `xx_ring_L6` to `L10` |
| `PeriodicSpinChainJ1J2` | J1-J2 Heisenberg ring | `spin` with two bond types on a chain with next-nearest-neighbour bonds | `physics.spin.ed.j1j2_L8_sz0_k0`, `k1`, `k3` |
| `TwoDimensionalHeisenbergModel` | XXZ model on the square lattice | `spin` on `square lattice` | `physics.spin.ed.heisenberg2d_3x3_jz0.7` |
| `SpinChainXYZ` (transverse-field Ising case) | Ising ring in a transverse field | `spin` with a transverse field (Gamma) | `physics.spin.tfim_ring_L8` |
| `HubbardSquareLatticeModel` | Hubbard model on the square lattice | `fermion Hubbard` on `square lattice` | `physics.hubbard.2x2.*`, `physics.hubbard.2x4.*`, `physics.hubbard.ed.square_2x3_n6_u2` |

The other spin programs stay checked by their own oracles and are not proposed for ALPS, because the ALPS model library has no ready counterpart: the AKLT and generalised AKLT chains (biquadratic terms), Haldane-Shastry (long-range couplings), the Potts, O'Brien-Fendley, double-triangle and pair-hopping chains, the Floquet, PEPS and DMRG programs, and the toric-code, Haah-code and X-cube entropies. Each could be written as a custom ALPS model later if the team wants a second check of it.

## Conventions

### Spin operators

Both codes use spin operators S, with eigenvalues of Sz from −S to S (±1/2 for spin 1/2), not Pauli matrices. The one exception on the list is `SpinChainXYZ`, whose couplings multiply Pauli matrices (below). Twice the spin is `-s` in DiagHam (`-s 1` is spin 1/2, `-s 2` is spin 1) and `S` in ALPS.

### XXZ chains and the square lattice

ALPS's spin model (paper, equation D.1):

H = −h Σᵢ Szᵢ + Σ⟨i,j⟩ [ (Jxy/2)(S⁺ᵢS⁻ⱼ + S⁻ᵢS⁺ⱼ) + Jz Szᵢ Szⱼ ]

The first bond term is Jxy (Sx Sx + Sy Sy). DiagHam's programs:

| ALPS | `GenericOpenSpinChain` | `GenericPeriodicSpinChain` | `TwoDimensionalHeisenbergModel` |
|---|---|---|---|
| Jxy | `-j Jxy` | `-j Jxy` | `--j-value Jxy` |
| Jz | `-z (Jz − Jxy)`: the option adds to the ZZ coupling | `--nn-coupling (Jz − Jxy)`: the same | `--jz-value Jz` |
| h | `--hz-value −h`: DiagHam adds +hz Sz | no field option | no field option |
| open or periodic | open; `--use-periodic` for a ring | periodic | periodic; `--open-boundaries` |

`GenericPeriodicSpinChain --nnn-coupling` is an Sz Sz term only, not a next-nearest Heisenberg coupling, so it is not a J1-J2 model. For J1-J2 use `PeriodicSpinChainJ1J2 -1 J1 -2 J2`, where both couplings are full Heisenberg S·S terms. In ALPS that is the spin model with bond-type couplings (J0 for nearest, J1 for next-nearest neighbours) on a chain lattice that has next-nearest bonds (to confirm: the lattice name and the bond-type parameter names in the installation used).

Set every ALPS coupling explicitly (Jxy, Jz, h, and per bond type). That way the comparison does not depend on ALPS's default values.

### Transverse-field Ising ring

`SpinChainXYZ` (established by matching, `tests/oracles/spin_hubbard.py`): H = Jx Σᵢ σˣᵢ σˣᵢ₊₁ + f Σᵢ Szᵢ, with Pauli matrices in the coupling and S in the field (`-x Jx -y 0 -z 0 -f f -b 1` for a ring). Since σˣσˣ = 4 SˣSˣ, this is 4Jx Σ SˣSˣ + f Σ Sz. Exchanging the x and z axes, which leaves the spectrum unchanged, gives the ALPS form Jz Σ SzSz − Γ Σ Sx with **Jz = 4Jx** and **Γ = −f**. The sign of Γ does not change the spectrum either. Checked numerically: on the L = 8 ring with Jx = 1 and f = 0.5, the ALPS form with Jz = 4 and Γ = ±0.5 gives −8.125499214731438, against the test reference −8.125499214731404 (4e-15 relative). ALPS's transverse-field parameter (Gamma) is in its model library, not in the paper (to confirm).

### Hubbard model on the square lattice

ALPS's `fermion Hubbard` model: H = −t Σ⟨i,j⟩,σ (c†ᵢσ cⱼσ + h.c.) + U Σᵢ nᵢ↑ nᵢ↓ − μ Σᵢ nᵢ (the form of the paper's bosonic Hubbard example, with fermions; to confirm against the model library). `HubbardSquareLatticeModel` has the same hopping sign and on-site U, established by the U = 0 tight-binding tests and `tests/oracles/lattice_fermions.py`:

| ALPS | `HubbardSquareLatticeModel` |
|---|---|
| t | `--nn-t t` |
| U | `--u-potential U` |
| number of electrons N | `-p N` |
| Sz | `--sz-value 2Sz` (twice Sz) |
| μ | no option: at fixed N, E(ALPS) = E(DiagHam) − μN |
| t' (next-nearest) | `--nnn-t`: its sign is not yet fixed by a test |

### Lattices two sites wide

On a periodic lattice two sites wide, DiagHam counts the bond from site 0 to site 1 twice: once directly and once across the boundary. This doubles the hopping or the coupling in that direction. The 2×2 Hubbard results (ground energy −8t at U = 0) and the 2×3 tests depend on it, and `TwoDimensionalTransverseFieldIsingModel` does the same. Whether ALPS's lattice generator does it too is to confirm. Until then, compare only lattices at least 3 sites wide in every periodic direction.

### Sectors and output

- **Sz:** DiagHam labels sectors by 2Sz (an integer) and, where the Sz ↔ −Sz symmetry is used, lists only 2Sz ≥ 0. ALPS labels them by Sz.
- **Momentum:** DiagHam labels a momentum sector by an integer k, meaning momentum 2πk/L. ALPS reports the total momentum as a value (to confirm: units and field name).
- **What to compare:** the safest comparison needs no sector mapping: the ground energy of an Sz sector, or the sorted full spectrum of a sector with DiagHam's other symmetries switched off (`--disable-momentum`, `--disable-szsymmetry`, `--disable-inversion` where the program has them). This is how the tests above compare DiagHam with its oracles.
- **Precision:** DiagHam writes 14 significant digits, the precision the team agreed. Compare at a relative 1e-12, the suite's default tolerance (`DIAGHAM_TEST_TOLERANCE`).

## Next steps

1. The team decides whether the spin-chain programs are a formal test priority or are kept mainly for ALPS. This list is the proposal for the second case.
2. On the first ALPS installation used, confirm the four marked items: the J1-J2 lattice and bond-type names, the transverse-field parameter, the fermion Hubbard model's form, and the bonds of lattices two sites wide.
3. Run one small case per program in both codes (for example L = 8 chains, 3×3 and 4×4 lattices), and record the pairs of results in a new crosscheck data file. Those become `crosscheck.alps.*` tests that run only where ALPS's output files are present, so the suite does not depend on an ALPS installation.
