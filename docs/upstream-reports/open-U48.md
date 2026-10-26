# Open finding U48: FTIGenericInteractionFromFileTwoBands applies the interaction file's band labels to the opposite bands unless --3band-convention is given

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U48), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The program's output is byte-identical on a pristine autotools build of r4493.

## What happens

Without spin or valley, the program reads one-body energies `kx ky m E` (`--singleparticle-file`) and two-body matrix elements `m1 kx1 ky1 ... m4 kx4 ky4 V` (`--interaction-file`), both with a band index m of 0 or 1, as on its wiki page. By default (no `--3band-convention`, no band caps) the Hilbert space is `FermionOnSquareLatticeWithSpinMomentumSpace`, which inherits from `FermionOnSphereWithSpin`. There the one-body operator `AdsigmaAsigma` reads band sigma at bit 2k + sigma, but the two-body operators `AsigmaAsigma` and `AdsigmaAdsigma` use bit 2k + 1 - sigma (`FQHE/src/HilbertSpace/FermionOnSphereWithSpin.h`, from line 1131). So the interaction written for band m acts on band 1 - m of the one-body file. With `--3band-convention` the program uses `FermionOnSquareLatticeWithSU2SpinMomentumSpace`, which overrides the two-body operators with bit 2k + sigma, and both files then mean the same bands. The option's help text calls it "the safest option"; the wiki's first example does not use it, and neither says what the default does.

Example, two particles in the smallest grid (`-x 2 -y 1`), one-body energies 0.0 and 0.1 (k = 0, bands 0 and 1), 3.0 and 3.7 (k = 1), and a single interaction line `0 0 0 0 1 0 0 0 0 0 1 0 1.0` (both particles in band 0, at k = 0 and k = 1): the default convention shifts the two-particle state with both particles in band 1 (energy 3.8 -> 2.8); `--3band-convention` shifts the state with both particles in band 0 (energy 3.0 -> 2.0).

With a generic interaction (`tests/oracles/fti_generic.py`), `--3band-convention` agrees with an independent exact diagonalisation of the Hamiltonian as written to 4e-14, and the default agrees to 4e-14 only once the interaction's band labels are exchanged.

## Fix direction

Make `FermionOnSphereWithSpin`'s one-body and two-body operators use the same bit for a given sigma, or make `--3band-convention` the default of this program, or document the default convention on the wiki and in the help text. For spin-1/2 problems whose terms are symmetric between up and down, the mismatch has no effect, which may be why it has gone unnoticed.

## Consequence for the tests

`physics.fti.generic_interaction_two_bands.*_3band_convention` (three cases) check the program against the Hamiltonian as written. `physics.fti.generic_interaction_two_bands.fermions_x3_y2_n2_default_convention` pins the default behaviour, with the label exchange stated in the oracle, so a change of the default is noticed.
