# Open defect U39: the one-body energy of the single-band many-body Hamiltonians is wrong

Purpose: reproducer for the DiagHam authors; no fix is carried here (`knownbug.fci.checkerboard.manybody.*`, `knownbug.fci.haldane.manybody.*` and `knownbug.fci.kagome.manybody.*` in ctest hold it).
Source: `docs/reference/known-defects.md` (U39), found 2026-09-29 while writing the many-body goldens; r4493 sources, `default` preset, GCC 15.

## What happens

Without an interaction the many-body spectrum must be a sum of band energies. It is not.

```
FCIHaldaneModel -p 1 -x 3 -y 3 --single-band --full-diag 5000 --full-momentum
FCIKagomeLatticeModel -p 1 -x 3 -y 3 --full-diag 5000 --full-momentum
FCICheckerboardLatticeModel -p 2 -x 3 -y 3 --single-band --u-potential 0 --full-diag 5000 --full-momentum
```

- Haldane: the nine energies are exactly half of the band energies of `--singleparticle-spectrum` (lowest −2.2675121777 against −4.5350243554; 1.3349 against 2.6697). N = 2 gives the sum of the halved energies.
- Kagome: −1.4 against −2.8 (also exactly half); the other momenta the same.
- Checkerboard (non-flat): −0.5277 (twice), 0 (six times), 0.5277 (twice) for N = 1 on 3×3; the band energies are −2.0000, −2.0792 and so on.

## What was expected

The band energy of the lowest band at each momentum (checked against numpy eigenvalues of the published Bloch Hamiltonians, `tests/oracles/fci_bands.py`, which the single-particle export of the same programs reproduces to 10⁻¹²). For the checkerboard at U > 0 the spectrum of the band-projected problem `sum_k eps_k n_k + U sum_bonds P n_A n_B P` (`tests/oracles/fci_manybody.py`, which reproduces the program's flat-band spectra, where eps = 0, to 4·10⁻¹² for N = 2 and 3).

## What is known

- Haldane and Kagome (and 24 sibling `ParticleOnLattice*SingleBandHamiltonian.cc` files, not run): `this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);` (`ParticleOnLatticeHaldaneModelSingleBandHamiltonian.cc:137`, `ParticleOnLatticeKagomeLatticeSingleBandHamiltonian.cc:135`). Whether the 0.5 is meant to be compensated elsewhere (the one-body term applied twice in the operator, or an N − 1 normalisation) was not found; N = 1 and N = 2 show it is not.
- Checkerboard: `ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian.cc:139` stores `GetEnergy(0, Index)`, and the block at lines 326 to 345 (executed when `FlatBand == false`) assigns `OneBodyInteractionFactors[Index1] = Sum` afterwards, where `Sum` is the real part of the interaction-derived exchange terms; it does not contain the band energy and does not scale with U.
- With `--flat-band` the one-body energy is zero and every check passes: the interaction matrix elements are right.
