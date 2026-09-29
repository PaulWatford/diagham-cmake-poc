# Open defect U35: TwoDimensionalTransverseFieldIsingModel ignores open boundaries for the Ising bonds

Purpose: reproducer for the DiagHam authors; no fix is carried here (`knownbug.spin.tfim2d.open_boundaries_ignored` in ctest holds it).
Source: `docs/reference/known-defects.md` (U35), found 2026-09-29 while writing the spin-model goldens; r4493 sources, `default` preset, GCC 15.

## What happens

```
TwoDimensionalTransverseFieldIsingModel -x 2 -y 3 --jz-value 1 --hx-value 0.7 --hz-value 0.2 --full-diag 5000
TwoDimensionalTransverseFieldIsingModel -x 2 -y 3 --jz-value 1 --hx-value 0.7 --hz-value 0.2 --full-diag 5000 --use-periodic --disable-momentum --disable-inversion
```

The first run writes `spin_1_2_ising_transversefield_open_n_6_x_2_y_3_...dat`, the second `..._closed_...dat`, and the two spectra are identical to 10⁻¹⁵ (64 levels). The same happens for 3×3 (512 levels).

## What was expected

With open boundaries a 2×3 cluster has 7 nearest-neighbour bonds and the spectrum of H = Jz Σ SᶻSᶻ + hx Σ Sˣ + hz Σ Sᶻ on them, from an independent dense diagonalisation (`tests/oracles/spin_models.py`, `tfim_2d(2, 3, 1, 0.7, 0.2, wrap=False)`), differs from the program's by up to 0.77. The program's numbers match the same oracle with the wrap-around bonds added, 12 bonds in all (for nx = 2 the wrap-around doubles the x bonds).

## What is known

In `Spin/src/Hamiltonian/TwoDimensionalTransverseFieldIsingHamiltonian.cc`, `EvaluateDiagonalMatrixElements` adds the bonds (0,k)-(nx−1,k) for every k and (j,0)-(j,ny−1) for every j after the bulk bonds, with no test of the periodic flag that the constructor receives from the program's `--use-periodic` option. The field terms are unaffected. Guarding those two loops with the flag is the obvious fix; it is not carried here because it changes computed numbers and the authors may prefer to decide what "open" should mean for nx = 2 (whether a doubled bond is intended when periodic).
