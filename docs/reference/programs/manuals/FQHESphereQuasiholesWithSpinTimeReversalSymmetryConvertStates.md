# FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates — manual

Source: DiagHam wiki page `FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates](../FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates.md)
[FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates](FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates.md) convert states written in the quasihole basis with a fixed number of particles to the quasihole basis with all particle numbers. To convert a single state, we just have to use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates --directory nphi_21_cylinder_perimeter_8.000000 --input-state fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_14_2s_21_sz_0_lz_0.0.vec*

nphi_21_cylinder_perimeter_8.000000 is the directory where the Hilbert space description is stored (created with e.g. FQHECylinderQuasiholeMatrixElements  *(wiki page FQHECylinderQuasiholeMatrixElements )*). The code will convert 
the vector fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_14_2s_21_sz_0_lz_0.0.vec defined in the quasihole subspace with 14 particles to the full space, storing the output in the new binary vector fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_fixedn_14_n_0_2s_21_sz_0_lz_0.0.vec. Note that the name is changed and includes a string that indicates that it comes from a state with a fixed number of particles (here fixedn_14).


If we want to convert several states belonging to the same quasihole subspace (same number of particles, angular momentum and number of flux quanta), we can use this alternate form

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates --directory nphi_21_cylinder_perimeter_8.000000 --degenerate-states states_n_14_2s_21.dat*

where states_n_14_2s_21.dat is a single column text file that contains the list of vectors to convert, like

```text
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_14_2s_21_sz_0_lz_0.0.vec
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_14_2s_21_sz_0_lz_0.1.vec
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_14_2s_21_sz_0_lz_0.2.vec
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_14_2s_21_sz_0_lz_0.3.vec
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_14_2s_21_sz_0_lz_0.4.vec
```
