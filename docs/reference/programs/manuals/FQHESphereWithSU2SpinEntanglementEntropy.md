# FQHESphereWithSU2SpinEntanglementEntropy — manual

Source: DiagHam wiki page `FQHESphereWithSU2SpinEntanglementEntropy`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereWithSU2SpinEntanglementEntropy](../FQHEOnSphere/FQHESphereWithSU2SpinEntanglementEntropy.md)
FQHESphereWithSU2SpinEntanglementEntropy computes the entanglement entropy and the entanglement spectrum for the orbital cut for spinful states. It is the analogue of [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md) for spinful bosonic or fermionic systems. Thus the usage of these two programs are similar. Assuming we have created the spin charge separated state fermions_sphere_su2_spinchargeseparated_n_8_2s_9_sz_0_lz_0.0.vec using [FQHESphereFermionsWithSpinThreeBodyGeneric](FQHESphereFermionsWithSpinThreeBodyGeneric.md), we just have to run 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU2SpinEntanglementEntropy fermions_sphere_su2_spinchargeseparated_n_8_2s_9_sz_0_lz_0.0.vec --use-lapack --density-matrix fermions_sphere_su2_spinchargeseparated_n_8_2s_9_sz_0_lz_0.0.full.ent*

We obtain two files fermions_sphere_su2_spinchargeseparated_n_8_2s_9_sz_0_lz_0.0.ent and fermions_sphere_su2_spinchargeseparated_n_8_2s_9_sz_0_lz_0.0.full.ent . The first one contains the entanglement entropy and should look like

```text
1 1.3381046410191 1
2 2.1586164719336 1
3 2.8690779459485 1
4 3.2941682135551 1
5 3.4359653908437 1
```

The first column is the number of kept orbitals, the second column is the entanglement entropy and the third column is the trace of the reduced density matrix (it should be equal to one). The second file is the entanglement spectrum and is only produced if we use the --density-matrix  option. Its contain is formatted as follow

```text
# l_a    N    Lz    Sz    lambda
1 0 0 0 0.39025648114054
1 1 0 -1 0.20974351885947
1 1 0 1 0.20974351885947
1 2 0 0 0.19025648114053
2 0 0 0 0.053661760300387
2 1 -1 -1 0.075250428649627
```

The file is column formatted (using a single space as separator) :
- the first column is the number of orbitals (l_a)
- the second column the number of particles (N)
- the third column is twice the projection of the angular momentum (Lz)
- the forth column is twice the projection of the spin (Sz)
- the fifth column is the eigenvalue (lambda)

For a better accuracy, we can use the --use-svd option.
