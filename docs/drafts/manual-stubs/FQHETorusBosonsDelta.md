# FQHETorusBosonsDelta (wiki page, not attached)

Source: DiagHam wiki page `FQHETorusBosonsDelta`, as of 2026-09-24. Reason: program renamed upstream; the current equivalent was not identified with certainty.

FQHETorusBosonsDelta provides exact diagonalization of two-body hardcore hamiltonian for bosons on the torus geometry. Only the translation along the y direction is taken into account. Many options are similar to  the ones of [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md). A typical usage will look like 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusBosonsDelta -p 3 -l 6 -n 1 --eigenstate --use-lapack*

This will compute the spectrum and the lowest eigenstate in each ky momentum sector for 3 particles and 6 flux quanta. The naming convention for the eigenstates is bosons_torus_kysym_delta_n_x_2s_y_ratio_r_ky_k.n.vec where 

- x is the number of particles
- y is the number of flux quanta
- r is the aspect ratio
- k is the momentum along ky
- n is an index to label different states belong to the same sector (defined through x,y and k)

A typical file name (for example, produced using the previous example) will look like bosons_torus_kysym_delta_n_3_2s_6_ratio_1.000000_ky_0.0.vec
