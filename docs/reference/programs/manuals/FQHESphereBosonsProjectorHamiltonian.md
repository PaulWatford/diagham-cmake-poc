# FQHESphereBosonsProjectorHamiltonian — manual

Source: DiagHam wiki page `FQHESphereBosonsProjectorHamiltonian`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereBosonsProjectorHamiltonian](../FQHEOnSphere/FQHESphereBosonsProjectorHamiltonian.md)
[FQHESphereBosonsProjectorHamiltonian](FQHESphereBosonsProjectorHamiltonian.md) allows to handle hamiltonians defined through a sum of k-body projectors (see [Simon et al.](http://arxiv.org/abs/0908.0947v2)). All basic options to tune the diagonalization are identical to those of [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md). The sum projector is defined through either a single state (for a single projector) or a list of states provided to a single column formatted text file. The Hamiltonian is then 

$H=\sum_i \left|\Psi_i(z_1,...,z_k)\right\rangle \left\langle \Psi_i(z_1,...,z_k)\right|$

All the states need to have the same number of particles and the same number of flux quanta. If there is a single state for the projector, it can be set *--projector-state*. Otherwise you can use the *--multiple-projectors* . It takes a single column formated text file that gives the list of vectors. Beware that [FQHESphereBosonsProjectorHamiltonian](FQHESphereBosonsProjectorHamiltonian.md) breaks the L symmetry (unless you build the projector to do so). A way to overcome this problem relies on using the *--use-hilbert* option and a subspace with a fixed total angular momentum.


## Generating a state with a given vanishing property

FQHESphereBosonsProjectorHamiltonian can be used to generate a state with a given vanishing property. For example, if you want to reproduce the results of [Simon et al.](http://arxiv.org/abs/0908.0947v2), here is the procedure. In this paper, the authors were interested in the densest $L=0$ states.

- You need to generate the subset of $L=0$ states with the root configuration corresponding to the vanishing propert. In the case of the $S_3$ model, this root configuration is 3000300030003... This can be done using [FQHESphereL2Diagonalize](FQHESphereL2Diagonalize.md). For N=6,9,12 and 15, there are respectively 2,3,6 and 10 such states. In the following, we will consider the N=12 case. So we generate 6 states with names *bosons_l2_n_12_2s_12_lz_0.0.vec*, *bosons_l2_n_12_2s_12_lz_0.1.vec* , ...

- There are two ways to vanish for 4 particles with a relative angular momentum of 4. These two directions have to be generated for the maximum $L$ and $L_z$. For this case, it corresponds to $L_z=4 N_{\Phi} - 8$.

- One needs to create the state that will be involved in the projector. In the case of $S_3$, we want to give some energy penalty to any state that does not have the right vanishing property. It means one has to create a state that is orthogonal to the vanishing direction we want to select. It is a linear combination of the two states built at previous. Let's name this state *bosons_tmp_n_4_2s_12_lz_40.0.vec*. 

- Know we have the direction we want to penalize,  the state that we look for can now be obtain this way :

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonsProjectorHamiltonian -p 12 -l 12 --nbr-lz 1 --projector-state bosons_tmp_n_4_2s_12_lz_40.0.vec -n 1 --eigenstate --use-hilbert basis_l2_n_12.dat --interaction-name myinteraction*

where *basis_l2_n_12.dat* is a text file that describes the basis built during the first step i.e.

```text
Basis=bosons_l2_n_12_2s_12_lz_0.0.vec bosons_l2_n_12_2s_12_lz_0.1.vec bosons_l2_n_12_2s_12_lz_0.2.vec bosons_l2_n_12_2s_12_lz_0.3.vec bosons_l2_n_12_2s_12_lz_0.4.vec bosons_l2_n_12_2s_12_lz_0.5.vec
```

- if everything went fine, you should obtain a spectrum with a single zero energy state

```text
# Lz E
0 2.4881430249479e-11
0 0.04049259429013
0 0.045933089627564
0 0.06278290682125
0 0.098049053486934
0 0.11741297257122
```
