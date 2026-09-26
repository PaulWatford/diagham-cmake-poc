# FQHETorusFermionsWithTranslationsNBodyHollowCore — manual

Source: DiagHam wiki page `FQHETorusFermionsWithTranslationsNBodyHollowCore`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusFermionsWithTranslationsNBodyHollowCore](../FQHEOnTorus/FQHETorusFermionsWithTranslationsNBodyHollowCore.md)
FQHETorusFermionsWithTranslationsNBodyHollowCore diagonalizes the shortest (i.e. hollowcore) n-body interaction for a fermionic system on the torus geometry. It is restricted to rectangular tori (for a generic n-body interaction, see [FQHETorusFermionsWithTranslationsGenericNBody](FQHETorusFermionsWithTranslationsGenericNBody.md)). For example, to generate the $Z_4$ Read-Rezayi state for 8 particles, we should use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusFermionsWithTranslationsNBodyHollowCore -p 8 -l 12 --nbr-nbody 5 --full-reducedbz --use-lapack*

The spectrum is stored in the file fermions_torus_5body_hollowcore_n_8_2s_12_ratio_1.000000.dat whose content should look like

```text
# Kx Ky E
0 0 6.981082378843e-13
0 0 1.3428147482841e-12
0 0 1.3770501376445
0 0 2.1259921043698
...
```


Using the option --full-reducedbz forces the program to look at all the momentum sectors of the reduced Brillouin zone, i.e. for $k_x,k_y$ going from 0 to ${\rm GCD}(N,N_{\Phi}) - 1$. For a square torus, not using the --full-reducedbz forces to look only at the sectors that are not related by the $C_4$ or inversion symmetries. As usual, all the parameters for the diagonalization (such the number for eigenstates, the threshold between full-diag and Lanczos,..) can be tuned through the options described in [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md) or in Lanczos *(wiki page Lanczos)*.

Computing the interaction matrix elements is time consuming. For that reason the code creates a file  5body_hollowcore_interactioncoefficient_2s_12_ratio_1.0000000000.dat. If this file is in the current directory, any diagonalization for 5-body hollow core with 12 flux quanta will re-read this file instead of generating the interaction matrix elements. If you want to regenerate the interaction file, you can either erase it or add the --regenerate-interactionelements option. 

The n-body hollowcore interaction is normalized such that the typical energy scale when considering n particles interacting on a large torus is one. Note that the normalization has NOT been obtained from an analytical calculation but rather by numerically solving the n particles spectrum on a very large square torus.

A two-body defined through pseudo-potentials can be added via the --interaction-file option. Its usage is similar to the one of [FQHEFermionsTorusWithTranslation](FQHETorusFermionsWithTranslations.md), including the input file syntax. It might be useful to tune the strength of the n-body interaction, this can be done through the option --nbody-strength. By default, this coefficient is set to one.

A linear interpolation between the n-body interaction and any two-body interaction can be achieved thanks to the --linear-interpolation option. Let's say we want to interpolate between the 5-body hollowcore interaction and the $V_1$ pseudo-potentials. First we need to write a pseudo-potential file pseudopotentials_torus_v1.dat

```text
Name = v1
Pseudopotentials = 0.0 1.0
```

Then, to diagonalize the hamiltonian $(1-\lambda) H_{\rm 5b} + \lambda V_1$ , we just use the command

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusFermionsWithTranslationsNBodyHollowCore -p 12 -l 18 --nbr-nbody 5 --use-lapack  --full-diag 2000 --interaction-file pseudopotentials_torus_v1.dat --linear-interpolation --linear-lambda 0.25*

The linear interpolation parameter $\lambda$ is set via --linear-lambda (here $\lambda =0.25$) . The spectrum is stored in fermions_torus_5body_hollowcore_v1_lambda_0.250000_n_12_2s_18_ratio_1.000000.dat , the output file name being created from the type of n-body interaction, the name provided in pseudopotentials_torus_v1.dat and the linear interpolation parameter $\lambda$.
