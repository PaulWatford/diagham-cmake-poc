# FQHETorusFermionsWithTranslations — manual

Source: DiagHam wiki page `FQHETorusFermionsWithTranslations`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusFermionsWithTranslations](../FQHEOnTorus/FQHETorusFermionsWithTranslations.md)
FQHETorusFermionsWithTranslations provides exact diagonalization for fermions on a torus geometry and any generic two body interaction. It fully supports magnetic translations in both directions. On many aspect, this program is similar to [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md). 

## Defining the interaction

The typical usage is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusFermionsWithTranslations -m 500 -p 4 -l 12 --use-lapack --interaction-file pseudopotentials.dat --ratio 1 -S --processors 4 --all-points*

pseudopotentials.dat is the two body pseudopotential file that looks like

```text
Name = laughlin
Pseudopotentials = 0.0 1.0
```


Notice that the pseudopotential file is a bit different from the one used on the sphere geometry. The number of pseudopotentials can be lower than the number of flux quanta (they are assumed to be zero if the sequence is truncated). The example shown above is thus valid to generate the Laughlin state for any number of flux quanta. The interaction name is directly provided through the *Name* parameter. If you want to use the Coulomb interaction, then *Name* and *Pseudopotentials* are not needed. Just use

```text
CoulombLandauLevel = 0
```

as the two body pseudopotential file, *ColombLandauLevel* could take any integer value (0 being the lowest Landau level).

The torus aspect ration is set through the --ratio option. The Wigner crystal energy contribution can be added using the --add-wigner option. --all-points allows to compute the spectrum in all momentum sectors and not only those in the first Brillouin zone.

The spectrum will be stored in a file *fermions_torus_laughlin_n_4_2s_12_ratio_1.000000.dat* using a similar naming convention (up to the additional aspect ratio information). It should look like

```text
# Kx Ky E
0 0 0.36097032062062 66
0 0 1.0484233056592 66
0 0 1.5109225013226 66
0 0 1.7828767452615 66
0 0 1.9911289164451 66
0 0 2.0038885272665 66
0 0 2.3791707154688 66
0 0 2.5280002238829 66
0 0 2.5434476154259 66
0 0 3.1279121773126 66
0 0 3.6171668147302 66
0 1 0.39267491023013 66
0 1 1.230957132042 66
```

The first and second columns give the momentum sector and the third one is the energy (the fourth column is for debugging purposes).

By default, the code uses a rectangular torus. This can be changed by setting the --angle option to a non-zero value. This latest option allows to set the angle (in radian units) between the two vectors 1 and $\tau$ that span the torus. Notice that the length of $\tau$ (or more precisely the ratio of lengths between $\tau$ and 1 is set via the --ratio option.

## Diagonalization options

The --eigenstate option requires to compute the eigenstates and the -n option indicates the number of wanted eigenstates. The eigenstates are binary vectors that can be managed with the corresponding tools *(wiki page Binary vectors)*. Their decomposition onto the n-body basis can be read using [FQHETorusShowBasis](FQHETorusShowBasis.md).

If more than one eigenvalue or eigenstate per momentum are needed and the Hilbert space requires to use the Lanczos algorithm, it is safer to use the options --block-lanczos --block-size x --force-reorthogonalize . Here x should match the number of eigenvalues set by the -n option.

## Mass anisotropy

A mass anisotropy term can be included within the interaction. This can be turn on using the two options --mass-anisotropy and --anisotropy where this latest is the anisotropy parameter $\alpha$ defined such that $q_g^2 = \alpha q_x^2 + q_y^2 / \alpha$.

## Other interactions

The double gated coulomb interaction, as defined in [FQHETorusFermionsWithSpinAndTranslations](FQHETorusFermionsWithSpinAndTranslations.md), is implemented. It can be activated either as command line option through --double-gate or within the interaction file using the same syntax as [FQHETorusFermionsWithSpinAndTranslations](FQHETorusFermionsWithSpinAndTranslations.md). Beware that taking the gate distance to infinity matches the bare coulomb up to a global shift due to the absence of zero mode in the bare coulomb.
