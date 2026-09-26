# FQHETorusFermionsWithTranslationsGenericNBody — manual

Source: DiagHam wiki page `FQHETorusFermionsWithTranslationsGenericNBody`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusFermionsWithTranslationsGenericNBody](../FQHEOnTorus/FQHETorusFermionsWithTranslationsGenericNBody.md)
FQHETorusFermionsWithTranslationsGenericNBody diagonalizes any n-body interaction defined by its Fourrier transform for a fermionic system on the torus geometry. The usage is similar to [FQHEFermionsTorusWithTranslation](FQHETorusFermionsWithTranslations.md). By default, the interaction is set to be the n-body hollowcore interaction. Notice that if you are interested in this specific interaction and you want to focus on rectangular tori, you should better use [FQHETorusFermionsWithTranslationsNBodyHollowCore](FQHETorusFermionsWithTranslationsNBodyHollowCore.md).

If we want to look at the Moore-Read state on the torus, we just have to use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusFermionsWithTranslationsGenericNBody -p 4 -l 8 --nbr-nbody 3 --full-reducedbz --use-lapack*

The generated spectrum fermions_torus_3body_hollowcore_n_4_2s_8_ratio_1.000000.dat will look like 

```text
# Kx Ky E
0 0 0.29992219194192
0 0 0.46812701274052
0 1 0.20197428837541
0 1 0.54404482703396
0 2 -2.2204460492503e-16
0 2 0.3961797599951
0 2 0.79602347618827
...
```

The code also store the interaction matrix elements. In this example, it creates  a file 3body_hollowcore_interactioncoefficient_2s_8_ratio_1.0000000000.dat. If the file is in the current directory, any diagonalization for 3-body hollow core with 8 flux quanta will re-read this file instead of generating the interaction matrix elements. If you want to regenerate the interaction file, you can either erase it or add the --regenerate-interactionelements option. Since generating the interaction matrix elements is highly time consuming (this part of the code is full parallelized, so you should use -S with the maximum number of processors), re-reading the interaction matrix elements from the hard drive spares a lot of time.

By default, the code uses a rectangular torus. This can be changed by setting the --angle option to a non-zero value. This latest option allows to set the angle (in pi units, contrary to [FQHEFermionsTorusWithTranslation](FQHETorusFermionsWithTranslations.md)) between the two vectors 1 and $\tau$ that span the torus. All output file names will contain the angle in radian units. Notice that the length of $\tau$ (or more precisely the ratio of lengths between $\tau$ and 1 is set via the --ratio option.
