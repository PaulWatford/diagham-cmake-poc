# FCIGenericModelWithSpin — manual

Source: DiagHam wiki page `FCIGenericModelWithSpin`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FCIGenericModelWithSpin](../FCI/FCIGenericModelWithSpin.md)
[FCIGenericModelWithSpin](FCIGenericModelWithSpin.md) handles interacting fermions with spin on a lattice, with a generic two-body interaction.
For now, the one-body terms of the Hamiltonian are set to those of a Checkerboard lattice for both spins, but the code can support a different tight-binding model for each spin. The typical usage is 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCIGenericModelWithSpin --use-lapack -p 4 -x 3 -y 3 --vupup-potential 0 --vdowndown-potential 0 --vupdown-potential 0 --u-potential 1 --real-space --checkerboard --fixed-sz --sz-value 0 --full-diag 3000*

The -p option sets the number of fermions, -x (resp. -y) the number of sites in the x (resp. y) direction. 

The strength of the interaction can be tuned with the options *--u-potential* (on-site interaction), *--vupdown-potential* (nearest neighbor interaction between fermions of opposite spins), *--vupup-potential* and *--vdowndown-potential* (nearest neighbor interaction between fermions with the same spin). The relative strength of the nearest neighbor terms is defined such that the interaction is SU(2) invariant when all three *v* terms take the same value.

The amplitude of the hopping terms can be tuned using the *--t1* (nearest neighbor), *--t2* (next nearest neighbor), *--tpp* (next next nearest neighbor). The default parameters are chosen to optimize the flatness ratio of the one-body spectrum.

The total magnetization Sz is a good quantum number. The *--fixed-sz* option has to be activated to take advantage of that symmetry. The value of Sz is set using the *--sz-value* option. The translation symmetries are also used (unless the *no-translation* option is activated). All momentum sectors are computed, unless otherwise specified (using the *--only-kx* and *--only-ky* options).

The gutzwiller projection (infinite on-site interaction) is activated with the option *--gutzwiller*. Activating this option significantly reduces the dimension of the Hilbert space.

An optional $S^2$ term can be added to the Hamiltonian. The prefactor in front of this term is set via the --s2-factor option. Its default value is 0 meaning the term is not added. Note that adding such a term might require additional memory for the precalculations (i.e. if you want to store the full Hamiltonian in RAM). The amount of memory dedicated for these additional precalculations can bet set in two different ways:
- it can be set using the --s2-memory option. This option is separated from the usual --memory option. In the standard output, the second reduction factor is related to the precalculations needed for the  $S^2$ term.
- it can define as a ratio of the amount of memory specified either via the --memory option or the MPI cluster description. In that case, the  amount of memory of the actual hamiltonian (i.e. the non $S^2$ part) is reduced by the amount of memory devoted to the $S^2$ part.
