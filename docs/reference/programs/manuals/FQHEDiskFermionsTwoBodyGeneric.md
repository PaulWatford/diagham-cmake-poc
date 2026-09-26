# FQHEDiskFermionsTwoBodyGeneric — manual

Source: DiagHam wiki page `FQHEDiskFermionsTwoBodyGeneric`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHEDiskFermionsTwoBodyGeneric](../FQHEOnDisk/FQHEDiskFermionsTwoBodyGeneric.md)
FQHEDiskFermionsTwoBodyGeneric provides exact diagonalization for fermions on a disk geometry and any generic two body interaction described by its pseudo-potentials. Its usage has many similarities with [FQHESphereFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md).

## Simple usage

By default the fermionic Hilbert space is only constraint by the number of particles and the system total angular momentum. The range of this latest can be set by the two options --maximum-momentum ( or -l) and --minimum-momentum. Unless specified, this latest is automatically set to the lowest possible one for the current number of particles. 

Let's consider the $V_1$ pseudopotential and a system with 5 electrons. We just have to run

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnDisk/FQHEDiskFermionsTwoBodyGeneric -p 5 --interaction-file pseudopotentials_v1.dat --interaction-name v1  --maximum-momentum 32 --use-lapack*

where the pseudo-potential file *pseudopotentials_v1.dat* is

```text
Pseudopotentials = 0.0 1.0
```

Note that only the pseudo-potentials from the 0 to the last non-zero one have to be provided. The energy spectrum is stored in *fermions_disk_v1_n_5_lz_32.dat* and should look like

```text
# Lz E
10 1.6164173907771
11 1.6164173907771
12 1.5293795312737
12 1.6164173907771
...
29 1.4503079006973
29 1.5293795312737
29 1.6164173907771
30 -3.441691376338e-15
30 0.13343225774049
30 0.13670167629347
30 0.13977592234891
...
```

The first zero energy state occurs at total angular momentum $L_z=3N(N-1)/2=30$ as expected for the Laughlin 1/3 state.

The number of orbitals can also be truncated by using the --force-maxmomentum option if it is set to a zero or positive (by default, it is set to -1).

A confining, i.e., one body, potential can be added by adding a line "Onebodypotentials = " in the interaction file. Note that the number of terms should match the number of orbitals either provided in the command line (--force-maxmomentum option) or from the maximal total angular momentum.
