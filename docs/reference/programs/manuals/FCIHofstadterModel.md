# FCIHofstadterModel — manual

Source: DiagHam wiki page `FCIHofstadterModel`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FCIHofstadterModel](../FCI/FCIHofstadterModel.md)
FCIHofstadterModel deals with interacting particles on the Hofstadter model. It provides support for both bosons and fermions, spinless or spinful, with or without Gutzwiller projection (i.e. hardcore particles).


## Basic usage

As a warm-up, we start with the Hofstadter model on a square lattice. For simplicity, we focus on hardcore bosons. A typical usage is 

*\$PATHTODIAGHAM//build/FTI/src/Programs/FCI/FCIHofstadterModel --real-space -p 2 -x 4 -y 2 -X 1 -Y 3 --use-lapack --hardcore --boson -q 1  --full-diag 5000*

Here we use two hardcore bosons having set the option --boson, --hardcore and -p 2. The lattice itself is defined through the following parameters :

- -x 4 is the number of magnetic unit cell along the x direction (here 4)
- -y 2 is the number of magnetic unit cell along the y direction (here 2)
- -X 3 is the number of sites along the x direction within the magnetic unit cell
- -Y 1 is the number of sites along the y direction within the magnetic unit cell

This means that the physical lattice has 3 * 4 * 2 = 24 sites, 12 along the x direction and 2 along the y direction. Unless explicitly specified through the -q option, the number of flux quanta through one magnetic unit cell is set to one. The flux density is equal to q / (X * Y). **Note that the code is only working  if -Y is set to one, another choice might lead to erroneous results.** 

If the hardcore constraint is removed, we can freely set the strength of the on-site density-density interaction using the --u-potential option.
