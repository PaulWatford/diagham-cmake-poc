# FQHETopInsulatorGetDimension (superseded wiki page)

Source: DiagHam wiki page `FQHETopInsulatorGetDimension`, as of 2026-09-24. Reason: the program is now `FTIGetDimension` and the wiki has a page under that name, which is the attached manual; this older page is kept for comparison.

[FTIGetDimension](../../reference/programs/manuals/FTIGetDimension.md) computes the dimensions of the Hilbert spaces that are involved in all codes related to fractional Chern insulators or fractional topological insulators. The typical usage is 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIGetDimension -p 8 -x 6 -y 4 --no-inversion*

The output shall look like

```text
(kx=0,ky=0) : 30709
(kx=0,ky=1) : 30624
(kx=0,ky=2) : 30704
(kx=0,ky=3) : 30624
(kx=1,ky=0) : 30624
(kx=1,ky=1) : 30624
```


[FTIGetDimension](../../reference/programs/manuals/FTIGetDimension.md) computes the Hilbert space dimension in each momentum sector.

Several options allows to select the system properties

- -s or --nbr-subbands sets the number of subbands (1,2 or 4)
- -p sets the number of particles
- -x, -y, -z, -t sets the number of sites in the x, y, z, and t directions
- --3d and --4d allows to select a 3d or 4d lattice instead of a 2d lattice
- --spin-conserved indicated that the spin is preserved in a 2d and two subband lattice model. In that case, the Hilbert space dimension is given per momentum and per Sz sector
- --no-inversion computes the dimensions, disregarding the inversion symmetry.
