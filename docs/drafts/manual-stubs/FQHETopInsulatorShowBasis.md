# FQHETopInsulatorShowBasis (superseded wiki page)

Source: DiagHam wiki page `FQHETopInsulatorShowBasis`, as of 2026-09-24. Reason: the program is now `FTIShowBasis` and the wiki has a page under that name, which is the attached manual; this older page is kept for comparison.

[FTIShowBasis](FTIShowBasis.md) displays the n-body basis that are involved in all codes related to fractional Chern insulators or fractional topological insulators. The typical usage is 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIShowBasis -p 8 -x 4 -y 3 --kx 0 --ky 0 -s 2 --sz 6 --spin-conserved*

```text
memory requested for Hilbert space = 9ko
memory requested for lookup table = 6Mo
792
[(0,0,-)(0,1,-)(2,1,-)(2,2,-)(3,0,-)(3,1,-)(3,2,+)(3,2,-)]
[(0,0,-)(0,2,-)(2,0,-)(2,2,-)(3,0,-)(3,1,-)(3,2,+)(3,2,-)]
[(0,2,-)(1,1,-)(1,2,-)(2,2,-)(3,0,-)(3,1,-)(3,2,+)(3,2,-)]
```

The basis description is given in such a way only occupied orbitals are written. So [(0,0,-)(0,1,-)(2,1,-)(2,2,-)(3,0,-)(3,1,-)(3,2,+)(3,2,-)] means one electron with momentum (0,0) and spin down, one electron with momentum (0,1) and spin down, ...

Several options allows to select the system properties

- -s or --nbr-subbands sets the number of subbands (1,2 or 4)
- -p sets the number of particles
- -x, -y, -z, -t sets the number of sites in the x, y, z, and t directions
- --kx, --ky, --kz --kt sets the momenta in the x, y, z, and t directions
- --3d and --4d allows to select a 3d or 4d lattice instead of a 2d lattice
- --spin-conserved indicated that the spin is preserved in a 2d and two subband lattice model.
- --sz allows to set  (twice the) Sz value when using the --spin-conserved option

One can combine the description of the Hilbert space with the components of a binary vector. 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIShowBasis -p 8 -x 4 -y 3 --kx 0 --ky 0 -s 2 --sz 6 --spin-conserved --no-autodetect --state myvector.vec*

will display the components of myvector.vec with in front of each of them, the corresponding n-body basis state.
