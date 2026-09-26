# FTIShowBasis — manual

Source: DiagHam wiki page `FTIShowBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIShowBasis](../FTI/FTIShowBasis.md)
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

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIShowBasis -p 4 -x 4 -y 3 --kx 2 --ky 0 --state fermions_singleband_kagomelattice_n_4_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_2_ky_0.0.vec*

will display the components of fermions_singleband_kagomelattice_n_4_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_2_ky_0.0.vec with in front of each of them, the corresponding n-body basis state.

```text
memory requested for Hilbert space = 516
memory requested for lookup table = 416ko
[(1,0)(3,0)(3,1)(3,2)] : (0.0053276704462768,0.18867981514898)
[(2,1)(2,2)(3,1)(3,2)] : (0.060096590070605,0.15131301227946)
[(0,1)(0,2)(3,1)(3,2)] : (-0.15366471611706,-0.0519641550345)
[(2,0)(2,1)(3,0)(3,2)] : (0.11899279225652,-0.042380798103326)
[(0,0)(0,1)(3,0)(3,2)] : (-0.18149761002157,3.9857609501266e-17)
[(0,0)(1,2)(2,2)(3,2)] : (-0.0040690944040187,-8.6811105224683e-05)
...
```

This example is for spinless fermions. For each line we first have on the left side of the semicolumn the configuration of occupied orbitals ( [(1,0)(3,0)(3,1)(3,2)] means one particle with momentum (kx,ky)=(1,0), ne particle with momentum (kx,ky)=(3,0) ...). On the right side of the semicolumn, we find the corresponding complex coefficient in the binary vector, complex numbers being displayed as (real part,imaginary part).

FTIShowBasis can also detect most of the quantum number sectors from the file name provided through --state. For example

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIShowBasis -s 3 --state  fermions_threeband_u_allowed_orbs_n_8_ns_48_x_4_y_6_kx_0_ky_0.0.vec --allowed-orbitals filter.txt*

will detect all the parameters of a three band calculation up to the list of allowed orbitals that should be provided with the --allowed-orbitals option.
