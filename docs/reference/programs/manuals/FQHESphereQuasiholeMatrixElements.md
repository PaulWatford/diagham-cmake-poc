# FQHESphereQuasiholeMatrixElements — manual

Source: DiagHam wiki page `FQHESphereQuasiholeMatrixElements`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereQuasiholeMatrixElements](../FQHEOnSphere/FQHESphereQuasiholeMatrixElements.md)
[FQHESphereQuasiholeMatrixElements](FQHESphereQuasiholeMatrixElements.md) computes the matrix elements of a simple operator between all quasihole states allowed for a given system size on the sphere geometry. For the time being, the code is limited to the fermionic Laughlin state $\nu = 1/3$ and to the annihilation operator $c_m$. The code computes all the matrix elements 

- $\langle \Psi_{N-1,N_\Phi,L_z-m,\alpha} | c_m |\Psi_{N,N_\Phi,L_z,\beta} \rangle$
- $\langle \Psi_{N,N_\Phi,L_z,\alpha} | c^\dagger_m c_m |\Psi_{N,N_\Phi,L_z,\beta} \rangle$

for any $N,L_z,m$ compatible with system size set by the number of flux quanta $N_\Phi$.

A typical usage of [FQHESphereQuasiholeMatrixElements](FQHESphereQuasiholeMatrixElements.md) is :

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholeMatrixElements -l 6*

where the option -l sets the number of flux quanta. If you also want to fix the number of particles (i.e.  $N$), you should add the -p option. The program displays many information on screen and stores many files on the harddrive. The screen output should look like

```text
memory requested for Hilbert space = 515ko
---------------------------------------
---------------------------------------
processing right states with N=1 LzMax=6 TotalLz=-6
found 1 quasihole states
admissible configurations : 
1 0 0 0 0 0 0 
memory requested for Hilbert space = 515ko
  computing <Psi_{N-1}|c_{-6}|Psi_{N}>
  found 1 quasihole states with N-1=0 LzMax=6 TotalLz=0
  admissible configurations : 
  0 0 0 0 0 0 0 
    <0|c_{-6}|0>=1
memory requested for Hilbert space = 515ko
---------------------------------------
---------------------------------------
processing right states with N=1 LzMax=6 TotalLz=-4
```

It gives for each set of $|\Psi_{N,N_\Phi,L_z,\beta}\rangle$ (right states), the number of quasihole states (the number of values for $\beta$), the corresponding admissible configurations. Then for each $m$, it finds all the $\langle \Psi_{N-1,N_\Phi,L_z-m,\alpha}|$ (left states) and it shows the matrix elements.

The code writes many files on the harddrive. All these files are automatically stored in a directory named nphi_xx_sphere, unless the --disable-directory option is used. Files like fermions_qh_states_k_1_r_2_n_3_nphi_6_lz_0.mat contains the orthonormal basis of the quasihole states. These binary encoded matrices ensure that the same gauge choice is used between different runs and they speeds up the calculations. 

Files like fermions_qh_k_1_r_2_n_3_nphi_9_lz_5_c_9.mat.txt are ASCII column-formatted text files where the matrix elements $\langle \Psi_{N-1,N_\Phi,L_z-m,\alpha} | c_m |\Psi_{N,N_\Phi,L_z,\beta} \rangle$ are stored. The file name provides all the useful information about the quantum numbers i.e.  $N$ (n_x_),  $N_\Phi$ (nphi_x_), $L_z$ (lz_2x) and $m$ (c_2x). The content of the file should look like

```text
0 0 -0.17456694261896
0 1 -0.036194827246396

1 0 0.56507787058479
1 1 -0.14969024371551

2 0 1.3877787807814e-17
2 1 0.80234732877373
```


The first column is the matrix row index (the value of $\alpha$), the second column is the matrix column index (the values of $\beta$) and the third column is the corresponding matrix element. Any of this file has a binary analogue using the same name without the .txt suffix. If you are only interested in the binary files, you can disable the ASCII file output with the option --disable-ascii.

Similarly, files such as fermions_qh_k_1_r_2_n_2_nphi_6_lz_0_cdc_2.mat.txt are ASCII column-formatted text files where the matrix elements $\langle \Psi_{N,N_\Phi,L_z,\alpha} | c^\dagger_m c_m |\Psi_{N,N_\Phi,L_z,\beta} \rangle$ are stored. The file content uses the same abovementioned convention.

In addition to these matrix element files, the description of the full quasihole subspace is stored in the  ASCII column-formatted text file fermions_qh_states_k_1_r_2_nphi_xx.dat where xx is the number of flux quanta. This file should look like

```text
# N_Phi N Lz dim
9 0 0 1
9 1 -9 1
9 1 -7 1
9 1 -5 1
9 1 -3 1
9 1 -1 1
9 1 1 1
9 1 3 1
9 1 5 1
9 1 7 1
9 1 9 1
9 2 -12 1
...
```


The first column is the number of flux quanta, the second column is the number of particles, the third column is twice the total angular momentum along z and the fourth column is the number of quasihole states with these quantum numbers.  Beware that this file and the quasihole state files are **not generated** if they are already present in the current directory.
