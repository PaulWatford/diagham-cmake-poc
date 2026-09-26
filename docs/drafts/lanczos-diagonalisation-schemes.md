# Diagonalisation schemes: full, Lanczos, Arnoldi (draft)

Purpose: explain how DiagHam chooses between full diagonalisation and Lanczos/Arnoldi, and the projector options for excited states.
Source: upstream wiki page "Lanczos" (as of 2026-09-24). Changed: wikitext converted to Markdown; every option named below checked against the generated program reference (r4493, LAPACK build): all present.
Status: **parked as incomplete** — the wiki page has two empty sections ("Real symmetric Hamiltonians" and "non-symmetric and non-hermitian Hamiltonians"); the "Complex hermitian" part (block Lanczos, projectors) is complete and usable. To complete: describe the real-symmetric path (`--full-diag` threshold, `--use-lapack`, plain Lanczos and its reorthogonalisation options) and the Arnoldi path for non-hermitian matrices; a DiagHam maintainer is the right author. Options in the text are listed under each program's page, e.g. [FQHESphereFermionsTwoBodyGeneric](../reference/programs/FQHEOnSphere/FQHESphereFermionsTwoBodyGeneric.md).

## Introduction

DiagHam uses several diagonalizations schemes. Below a threshold usually set by the option --full-diag, a full diagonalization (using the LAPACK library if the --use-lapack is there, which is HIGHLY recommended) is performed. Above this threshold, Lanczos or Arnoldi algorithms are used. Depending on the Hamiltonian (or more generally matrix) that is considered, we list here all the available options.

## Real symmetric Hamiltonians

## Complex hermitian Hamiltonians

If more than one eigenvalue are needed and the Hilbert space requires to use the Lanczos algorithm.

They are several options :

### Block Lanczos

The safest way is to use the options --block-lanczos --block-size x --force-reorthogonalize . Here x should match the number of eigenvalues set by the -n option. This technique guarantees that if there is less than x degenerate eigenstates, no spurious degeneracy will appear.

### Lanczos with projectors

$H \longrightarrow H + \alpha \sum_{i=1}^N | \Psi_i \rangle \langle \Psi_i |$

The value of $\alpha$ is set via the --addprojector-factor option (its default value is 10.0). 

The first run can be done by requiring a single state (with or without --fast-disk) in order to get the groundstate. Let denote eigenvector.0.vec the groundstate. To get the first excited state, you need to run the code with the following Lanczos options

- -n 1 --eigenstate --fast-disk  --add-projector projector.dat*

where projector.dat is a single column ASCII file that contains the list of the $| \Psi_i \rangle$ vectors. Here this file should be 

```text
eigenvector.0.vec
```

The output file name for the new eigenstate won't collide with the previously computed eigenstates : its index is deduced from the number of vectors provided in projector.dat (this feature can be disabled using the --addprojector-noshift option). Once the first excited state eigenvector.1.vec has been obtained, we can get the second excited states by modifying the projector.dat file to look like

```text
eigenvector.0.vec
eigenvector.1.vec
```

and rerunning the same command.


The option --auto-addprojector automatically builds the projector subspace. It can be used with or without the --fast-disk option. It requires to ask more than one eigenvalue. Indeed, it requires the typical following Lanczos options (for the three lowest eigenstates)

- -n 3 --auto-addprojector*

A typical run will look like 

```text
-19.358395052191 -18.589327180734 7.4967711046297e-09 
-19.358395052192 -18.589327283999 5.5550781533772e-09 
-19.358395052192 -18.589327343052 3.1767237258694e-09 
found eigenstate 0
computing eigenstate 0
69/70           
reinitializing Lanczos 
starting Lanczos 
-19.358395052192 8.7386994940181 1.4700922918163 
-19.358395052192 -1.0216147898183 1.1169069597276 
-19.358395052192 -4.175599016053 3.0872538824497 
... 
-19.358395052192 -18.58932733565 3.8959799201819e-09 
-19.358395052192 -18.589327370544 1.8770670799477e-09 
-19.358395052192 -18.589327386414 8.5375352996193e-10 
found eigenstate 1
computing eigenstate 1
66/67           
checking orthogonality with previously found eigenstates
< 0 | 1 > = (7.5059372186279e-15,3.4867320489148e-14)
reinitializing Lanczos 
starting Lanczos 
-19.358395052192 -2.7151489474257 0.85394044168537 
-19.358395052192 -9.7166469138983 2.5786791450653 
-19.358395052192 -12.151705731806 0.25060690580667 
...
-19.358395052192 -18.589327398957 2.6565092472415e-14 
-19.358395052192 -18.589327398957 1.1658062164153e-14 
-19.358395052192 -18.589327398957 6.1157047418507e-15 

-19.358395052192 -18.589327398957 6.1157047418507e-15  Nbr of iterations = 209
```


Each time an eigenvalue is found, the eigenstate is added to the projector subspace. A sanity check is performed by computing the overlap which each state of the subspace. Note that the two options --auto-addprojector and --add-projector can be combined in case some eigenstates were already obtained.

## non-symmetric and non-hermitian Hamiltonians
