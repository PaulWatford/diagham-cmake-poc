# SpinChainMultipleEntanglementSpectra — manual

Source: DiagHam wiki page `SpinChainMultipleEntanglementSpectra`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinChainMultipleEntanglementSpectra](../Spin/SpinChainMultipleEntanglementSpectra.md)

> **Options in this manual that the current program does not list**: `--all-eigenstates`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[SpinChainMultipleEntanglementSpectra](SpinChainMultipleEntanglementSpectra.md) computes the entanglement spectra of many eigenstates stored in a matrix. It takes as an input the matrices produced the spin exact diagonalization codes such as [PeriodicSpinChainJ1J2](PeriodicSpinChainJ1J2.md) when using the --all-eigenstates options. A typical usage is 

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainMultipleEntanglementSpectra --multiple-states spin_1_periodicaklt_n_14_sz_0_invsym_1_szsym_-1_k_0.eigenvec.mat -S --processors 4*

By default, [SpinChainMultipleEntanglementSpectra](SpinChainMultipleEntanglementSpectra.md) computes the entanglement spectrum of all the eigenstates stored in the matrix whose file name is provided by the  --multiple-states  option. The range of computed entanglement spectra can be modified by the --min-multiplestates and --max-multiplestates options. The entanglement spectra are stored in a binary file whose name is deduced from the input file. In the previous example, it will be named spin_1_periodicaklt_n_14_sz_0_invsym_1_szsym_-1_k_0.eigenvec_la_7.full.ent . This file can be processed by a tool such as [SpinChainMultipleEntanglementSpectrumLevelStatistics](SpinChainMultipleEntanglementSpectrumLevelStatistics.md).


## Entanglement spectrum counting


[SpinChainMultipleEntanglementSpectra](SpinChainMultipleEntanglementSpectra.md) can also output the number of non-zero eigenvalues in the reduced density matrix (i.e. the number of levels in the entanglement spectrum). For that purpose, we just have to add the --counting options. The --counting-error allows to set the threshold below which a reduced density matrix eigenvalue is considered to be zero. The number of non-zero eigenvalues will be stored in an ASCII file (here in the above example, named spin_1_periodicaklt_n_14_sz_0_invsym_1_szsym_-1_k_0.eigenvec_la_7.full.ent.counting) That should look like

```text
# index total_counting total_nbr counting(Sza=-14) counting(Sza=-12) counting(Sza=-10) counting(Sza=-8) counting(Sza=-6) counting(Sza=-4) counting(Sza=-2) counting(Sza=0) counting(Sza=2) counting(Sza=4) counting(Sza=6) counting(Sza=8) counting(Sza=10) counting(Sza=12) counting(Sza=14)
0 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
1 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
2 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
3 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
4 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
5 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
6 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
7 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
8 2150 2187 0 6 28 76 160 266 356 366 356 266 160 76 28 6 0
9 2146 2187 0 6 26 76 160 266 356 366 356 266 160 76 26 6 0
...
```

- The first column is the state index.
- The second column is the number of non-zero eigenvalues.
- The third column is the maximum number of eigenvalues for the reduced density matrix at this system size.
- The remaining columns are the number of non-zero eigenvalues per diagonal block (with respect to the quantum numbers, here the spin projection in the region A) of the reduced density matrix.\



## Entanglement spectrum binary file format

The binary file format used to stored the entanglement is rather simple. As usual for DiagHam binary files, little endian encoding is assumed. First the header has the following structure

- 4-bytes int : total number of entanglement spectra.
- 4-bytes int : number of Sza sectors per entanglement spectrum.
- (4-bytes int) * (number of Sza sectors) : list of all twice all Sza values.

Then each entanglement spectrum is stored in the following way

- 4-bytes int : number of levels (nbr_levels) for the current spectrum and the current Sza sector
- (8-bytes double) * (nbr_levels)  : the sorted reduced density matrix eigenvalues (down ordering) for the current spectrum and the current Sza sector

We give below a simple python code to extract a single spectrum whose index is set by EntSpectrumIndex out of the binary file spin_1_periodicaklt_n_8_sz_0_invsym_1_szsym_1_k_0.eigenvec_la_4.full.ent

```text
import numpy as np
import os

InputFile = open('spin_1_periodicaklt_n_8_sz_0_invsym_1_szsym_1_k_0.eigenvec_la_4.full.ent','rb')
EntSpectrumIndex = 10

NbrEntSpectra = np.fromfile(InputFile, '<i4', count=1)[0]
print ('Nbr of entanglement spectra = '+str(NbrEntSpectra))

NbrSzaSectors = np.fromfile(InputFile, '<i4', count=1)[0]
print ('Nbr of Sza sectors = ' + str(NbrSzaSectors))

SzaSectors = np.fromfile(InputFile, '<i4', count=NbrSzaSectors)
print ('Sza sectors = ' + str(SzaSectors))

i = 0
while i < NbrEntSpectra:
    j = 0
    while j < NbrSzaSectors:
        NbrLevels = np.fromfile(InputFile, '<i4', count=1)[0]
        Spectrum = np.fromfile(InputFile, '<d', count=NbrLevels)
        if i == EntSpectrumIndex:
            for Xi in Spectrum:
                print (str(SzaSectors[j])+' '+str(Xi))
        j += 1
    i += 1
InputFile.close()
```
