# HaahCodeEntropy — manual

Source: DiagHam wiki page `HaahCodeEntropy`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [HaahCodeEntropy](../Spin/HaahCodeEntropy.md)
[HaahCodeEntropy](HaahCodeEntropy.md) evaluates the entanglement entropy fro one ground state of the Haah code. Its usage is 

*\$PATHTODIAGHAM/build/Spin/src/Programs/HaahCodeEntropy -x 3 -y 3 -z 3 --nbra-sitex 2 --nbra-sitey 2 --nbra-sitez 1*

The program output should look like

```text
67108864 components generated for the groundstate
Hilbert space dimension for the A region = 256
Building the entanglement matrix for the parity sector 0
Building the reduced density matrix
Diagonalizing the reduced density matrix (64x64)
Building the entanglement matrix for the parity sector 1
Building the reduced density matrix
Diagonalizing the reduced density matrix (64x64)
Building the entanglement matrix for the parity sector 2
Building the reduced density matrix
Diagonalizing the reduced density matrix (64x64)
Building the entanglement matrix for the parity sector 3
Building the reduced density matrix
Diagonalizing the reduced density matrix (64x64)
Trace of the reduced density matrix = 1
Number of non zero eigenvalues for the reduced density matrix = 256
Entangement entropy = 8 * log 2
```


The system size (assuming periodic boundary conditions) is set the options -x, -y and -z. The size subregion for the  entanglement entropy, which is assumed to be a cube, is set by the options --nbra-sitex, --nbra-sitey, and --nbra-sitez.

Instead of a cube for the subregion A, we can define it through a text file that should look like following

```text
0 0 0 1
1 0 0 0
1 1 0 1
0 0 1 0
1 0 1 0
1 0 1 1
0 1 1 1
1 1 1 0
```


Each line describes a spin belonging to the subregion A. The first column is the x coordinate, the second column is the y coordinate, the third column is the z coordinate and the fourth column is the spin index (either 0 or 1). Naming this text file haah_n_8.dat, we just have to run 

*\$PATHTODIAGHAM/build/Spin/src/Programs/HaahCodeEntropy -x 3 -y 3 -z 3 --kept-sites haah_n_8.dat*

[HaahCodeEntropy](HaahCodeEntropy.md) has a specific mode with a lower memory imprint but larger CPU usage. For that purpose, we just have to append the option --low-memory

*\$PATHTODIAGHAM/build/Spin/src/Programs/HaahCodeEntropy -x 3 -y 3 -z 3 --nbra-sitex 2 --nbra-sitey 2 --nbra-sitez 2  --low-memory*

The spectrum of the reduced density can also be exported in an ASCII file by adding the --export-entspectrum option. For example

*\$PATHTODIAGHAM/build/Spin/src/Programs/HaahCodeEntropy -x 3 -y 3 -z 3 --nbra-sitex 2 --nbra-sitey 2 --nbra-sitez 2 --export-entspectrum*

generates a file haahcode_entspectrum_x_3_y_3_z_3_xa_2_ya_2_za_2.dat that should contain

```text
# Trace of the reduced density matrix before normalization = 1
# Number of non zero eigenvalues for the reduced density matrix = 16384
# Entangement entropy = 14.000000000005 * log 2
6.103515625e-05
6.103515625e-05
6.103515625e-05
```
