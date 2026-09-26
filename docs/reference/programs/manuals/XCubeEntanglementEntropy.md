# XCubeEntanglementEntropy — manual

Source: DiagHam wiki page `XCubeEntanglementEntropy`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [XCubeEntanglementEntropy](../Spin/XCubeEntanglementEntropy.md)
[XCubeEntanglementEntropy](XCubeEntanglementEntropy.md) is the analogue of [HaahCodeEntropy](HaahCodeEntropy.md) for the X-cube model. It has the same usage i.e.

*\$PATHTODIAGHAM/build/Spin/src/Programs/XCubeEntanglementEntropy -x 3 -y 2 -z 2 --nbra-sitex 2 --nbra-sitey 2 --nbra-sitez 1*

which should produce this type of output 
```text
         
product of all the type 1 z term is equal to the identity
8388608 intermediate components generated for the groundstate (done in 0.020914s)
262144 distinct components generated for the groundstate (done in 0.499167s)
Hilbert space dimension for the A region = 1024 (done in 0.00495s)
Hilbert space dimension for the region A in  the parity sector = 0 (done in 1.9e-05s)
Building the entanglement matrix
done in 0.010192s
Building the reduced density matrix
done in 0.002265s
Diagonalizing the reduced density matrix (512x512)
done in 0.043355s
Hilbert space dimension for the region A in  the parity sector = 1 (done in 3e-06s)
Building the entanglement matrix
done in 0.009706s
Building the reduced density matrix
done in 0.002181s
Diagonalizing the reduced density matrix (512x512)
done in 0.042574s
Trace of the reduced density matrix before normalization = 1
Number of non zero eigenvalues for the reduced density matrix = 636
Entangement entropy = 8 * log 2
```
