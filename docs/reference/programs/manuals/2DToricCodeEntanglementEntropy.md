# 2DToricCodeEntanglementEntropy — manual

Source: DiagHam wiki page `2DToricCodeEntanglementEntropy`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [2DToricCodeEntanglementEntropy](../Spin/2DToricCodeEntanglementEntropy.md)
[2DToricCodeEntanglementEntropy](2DToricCodeEntanglementEntropy.md) is the for the two dimensional toric code using the same conventions than [XCubeEntanglementEntropy](XCubeEntanglementEntropy.md). It has the same usage i.e.

\$PATHTODIAGHAM/build/Spin/src/Programs/2DToricCodeEntanglementEntropy -x 3 -y 3 --nbra-sitex 1 --nbra-sitey 1

A typical output should look like

```text
256 intermediate components generated for the groundstate (done in 1e-06s)
256 distinct components generated for the groundstate (done in 1.2e-05s)
Hilbert space dimension for the A region = 4 (done in 2e-06s)
Hilbert space dimension for the region A in  the parity sector = 0 (done in 0s)
Building the entanglement matrix
done in 4e-06s
Building the reduced density matrix
done in 6e-06s
Diagonalizing the reduced density matrix (2x2)
done in 0.000198s
Hilbert space dimension for the region A in  the parity sector = 1 (done in 0s)
Building the entanglement matrix
done in 1.3e-05s
Building the reduced density matrix
done in 2e-06s
Diagonalizing the reduced density matrix (2x2)
done in 3e-06s
Trace of the reduced density matrix before normalization = 1
Number of non zero eigenvalues for the reduced density matrix = 4
Entangement entropy = 2 * log 2
```
