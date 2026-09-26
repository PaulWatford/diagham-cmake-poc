# SpinSystemSzToS — manual

Source: DiagHam wiki page `SpinSystemSzToS`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinSystemSzToS](../Spin/SpinSystemSzToS.md)
SpinSystemSzToS converts a SU(2) invariant spectrum from a $S_z$ basis to a $S$. For example, if you compute a spectrum of the $J_1-J_2$ model using [PeriodicSpinChainJ1J2](PeriodicSpinChainJ1J2.md), we convert the spectrum spin_1_2_periodicj1j2_j1_1.000000_j2_0.500000_n_12.dat using

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinSystemSzToS -s spin_1_2_periodicj1j2_j1_1.000000_j2_0.500000_n_12.dat -c 2 -k 1 -z 0*


This creates an output file spin_1_2_periodicj1j2_j1_1.000000_j2_0.500000_n_12_s.dat that should look like

```text
# 2S K Energy deg.
0 0 -4.5 1
0 0 -4.0372851831344 1
0 0 -2.885900915068 1
0 0 -2.6363155510329 1
0 0 -2.0402880262247 1
0 0 -1.8940110920444 1
```


- The first column is twice the $S$ value.
- The second column is the momentum.
- The third column is the energy.
- The fourth column is the SU(2) degeneracy.
