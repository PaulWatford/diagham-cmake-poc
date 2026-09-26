# HubbardDensity — manual

Source: DiagHam wiki page `HubbardDensity`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [HubbardDensity](../HubbardModels/HubbardDensity.md)
HubbardDensity computes the density matrix $\langle\Psi_{left}|c^\dagger_{i,\sigma} c_{j,\sigma'}|\Psi_{right}\rangle$. Its usage is identical to the one of [HubbardSuperconductorOrderParameter](HubbardSuperconductorOrderParameter.md). If more than one state are provided, HubbardDensity will then compute $\hat{P}_{N+2}c^\dagger_{i,\sigma} c_{j,\sigma'} \hat{P}_{N} c^\dagger_{j,\sigma'} c_{i,\sigma}$ as defined in [HubbardSuperconductorOrderParameter](HubbardSuperconductorOrderParameter.md).
