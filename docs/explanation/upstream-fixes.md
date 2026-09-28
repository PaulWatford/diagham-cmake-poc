# Upstream fixes carried on `main`

Purpose: the index of every change this repository makes to DiagHam's own sources, why, and its status with respect to upstream.
Source: the commit trailers `Upstream-Patch:` on `main` (2026-09-26) and `patches/PATCHES.md`, the audit trail kept verbatim from the proof-of-concept era. Changed: this is a new index; the patch files in `patches/` are no longer applied at configure time — each is an ordinary commit.

To send a fix upstream: `docs/upstream-reports/` holds one report and `patches/upstream/` one Subversion-applicable diff per row below, generated from the commits by `scripts_cmake/upstream_pack.py` (`docs/how-to/develop/send-a-fix-upstream.md`).

Rules (from `CONTRIBUTING.md` and `AGENTS.md`): one fix per commit with its
root cause and verification in the message; fixes that change computed
physics are marked **for maintainer review** and go upstream as bug reports;
fixes upstream has since made itself are not carried.

| Class | Commit subject (abridged) | Files | Kind | Upstream status at r4493 |
|---|---|---|---|---|
| A | `FQHECylinderDensity`: gate `LapackDiagonalize` behind `__LAPACK__` (01) | 1 | build failure without LAPACK | still needed |
| A | `FQHECylinderWithSU2SpinDensity` (02) | 1 | same | still needed |
| A | `FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity` (03) | 1 | same, 23 sites | still needed |
| A | `FCIHofstadterModelCompositeFermions` (04) | 1 | same | still needed |
| A | `HubbardSquareLatticeModelJ2S2` (05) | 1 | same | still needed |
| C | `FCIDiceLatticeModel`: use the existing header, not `*New.h` (06) | 1 | include of a never-existing header | still needed; program was not in upstream's build list |
| D | `FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition`: remove undeclared `SubsystemSize` (07) | 1 | copy-paste from the orbital-partition sibling | still needed; not in upstream's build list |
| E | results written at full double precision: `precision(14)` → `max_digits10` (08, regenerated on r4493) | 544 | output precision; 978 sites; deliberate display widths untouched | design change, not a bug fix; proposed upstream |
| E | `FCIHofstadterCorrelation`: same (09; the file is CRLF upstream) | 1 | as 08 | as 08 |
| F | `FermionOnSphere`, `FermionOnTorusWithMagneticTranslations`: public getters (10) | 2 | enables `FCIWannierConstruction` | still needed; not in upstream's build list |
| G | `FCIDiceLatticeModel`: wire `t1/t2/l1/l2`, use the completed SU(2) boson class (11) | 1 | enables the program; physics defaults from the Kagome sibling | **for maintainer review** (defaults; `// IS IT CORRECT?` mapping) |
| H | `ThreeDTwoParticles::PrintState`: chain the prints (12) | 1 | invalid C++11 that broke every QuantumDots program | still needed |
| A | Spin programs: gate 21 `LapackDiagonalize` calls (13, regenerated) + one in `PairHoppingModelFSA` | 9 | same as class A | upstream guarded one file after 2020; the rest still needed |
| I | spinful torus Coulomb: basis order vs `FindStateIndex`; missing up-down interaction (14) | 3 | **physics**: results wrong for N ≥ 3 and for every unpolarised case | **for maintainer review**; applies cleanly to r4493 (files unchanged since 2020) — see `torus-su2-coulomb-defect.md` |
| J | plain Lanczos (`BasicLanczosAlgorithm`, `ComplexBasicLanczosAlgorithm`): detect the closure of the Krylov space and stop instead of iterating on round-off (16) | 4 | **physics**: `-n 1` returned energies below the ground state on degenerate spectra (U30) | still needed |
| J | `ParticleOnCylinderCoulombHamiltonian`: fail clearly without GSL instead of undefined behaviour (17) | 1 | same as U28: a segfault in every no-GSL run of `FQHECylinderFermionsCoulomb` | still needed |
| — | `SpinChainHamiltonianWithTranslations` uninitialised members (was 15) | — | dropped | **fixed upstream after 2020** (r4493 initialises them); not carried |

Also on `main`, not upstream changes: two FQHEOnDisk programs renamed to
resolve install-name collisions (`FQHEDiskBosonsDelta`,
`FQHEDiskLaughlinMonteCarloOverlap`, with their `Makefile.am` entries), and
one `--help` smoke-test exclusion (`EvaluateBroadening`, U26).

How each was found and verified, sibling comparisons and the numbers, are
in `patches/PATCHES.md` (classes A–I), which is kept as the audit trail;
the defect register with provenance is
[../reference/known-defects.md](../reference/known-defects.md).
