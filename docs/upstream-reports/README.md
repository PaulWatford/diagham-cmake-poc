# Reports for the DiagHam authors

Purpose: everything needed to offer this repository's fixes and findings upstream — one report and one `svn patch`-applicable diff per fix, and a report per open defect.
Source: generated from the commits on `main` that carry an `Upstream-Patch:` trailer (`scripts_cmake/upstream_pack.py`; `--check` applies the diffs in series to a pristine `upstream` tree); the open-defect reports are hand-written from `docs/reference/known-defects.md`. Nothing here is sent automatically: `docs/how-to/develop/send-a-fix-upstream.md` says how.

## Fixes carried on `main` (diffs in `patches/upstream/`, a series to apply in order)

| # | Fix | Files | Diff |
|---|---|---|---|
| 01 | [FQHECylinderDensity: gate LapackDiagonalize behind __LAPACK__ (upstream build failure with](01-FQHECylinderDensity.md) | 1 | `patches/upstream/01-FQHECylinderDensity.diff` |
| 02 | [FQHECylinderWithSU2SpinDensity: gate LapackDiagonalize behind __LAPACK__](02-FQHECylinderWithSU2SpinDensity.md) | 1 | `patches/upstream/02-FQHECylinderWithSU2SpinDensity.diff` |
| 03 | [FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity: gate LapackDiagonalize behind __L](03-FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity.md) | 1 | `patches/upstream/03-FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity.diff` |
| 04 | [FCIHofstadterModelCompositeFermions: gate LapackDiagonalize behind __LAPACK__](04-FCIHofstadterModelCompositeFermions.md) | 1 | `patches/upstream/04-FCIHofstadterModelCompositeFermions.diff` |
| 05 | [HubbardSquareLatticeModelJ2S2: gate LapackDiagonalize behind __LAPACK__](05-HubbardSquareLatticeModelJ2S2.md) | 1 | `patches/upstream/05-HubbardSquareLatticeModelJ2S2.diff` |
| 06 | [FCIDiceLatticeModel: use the existing BosonOnSquareLatticeWithSpinMomentumSpace header (ne](06-FCIDiceLatticeModel.md) | 1 | `patches/upstream/06-FCIDiceLatticeModel.diff` |
| 07 | [FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition: remove undeclared Subsyste](07-FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition.md) | 1 | `patches/upstream/07-FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition.diff` |
| 08 | [FermionOnSphere, FermionOnTorusWithMagneticTranslations: public getters needed by FCIWanni](08-FCIWannierConstructionEnablement.md) | 2 | `patches/upstream/08-FCIWannierConstructionEnablement.diff` |
| 09 | [FCIDiceLatticeModel: wire t1/t2/l1/l2 tight-binding options and use the completed SU2 boso](09-FCIDiceLatticeModelEnablement.md) | 1 | `patches/upstream/09-FCIDiceLatticeModelEnablement.diff` |
| 10 | [FCIHofstadterCorrelation: write results at full double precision (max_digits10); file is C](10-FCIHofstadterCorrelationPrecision.md) | 1 | `patches/upstream/10-FCIHofstadterCorrelationPrecision.diff` |
| 11 | [ThreeDTwoParticles::PrintState: chain the sub-state prints (invalid ostream << ostream& br](11-QuantumDotsThreeDTwoParticlesPrintState.md) | 1 | `patches/upstream/11-QuantumDotsThreeDTwoParticlesPrintState.diff` |
| 12 | [Spinful torus Coulomb: generate the basis in the order FindStateIndex searches, and add th](12-TorusSU2CoulombBasisOrderAndInterSpinInteraction.md) | 3 | `patches/upstream/12-TorusSU2CoulombBasisOrderAndInterSpinInteraction.diff` |
| 13 | [Spin programs: gate LapackDiagonalize calls behind __LAPACK__ (patch 13 regenerated agains](13-SpinLapackDiagonalizeGating.md) | 8 | `patches/upstream/13-SpinLapackDiagonalizeGating.diff` |
| 14 | [Write computed results at full double precision: precision(14) -> numeric_limits<double>::](14-IEEE754PrecisionForEigenvalueOutput.md) | 544 | `patches/upstream/14-IEEE754PrecisionForEigenvalueOutput.diff` |
| 15 | [ParticleOnCylinderPseudopotentialHamiltonian: fail clearly without GSL instead of undefine](15-cylinder-pseudopotential-needs-gsl.md) | 1 | `patches/upstream/15-cylinder-pseudopotential-needs-gsl.diff` |
| 16 | [Plain Lanczos: detect the closure of the Krylov space instead of iterating on round-off (U](16-LanczosKrylovSpaceClosure.md) | 4 | `patches/upstream/16-LanczosKrylovSpaceClosure.diff` |
| 17 | [ParticleOnCylinderCoulombHamiltonian: fail clearly without GSL instead of undefined behavi](17-cylinder-coulomb-needs-gsl.md) | 1 | `patches/upstream/17-cylinder-coulomb-needs-gsl.diff` |

## Registered defects without a fix (reproducers)

- [Open defect U29: FQHEDiskFermionsTwoBodyGeneric hangs after printing "start"](open-U29.md)
- [Open defect U31: FCICheckerboardLatticeModel segfaults without --single-band](open-U31.md)
- [Open defect U34: the disk pseudopotential Hamiltonian gives wrong energies for N ≥ 3](open-U34.md)
- [Open defect U35: TwoDimensionalTransverseFieldIsingModel ignores open boundaries for the Ising bonds](open-U35.md)
- [Open defect U36: FTI3DHopf segfaults after writing its band file](open-U36.md)
- [Open defect U37: 3D --export-onebodytext writes only the kz = 0 slice](open-U37.md)
- [Open defect U39: the one-body energy of the single-band many-body Hamiltonians is wrong](open-U39.md)
- [Open defect U40: SpinChainEntanglementEntropy aborts in a build without LAPACK](open-U40.md)
- [Open defect U41: 2DToricCodeEntanglementEntropy loses its spectrum when --kept-sites has a path](open-U41.md)

Generated from 17 fix commits; regenerate with `python3 scripts_cmake/upstream_pack.py`.
