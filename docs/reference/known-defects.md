# Known defects in the DiagHam sources — register

Purpose: one row per defect found in the upstream code during the migration: where, what, provenance from the history, what was done, whether upstream has fixed it.
Source: the migration's defect register (2026-09-24/26), built from the r4114 mirror's full history (`git log --follow`), the SVN log to r4493, the pristine autotools build and the CMake build. "First/last touch" are the file's first and last commits in the git history; many first dates are 2006-08-20, the day the tree was split into packages, so real origins can be earlier.

| # | File | Defect | First / last touch, author | Done here | Upstream (r4493) |
|---|---|---|---|---|---|
| U01 | `FQHE/src/Programs/FQHEOnCylinder/FQHECylinderDensity.cc` | `LapackDiagonalize` not gated by `#ifdef __LAPACK__` → fails in the default (no-LAPACK) build | 2016-06-22, regnault | gated (class A) | genuine autotools failure; still present |
| U02 | `…/FQHECylinderWithSU2SpinDensity.cc` | same | 2018-03-22 / 2018-10-23, regnault | gated | still present |
| U03 | `…/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity.cc` | same, 23 sites | 2016-06-19 / 2017-05-26, regnault | gated | still present |
| U04 | `FTI/src/Programs/FCI/FCIHofstadterModelCompositeFermions.cc` | same | 2014-09-18 / 2016-10-27, moller | gated | still present |
| U05 | `FTI/src/Programs/HubbardModels/HubbardSquareLatticeModelJ2S2.cc` | same | 2016-02-12, regnault | gated | still present |
| U06 | `FTI/src/Programs/FCI/FCIDiceLatticeModel.cc` | includes a non-existent `…MomentumSpaceNew.h`; tight-binding parameters unwired; abandoned stub Hilbert class; `// IS IT CORRECT?` | 2019-06-25, moller | fixed (06, 11); **defaults need a maintainer** | never in upstream's `bin_PROGRAMS` |
| U07 | `…/FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition.cc` | undeclared `SubsystemSize` (copy-paste from the orbital-partition sibling) | 2010-06-28 / 2013-02-28, sterdyni | fixed (07) | never in `bin_PROGRAMS` |
| U08 | 544 files at r4493 | `precision(14)` truncates written eigenvalues below double precision | e.g. `src/MainTask/GenericRealMainTask.cc` 2009 / 2020, regnault | 14-digit result sites → `max_digits10` (08, 09); deliberate 4–12-digit display sites left alone | design choice upstream; proposed |
| U09 | `FTI/src/Programs/FCI/FCIWannierConstruction.cc` (+ two Hilbert-space headers) | calls protected members with no public getters | 2012-09-22, regnault | getters added (10) | never in `bin_PROGRAMS` |
| U10 | `FQHE/src/Programs/FQHEOnTorus/QHEFermionsTorusWithSpin.cc` | `#include …;` syntax error, wrong constructor arity, class hierarchy mismatch — dead duplicate of `FQHETorusFermionsWithSpin` | 2006-08-20 / 2009-11-27, moller | excluded from the build (`deferred-code.md`) | never in `bin_PROGRAMS`; recommend deletion |
| U11 | `FQHE/src/HilbertSpace/FermionOnTorusWithSpinNew.cc`, `FQHE/src/Hamiltonian/ParticleOnTorusCoulombWithSpinHamiltonian.cc` | **physics**: spinful torus Coulomb wrong for N ≥ 3 and for every unpolarised case (basis generated in an order `FindStateIndex` cannot search; same-orbital up-down terms dropped) | 2009-12-01 / 2015-11-01, regnault | fixed (14), verified against two independent programs to 1e-12; **for maintainer review** | files unchanged since 2020; upstream develops the `…AndMagneticTranslations` path instead — see `torus-su2-coulomb-defect.md` |
| U12 | `FTI/src/Tools/FTITightBinding/TightBindingModelDiceLattice.cc`, `FTI/src/Hamiltonian/ParticleOnLatticeDiceLattice{Single,Two}BandHamiltonian.cc` | complete classes never listed in any `Makefile.am` | 2017-12-07 / 2019-06-25, moller | recovered by the generator's allow-list | orphaned |
| U13 | `FTI/src/Tools/FTITightBinding/TightBindingModelKapitMueller.cc` | orphaned *because* broken (references a removed `HofstadterSquare` API) | 2017-03-02, moller | not recovered | orphaned |
| U14 | `FQHE/src/HilbertSpace/BosonOnSquareLatticeWithSpinMomentumSpace.cc` | abandoned 2012 draft: `GenerateStates` mostly commented out | 2012-01-13 / 2015-04-28, sterdyni | Dice program rerouted to the completed SU(2) class | dead class still in tree |
| U15 | `FQHE/src/HilbertSpace/Makefile.am:5` → `Fermions.cc` | listed as first source of `libQHEHilbertSpace` but no such file exists | — | generator skips it | harmless Makefile.am entry |
| U16 | `src/Vector/DelocalizedRealVector.cc` | body entirely under `#ifdef USE_CLUSTER_ARCHITECTURE` without including `config.h` → empty object | 2004-06-15 / 2004-12-01, regnault | untouched | dead code |
| U17 | `QuantumDots/src/HilbertSpace/ThreeDTwoParticles.cc:144` | streams an `ostream&` — invalid C++11; blocked every QuantumDots program | 2006-08-20 / 2011-05-23, regnault | fixed (12) | still present |
| U18 | `QuantumDots/src/Programs/Makefile.am` | 2 sources (`ExplicitPeriodic3DQuantumDots`, `VisualPeriodic2D`) unlisted | — | not built | as upstream |
| U19 | 9 Spin programs | ungated `LapackDiagonalize` (as U01) | 2015–2019, regnault/sterdyni | gated (13) | one file guarded upstream after 2020 |
| U20 | `Spin/src/Programs/SUNSpinsOnLatticeCorrelations.cc:18` | `#include "MathTools/h"` typo | 2009-10-06, moller | not built | never in `bin_PROGRAMS` |
| U21 | `FQHE/src/Programs/FQHEOnSphere/Makefile.am` | one program listed twice in `bin_PROGRAMS` | — | generator collapses duplicates | harmless |
| U22 | `configure.ac` | `--enable-anyons` refers to an `Anyons/` package that does not exist | — | not mapped | dead flag |
| U23 | `configure.ac` (r4114) | `--lpthread` typo in the MKL link test | — | — | fixed upstream after 2020 |
| U24 | `src/Programs/TestDiagHamConf.cc` | labels the `__LP64__` branch "Pathscale compiler found" | — | — | cosmetic |
| U25 | documentation | upstream `ChangeLog` stops in 2005; `INSTALL` is the generic GNU text; `README` empty | — | `docs/history/`, this knowledge base | — |
| U26 | `QuantumDots/src/Tools/Analysis/EvaluateBroadening.cc` | tests its required `--input` before honouring `--help`; help unreachable, exit −1 | — | excluded from the `--help` smoke test | usability quirk |
| U27 | `Spin/src/Programs/PairHoppingModelFSA.cc:206` | one `LapackDiagonalize` call still ungated at r4493 (upstream guarded a different site) | 2019, regnault | gated (with 13) | present |
| U28 | `FQHE/src/Hamiltonian/ParticleOnCylinderPseudopotentialHamiltonian.cc` | `PseudopotentialMatrixElement` and `LineChargeMatrixElement` are entirely inside `#ifdef HAVE_GSL`; without GSL the non-void functions fall off the end (undefined behaviour) — every `FQHECylinderFermionsTwoBodyGeneric` run of a build without GSL aborted with "stack smashing detected" | — | `#else` branch: prints that GSL is needed and exits 1 (commit `cylinder-pseudopotential-needs-gsl`) | present at r4493 |
| U30 | `src/LanczosAlgorithm/BasicLanczosAlgorithm*` (plain Lanczos, `-n 1`, no reorthogonalisation) | on a degenerate spectrum the Krylov recurrence breaks down undetected and the returned "ground energy" can lie **below** the true one: `HubbardSquareLatticeModel -p 8 -x 2 -y 4 --u-potential 0 -n 1 --full-diag 100` gives −18.41 (sector 1,2) and −17.87 (0,0); exact is −16.000 (full diagonalisation, `-n 4`, or `--force-reorthogonalize` all agree) | — | open; `known-bug` test `knownbug.hubbard.2x4.U0.plain_lanczos_below_ground_state` (WILL_FAIL) | present at r4493 |
| U29 | `FQHE/src/Programs/FQHEOnDisk/FQHEDiskFermionsTwoBodyGeneric.cc` | hangs after printing `start` for every case tried (`-p 4 --interaction-file v1 --minimum-momentum 18 --maximum-momentum 18`, also Lz 20, 22, with and without `--force-maxmomentum`; the bosonic twin with the same options runs in milliseconds) | — | open; no fermionic disk golden until understood | present at r4493 |
| — | `Spin/src/Hamiltonian/SpinChainHamiltonianWithTranslations.cc` | data constructor left `FastMultiplicationFlag` etc. uninitialised (destructor used them) | 2003 / 2015 | **not carried**: upstream fixed it | **fixed at r4493** |

Physics-affecting entries (U06 defaults, U11) are explained in the
explanation section; everything else is a build or hygiene matter.
