# Drafts: incomplete documentation, and where retired pages went

Purpose: one place that lists every page that is not finished (what is missing and who can complete it) and every page that was removed, merged or moved during the 2026-09 documentation reorganisation, so nothing is lost silently.
Source: new (2026-09-26). Rule (from `CONTRIBUTING.md`): a page leaves this folder when its provenance line no longer says "parked" and it has been checked against the r4493 build.

## Incomplete pages parked here

| Page | What is missing | Who can complete it | Parked |
|---|---|---|---|
| [kent-cluster-build.md](kent-cluster-build.md) | the Kent cluster's MKL/compiler/MPI paths and module names | a Kent cluster user with the administrators | 2026-09-26 |
| [lanczos-diagonalisation-schemes.md](lanczos-diagonalisation-schemes.md) | the real-symmetric and non-hermitian sections (empty on the wiki) | a DiagHam maintainer | 2026-09-26 |
| [manual-stubs/](manual-stubs/INDEX.md) (25 wiki pages: 14 stubs under 200 bytes, 1 with no matching program or an uncertain rename, 10 superseded old-name pages) | purpose, inputs, outputs, an example | any contributor who has run the program; the rename questions need a maintainer | 2026-09-25 |
| [programs-without-manual.md](programs-without-manual.md) (471 of 603 programs) | a manual; the generated `--help` page exists for each | any contributor who has run the program | 2026-09-26 |

## Upstream wiki pages: where each went

Every non-program page of the DiagHam wiki (archived 2026-09-24); program pages are covered by the manuals and the stubs index above.

| Wiki page | Action | Where |
|---|---|---|
| Main_Page | retired | hand-kept index; superseded by `docs/README.md` (the documentation index) and the generated [program reference](../reference/programs/README.md) |
| FQHE_programs, FCI_programs, Spin_chain_programs | retired | hand-kept program lists; superseded by the generated [program reference](../reference/programs/README.md); entries on the wiki lists with no program in the r4493 build: FQHETorusBosonsDelta |
| Install | split | [how-to/build/](../how-to/build/README.md): linux-gcc, linux-clang, intel-mkl, macos (svn checkout → git clone; configure → presets) |
| MPI | updated | [how-to/run-mpi.md](../how-to/run-mpi.md) |
| Scalapack | updated | [how-to/run-scalapack.md](../how-to/run-scalapack.md) |
| Binary_vectors, Binary_matrices | merged | [reference/data-formats.md](../reference/data-formats.md), verified against the r4493 source |
| Laughlin_state | merged | [tutorials/first-fqhe-run.md](../tutorials/first-fqhe-run.md), re-run at r4493 |
| Lanczos | parked | [lanczos-diagonalisation-schemes.md](lanczos-diagonalisation-schemes.md) — two sections empty on the wiki |
| Articles | updated | [reference/publications.md](../reference/publications.md); also the source of `CITATION.cff` |
| People | merged | `AUTHORS` at the root (names only) |
| Create_new_FQHE_code | updated | [how-to/develop/create-a-program.md](../how-to/develop/create-a-program.md) |
| Create_a_new_spin_Hilbert_space_with_a_discrete_symmetry | updated | [how-to/develop/new-spin-hilbert-space-with-symmetry.md](../how-to/develop/new-spin-hilbert-space-with-symmetry.md) |
| Long_double | updated | [how-to/develop/long-double.md](../how-to/develop/long-double.md); its two attachments were not archived |
| Logo | retired | 96-byte image stub; nothing to carry |

## Repository documents: removed, merged or moved

| Was | Action | Now |
|---|---|---|
| upstream `README` (0 bytes) | retired | empty file; `README.md` is the README (prompt 2) |
| `INSTALL`, `ChangeLog`, `TODO`, `NEWS` (upstream root) | history | [docs/history/](../history/README.md), verbatim with a note on what superseded each |
| `HPC.md` | merged | [how-to/build/hpc-cluster.md](../how-to/build/hpc-cluster.md) and the MPI/ScaLAPACK how-tos |
| `CONFIGURE_FLAGS.md` | moved | [reference/configure-flag-map.md](../reference/configure-flag-map.md) |
| `MODULE_MAP.md` | merged | [reference/build-system.md](../reference/build-system.md) |
| `HUBBARD_BENCHMARK.md` | merged | [tutorials/hubbard-walkthrough.md](../tutorials/hubbard-walkthrough.md); the full log stays at `benchmarks/BENCHMARK.md` |
| `TESTING.md` | moved | [reference/tests.md](../reference/tests.md); its "Adding a golden" section is [how-to/develop/add-a-golden-test.md](../how-to/develop/add-a-golden-test.md) |
| `MIGRATION_ROADMAP.md`, `DEFERRED.md`, `BUG_torus_su2_coulomb.md` | moved | [explanation/](../explanation/migration-history.md): migration-history, deferred-code, torus-su2-coulomb-defect (+ [reference/known-defects.md](../reference/known-defects.md)) |
| `patches/PATCHES.md` | kept, frozen | audit trail of the proof-of-concept fixes; the live index is [explanation/upstream-fixes.md](../explanation/upstream-fixes.md) |
| `docs/how-to/build/kent.md` | parked | [kent-cluster-build.md](kent-cluster-build.md) — site values unknown |
| `COPYING` and `LICENSE` | renamed | `COPYING` (GPL-2 text) → `LICENSE`; the contribution notice → `NOTICE` (2026-09-27, so the licence is machine-detected) |
| `install-sh`, `mkinstalldirs` | kept | autotools helper scripts used by the upstream `Makefile.in` files; not documentation |
