# Build on the University of Kent cluster (draft)

Purpose: the site-specific recipe for Kent, the planned production home of this repository.
Source: `cmake/KentDefaults.cmake` (a scaffold modelled on Gunnar Möller's `TCMDefaults.cmake` from BDMC_UFL) and `docs/explanation/migration-history.md`. Status: **parked as incomplete (2026-09-26) — the site values are not yet known.** Who can complete it: someone with an account on the Kent cluster, with the administrators for module names and paths. This page states what exists and what is missing so the gaps are visible rather than guessed.

## What exists

`cmake/KentDefaults.cmake` is included when `DIAGHAM_KENT_DEFAULTS=ON`. It is
the one place Kent's paths belong, following the same pattern as the TCM
group's file at Cambridge: a cache flag, off by default, that sets the
site's MKL root, compiler root, MPI wrapper and HDF5 preference so users
need no per-machine flags. Today every path in it is a placeholder and the
file stops with a `FATAL_ERROR` on purpose, so nobody can silently build
against empty paths.

## What is missing (to be filled in with the cluster's administrators)

| Value | Cache variable | TCM equivalent |
|---|---|---|
| MKL installation root | `DIAGHAM_KENT_MKL_ROOT` | `MKL_ROOT=/misc/shared/mkl` |
| Compiler root (if Intel) | `DIAGHAM_KENT_COMPILER_ROOT` | `INTEL_ROOT` |
| MPI C++ wrapper path | `DIAGHAM_KENT_MPI_CXX_COMPILER` | `MPI_CXX_COMPILER=/usr/local/shared/MPI/OpenMPI-…/bin/mpiCC` |
| Module names to load | (documented here) | — |
| Scheduler and `mpirun` conventions | (documented in [../run-mpi.md](../how-to/run-mpi.md)) | — |
| LP64 vs ILP64 MKL interface | decided in `FindMKL`/[intel-mkl.md](../how-to/build/intel-mkl.md) | — |

## Intended use once filled in

```
module load <kent toolchain>
cmake --preset hpc -DDIAGHAM_KENT_DEFAULTS=ON      # or --preset mkl -DDIAGHAM_USE_MPI=ON for the Intel stack
cmake --build --preset hpc -j
ctest --preset hpc
```

Until then, build on Kent as on any cluster: [hpc-cluster.md](../how-to/build/hpc-cluster.md).
