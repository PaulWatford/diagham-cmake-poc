# Building DiagHam — which guide?

Purpose: send you to the one build guide for your platform.
Source: this index is new (2026-09-25); the guides it points to are adapted from the upstream wiki "Install" page (as of 2026-09-24) and the earlier `HPC.md`.

| You have | Read |
|---|---|
| Linux, GCC (the default and the configuration CI tests) | [linux-gcc.md](linux-gcc.md) |
| Linux, Clang / LLVM | [linux-clang.md](linux-clang.md) |
| Intel compilers + MKL (+ Intel MPI) | [intel-mkl.md](intel-mkl.md) |
| Apple silicon Mac (M1/M2/…) | [macos.md](macos.md) |
| An HPC cluster with a module system (MPI, ScaLAPACK) | [hpc-cluster.md](hpc-cluster.md) |
| The University of Kent cluster | as any cluster, [hpc-cluster.md](hpc-cluster.md); the site recipe is a [draft](../../drafts/kent-cluster-build.md) until the site values are known |

Every guide uses the same three commands with a different **preset**; the
presets and what each turns on are in [presets.md](presets.md). The
one-to-one map from the old `configure` flags to CMake options is in
[../../reference/configure-flag-map.md](../../reference/configure-flag-map.md).

Common to all platforms:

- CMake ≥ 3.21, a C++11 compiler, pthreads, Python 3 (the generator and the
  Python cross-check test; numpy for the latter).
- Build out of source: the presets put the build tree in `build/<preset>`;
  an in-source configure is refused.
- After building, run `ctest --preset <preset>` — the physics goldens take
  seconds and tell you the binaries are right, not just linked.
- To check what the build actually enabled, run `build/<preset>/src/Programs/TestDiagHamConf`
  (see [../../tutorials/first-build-and-verify.md](../../tutorials/first-build-and-verify.md)).
