# Full diagonalisation with ScaLAPACK

Purpose: use ScaLAPACK for full diagonalisation on a parallel machine.
Source: upstream wiki "Scalapack" page (as of 2026-09-24). Changed: build instructions point to the CMake guides; the example's program path was updated (the wiki's `FQHETopInsulator/FQHECheckerboardLatticeModel` is now `FTI/src/Programs/FCI/FCICheckerboardLatticeModel`); `--use-scalapack` verified present in 113 program sources at r4493.

DiagHam can rely on [ScaLAPACK](http://www.netlib.org/scalapack) to perform
full diagonalisation on a parallel machine. It is a build-time choice that
requires MPI and LAPACK (`--enable-mpi --enable-scalapack --enable-lapack`):
the `hpc` preset turns all three on, or add `-DDIAGHAM_USE_SCALAPACK=ON` to
an MPI build — see [build/hpc-cluster.md](build/hpc-cluster.md). If the
library is not found automatically, give the link line with
`-DDIAGHAM_SCALAPACK_LIBS="..."` (the old `--with-scalapack-libs`).

Programs that support it show this line in their `--help`:

```
--use-scalapack : use SCALAPACK libraries instead of DiagHam or LAPACK libraries
```

Not every program does: if `--use-scalapack` is absent from a program's
help, that program does not support it (113 of the 603 programs do at
r4493).

Typical use — `--mpi` or `--mpi-smp` activates MPI, `--use-scalapack` selects
ScaLAPACK for the full diagonalisation:

```
mpirun --hostfile hostfile \
  $DIAGHAM/build/hpc/FTI/src/Programs/FCI/FCICheckerboardLatticeModel \
  -p 6 -x 6 -y 3 --single-band --memory 0 --flat-band --use-scalapack \
  --only-ky 0 --only-kx 0 --mpi-smp cluster.desc --cluster-profil cluster.log --full-diag 5000
```

The cluster-description file is the one described in [run-mpi.md](run-mpi.md).

## Common issues and limitations

Some MPI implementations limit the size of a single transfer, which can
crash a ScaLAPACK run. The wiki's experience: without eigenstates, real
symmetric matrices up to 65536×65536 and complex Hermitian matrices up to
32768×32768 gave no trouble; with eigenstates, make sure the eigenstate
matrix per slot stays under 2 GB by using enough slots/nodes. Sizes reached:
81,828 (real symmetric, with eigenstates, 32 MPI processes, 105 h wall /
8,470 h CPU on 16-core Xeon E5-2630 nodes) and 75,910 (complex Hermitian,
no eigenstates, 64 processes over 3 nodes, 11 h).

With Intel MPI the ScaLAPACK link can fail with

```
libscalapack.a(pzhetd2.o): undefined reference to `zhemv_' / `zher2_'
```

The fix was `--with-scalapack-libs="-lscalapack -llapack -lblas"`, i.e.
`-DDIAGHAM_SCALAPACK_LIBS="-lscalapack -llapack -lblas"` here, on top of the
Intel-MPI `-DMPICH_IGNORE_CXX_SEEK` flag.
