# Run DiagHam under MPI

Purpose: run a DiagHam program on a cluster or any distributed-memory machine, in MPI-only or mixed MPI/SMP mode, and profile it.
Source: upstream wiki "MPI" page (as of 2026-09-24). Changed: the build section now points to the CMake guides; option names (`--mpi`, `--mpi-smp`, `--cluster-profil`, `--auto-loadbalancing`) and the cluster-file format were checked against `src/Architecture/MixedMPISMPArchitecture.cc` and `ArchitectureManager.cc` at r4493 and are unchanged; the example command's program path was verified.

DiagHam relies on MPI for distributed memory. MPI support is a build-time
choice — the `hpc` preset, or `-DDIAGHAM_USE_MPI=ON` (the equivalent of
`--enable-mpi`); see [build/hpc-cluster.md](build/hpc-cluster.md). To use a
particular MPI C++ wrapper (the old `--with-mpi-cxx=...`), set
`-DMPI_CXX_COMPILER=mpicxx` or configure with `CXX=mpicxx`. DiagHam has been
extensively tested with OpenMPI; Intel MPI needs the
`-DMPICH_IGNORE_CXX_SEEK` flag described in [build/intel-mkl.md](build/intel-mkl.md).

Two modes exist: MPI-only, and mixed MPI/SMP. Running on a single machine
with many cores does **not** need MPI: the SMP mode (`-S --processors N`,
pthreads) covers that and is available in every build.

## MPI-only mode

Add `--mpi` to the program's command line and launch it with `mpirun`
(one MPI process per core).

## MPI/SMP mode

`--mpi-smp <cluster file>` runs one MPI process per node and uses threads
within the node. It needs a text file describing the cluster, readable at
least by the master node.

### Full cluster description

One line per node: hostname, number of cores to use, then two optional
columns — a relative performance index for static load balancing, and the
memory in MB the node may use (a per-node `--memory` value):

```
nostromo1 6 0.11 60000
nostromo2 6 0.15 50000
nostromo3 4 0.14 60000
```

Because MPI/SMP assumes one MPI process per node, the MPI hostfile must give
each node one slot:

```
nostromo1 slots=1
nostromo2 slots=1
nostromo3 slots=1
```

### Generic descriptions

With a batch queue you rarely know the node names in advance. Two other forms
are accepted (both implemented in `MixedMPISMPArchitecture`):

- `default` applies to every node not otherwise listed, and `master` may
  override the master node:

  ```
  master 1 0.1 100
  default 2 0.15 400
  ```

- a trailing `*` matches every hostname with that prefix, so two machine
  types can be described in two lines, and these can be mixed with exact
  hostnames:

  ```
  nostromo* 4 0.2 10000
  discovery* 8 0.15 60000
  ```

## Threads per MPI process under a batch system

Under SLURM or similar, look for the system's MPI+OpenMP ("hybrid") section
to request cores per task. With OpenMPI directly, process binding can cap the
threads a process may use; `mpirun --bind-to none` lifts that. A typical
command (the wiki's, program path verified at r4493):

```
mpirun --bind-to none -hostfile hostfile \
  $DIAGHAM/build/hpc/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsTwoBodyGeneric \
  --symmetrized-basis -p 14 -l 34 --interaction-file pseudopotential_coulomb_l_0_2s_34.dat \
  --interaction-name coulomb_l_0 --nbr-lz 1 -n 1 --fast-disk \
  --mpi-smp cluster.desc --cluster-profil cluster.log --auto-loadbalancing \
  --eigenstate --lanczos-precision 1e-13 --show-itertime --memory 0
```

## Profiling and static load balancing

`--cluster-profil <log>` writes how much time each node spends in every
MPI/SMP operation (`VectorHamiltonianMultiply core operation done in 3.002
seconds`, …). Dynamic load balancing is not available; the log is used to
tune the static one. `scripts/ClusterProfiling.pl cluster.log cluster.desc`
prints per-operation totals, per-node means and, at the end, an "optimized
performance index" block in cluster-file format that can replace the
original description. `--auto-loadbalancing` asks the run to balance from
the description it was given.

## Not yet covered by the test suite

The `hpc` preset compiles and links the MPI code paths and passes the serial
tests, but no ctest launches a program under `mpirun`; see `TESTING.md`.
