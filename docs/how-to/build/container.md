# Build and run DiagHam in a container

Purpose: get a tested DiagHam without installing a toolchain, and capture exactly what a result was computed with.
Source: new (2026-09-28); `container/Dockerfile`, `container/diagham.def`, `scripts_cmake/environment.sh`, the `container` job in `.github/workflows/ci.yml`. Status: the Docker recipe is built and its full test suite run by CI on every push to `main`; neither Docker nor Apptainer is available on the machine this was written on, so the Apptainer definition mirrors the Dockerfile step for step but has not been built here — the first person with Apptainer should run `apptainer build` and report.

## Docker

```
docker build -t diagham -f container/Dockerfile .
docker run --rm diagham FQHESphereFermionsTwoBodyGeneric --help
docker run --rm -v "$PWD":/work -w /work diagham FQHESphereJackGenerator -a -2 ...
```

The build stage compiles the `full` preset (LAPACK, GSL, GMP, FFTW3, bzip2), runs the **whole ctest suite** — the image is not produced if a test fails — and installs to `/opt/diagham`. The runtime image carries the installed programs, libraries and headers plus the shared libraries they need — no sources, objects or build tree; the installed tree is still about 1.3 GB because the 61 libraries are static and each of the 603 programs links what it needs. `docker build --build-arg PRESET=default` builds without optional libraries; `JOBS` sets the parallelism.

## Apptainer / Singularity (clusters)

```
apptainer build diagham.sif container/diagham.def    # from the repository root
apptainer exec diagham.sif FQHESphereFermionsTwoBodyGeneric --help
apptainer test diagham.sif                            # runs a golden inside the image
```

Same recipe; programs are on `PATH` inside the image and your working directory is bind-mounted as usual. MPI runs need the host's MPI matched with the image's — the definition builds without MPI; for a cluster build follow [hpc-cluster.md](hpc-cluster.md) natively.

## What a result was computed with

Every build can record its environment:

```
scripts_cmake/environment.sh build/full > environment.txt
```

writes the repository commit and DiagHam revision, the CMake options that were on, the compiler, CMake, Python/numpy and MPI versions, the LAPACK/GSL/GMP/FFTW/bzip2 package versions, loaded modules if any, and the output of `TestDiagHamConf` (what the compiled code itself reports). Keep that file next to the spectra of a calculation you publish. CI stores it (with the ctest log) as an artefact of every run; the container images carry it at `/opt/diagham/environment.txt`.

## Limits

- Bit-for-bit reproducibility across machines is not promised: different BLAS/LAPACK builds and CPUs change results at the 10⁻¹⁴–10⁻¹² level, which the goldens' tolerances absorb; a change of eigensolver can reorder degenerate states (see the regression-case rule in [../../reference/tests.md](../../reference/tests.md)).
- The images pin Ubuntu 24.04 but not package versions: a rebuild months later may carry newer libraries — the environment file says which.
