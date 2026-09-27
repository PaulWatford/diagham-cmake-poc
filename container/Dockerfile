# DiagHam (this repository's CMake build) in a container.
#
#   docker build -t diagham -f container/Dockerfile .
#   docker run --rm diagham FQHESphereFermionsTwoBodyGeneric --help
#   docker run --rm -v "$PWD":/work -w /work diagham FQHESphereJackGenerator ...
#
# Stage "build" compiles the `full` preset (LAPACK, GSL, GMP, FFTW3, bzip2),
# runs the whole ctest suite -- the image is not produced if a test fails --
# and installs to /opt/diagham. Stage "runtime" carries only the installed
# programs, libraries, headers and the runtime shared libraries.
# Build arguments: PRESET (default full), JOBS.
ARG UBUNTU=24.04
FROM ubuntu:${UBUNTU} AS build
ARG PRESET=full
ARG JOBS=4
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates git cmake ninja-build g++ python3 python3-numpy \
        liblapack-dev libblas-dev libgsl-dev libgmp-dev libfftw3-dev libbz2-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . /src
RUN cmake --preset ${PRESET} -G Ninja \
    && cmake --build --preset ${PRESET} --parallel ${JOBS} \
    && ctest --preset ${PRESET} --parallel ${JOBS} \
    && cmake --install build/${PRESET} --prefix /opt/diagham \
    && scripts_cmake/environment.sh build/${PRESET} > /opt/diagham/environment.txt

FROM ubuntu:${UBUNTU} AS runtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        liblapack3 libblas3 libgsl27 libgmp10 libfftw3-double3 libbz2-1.0 libgomp1 python3 python3-numpy \
    && rm -rf /var/lib/apt/lists/* \
    && (apt-get install -y --no-install-recommends libgsl28 2>/dev/null || true)
COPY --from=build /opt/diagham /opt/diagham
ENV PATH=/opt/diagham/bin:$PATH
LABEL org.opencontainers.image.title="DiagHam (CMake migration)" \
      org.opencontainers.image.source="https://github.com/PaulWatford/diagham-cmake-poc" \
      org.opencontainers.image.description="DiagHam exact-diagonalization toolkit, built and tested with this repository's CMake build; /opt/diagham/environment.txt says what it was built with"
WORKDIR /work
CMD ["bash"]
