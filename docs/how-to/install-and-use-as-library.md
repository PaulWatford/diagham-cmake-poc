# Install DiagHam and use it as a library

Purpose: install the programs and libraries to a prefix, and build your own code against them with `find_package(DiagHam)`.
Source: `cmake/DiagHamInstall.cmake` and `tests/consumer/` (2026-09-24); the steps below were run on 2026-09-26 against the r4493 `lapack` build (603 programs, 61 libraries installed; the example program configured, built and ran).

## Install

```
cmake --install build/lapack --prefix $HOME/opt/diagham
```

gives

```
bin/                      the programs                          (component Runtime)
lib/diagham/lib<NAME>.a   the 61 static libraries               (component Development)
include/diagham/...       the headers, one merged tree          (component Development)
lib/cmake/DiagHam/        DiagHamConfig.cmake, exported targets, Find modules
```

`--component Runtime` installs only the programs; `--component Development`
only the libraries, headers and CMake package. On a shared cluster
filesystem install once per toolchain and point users at `bin/`.

## Use from another CMake project

`CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.16)
project(MyDiagHamCode CXX)

find_package(DiagHam REQUIRED)          # CMAKE_PREFIX_PATH=<install prefix>

add_executable(FQHESphereMyCode FQHESphereMyCode.cc)
target_link_libraries(FQHESphereMyCode PRIVATE DiagHam::DiagHam)
```

```
cmake -S . -B build -DCMAKE_PREFIX_PATH=$HOME/opt/diagham
cmake --build build
./build/FQHESphereMyCode --help
```

`DiagHam::DiagHam` carries the include path, the compile definitions the
headers need (`HAVE_CONFIG_H`, `MACHINE_PRECISION`), every static library
(grouped for the mutually recursive link) and the external dependencies
the installed build was made with (pthreads, LAPACK, …). Include headers
exactly as DiagHam does (`#include "Matrix/RealSymmetricMatrix.h"`,
`#include "config.h"`). The package also sets `DiagHam_HAS_LAPACK`,
`DiagHam_HAS_MPI`, `DiagHam_HAS_GSL`, … so your code can adapt.

A complete example program is in
[develop/create-a-program.md](develop/create-a-program.md); the same test
(`install.find_package_consumer`) runs on every build, so the package is
known to work.

## Caveats

- The exported dependencies record the absolute paths of the external
  libraries found at build time; an installed tree copied to another
  machine or module path may need `CMAKE_PREFIX_PATH` pointing at
  equivalent libraries there.
- Two FQHEOnDisk programs were renamed in this repository because they
  collided with FQHEOnSphere ones at install (`FQHEDiskBosonsDelta`,
  `FQHEDiskLaughlinMonteCarloOverlap`); scripts written for the old names
  need updating.
