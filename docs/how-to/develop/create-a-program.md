# Add a new program to DiagHam

Purpose: add your own program, either inside the DiagHam tree (so it is built, installed and smoke-tested with everything else) or as a separate project that links the installed library.
Source: upstream wiki "Create_new_FQHE_code" (sections "Adding a new single program" and "Adding a new package", as of 2026-09-24). Changed: the `Makefile.am` dependency/link-line steps and `configure.in`/`bootstrap.sh` are replaced by the generator; the example program was rebuilt and run against the r4493 package on 2026-09-26 (it now gates `LapackDiagonalize` so it also works without LAPACK, and prints at full precision).

## Option A: inside the tree

Say you want `FQHESphereMyCode`, built from `FQHESphereMyCode.cc`. Since it
is FQHE on the sphere, it belongs in `FQHE/src/Programs/FQHEOnSphere/`.

1. Put `FQHESphereMyCode.cc` in that directory.
2. Register it in the directory's `Makefile.am`, exactly as upstream does —
   add the name to `bin_PROGRAMS` and one `_SOURCES` line:

   ```
   bin_PROGRAMS=QHEFermions QHEBosons ... FQHESphereQuasiholeMatrixElements FQHESphereMyCode

   FQHESphereMyCode_SOURCES=FQHESphereMyCode.cc
   ```

   The long `_DEPENDENCIES` and `_LDADD` lines the wiki asked for are **not
   needed**: every program links every DiagHam library and the external
   dependencies through one helper. (Keeping `Makefile.am` in step means
   the autotools build and the parity check still see your program, and
   the change can go upstream unchanged.)
3. Regenerate and build:

   ```
   python3 scripts_cmake/extract_autotools.py .
   cmake --build --preset default -j
   ```

   The generator emits `diagham_add_programs(... FQHESphereMyCode FQHESphereMyCode.cc ...)`
   in that directory's `CMakeLists.txt`; the binary is
   `build/default/FQHE/src/Programs/FQHEOnSphere/FQHESphereMyCode`, the
   target `FQHEOnSphere_FQHESphereMyCode`, and `ctest` gains
   `smoke.FQHEOnSphere_FQHESphereMyCode` (your program must exit 0 on
   `--help`). Commit the source, the `Makefile.am` change and the
   regenerated `CMakeLists.txt` together.

A new *directory* of programs (the wiki's "MyCodes" package) is the same
recipe one level up: create `FQHE/src/Programs/MyCodes/` with a
`Makefile.am` containing `SUBDIRS=`, `bin_PROGRAMS=…` and the `_SOURCES`
lines, add `MyCodes` to `SUBDIRS` in `FQHE/src/Programs/Makefile.am`, and
regenerate. No `configure.in` edit, no `bootstrap.sh`.

## Option B: a separate project against the installed library

If you would rather not touch the DiagHam tree, install DiagHam once
(`cmake --install build/lapack --prefix $HOME/opt/diagham`) and build your
program with three lines of CMake — see
[../install-and-use-as-library.md](../install-and-use-as-library.md). The
example below builds either way.

## The example program

`FQHESphereMyCode.cc` reads a real symmetric matrix from a DiagHam binary
matrix file and prints its eigenvalues. It shows the option-manager
pattern every DiagHam program uses:

```cpp
#include "config.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Options/Options.h"

#include <iostream>
#include <limits>

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(std::numeric_limits<double>::max_digits10);

  // this part of code takes care of the command line options
  OptionManager Manager ("FQHESphereMyCode" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "sym-matrix", "name of the file containing the real symmetric matrix to diagonalize");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMyCode -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("sym-matrix") == 0)
    {
      cout << "no real symmetric matrix has been provided" << endl;
      return 0;
    }

  RealSymmetricMatrix Matrix2;
  if (Matrix2.ReadMatrix(Manager.GetString("sym-matrix")) == false)
    {
      cout << "can't read " << Manager.GetString("sym-matrix") << endl;
      return 0;
    }
  RealDiagonalMatrix TmpDiagonalMatrix (Matrix2.GetNbrRow());
#ifdef __LAPACK__
  Matrix2.LapackDiagonalize(TmpDiagonalMatrix);
#else
  Matrix2.Diagonalize(TmpDiagonalMatrix);
#endif
  cout << "eigenvalues of " << Manager.GetString("sym-matrix") << " : " << endl;
  for (int i = 0; i < TmpDiagonalMatrix.GetNbrRow(); ++i)
    cout << TmpDiagonalMatrix[i] << endl;
  return 0;
}
```

Three habits worth copying from the existing programs: check `--help`
*before* checking required arguments (or `--help` never works — see
`EvaluateBroadening` in the defect register); gate `LapackDiagonalize`
behind `#ifdef __LAPACK__` with the built-in `Diagonalize` as the fallback
(the most common upstream build failure); and print results with
`max_digits10`, not `precision(14)`.
