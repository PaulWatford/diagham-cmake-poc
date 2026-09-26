# Data formats and file tools

Purpose: the on-disk formats DiagHam reads and writes (binary vectors, binary matrices, spectrum files, interaction files) and the programs that manipulate them.
Source: upstream wiki pages "Binary_vectors" and "Binary_matrices" (as of 2026-09-24), merged; formats checked against `src/Vector/RealVector.cc`, `src/Vector/ComplexVector.cc` and the `MainTask` writers at r4493. Changed: the large-vector header (which the wiki omits) and the complex-vector layout added from the source; tool names linked to the program reference; spectrum, interaction and test-checker formats added.

All binary data is little-endian (the native order on x86 and ARM), so the
files are readable from any language, as the examples below show.

## Binary vectors (`.vec`)

Written by `RealVector::WriteVector` / `ComplexVector::WriteVector`:

| Field | Type | Meaning |
|---|---|---|
| dimension | `int` (4 bytes) | number of components, or **−1** if the vector has 2³¹ components or more |
| large dimension | `long` (8 bytes), only if dimension = −1 | the number of components |
| components | `double` × dimension | for a real vector |
| components | (`double` re, `double` im) × dimension | for a complex vector: real and imaginary parts interleaved |

A real vector `test.vec` read with numpy (the wiki's example, valid for
dimension < 2³¹):

```python
import numpy as np
with open('test.vec', 'rb') as f:
    dimension = np.fromfile(f, '<i4', count=1)[0]
    vector = np.fromfile(f, '<d', count=dimension)
print('dimension', dimension, 'square norm', float(vector @ vector))
```

For a complex vector read `2*dimension` doubles and take even/odd entries as
real/imaginary parts. For the large format, read one `<i4` (−1), then one
`<i8` for the dimension.

## Binary matrices (`.mat`)

| Field | Type | Meaning |
|---|---|---|
| type | `int` | bit 0 set: real matrix; bit 1 set: complex matrix |
| rows m | `int` | |
| columns n | `int` | |
| elements | `double` × m·n (real) or 2·m·n (complex) | row-major: (0,0), (0,1), …, (1,0), …; complex as Re, Im consecutively |

The wiki's C++ reader, unchanged:

```cpp
ifstream File;
File.open("mymatrix.mat", ios::binary | ios::in);
int TmpType, TmpNbrRow, TmpNbrColumn;
File.read((char*) &TmpType, sizeof(int));
File.read((char*) &TmpNbrRow, sizeof(int));
File.read((char*) &TmpNbrColumn, sizeof(int));
if ((TmpType & 1) != 0)
  {
    double Tmp;
    for (int i = 0; i < TmpNbrRow; ++i)
      for (int j = 0; j < TmpNbrColumn; ++j)
        { File.read((char*) &Tmp, sizeof(double)); cout << i << " " << j << " " << Tmp << endl; }
  }
else
  {
    double TmpRe, TmpIm;
    for (int i = 0; i < TmpNbrRow; ++i)
      for (int j = 0; j < TmpNbrColumn; ++j)
        {
          File.read((char*) &TmpRe, sizeof(double));
          File.read((char*) &TmpIm, sizeof(double));
          cout << i << " " << j << " " << TmpRe << " " << TmpIm << endl;
        }
  }
```

## Spectrum files (`.dat`)

Diagonalisation programs write one text file per run into the current
directory, named after the system and its parameters (for example
`fermions_hubbard_square_x_2_y_2_n_4_ns_4_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat`;
at U=0 the `_u_…` token is omitted). The first line is a `#` header naming
the columns; every following line is one eigenvalue preceded by its quantum
numbers. The headers written by the `MainTask` classes at r4493:

| Header | Geometry / program family | Columns |
|---|---|---|
| `# Lz E` | sphere | total Lz (or 2Lz), energy |
| `# Ky E`, `# Kx Ky E` | torus | momentum sector(s), energy |
| `# kx ky sz E` | Hubbard / lattice models | momenta, Sz, energy |
| `# E |E| arg(E)/pi` | non-Hermitian problems | complex eigenvalue as modulus and phase |

Energies are written with 17 significant digits (`max_digits10`) since the
precision fix on `main`; the wiki-era output had 14. A Lanczos run appends
`# Lanczos terminated at step N` when it stops.

## Interaction files

Two-body programs (`FQHESphereFermionsTwoBodyGeneric`,
`FQHETorusFermionsTwoBodyGeneric`, …) take `--interaction-file` with the
pseudopotentials as a space-separated list:

```
Pseudopotentials = 0 1
```

(`tests/data/v1_pseudopotential.dat`: the pure V₁ Haldane pseudopotential,
whose zero-energy ground state on the sphere at 2S = 3(N−1) is the Laughlin
state.) Programs such as `CoulombPseudopotentials` generate these files.

## Programs for vectors and matrices

Each name links to its generated option page in the program reference.

| Program | Does |
|---|---|
| [VectorBinary2Ascii](programs/core/VectorBinary2Ascii.md) | binary vector → text column |
| [VectorAscii2Binary](programs/core/VectorAscii2Binary.md) | text column → binary vector |
| [VectorComplex2RealImaginary](programs/core/VectorComplex2RealImaginary.md), [VectorRealImaginary2Complex](programs/core/VectorRealImaginary2Complex.md) | split / join complex vectors |
| [VectorRational2Double](programs/core/VectorRational2Double.md) | exact-rational vector → double |
| [VectorPhaseMultiply](programs/core/VectorPhaseMultiply.md) | multiply by a complex phase |
| [NormalizeVector](programs/core/NormalizeVector.md), [ZeroingVector](programs/core/ZeroingVector.md), [CountingZero](programs/core/CountingZero.md) | normalise; zero components; count zero components |
| [GenericOverlap](programs/core/GenericOverlap.md), [DiffBinaryVectors](programs/core/DiffBinaryVectors.md) | overlaps between vectors; compare two vectors |
| [BuildSuperPosition](programs/core/BuildSuperPosition.md), [MergeVectors](programs/core/MergeVectors.md) | linear superposition; merge |
| [ExtractLinearlyIndependentVectors](programs/core/ExtractLinearlyIndependentVectors.md), [ReorthogonalizeVectorSet](programs/core/ReorthogonalizeVectorSet.md) | independent subset; re-orthogonalise a set |
| [MatrixBinary2Ascii](programs/core/MatrixBinary2Ascii.md) | binary matrix → text |
| [MatrixExtractColumns](programs/core/MatrixExtractColumns.md), [MatrixElement](programs/core/MatrixElement.md) | columns as vectors; one element |
| [GenericHamiltonianDiagonalization](programs/core/GenericHamiltonianDiagonalization.md), [GenericMatrixMultiplication](programs/core/GenericMatrixMultiplication.md) | diagonalise a stored matrix (or a combination); multiply |
| [ReplayLanczos](programs/core/ReplayLanczos.md), [ReplayFastLanczos](programs/core/ReplayFastLanczos.md), [ResizeLanczos](programs/core/ResizeLanczos.md) | work with saved Lanczos runs |

## Checking a spectrum file (test suite)

`tests/check_spectrum.cc` is the comparison tool behind the ctest goldens;
it needs no Python and no DiagHam library:

```
check_spectrum min      FILE COLUMN EXPECTED MAX_ULP   lowest value within MAX_ULP of EXPECTED
check_spectrum min-abs  FILE COLUMN EXPECTED TOL       lowest value within |TOL|
check_spectrum spectrum FILE COLUMN REFERENCE REF_COLUMN TOL   sorted columns agree element-wise
check_spectrum count    FILE COLUMN VALUE TOL N        exactly N values within TOL of VALUE
```

`COLUMN` is 0-based; −1 is the last column (the energy in every header above).
