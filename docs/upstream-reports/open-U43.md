# Open defect U43: HaahCodeEntropy sizes its state array 2^(number of Z terms)

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/known-defects.md` (U43), found 2026-09-30; r4493 sources, GCC 15.

## What happens

`Spin/src/Programs/HaahCodeEntropy.cc`, line 168:

```
long GroundStateDimension = 1l << MaximumNumberZTerms;
```

The ground state is the uniform superposition over the GF(2) span V of the Z terms (translates of an eight-spin mask), so it has 2^rank components, but the array is allocated for every subset of the terms and duplicates are removed afterwards. The number of terms is L^3; the rank is smaller:

| L | terms | rank of V | array requested |
|---|-------|-----------|-----------------|
| 2 | 8 | 5 | 2^8 |
| 3 | 27 | 26 | 2^27 (1 GB of `long`) |
| 4 | 64 | 57 | 2^63 -> `std::bad_array_new_length` |

`-x 4 -y 4 -z 4 --nbra-sitex 2 --nbra-sitey 2 --nbra-sitez 2 --export-entspectrum` therefore throws, although the entropy is a rank computation: for an a x a x a cube on the L^3 torus S / ln 2 = 6a^2 - 6a + 2, i.e. 14 for a = 2 (`tests/oracles/haah.py`; rank(V) for L = 2..7: 5, 26, 57, 124, 213, 342).

## Fix direction

Size the array from the rank of the terms (reduce them first), or generate the span by Gray-code enumeration of a basis.

## Consequence for the tests

`knownbug.haah_code.ground_state_array_2_pow_nterms.4x4x4_region_2x2x2` (WILL_FAIL) expects the oracle value 14. The 65-spin `ULONGLONG` branch of the program cannot be exercised for the same reason. The XCube variant with 3x3x3 also fails with `std::bad_alloc`; not investigated.
