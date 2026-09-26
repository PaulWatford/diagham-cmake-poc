# Tutorial: the Hubbard model end to end, at machine precision

Purpose: run DiagHam's Hubbard-model program on cases with known answers, read its output, and see how the results are checked against an independent calculation.
Source: `HUBBARD_BENCHMARK.md` (proof-of-concept demonstration, 2026, corrected 22/09) and `benchmarks/BENCHMARK.md` (the verification log), merged; every number below is reproduced by the ctest goldens on the r4493 build. Changed: the two overlapping documents become one walkthrough; the log itself stays in `benchmarks/BENCHMARK.md`.

## The model and the program

The 2D Hubbard model on an L×L square lattice with periodic boundary
conditions:

```
H = -t Σ_{<ij>,σ} (c†_{iσ} c_{jσ} + h.c.) + U Σ_i n_{i↑} n_{i↓}
```

`HubbardSquareLatticeModel` (in `FTI/src/Programs/HubbardModels/`)
diagonalises it in every momentum sector (kx, ky) separately, with full
diagonalisation for small Hilbert spaces or Lanczos for larger ones. Half
filling means N = L² electrons, split evenly between the spins (Sz = 0).

One convention to know: `--nn-t T` is the coefficient of c†c directly, i.e.
H = +T Σ c†c. On 2×L lattices the spectrum is symmetric under t → −t so it
does not matter; on odd-length lattices (3×3) pass `--nn-t -1` to match the
textbook sign.

## Test 1: 2×2 at U = 4 — an exact answer

```
mkdir run && cd run
../build/default/FTI/src/Programs/HubbardModels/HubbardSquareLatticeModel \
    -p 4 -x 2 -y 2 --u-potential 4 --nn-t 1.0 --full-diag 1000
```

The Hilbert space at Sz = 0 has 36 states over four momentum sectors
(dimensions 12, 8, 8, 8). The output file
`fermions_hubbard_square_x_2_y_2_n_4_ns_4_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat`
has one eigenvalue per line, `# kx ky sz E`; the lowest is in the (0,0)
sector:

```text
0 0 0 -5.6568542494923806
```

The 2×2 model at U = 4 is one of the canonical exactly solvable many-body
problems: group theory gives the ground state in closed form, −4√2. That is
−5.65685424949238019520… and its correctly rounded double is
`-5.6568542494923806` — DiagHam's answer is bit-identical (0 ULP). Getting
this right requires the Hamiltonian construction, the basis enumeration,
the fermion sign tracking and the diagonaliser all to be correct to every
digit.

(Before the output-precision fix on `main` the file showed
`-5.6568542494924`, 14 digits; the accuracy was the same, the print was not.)

## Test 2: U = 0 — free fermions

Same command with `--u-potential 0` (the file name then omits the `_u_…`
token). The single-particle energies on the 2×2 periodic lattice are
{−4, 0, 0, +4}; filling two per spin gives exactly −8:

```text
0 0 0 -8.0000000000000053
```

found twice (sectors (0,0) and (1,1)), the expected two-fold degeneracy from
the zero-energy single-particle states.

## Test 3: a U sweep

```
for U in 0 1 2 3 4 5 6 7 8; do
  ../build/default/FTI/src/Programs/HubbardModels/HubbardSquareLatticeModel \
      -p 4 -x 2 -y 2 --u-potential $U --nn-t 1.0 -S --processors 2 -n 1
done
```

| U/t | E₀ | E₀/N_site |
|---:|---:|---:|
| 0 | −8.0000000000000053 | −2.0000 |
| 1 | −7.2979734435359536 | −1.8245 |
| 2 | −6.6816952344966900 | −1.6704 |
| 3 | −6.1381306043250996 | −1.5345 |
| 4 | −5.6568542494923806 | −1.4142 (= −√2) |
| 5 | −5.2294257740670140 | −1.3074 |
| 6 | −4.8488578017961066 | −1.2122 |
| 7 | −4.5092407162187769 | −1.1273 |
| 8 | −4.2054969669241498 | −1.0514 |

Every value agrees with the independent Python diagonaliser below to better
than 10⁻¹³ (max 1.7×10⁻¹⁴).

## Test 4: strong coupling

At U ≫ t the half-filled model maps onto a Heisenberg antiferromagnet. For
the 2×2 periodic lattice the effective Hamiltonian gives E₀·U → −48. Running
U ∈ {50, 100, 200, 500, 1000, 2000} and fitting E₀·U = a + b/U + c/U² + d/U³
gives a → −47.9999 (cubic fit), 10⁻⁴ from theory — a check that the bond
multiplicity of the 2×2 periodic topology, the full Hubbard dynamics and the
t–J limit are all handled correctly.

## Test 5: 2×4 — Lanczos

```
../build/default/FTI/src/Programs/HubbardModels/HubbardSquareLatticeModel \
    -p 8 -x 2 -y 4 --u-potential 4 --nn-t 1.0 -S --processors 4 -n 1
```

4,900 states per momentum sector; Lanczos converges in about 40
iterations. Lowest eigenvalues by sector:

```text
0 0 0 -10.252952955263581    <- ground state
1 2 0  -9.8328986753265575
1 1 0  -8.6073071915857
1 3 0  -8.6073071915857
0 1 0  -8.2723788881326
```

The ground state sits in (0,0), as expected for a half-filled bipartite
lattice; E₀/N_site = −1.28 t agrees with published 2×4 benchmarks at U/t = 4.
With `--full-diag 5000` all 4,900 eigenvalues can be compared with the
Python code: they agree to 10⁻¹² for 99.9% of levels, and the degeneracy
multiplicities at the integer-valued levels (37, 40, 90, 37) match exactly.

## The independent check: a diagonaliser from scratch

`benchmarks/hubbard_ed.py` is a 100-line Python program that builds the
same Hamiltonian in the second-quantised Fock basis — Jordan–Wigner signs
by hand, dense numpy diagonalisation — sharing nothing with DiagHam but the
model's definition. `python3 benchmarks/hubbard_ed.py` runs the 2×2 and 2×4
cases (about ten seconds). The ctest `python.hubbard_ed_cross_check` runs
both codes on four cases and requires agreement to 10⁻¹⁰; observed
differences are 10⁻¹⁵–10⁻¹⁴.

## What this establishes, and what it does not

It establishes that the CMake-built `HubbardSquareLatticeModel` reproduces
exact and independently computed results to machine precision on 2×2, 2×4
and 3×3 cases, with the strong-coupling limit recovered. It does not
establish the correctness of the other 602 programs, of larger lattices
where Lanczos convergence parameters matter, or of FQHE programs — those
have their own goldens (`docs/reference/tests.md`) and, for FQHE,
[your first FQHE run](first-fqhe-run.md).

Full log with every table: `benchmarks/BENCHMARK.md`.
