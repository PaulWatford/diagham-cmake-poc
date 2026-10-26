# Open finding U53: every fermionic Laughlin MPS on the sphere or cylinder exits 0 without computing anything

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U53), found 2026-10-26 while writing the first independent tests of `FQHESphereMPSCreateState`; r4493 sources, GCC 15. The same on a pristine autotools build of r4493.

## What happens

The sphere and cylinder constructor of `FQHEMPSLaughlinMatrix` (`FQHE/src/Tools/FQHEMPS/FQHEMPSLaughlinMatrix.cc`, line 107) ends with an unconditional call to `ComputeSiteDependentMatrices(-8, 3)`. That function (line 1281) starts with

    if (this->BosonicVersion == false)
      {
        cout << "FQHEMPSLaughlinMatrix::ComputeSiteDependentMatrices only works using the MPS bosonic convention" << endl;
        exit(0);
      }

and `BosonicVersion` is the program's `--boson` flag (`FQHEMPSMatrixManager.cc`, around line 563). So any run with a fermionic Laughlin state stops while the B matrices are being built, prints that one line and exits with status 0, whether or not `--site-dependent` was asked for. No state, spectrum or density is written.

Both lines came in at r4149 (2022-01-05, "start improving MPS RSES code").

Programs checked, each with `--laughlin-index 3` and no `--boson`: `FQHESphereMPSCreateState`, `FQHESphereMPSEntanglementSpectrumParticlePartition`, `FQHESphereMPSEntanglementSpectrum`, `FQHESphereMPSEntanglementSpectrumUsingSVD`, `FQHESphereMPSDensity`, `FQHESphereMPSCorrelation`, `FQHEMPSEMatrix`, `FQHEMPSBMatrix`. All eight print the line and exit 0. The wiki's own first examples for `FQHESphereMPSCreateState` (`--laughlin-index 3 --reference-file laughlin3_n_4_2s_11.dat --use-padding --p-truncation 1`) and `FQHESphereMPSEntanglementSpectrumParticlePartition` (the nu = 1/3 cylinder orbital entanglement spectrum) are among them. The torus constructor does not make the call, and the bosonic Laughlin state is built correctly: its coefficients match the exact expansion of prod (z_i - z_j)^m (see below).

## Fix direction

Compute the site-dependent matrices only when they are asked for (`--site-dependent` or the real-space entanglement code that needs them), or skip the call for the fermionic convention instead of exiting.

## Consequence for the tests

`openfinding.fqhe.sphere.mps_create_state.fermion_laughlin_exits_without_state.m3_n4` (WILL_FAIL) asks `FQHESphereMPSCreateState` for the fermionic nu = 1/3 Laughlin state of 4 particles and compares it with the exact Slater-determinant expansion (`tests/oracles/mps_laughlin.py`); it fails while nothing is written. The bosonic Laughlin states (m = 2, 4, 6) are physics tests, in the monomial and the normalised sphere bases.
