# Open finding U59: HubbardSuperconductorOrderParameter -o writes to a file with an uninitialised name

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U59), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The same on a pristine autotools build of r4493.

## What happens

`FTI/src/Programs/HubbardModels/HubbardSuperconductorOrderParameter.cc`, lines 400 to 404:

    if (Manager.GetString("output-file") != 0)
      {
        OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
        File.open(OutputFileName, ios::binary | ios::out);
      }

The buffer is allocated but the option's value is never copied into it, so the order parameter is written to a file whose name is whatever the uninitialised buffer holds (for example a few non-printable bytes in the working directory), and the file asked for is not created. The run exits 0. Without `-o` the default name (the left state's name with `.vec` replaced by `.orderparam.dat`) works.

A static scan of the tree finds the same three lines (allocation, no copy, open) in `Spin/src/Programs/SpinComputeBondEnergyTriangleLatticeProjectedFromKagome.cc` (lines 189 and 375), `Spin/src/Programs/SpinKagomeComputeSpinSpinCorrelations.cc` (line 202) and `Spin/src/Programs/SpinApplyMirrorSymmetry.cc` (line 172); those were not run.

## Fix direction

Add `strcpy (OutputFileName, Manager.GetString("output-file"));` after each allocation.

## Consequence for the tests

`openfinding.hubbard.superconductor_order_parameter.output_option_unnamed_file.ns4_n2_to_n4` (WILL_FAIL) runs a physics case with `-o order.dat` and fails while that file is not written. The physics tests use the default name.
