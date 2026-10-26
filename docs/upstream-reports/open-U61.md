# Open finding U61: FQHETorusMPSCreateState --topological-sector selects a different state for each P truncation mod m

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U61), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15.

## What happens

`--topological-sector s` (0 to m - 1) should pick one of the m Laughlin states on the torus, whatever the truncation. `FQHEMPSLaughlinMatrix::GetTopologicalSectorIndices` (`FQHE/src/Tools/FQHEMPS/FQHEMPSLaughlinMatrix.cc`, line 1052) keeps the bond indices whose charge index satisfies `CurrentQValue % LaughlinIndex == topologicalSector`, but the charge index is an internal one, shifted by `NValueGlobalShift`, which the constructor sets to the P truncation. So the selected sector moves with P mod m.

Example, nu = 1/3, N_phi = 6, aspect ratio 2 (the wiki's first example, with `--trim-qsector`): the Ky of the state written for `--topological-sector 0, 1, 2` is

| P | sector 0 | sector 1 | sector 2 |
|---|---|---|---|
| 3 | 1 | 5 | 3 |
| 4 | 3 | 1 | 5 |
| 5 | 5 | 3 | 1 |
| 6 | 1 | 5 | 3 |
| 7 | 3 | 1 | 5 |
| 8 | 5 | 3 | 1 |

Each state is itself a correct Laughlin state: it is the V1 zero mode of its Ky sector to 5e-15 (`tests/oracles/torus_mps.py`). But a convergence study in P at fixed `--topological-sector` compares different states, and the sector numbering of the wiki ("the angular momentum along the y axis is 1 in the above mentioned example") holds only for P = 0 mod 3.

## Fix direction

Select on the physical charge, `(CurrentQValue - NValueGlobalShift)` modulo m (taken positive), or document that the label is relative to the truncation.

## Consequence for the tests

`openfinding.fqhe.torus.mps_create_state.topological_sector_moves_with_truncation.nphi6_r2_p3_p4` (WILL_FAIL) asks for the same Ky at P = 3 and P = 4 for each sector. The physics tests compare each written state with the zero mode of the Ky sector in its file name.
