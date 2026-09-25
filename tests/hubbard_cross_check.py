#!/usr/bin/env python3
"""
Cross-check DiagHam's HubbardSquareLatticeModel against the independent
Python exact diagonalisation in benchmarks/hubbard_ed.py.

For each lattice, DiagHam's global ground state (lowest eigenvalue over all
momentum sectors, Sz = 0) must agree with the Python ED ground state to
1e-10. The two codes share nothing but the Hamiltonian's definition.

Usage: hubbard_cross_check.py HUBBARD_PROGRAM WORK_DIR
"""
import glob
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks'))
from hubbard_ed import ground_state  # noqa: E402

CASES = [  # (Lx, Ly, U)
    (2, 2, 0.0),
    (2, 2, 4.0),
    (2, 2, 100.0),
    (2, 4, 4.0),
]
TOLERANCE = 1e-10


def diagham_ground_state(program, work_dir, lx, ly, u):
    run_dir = os.path.join(work_dir, f'{lx}x{ly}_U{u:g}')
    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(run_dir)
    n = lx * ly
    subprocess.run([program, '-p', str(n), '-x', str(lx), '-y', str(ly),
                    '--u-potential', repr(u), '--nn-t', '1.0', '-n', '1'],
                   cwd=run_dir, check=True, stdout=subprocess.DEVNULL)
    (spectrum,) = glob.glob(os.path.join(run_dir, '*_sz_0.dat'))
    energies = []
    with open(spectrum) as f:
        for line in f:
            if line.strip() and not line.lstrip().startswith('#'):
                energies.append(float(line.split()[-1]))
    return min(energies)


def main():
    program, work_dir = sys.argv[1], sys.argv[2]
    failed = False
    for lx, ly, u in CASES:
        reference = ground_state(lx, ly, u)
        diagham = diagham_ground_state(program, work_dir, lx, ly, u)
        diff = abs(diagham - reference)
        ok = diff <= TOLERANCE
        failed |= not ok
        print(f'{lx}x{ly} U={u:<5g} DiagHam {diagham:.16g}  Python ED {reference:.16g}  '
              f'|diff| {diff:.2e}  {"PASS" if ok else "FAIL"}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
