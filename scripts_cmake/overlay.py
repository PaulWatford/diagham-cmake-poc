#!/usr/bin/env python3
"""
overlay.py: lay this repository's CMake build over an upstream DiagHam tree.

This repository is build scaffolding, not the DiagHam sources. To build,
its files are copied into an upstream checkout (the guysoft/DiagHam mirror,
or a Kent GitLab / SVN working copy), and the per-directory CMakeLists.txt
files are generated from that checkout's own Makefile.am files:

    python3 scripts_cmake/overlay.py /path/to/DiagHam
    cd /path/to/DiagHam
    cmake --preset default
    cmake --build --preset default
    ctest --preset default

Copied: CMakeLists.txt, CMakePresets.json, cmake/, patches/, scripts_cmake/,
tests/ and benchmarks/. Nothing in the upstream tree is overwritten except
files this script itself wrote on an earlier run (the copied files and
the generated CMakeLists.txt), so re-running after an upstream update or a
change here is safe.

Patches in patches/ are applied by CMake at configure time, not by this
script (see cmake/ApplyUpstreamPatches.cmake).
"""
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
FILES = ('CMakeLists.txt', 'CMakePresets.json')
DIRECTORIES = ('cmake', 'patches', 'scripts_cmake', 'tests', 'benchmarks')


def main(argv):
    if len(argv) != 2 or argv[1] in ('-h', '--help'):
        print(__doc__.strip())
        return 0 if len(argv) == 2 else 2
    target = Path(argv[1]).resolve()
    if not (target / 'configure.ac').is_file() or not (target / 'src' / 'config.h').is_file():
        print(f'overlay.py: {target} does not look like a DiagHam tree '
              '(expected configure.ac and src/config.h)', file=sys.stderr)
        return 1
    if target == REPO_ROOT:
        print('overlay.py: target is this repository itself', file=sys.stderr)
        return 1

    for name in FILES:
        shutil.copy2(REPO_ROOT / name, target / name)
    for name in DIRECTORIES:
        destination = target / name
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(REPO_ROOT / name, destination,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))

    subprocess.run([sys.executable, str(target / 'scripts_cmake' / 'extract_autotools.py'),
                    str(target)], check=True)
    print(f'Overlay complete. Next: cd {target} && cmake --preset default')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
