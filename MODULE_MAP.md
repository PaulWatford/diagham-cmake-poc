# Task-oriented map

This repo's files are organised the way autotools left them (a
top-level build description, a `cmake/` support directory, a flat
`patches/` directory). That's fine for building, but it's not how a
contributor thinks about a task. This map goes the other way: start
from what you're trying to do, and it points at the file(s) that matter.
It supersedes navigating by directory listing for anything below "I want
to build it" — the plain directory tree is still in `README.md` for
that.

## "I want to configure or build DiagHam with this CMake port"

- `CMakeLists.txt` — the top-level build: options
  (`DIAGHAM_BUILD_FQHE`, `DIAGHAM_USE_LAPACK`, `DIAGHAM_USE_MPI`,
  `DIAGHAM_USE_GSL`, `DIAGHAM_KENT_DEFAULTS`, ...), `find_package`
  blocks, and the configuration summary printed at the end of
  `cmake -S . -B build`.
- `cmake/DiagHamHelpers.cmake` — the two functions every
  per-subdirectory `CMakeLists.txt` actually calls:
  `diagham_add_library()` and `diagham_add_program()`. Read this if a
  build target isn't showing up or is linking against the wrong
  libraries.
- `cmake/verify_build.sh` — the 79-check build-equivalence suite. Run
  this after any build-system change; see `AGENTS.md`'s "Locked facts"
  before touching anything it checks.

## "I want to add or change a dependency (LAPACK, MPI, GSL, MKL, ...)"

- `CMakeLists.txt` — where each dependency's `option()` and
  `find_package()`/module-path wiring lives (search for the
  dependency's name).
- `cmake/Find*.cmake` — one file per optional dependency not covered by
  CMake's own bundled Find modules (`FindGSL.cmake`, `FindMKL.cmake`,
  `FindNAG.cmake`, `FindMPC.cmake`, `FindMPFR.cmake`, `FindGMP.cmake`,
  `FindFFTW3.cmake`, `FindDLR.cmake`). Only `FindGSL.cmake` and
  `FindMKL.cmake` are currently wired into `CMakeLists.txt` — the rest
  are staged for when their dependency is actually enabled (see
  `MIGRATION_ROADMAP.md` for why they're present but unused).
- `cmake/config_ac.h.in` — the generated config header template. A new
  `HAVE_<X>` needs a `#cmakedefine HAVE_<X>` line here, matching the
  `set(HAVE_<X> 1)` pattern in `CMakeLists.txt` (never inject the
  downstream `__X__` macro directly — see the LAPACK/MPI root-cause
  writeup in `patches/PATCHES.md` for why that's a bug, not a style
  choice).

## "I want to add a new source file / new program"

Nothing to edit for a new `.cc` in an existing `Programs/` directory —
`diagham_add_programs_in_directory()` in
`cmake/DiagHamHelpers.cmake` globs `*.cc` automatically. For a new
*library* source file or a new subdirectory:

- `scripts_cmake/extract_autotools.py` — regenerate the per-subdirectory
  `CMakeLists.txt` files from the upstream `Makefile.am`s rather than
  hand-editing them. It's the source of truth for what's a library
  target, not this repo's own tree.

## "I want to understand or extend an upstream-bug patch"

- `patches/PATCHES.md` — the audit trail: every patch's root cause,
  what was tried, what shipped, and how it was verified. Read the
  relevant dated section before touching a `.patch` file.
- `patches/NN-*.patch` — the patches themselves, applied at configure
  time by `cmake/ApplyUpstreamPatches.cmake`.
- `DEFERRED.md` — files explicitly *not* patched, with the evidence for
  why, and what would need to change for that to be revisited.

## "I want to verify the physics is still correct"

- `benchmarks/BENCHMARK.md` — the verification log for this PoC
  iteration, including the independent Python cross-check.
- `benchmarks/hubbard_ed.py` — the from-scratch 100-line Python exact
  diagonalisation used as an independent reference, not derived from
  DiagHam's own code.
- `benchmarks/*.dat.txt` — saved DiagHam output the benchmark compares
  against; these are the Hubbard golden files referenced in
  `AGENTS.md`'s "Locked artefacts."
- `HUBBARD_BENCHMARK.md` — the original U-sweep demonstration at the
  repo root (predates the `benchmarks/` directory's own log).
- `BUG_torus_su2_coulomb.md` — a standalone writeup of one upstream
  physics bug (spinful torus Coulomb, N >= 3) surfaced but not fixed by
  this migration; not something this repo's patches touch.

## "I want to know what's decided vs. still open for the real migration"

- `MIGRATION_ROADMAP.md` — reconciles this PoC's scope and reference
  material (a GitHub scratch repo, a single-commit git mirror) against
  what was actually agreed for production (Kent GitLab, full SVN
  history, single repo, static linking).
- README.md's "Beyond this iteration" section — build-coverage gaps
  specifically (Spin/QuantumDots modules, remaining optional
  dependencies, packaging, CI, docs generation).

## "I'm an AI agent, or reviewing an AI-suggested change"

- `AGENTS.md` — the whole rulebook: scope, authority, the gate/damage
  pipeline, what must never happen, and the approval process. Read
  this first, before editing anything in this repo as an agent.

## What's still organised the old way (and why that's fine)

`src/`, `Base/src/`, `FQHE/src/`, `FTI/src/` (in the upstream DiagHam
tree this CMake build targets, not in this repo) keep their autotools-era
layout on purpose — `scripts_cmake/extract_autotools.py` derives CMake
from that structure precisely so this migration doesn't also have to be
a source-tree reorganisation. Modularisation here means *this repo's*
files are easy to navigate by task; it does not mean restructuring
DiagHam's own source tree, which is explicitly out of scope (see
`AGENTS.md`'s Scope section).
