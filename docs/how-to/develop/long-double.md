# Build DiagHam with long-double precision

Purpose: switch the whole code base from `double` to `long double`.
Source: upstream wiki "Long_double" (as of 2026-09-24). Changed: `svn checkout`/`bootstrap`/`configure` replaced by the CMake steps; the two vector-conversion programs the wiki attached as media files (`LongDoubleToDouble.cc`, `DoubleToLongDouble.cc`) were not archived and are noted as missing.

DiagHam uses double precision. If that is not enough, the tree can be
switched to `long double`. Beware: **LAPACK support is then unavailable**
(LAPACK is double precision), so use a preset without LAPACK. The proper
way, as on the wiki, is to do it on a fresh copy of the tree, because the
script rewrites the sources in place:

1. Make a fresh clone (or worktree) of this repository for the long-double
   build — do not run the script in your working tree.
2. From the top directory run upstream's script:
   ```
   ./scripts/switch2longdouble.pl
   ```
3. Configure and build with a preset that has no LAPACK, e.g.
   ```
   cmake --preset default
   cmake --build --preset default -j
   ```
4. Vectors written by a long-double build are not readable by a double
   build and vice versa (16-byte components on x86-64 instead of 8). The
   wiki offered two small converters, `LongDoubleToDouble.cc` and
   `DoubleToLongDouble.cc`, as attachments; they are not in the tree and
   were not archived. Writing one from the
   [binary vector format](../../reference/data-formats.md) is a few lines.

Note that the precision fix on `main` writes results with
`std::numeric_limits<double>::max_digits10`; in a long-double tree that
expression should read `<long double>` to print the extra digits.
