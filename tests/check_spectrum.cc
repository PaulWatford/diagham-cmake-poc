// check_spectrum: compare a DiagHam spectrum file against a golden value.
//
// DiagHam writes spectra as whitespace-separated text, one eigenvalue per
// line, with quantum numbers in the leading columns and '#' comment lines,
// e.g. HubbardSquareLatticeModel's "kx ky sz E". This tool reads one
// column of such a file and checks it, so the ctest suite needs neither
// Python nor any DiagHam library.
//
// Usage (COLUMN is 0-based; a negative COLUMN counts from the end, so -1
// is the last column):
//
//   check_spectrum min      FILE COLUMN EXPECTED MAX_ULP
//       lowest value in COLUMN equals EXPECTED to within MAX_ULP units in
//       the last place (0 means bit-identical)
//
// Every TOL below is either a number (absolute: |a - b| <= TOL) or "rel:X" (relative:
// |a - b| <= X * max(|a|, |b|, 1)). Every check prints the largest absolute and relative difference
// it found, whichever form its tolerance takes.
//
//   check_spectrum min-abs  FILE COLUMN EXPECTED TOL
//       lowest value in COLUMN equals EXPECTED to within TOL
//
//   check_spectrum spectrum FILE COLUMN REFERENCE REF_COLUMN TOL
//       the sorted COLUMN of FILE matches the sorted REF_COLUMN of
//       REFERENCE element by element to within |TOL|, same length
//
//   check_spectrum count    FILE COLUMN VALUE TOL N
//       exactly N values in COLUMN lie within |TOL| of VALUE (degeneracy
//       checks, e.g. the torus Laughlin ground-state multiplet)
//
//   check_spectrum line     FILE KEY REFERENCE TOL
//       the numbers following "KEY =" on the first line of FILE that
//       starts with KEY match those on the same line of REFERENCE, one by
//       one, to within |TOL| (pseudopotential files)
//
//   check_spectrum lowest   FILE COLUMN REFERENCE REF_COLUMN K TOL
//       the K lowest values of COLUMN match the K lowest of REF_COLUMN to
//       within |TOL| (a Lanczos run against a full diagonalisation)
//
//   check_spectrum numbers  FILE REFERENCE TOL
//       every line of FILE made only of numbers matches the corresponding
//       line of REFERENCE number by number to within |TOL|; lines with any
//       non-numeric token (text, timings) are ignored on both sides
//       (regression checks on a program's standard output)
//
//   check_spectrum nonzero  FILE COLUMN THRESHOLD N [FILTER_COLUMN FILTER_VALUE]...
//       exactly N rows have COLUMN > THRESHOLD, among the rows whose
//       FILTER_COLUMNs equal the FILTER_VALUEs (entanglement-spectrum
//       level counts per sector)
//
// Exit status: 0 on pass, 1 on a failed check, 2 on bad usage or input.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool ReadColumn(const char* path, int column, std::vector<double>& values)
{
  std::ifstream in(path);
  if (!in)
    {
      std::cerr << "check_spectrum: cannot open " << path << std::endl;
      return false;
    }
  std::string line;
  int lineNumber = 0;
  while (std::getline(in, line))
    {
      ++lineNumber;
      std::string::size_type first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos || line[first] == '#')
        continue;
      std::istringstream fields(line);
      std::vector<std::string> tokens;
      std::string token;
      while (fields >> token)
        tokens.push_back(token);
      int index = column < 0 ? (int) tokens.size() + column : column;
      if (index < 0 || index >= (int) tokens.size())
        {
          std::cerr << "check_spectrum: " << path << ":" << lineNumber
                    << " has no column " << column << std::endl;
          return false;
        }
      char* end = 0;
      double value = std::strtod(tokens[index].c_str(), &end);
      if (end == tokens[index].c_str() || *end != '\0')
        {
          std::cerr << "check_spectrum: " << path << ":" << lineNumber
                    << " column " << column << " is not a number: '"
                    << tokens[index] << "'" << std::endl;
          return false;
        }
      values.push_back(value);
    }
  if (values.empty())
    {
      std::cerr << "check_spectrum: no data rows in " << path << std::endl;
      return false;
    }
  return true;
}

// Distance between two finite doubles in units in the last place.
int64_t UlpDistance(double a, double b)
{
  int64_t ia, ib;
  std::memcpy(&ia, &a, sizeof(double));
  std::memcpy(&ib, &b, sizeof(double));
  // Map the sign-magnitude IEEE-754 ordering onto a monotonic integer line.
  if (ia < 0)
    ia = INT64_MIN - ia;
  if (ib < 0)
    ib = INT64_MIN - ib;
  return ia > ib ? ia - ib : ib - ia;
}

double ParseDouble(const char* text, const char* what)
{
  char* end = 0;
  double value = std::strtod(text, &end);
  if (end == text || *end != '\0')
    {
      std::cerr << "check_spectrum: " << what << " is not a number: '" << text << "'" << std::endl;
      std::exit(2);
    }
  return value;
}

// A tolerance: "1e-9" is absolute (|a - b| <= 1e-9); "rel:1e-12" is relative to the size of the numbers,
// |a - b| <= 1e-12 * max(|a|, |b|, 1), so values of order one or less are held to 1e-12 absolutely and a
// zero mode stays well defined.
struct Tolerance
{
  double value;
  bool relative;
};

Tolerance ParseTolerance(const char* text)
{
  Tolerance t;
  t.relative = std::strncmp(text, "rel:", 4) == 0;
  t.value = std::fabs(ParseDouble(t.relative ? text + 4 : text, "TOL"));
  return t;
}

double Scale(double a, double b)
{
  return std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

// the measure the tolerance is compared with: |a - b|, or |a - b| / max(|a|, |b|, 1)
double Deviation(double a, double b, const Tolerance& t)
{
  double d = std::fabs(a - b);
  return t.relative ? d / Scale(a, b) : d;
}

std::ostream& operator<<(std::ostream& out, const Tolerance& t)
{
  return out << (t.relative ? "relative " : "absolute ") << t.value;
}

// element-by-element comparison of the first n values; reports the worst absolute and relative differences
struct Worst
{
  double deviation, absolute, relative;
  size_t index;
};

Worst Compare(const std::vector<double>& values, const std::vector<double>& reference, size_t n, const Tolerance& t)
{
  Worst w = {0.0, 0.0, 0.0, 0};
  for (size_t i = 0; i < n; ++i)
    {
      double d = std::fabs(values[i] - reference[i]);
      double dev = Deviation(values[i], reference[i], t);
      w.absolute = std::max(w.absolute, d);
      w.relative = std::max(w.relative, d / Scale(values[i], reference[i]));
      if (dev > w.deviation)
        {
          w.deviation = dev;
          w.index = i;
        }
    }
  return w;
}

void Report(const Worst& w, const std::vector<double>& values, const std::vector<double>& reference,
            const Tolerance& t)
{
  std::cout << "max |diff| = " << w.absolute << ", max relative diff = " << w.relative << ", worst at index "
            << w.index << " (" << values[w.index] << " vs " << reference[w.index] << "), tolerance " << t
            << std::endl;
}

int Usage()
{
  std::cerr << "usage: check_spectrum min      FILE COLUMN EXPECTED MAX_ULP\n"
               "       check_spectrum min-abs  FILE COLUMN EXPECTED TOL\n"
               "       check_spectrum spectrum FILE COLUMN REFERENCE REF_COLUMN TOL\n"
               "       check_spectrum count    FILE COLUMN VALUE TOL N\n"
               "       check_spectrum line     FILE KEY REFERENCE TOL\n"
               "       check_spectrum lowest   FILE COLUMN REFERENCE REF_COLUMN K TOL\n"
               "       check_spectrum numbers  FILE REFERENCE TOL\n"
               "       check_spectrum nonzero  FILE COLUMN THRESHOLD N [FILTER_COLUMN FILTER_VALUE]...\n";
  return 2;
}

// numbers after "KEY =" on the first line starting with KEY
bool ReadKeyedLine(const char* path, const std::string& key, std::vector<double>& values)
{
  std::ifstream in(path);
  if (!in)
    {
      std::cerr << "check_spectrum: cannot open " << path << std::endl;
      return false;
    }
  std::string line;
  while (std::getline(in, line))
    {
      std::string::size_type first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos || line.compare(first, key.size(), key) != 0)
        continue;
      std::string::size_type eq = line.find('=', first + key.size());
      if (eq == std::string::npos)
        continue;
      std::istringstream fields(line.substr(eq + 1));
      double x;
      while (fields >> x)
        values.push_back(x);
      return true;
    }
  std::cerr << "check_spectrum: no line starting with '" << key << " =' in " << path << std::endl;
  return false;
}

// the numeric lines of a file (lines whose every token parses as a double), flattened
bool ReadNumericLines(const char* path, std::vector<double>& values, size_t& lines)
{
  std::ifstream in(path);
  if (!in)
    {
      std::cerr << "check_spectrum: cannot open " << path << std::endl;
      return false;
    }
  std::string line;
  lines = 0;
  while (std::getline(in, line))
    {
      std::istringstream fields(line);
      std::vector<double> row;
      std::string token;
      bool numeric = true;
      while (fields >> token)
        {
          char* end = 0;
          double x = std::strtod(token.c_str(), &end);
          if (end == token.c_str() || *end != '\0')
            {
              numeric = false;
              break;
            }
          row.push_back(x);
        }
      if (numeric && !row.empty())
        {
          values.insert(values.end(), row.begin(), row.end());
          ++lines;
        }
    }
  return true;
}

// all rows of a whitespace table as doubles ('#' lines skipped)
bool ReadRows(const char* path, std::vector<std::vector<double> >& rows)
{
  std::ifstream in(path);
  if (!in)
    {
      std::cerr << "check_spectrum: cannot open " << path << std::endl;
      return false;
    }
  std::string line;
  while (std::getline(in, line))
    {
      std::string::size_type first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos || line[first] == '#')
        continue;
      std::istringstream fields(line);
      std::vector<double> row;
      double x;
      while (fields >> x)
        row.push_back(x);
      if (!row.empty())
        rows.push_back(row);
    }
  return true;
}

}  // namespace

int main(int argc, char** argv)
{
  if (argc < 4)
    return Usage();
  std::string mode = argv[1];
  const char* file = argv[2];
  int column = std::atoi(argv[3]);
  std::vector<double> values;
  std::cout << std::setprecision(17);

  if ((mode == "min" || mode == "min-abs") && argc == 6)
    {
      if (!ReadColumn(file, column, values))
        return 2;
      double expected = ParseDouble(argv[4], "EXPECTED");
      double lowest = *std::min_element(values.begin(), values.end());
      int64_t ulp = UlpDistance(lowest, expected);
      std::cout << "lowest   = " << lowest << "\nexpected = " << expected
                << "\ndiff     = " << (lowest - expected) << " (" << ulp << " ulp)"
                << "\nmax |diff| = " << std::fabs(lowest - expected) << ", max relative diff = "
                << std::fabs(lowest - expected) / Scale(lowest, expected) << std::endl;
      bool pass;
      if (mode == "min")
        pass = ulp <= std::atoll(argv[5]);
      else
        pass = Deviation(lowest, expected, ParseTolerance(argv[5])) <= ParseTolerance(argv[5]).value;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "spectrum" && argc == 7)
    {
      std::vector<double> reference;
      if (!ReadColumn(file, column, values) || !ReadColumn(argv[4], std::atoi(argv[5]), reference))
        return 2;
      Tolerance tolerance = ParseTolerance(argv[6]);
      if (values.size() != reference.size())
        {
          std::cout << "FAIL: " << values.size() << " values, reference has "
                    << reference.size() << std::endl;
          return 1;
        }
      std::sort(values.begin(), values.end());
      std::sort(reference.begin(), reference.end());
      Worst w = Compare(values, reference, values.size(), tolerance);
      std::cout << values.size() << " eigenvalues (sorted), ";
      Report(w, values, reference, tolerance);
      bool pass = w.deviation <= tolerance.value;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "count" && argc == 7)
    {
      if (!ReadColumn(file, column, values))
        return 2;
      double target = ParseDouble(argv[4], "VALUE");
      Tolerance tolerance = ParseTolerance(argv[5]);
      long expectedCount = std::atol(argv[6]);
      long count = 0;
      for (size_t i = 0; i < values.size(); ++i)
        if (Deviation(values[i], target, tolerance) <= tolerance.value)
          ++count;
      std::cout << count << " of " << values.size() << " values within " << tolerance
                << " of " << target << ", expected " << expectedCount << std::endl;
      bool pass = count == expectedCount;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "lowest" && argc == 8)
    {
      std::vector<double> reference;
      if (!ReadColumn(file, column, values) || !ReadColumn(argv[4], std::atoi(argv[5]), reference))
        return 2;
      size_t k = std::atol(argv[6]);
      Tolerance tolerance = ParseTolerance(argv[7]);
      if (values.size() < k || reference.size() < k)
        {
          std::cout << "FAIL: need " << k << " values, have " << values.size() << " and " << reference.size() << std::endl;
          return 1;
        }
      std::sort(values.begin(), values.end());
      std::sort(reference.begin(), reference.end());
      Worst w = Compare(values, reference, k, tolerance);
      std::cout << k << " lowest eigenvalues, ";
      Report(w, values, reference, tolerance);
      bool pass = w.deviation <= tolerance.value;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "numbers" && argc == 5)
    {
      std::vector<double> reference;
      size_t n1 = 0, n2 = 0;
      if (!ReadNumericLines(file, values, n1) || !ReadNumericLines(argv[3], reference, n2))
        return 2;
      Tolerance tolerance = ParseTolerance(argv[4]);
      if (values.size() != reference.size() || n1 != n2)
        {
          std::cout << "FAIL: " << n1 << " numeric lines / " << values.size() << " numbers, reference has "
                    << n2 << " / " << reference.size() << std::endl;
          return 1;
        }
      Worst w = Compare(values, reference, values.size(), tolerance);
      std::cout << n1 << " numeric lines, " << values.size() << " numbers, ";
      Report(w, values, reference, tolerance);
      bool pass = w.deviation <= tolerance.value;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "line" && argc == 6)
    {
      std::vector<double> reference;
      if (!ReadKeyedLine(file, argv[3], values) || !ReadKeyedLine(argv[4], argv[3], reference))
        return 2;
      Tolerance tolerance = ParseTolerance(argv[5]);
      if (values.size() != reference.size())
        {
          std::cout << "FAIL: " << values.size() << " numbers after '" << argv[3] << " =', reference has "
                    << reference.size() << std::endl;
          return 1;
        }
      Worst w = Compare(values, reference, values.size(), tolerance);
      std::cout << values.size() << " numbers, ";
      Report(w, values, reference, tolerance);
      bool pass = w.deviation <= tolerance.value;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "nonzero" && argc >= 6 && (argc - 6) % 2 == 0)
    {
      std::vector<std::vector<double> > rows;
      if (!ReadRows(file, rows))
        return 2;
      double threshold = ParseDouble(argv[4], "THRESHOLD");
      long expectedCount = std::atol(argv[5]);
      long count = 0, considered = 0;
      for (size_t r = 0; r < rows.size(); ++r)
        {
          const std::vector<double>& row = rows[r];
          bool selected = true;
          for (int a = 6; a + 1 < argc; a += 2)
            {
              int fc = std::atoi(argv[a]);
              int fi = fc < 0 ? (int) row.size() + fc : fc;
              if (fi < 0 || fi >= (int) row.size() || row[fi] != ParseDouble(argv[a + 1], "FILTER_VALUE"))
                selected = false;
            }
          if (!selected)
            continue;
          int index = column < 0 ? (int) row.size() + column : column;
          if (index < 0 || index >= (int) row.size())
            continue;
          ++considered;
          if (row[index] > threshold)
            ++count;
        }
      std::cout << count << " of " << considered << " selected rows have column " << column << " > "
                << threshold << ", expected " << expectedCount << std::endl;
      bool pass = count == expectedCount;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  return Usage();
}
