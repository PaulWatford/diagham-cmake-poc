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
//   check_spectrum min-abs  FILE COLUMN EXPECTED TOL
//       lowest value in COLUMN equals EXPECTED to within |TOL|
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
                << std::endl;
      bool pass;
      if (mode == "min")
        pass = ulp <= std::atoll(argv[5]);
      else
        pass = std::fabs(lowest - expected) <= std::fabs(ParseDouble(argv[5], "TOL"));
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "spectrum" && argc == 7)
    {
      std::vector<double> reference;
      if (!ReadColumn(file, column, values) || !ReadColumn(argv[4], std::atoi(argv[5]), reference))
        return 2;
      double tolerance = std::fabs(ParseDouble(argv[6], "TOL"));
      if (values.size() != reference.size())
        {
          std::cout << "FAIL: " << values.size() << " values, reference has "
                    << reference.size() << std::endl;
          return 1;
        }
      std::sort(values.begin(), values.end());
      std::sort(reference.begin(), reference.end());
      double worst = 0.0;
      size_t worstIndex = 0;
      for (size_t i = 0; i < values.size(); ++i)
        {
          double diff = std::fabs(values[i] - reference[i]);
          if (diff > worst)
            {
              worst = diff;
              worstIndex = i;
            }
        }
      std::cout << values.size() << " eigenvalues, max |diff| = " << worst
                << " at sorted index " << worstIndex << " (" << values[worstIndex]
                << " vs " << reference[worstIndex] << "), tolerance " << tolerance
                << std::endl;
      bool pass = worst <= tolerance;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "count" && argc == 7)
    {
      if (!ReadColumn(file, column, values))
        return 2;
      double target = ParseDouble(argv[4], "VALUE");
      double tolerance = std::fabs(ParseDouble(argv[5], "TOL"));
      long expectedCount = std::atol(argv[6]);
      long count = 0;
      for (size_t i = 0; i < values.size(); ++i)
        if (std::fabs(values[i] - target) <= tolerance)
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
      double tolerance = std::fabs(ParseDouble(argv[7], "TOL"));
      if (values.size() < k || reference.size() < k)
        {
          std::cout << "FAIL: need " << k << " values, have " << values.size() << " and " << reference.size() << std::endl;
          return 1;
        }
      std::sort(values.begin(), values.end());
      std::sort(reference.begin(), reference.end());
      double worst = 0.0;
      size_t worstIndex = 0;
      for (size_t i = 0; i < k; ++i)
        {
          double diff = std::fabs(values[i] - reference[i]);
          if (diff > worst)
            {
              worst = diff;
              worstIndex = i;
            }
        }
      std::cout << k << " lowest eigenvalues, max |diff| = " << worst << " at index " << worstIndex
                << " (" << values[worstIndex] << " vs " << reference[worstIndex] << "), tolerance " << tolerance << std::endl;
      bool pass = worst <= tolerance;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "numbers" && argc == 5)
    {
      std::vector<double> reference;
      size_t n1 = 0, n2 = 0;
      if (!ReadNumericLines(file, values, n1) || !ReadNumericLines(argv[3], reference, n2))
        return 2;
      double tolerance = std::fabs(ParseDouble(argv[4], "TOL"));
      if (values.size() != reference.size() || n1 != n2)
        {
          std::cout << "FAIL: " << n1 << " numeric lines / " << values.size() << " numbers, reference has "
                    << n2 << " / " << reference.size() << std::endl;
          return 1;
        }
      double worst = 0.0;
      size_t worstIndex = 0;
      for (size_t i = 0; i < values.size(); ++i)
        {
          double diff = std::fabs(values[i] - reference[i]);
          if (diff > worst)
            {
              worst = diff;
              worstIndex = i;
            }
        }
      std::cout << n1 << " numeric lines, " << values.size() << " numbers, max |diff| = " << worst << " at index "
                << worstIndex << " (" << values[worstIndex] << " vs " << reference[worstIndex] << "), tolerance "
                << tolerance << std::endl;
      bool pass = worst <= tolerance;
      std::cout << (pass ? "PASS" : "FAIL") << std::endl;
      return pass ? 0 : 1;
    }

  if (mode == "line" && argc == 6)
    {
      std::vector<double> reference;
      if (!ReadKeyedLine(file, argv[3], values) || !ReadKeyedLine(argv[4], argv[3], reference))
        return 2;
      double tolerance = std::fabs(ParseDouble(argv[5], "TOL"));
      if (values.size() != reference.size())
        {
          std::cout << "FAIL: " << values.size() << " numbers after '" << argv[3] << " =', reference has "
                    << reference.size() << std::endl;
          return 1;
        }
      double worst = 0.0;
      size_t worstIndex = 0;
      for (size_t i = 0; i < values.size(); ++i)
        {
          double diff = std::fabs(values[i] - reference[i]);
          if (diff > worst)
            {
              worst = diff;
              worstIndex = i;
            }
        }
      std::cout << values.size() << " numbers, max |diff| = " << worst << " at index " << worstIndex
                << " (" << values[worstIndex] << " vs " << reference[worstIndex] << "), tolerance "
                << tolerance << std::endl;
      bool pass = worst <= tolerance;
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
