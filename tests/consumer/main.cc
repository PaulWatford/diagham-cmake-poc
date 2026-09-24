// Downstream use of the installed DiagHam library API.
//
// Diagonalises the tight-binding ring H = -sum_i (|i><i+1| + h.c.) on L
// sites with DiagHam's own Householder + QL routines and compares with the
// exact spectrum -2 cos(2 pi k / L), k = 0 .. L-1.

#include "config.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
  const int L = 6;
  RealSymmetricMatrix H(L, true);
  for (int i = 0; i < L; ++i)
    H.SetMatrixElement(i, (i + 1) % L, -1.0);

  RealTriDiagonalSymmetricMatrix T(L);
  H.Householder(T, 1e-14);
  T.Diagonalize();

  std::vector<double> computed, exact;
  for (int i = 0; i < L; ++i)
    {
      computed.push_back(T.DiagonalElement(i));
      exact.push_back(-2.0 * std::cos(2.0 * M_PI * i / L));
    }
  std::sort(computed.begin(), computed.end());
  std::sort(exact.begin(), exact.end());

  double worst = 0.0;
  for (int i = 0; i < L; ++i)
    worst = std::max(worst, std::fabs(computed[i] - exact[i]));
  std::cout << "tight-binding ring, L = " << L << ", max |E - E_exact| = " << worst << std::endl;
  if (worst > 1e-12)
    {
      std::cout << "FAIL" << std::endl;
      return 1;
    }
  std::cout << "PASS" << std::endl;
  return 0;
}
