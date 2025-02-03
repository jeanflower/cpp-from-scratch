#pragma once

#include <vector>

namespace nr_explore {
  void gaussj(
  // Linear equation solution by Gauss-Jordan elimination, 
  // equation (2.1.1) above. 
  // The input matrix is a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is input 
  // containing the m right-hand side vectors. On output, a is replaced by its 
  // matrix inverse, and b is replaced by the corresponding set of solution vectors.
    std::vector<std::vector<double>> &a, 
    std::vector<std::vector<double>> &b);

  void testGaussj();
}