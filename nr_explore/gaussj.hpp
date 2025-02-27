#pragma once

#include <vector>

namespace nr_explore {
  // Linear equation solution by Gauss-Jordan elimination, 
  // equation (2.1.1) above. 
  // The input matrix is a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is input 
  // containing the m right-hand side vectors. On output, a is replaced by its 
  // matrix inverse, and b is replaced by the corresponding set of solution vectors.
  template<class T>
  void gaussj(
    std::vector<std::vector<T>> &a, 
    std::vector<std::vector<T>> &b
  );

  template<class T>
  void make2dIdentity(    
    std::vector<std::vector<T>> &b
  );

  template<class T>
  void make3dIdentity(
    std::vector<std::vector<T>> &a
  );

  template<class T>
  void printMatrix(
    const std::vector<std::vector<T>> &a
  );

  template<class T>
  void printMatrices(
    const std::vector<std::vector<T>> &a, 
    const std::vector<std::vector<T>> &b
  );

  void testFunc();

  void testGaussj();
}