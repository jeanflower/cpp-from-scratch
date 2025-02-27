#include <iostream>
#include <iomanip> // for setprecision

#include "gaussj.hpp"

namespace nr_explore{

  template<class T>
  void SWAP(T& a, T& b) {
      const T c = a;
      a = b;
      b = c;
  }

  // inputs should be same-size square matrices
  template<class T>
  void printMatrix(
    const std::vector<std::vector<T>> &a
   ) {
    std::cout << std::fixed << std::setprecision(6) << std::setw(5);
    for (int i = 0; i < a.size(); i++) {
      for (int j = 0; j < a.size(); j++) {
        std::cout << " \t" << a[i][j];
      }
      std::cout << "\n";
    }
    std::cout << std::defaultfloat;
  }
  template void printMatrix<double>(const std::vector<std::vector<double>> &);
  template void printMatrix<float>(const std::vector<std::vector<float>> &);

  // inputs should be same-size square matrices
  template<class T>
  void printMatrices(
    const std::vector<std::vector<T>> &a, 
    const std::vector<std::vector<T>> &b
  ) {
    std::cout << std::fixed << std::setprecision(6) << std::setw(5);
    std::cout << "\ta\t";
    for (int j = 0; j < a[0].size(); j++) {
      std::cout << "\t\t";
    }
    std::cout << "b\n";

    for (int i = 0; i < a.size(); i++) {
      for (int j = 0; j < a[0].size(); j++) {
        std::cout << " \t" << a[i][j];
      }
      std::cout << " \t";
      for (int j = 0; j < b[0].size(); j++) {
        std::cout << " \t" << b[i][j];
      }
      std::cout << "\n";
    }
    std::cout << std::defaultfloat;
  }
  template void printMatrices<double>(
    const std::vector<std::vector<double>> &,
    const std::vector<std::vector<double>> &
  );
  template void printMatrices<float>(
    const std::vector<std::vector<float>> &,
    const std::vector<std::vector<float>> &
  );

  void testFunc(){
    std::cout << "hello\n";
  }

  // a is n rows by m cols, b is m rows by k cols, output is resized accordingly
  void multiply(
    const std::vector<std::vector<double>> &a, 
    const std::vector<std::vector<double>> &b,
    std::vector<std::vector<double>> &out
  ){
    const int n = a.size();
    const int m = a[0].size();
    const int k = b[0].size();

    out.resize(n);
    for (int i = 0; i < n; i++) {
      out[i].resize(k);
      for (int j = 0; j < k; j++) {
        out[i][j] = 0.0;
        for (int l = 0; l < m; l++) {
          out[i][j] += a[i][l] * b[l][j];
        }
      }
    }
  }

  void makeCopy(
    const std::vector<std::vector<double>> &a, 
    std::vector<std::vector<double>> &out
  ){
    const int n = a.size();
    const int m = a[0].size();

    out.resize(n);
    for (int i = 0; i < n; i++) {
      out[i].resize(m);
      for (int j = 0; j < m; j++) {
        out[i][j] = a[i][j];
      }
    }
  }

  template<class T>
  void gaussj( // on p44 of my NR book
  // Linear equation solution by Gauss-Jordan elimination, 
  // The input matrix is a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is input 
  // containing the m right-hand side vectors. On output, a is replaced by its 
  // matrix inverse, and b is replaced by the corresponding set of solution vectors.
    std::vector<std::vector<T>> &a, 
    std::vector<std::vector<T>> &b
  ) {
    const bool printDebug = false;
    int icol, irow, ll, n = a.size(), m = b[0].size();
    T dum, pivinv;
    std::vector<int> indxc(n), indxr(n), ipiv(n);  // These integer arrays are used for bookkeeping on the pivoting.
    // clear the ipiv array
    for (int j = 0; j < n; j++) {
      ipiv[j] = 0; // no rows have been pivoted
    }
    for (int i = 0; i < n; i++) { // This is the main loop over the columns to be reduced
      if (printDebug) {
        std::cout << "start of main loop ----- !!!!!!!!!!!!!!! \n";
        std::cout << "pivot data is " << "\n";
        for (int j = 0; j < ipiv.size(); j++) {
          std::cout << ipiv[j];
          if (j < ipiv.size() - 1) {
            std::cout << ", ";
          }
        }
        std::cout << "\n";
      }
      T big = 0.0;
      for (int j = 0; j < n; j++) { // This is the outer loop of the search for a pivot element.
        if (ipiv[j] != 1) {
          for (int k=0; k<n; k++) {
            if (ipiv[k] == 0) {
              if (abs(a[j][k]) >= big) {
                big = abs(a[j][k]);
                if (printDebug) {
                  std::cout << "big = " << big << " from [" << j << "][" << k << "]\n";
                }
                irow = j;
                icol = k;
              }
            }
          }
        }
      }
      ++(ipiv[icol]); // this line was, wrongly, in the loop above but it should be set once icol is known

      if (printDebug) {
        std::cout << "biggest value is " << big << " found at irow = " << irow << ", and icol = " << icol <<"\n";
      }

      // We now have the pivot element, so we interchange rows,
      // if needed, to put the pivot element on the diagonal.
      // The columns are not physically interchanged, only relabeled: 
      // indxc[i], the column of the .i C 1/th pivot element,
      // is the .i C 1/th column that is reduced, while
      // indxr[i] is the row in which that pivot element was 
      // originally located. If indxr[i] ¤ indxc[i], there is an 
      // implied column interchange. With this form of bookkeeping, 
      // the solution b’s will end up in the correct order, 
      // and the inverse matrix will be scrambled by columns.
      if (irow != icol) {
        for (int l=0;l<n;l++) 
          SWAP(a[irow][l], a[icol][l]);
        for (int l=0;l<m;l++) 
          SWAP(b[irow][l], b[icol][l]);
      }

      if (printDebug) {
        std::cout << "after swapping  " << big << " is on the diagonal ------------------- \n";
        printMatrices(a, b);
      }

      indxr[i] = irow; // We are now ready to divide the pivot row by the
      indxc[i] = icol; // pivot element, located at irow and icol.
      
      // look at the diagonal element
      if (a[icol][icol] == 0.0) {
        throw("gaussj: Singular Matrix"); 
      }

      pivinv = 1.0/a[icol][icol];
      a[icol][icol] = 1.0;
      for (int l = 0; l<n; l++) {
        if (l == icol) continue; // this important line was missing from my NR book
        a[icol][l] *= pivinv;
      }
      for (int l=0;l<m;l++) {
        b[icol][l] *= pivinv;
      }

      if (printDebug) {
        std::cout << "after dividing ------------------- \n";
        printMatrices(a, b);
      }

      for (int ll=0; ll<n; ll++) { //Next, we reduce the rows...
        if (ll != icol) {  // ...except for the pivot one, of course. 
          //std::cout << "Reduce " << ll << "th row\n";
          dum = a[ll][icol];
          //std::cout << "subtract " << dum << " * pivot row\n";
          a[ll][icol] = 0.0;

          for (int l=0; l<n;l++) {
            if (l == icol) continue; // this important line was missing from my NR book
            a[ll][l] -= a[icol][l] * dum;
          }
          for (int l=0; l<m;l++) {
            b[ll][l] -= b[icol][l] * dum; 
          }
        }
      }
      if (printDebug) {
        std::cout << "after reducing ------------------- \n";
        printMatrices(a, b);
      }
    }
    // This is the end of the main loop over columns of the reduction.
    // It only remains to unscramble the solution in view of the column
    // interchanges. We do this by interchanging pairs of columns in the
    // reverse order that the permutation was built up.
    for (int l = n-1; l >= 0; l--) {
      if (indxr[l] != indxc[l]) {
        for (int k=0; k<n; k++) {
          SWAP(a[k][indxr[l]], a[k][indxc[l]]);
        }
      }
    }
    if (printDebug) {
      std::cout << "after unscrambling ------------------- \n";
      printMatrices(a, b);
    }
    // And we are done.
  }
  template void gaussj<double>(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
  template void gaussj<float>(std::vector<std::vector<float>> &, std::vector<std::vector<float>> &);


  template<class T>
  void make1dIdentity(    
    std::vector<std::vector<T>> &b
  ) {
    std::vector<T> b0 = {1};
    b.resize(1);
    b[0] = b0;
  }

  template<class T>
  void make2dIdentity(    
    std::vector<std::vector<T>> &b
  ) {
    std::vector<T> b0 = {1, 0};
    std::vector<T> b1 = {0, 1};
    b.resize(2);
    b[0] = b0;
    b[1] = b1;
  }
  template void make2dIdentity<double>(std::vector<std::vector<double>> &);
  template void make2dIdentity<float>(std::vector<std::vector<float>> &);

  template<class T>
  void make3dIdentity(
    std::vector<std::vector<T>> &a
  ) {
    std::vector<T> a0 = {1, 0, 0};
    std::vector<T> a1 = {0, 1, 0};
    std::vector<T> a2 = {0, 0, 1};
    a.resize(3);
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
  }
  template void make3dIdentity<double>(std::vector<std::vector<double>> &);
  template void make3dIdentity<float>(std::vector<std::vector<float>> &);

  void make1dExample(
    std::vector<std::vector<double>> &a
  ) {
    std::vector<double> a0 = {3};
    a.resize(1);
    a[0] = a0;
  }

  void make2dExample1(
    std::vector<std::vector<double>> &a
  ) {
    std::vector<double> a0 = {0, 3};
    std::vector<double> a1 = {2, 0};
    a.resize(2);
    a[0] = a0;
    a[1] = a1;
  }

  void make2dExample2(
    std::vector<std::vector<double>> &a
  ) {
    std::vector<double> a0 = {1, 3};
    std::vector<double> a1 = {2, 1};
    a.resize(2);
    a[0] = a0;
    a[1] = a1;
  }

  void make3dExample1(
    std::vector<std::vector<double>> &a
  ) {
    std::vector<double> a0 = {0, 3, 0};
    std::vector<double> a1 = {2, 0, 0};
    std::vector<double> a2 = {0, 0, 4};
    a.resize(3);
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
  }

  void make3dExample2(
    std::vector<std::vector<double>> &a
  ) {
    std::vector<double> a0 = {1, 3, 3};
    std::vector<double> a1 = {2, 2, 2};
    std::vector<double> a2 = {2, 2, 4};
    a.resize(3);
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
  }

  ///////////////////////// ENTER TEST CODE HERE

  void testGaussj() {
    std::vector<std::vector<double>> a;
    std::vector<std::vector<double>> b;
    std::vector<std::vector<double>> acopy;
    std::vector<std::vector<double>> product = {};

    //make_1d_example(a);
    //make_1d_identity(b);

    //make_2d_example1(a); // zeros on diagonal, non singular matrix
    //make_2d_example2(a); // all non zero, non singular matrix
    //make_2d_identity(b);

    //make_3d_example1(a); // three non-zero entries, non singular matrix
    make3dExample2(a); // all non-zero entries, non singular matrix
    make3dIdentity(b);

    makeCopy(a, acopy);
    /*
    multiply(b, b, product);
    std::cout << "Multiply: identity is \n";
    printMatrix(product);

    std::cout << "Multiple: this times identity\n";
    printMatrix(a);
    multiply(a, b, product);
    std::cout << "is\n";
    printMatrix(product);
    */
    
    std::cout << "Task: find the inverse of\n";
    printMatrix(a);

    gaussj(a, b);

    // this is the inverse
    std::cout << "Outcome: the inverse is\n";
    printMatrix(b);

    multiply(acopy, b, product);

    std::cout << "Check: the product is\n";
    printMatrix(product);
  }

}