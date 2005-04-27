#include <iostream>
#include "CMatrix.h"
using namespace std;

int main()
{
  int nrows = 3;
  int ncols = 3;
  double a[9] = {1.1, -1, 3, 4, 5.4, 6, 7, 8.2, 9};
  CMatrix A(nrows, ncols, a);
  A.inv();
  cout << "A = " << endl;
  cout << A;
  double b[9] = {18.0455,    4.9297,    3.6417,    4.9297,    5.9624,   -4.7233,    3.6417,   -4.7233,    7.8203};
  CMatrix B(nrows, ncols, b);
  B.setSymmetric(true);
  CMatrix C = chol(B);
  C.gemm(B, A, 1, 1, "n", "n");
  CMatrix D = inv(A);
  CMatrix E = multiply(A, D);
  cout << E;
  
  double x = 1.03;
  cout << ndlutil::digamma(x) << endl;
  cout << ndlutil::gammaln(x) << endl;
  cout << ndlutil::gamma(x) << endl;

}

