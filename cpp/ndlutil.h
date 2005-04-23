#ifndef NDLUTIL_H
#define NDLUTIL_H
#include <cmath>
#include <climits>
#include <cfloat>
#include "xpose.h"

using namespace std;
namespace ndlutil {
  const double MATCHTOL = 1e-12;
  const double EPS=DBL_EPSILON;
  const double DISPEPS = 1e-14;

  // Probability of a standard Gaussian (mean 0 and variance 1).
  double ngaussian(double x);
  // cumulative Gaussian distribution function.
  double cumGaussian(double x);
  double invCumGaussian(double x);
  // Gradient of the log fo the cumulative Gaussian distribution function.
  double gradLnCumGaussian(double x);
  double lnCumGaussian(double x);
  // inverse of error function in double precision from http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html
  double erfcinv(double y);
}
#endif
