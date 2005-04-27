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
  const double SMALLVAL=1e-6;
  const double DISPEPS = 1e-14;
  const double LOGTWOPI = log(2*M_PI);
  const double HALFLOGTWOPI = 0.5*LOGTWOPI;
  const double HALFSQRTTWO = 0.7071067811865476;
  const double SQRTTWOPI = sqrt(2.0*M_PI);

  // Probability of a standard Gaussian (mean 0 and variance 1).
  double ngaussian(double x);
  // cumulative Gaussian distribution function.
  double cumGaussian(double x);
  double invCumGaussian(double x);
  // Gradient of the log fo the cumulative Gaussian distribution function.
  double gradLnCumGaussian(double x);
  double lnCumGaussian(double x);
  // The log of the weighted sum of two cumulative Gaussians.
  double lnCumGaussSum(double u1, double u2, double w1, double w2);
  // the sigmoid (logistic) function.
  double sigmoid(double x);
  // inverse of the sigmoid (logistic) function.
  double invSigmoid(double x);
  // inverse of error function in double precision from http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html
  double erfcinv(double x);

  // variations of the gamma function.
  double gamma(double x);
  double gammaln(double x);
  double digamma(double x);

  double xlogy(double x, double y);

}
#endif
