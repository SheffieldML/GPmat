#ifndef NDLUTIL_H
#define NDLUTIL_H
#include <cmath>
#include "xpose.h"

using namespace std;

// Probability of a standard Gaussian (mean 0 and variance 1).
double ngaussian(double x);
// cumulative Gaussian distribution function.
double cumGaussian(double x);
// Gradient of the log fo the cumulative Gaussian distribution function.
double gradLnCumGaussian(double x);
double lnCumGaussian(double x);


#endif
