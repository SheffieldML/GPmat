#ifndef NDLUTIL_H
#define NDLUTIL_H
#include <cmath>
#include <climits>
#include <cfloat>
#include "ndlfortran.h"



using namespace std;
namespace ndlutil {
  const double MATCHTOL = 1e-12;
  const double EPS=DBL_EPSILON;
  const double SMALLVAL=1e-6;
  const double GRADCHANGE=1e-6;
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
  double randn(); // normal deviates from N(0, 1.0)
  double rand(); // uniform deviates from (0.0, 1.0)

  /******** Mersenne Twister Random number generator ********/
  /* Period parameters */  
  const unsigned long MATRIX_A=0x9908b0dfUL;   /* constant vector a */
  const unsigned long UPPER_MASK=0x80000000UL; /* most significant w-r bits */
  const unsigned long LOWER_MASK=0x7fffffffUL; /* least significant r bits */
  const unsigned long N=624;
  const unsigned long M=397;
    
  static unsigned long mt[N]; /* the array for the state vector  */
  static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */
  void init_genrand(unsigned long s);
  void init_by_array(unsigned long init_key[], int key_length);
  unsigned long genrand_int32(void);
  long genrand_int31(void);
  double genrand_real1(void);
  double genrand_real2(void);
  double genrand_real3(void);
  double genrand_res53(void);
  
}
#endif
