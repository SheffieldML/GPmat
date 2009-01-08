#ifndef NDLUTIL_H
#define NDLUTIL_H
#include <cmath>
#include <climits>
#include <cfloat>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530941723212145818
#endif
#ifdef _MSC_VER
#define __DBL_MIN__ 2.2250738585072014e-308
#define isnan(x) _isnan(x)
#define isinf(x) !_finite(x)
#pragma warning(disable:4018) // < signed/unsigned mismatch
#pragma warning(disable:4267) // conversion from size_t to int
#pragma warning(disable:4800) // forcing double to bool (performance warning)
#pragma warning(disable:4244) // conversion from difference_type to int
#endif
#include "ndlassert.h"
#include "ndlfortran.h"
#include "ndlexceptions.h"
using namespace std;
namespace ndlutil {
  const double MATCHTOL = 1e-10;
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
  // cumulative Gaussian distribution function and its inverse.
  double cumGaussian(double x);
  double invCumGaussian(double x);
  // Gradient of the log of the cumulative Gaussian distribution function.
  double gradLnCumGaussian(double x);
  // The log of a cumulative Gaussian.
  double lnCumGaussian(double x);
  // The log of the weighted sum of two cumulative Gaussians.
  double lnCumGaussSum(double u1, double u2, double w1, double w2);
  // Computes the log of the difference between two cumulative Gaussians.
  double lnDiffCumGaussian(double u, double uprime);
  // A Gaussian over the difference between two cumulative Gaussians.
  double gaussOverDiffCumGaussian(double x, double xp, int order);
  // the sigmoid (logistic) function and its inverse.
  double sigmoid(double x);
  double invSigmoid(double x);
  // inverse of error function in double precision (from http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html)
  double erfcinv(double x);

  // variations of the gamma function.
  double gamma(double x);
  double gammaln(double x);
  double digamma(double x);

  // xlogy returns 0.0 if x is 0.0
  double xlogy(double x, double y);
  // normal deviates from N(0, 1.0)
  double randn();
  // uniform deviates from (0.0, 1.0)
  double rand(); 

  // Give a random perumation of numbers.
  vector<unsigned long> randperm(unsigned long maxVal);
  // give a truncation of a random permutation of numbers
  vector<unsigned long> randpermTrunc(unsigned long maxVal, unsigned long length);

  /******** Mersenne Twister Random number generator ********/
  /* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
  */

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
