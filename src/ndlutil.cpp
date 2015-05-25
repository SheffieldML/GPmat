#include "ndlutil.h"

#ifdef _MSC_VER 
extern "C" double erf(double);
#endif /* _MSC_VER */

namespace ndlutil {

  const double ROBUSTADD=1e-300; // added in various points to prevent log of zero.
  double ngaussian(double x)
  {
    x *= x;
    x = exp(-.5*x);
    x = x/SQRTTWOPI;
    return x;
  }
  double cumGaussian(double x)
  {
    x *= HALFSQRTTWO;
    x = erf(x);
    x++;
    x*=0.5;
    return x;
  }
  double invCumGaussian(double x)
  {
    return -sqrt(2.0)*erfcinv(2.0*x);
  }
  double gradLnCumGaussian(double x)
  {
    if(x>0)
      x = ngaussian(x)/cumGaussian(x);
    else
      x = 1.0/(SQRTTWOPI*0.5*derfcx_(-HALFSQRTTWO*x));
    return x;
  }
  double lnCumGaussian(double x)
  {
    if(x<0)
      x = -.5*x*x + log(0.5) + log(derfcx_(-HALFSQRTTWO*x));
    else
      x = log(cumGaussian(x));
    return x;
  }

  double lnCumGaussSum(double u1, double u2, double w1, double w2)
  {
    
    if(u1 > 0 && u2 > 0)
      return log(w1*cumGaussian(u1)+ w2*cumGaussian(u2));
    
    else if(u1>u2)
      return log(w1) + lnCumGaussian(u1)
      + log(1.0 + w2/w1*exp(lnCumGaussian(u2)
			    -lnCumGaussian(u1)));
    else
      return log(w2) + lnCumGaussian(u2)
      + log(1.0 + w1/w2*exp(lnCumGaussian(u1)
			    -lnCumGaussian(u2)));
  }
  // f = log(\phi(u) - phi(uprime))
  double lnDiffCumGaussian(double u, double uprime)
  {
    double arg = gaussOverDiffCumGaussian(u, uprime,  1) + ROBUSTADD;
    double retVal = -log(arg) - .5*u*u - HALFLOGTWOPI;
    return retVal;
  }
  double gaussOverDiffCumGaussian(double x, double xp, int order)
  {
    double xp2 = xp*xp;
    double x2 = x*x;
    double expRatio = 0.0;
    switch(order)
    {
    case 1:
      expRatio = exp(.5*(x2-xp2));
      if(x<=0)
	return 2/(SQRTTWOPI*(derfcx_(-HALFSQRTTWO*x)-expRatio*derfcx_(-HALFSQRTTWO*xp)+ROBUSTADD));
      else
	return 2/(SQRTTWOPI*(expRatio*derfcx_(HALFSQRTTWO*xp)-derfcx_(HALFSQRTTWO*x)+ROBUSTADD));
      break;
    case 2:
      expRatio = exp(.5*(xp2-x2));
      if(x<=0)
	return 2/(SQRTTWOPI*(expRatio*derfcx_(-HALFSQRTTWO*x)-derfcx_(-HALFSQRTTWO*xp)+ROBUSTADD));
      else
	return 2/(SQRTTWOPI*(derfcx_(HALFSQRTTWO*xp)-expRatio*derfcx_(HALFSQRTTWO*x)+ROBUSTADD));
      break;
    default:
      throw ndlexceptions::Error("Incorrect order parameter in gaussOverDiffCumGaussian");
    }
  }
  double invSigmoid(double x)
  {
    return log(x/(1.0-x));
  }
  double sigmoid(double x)
  {
    return 1.0/(1.0+exp(-x));
  }
  
  double erfcinv(double x)
  {
    double s, t, u, w, y, z;
    
    z = x;
    if (x > 1.0) 
    {
      z = 2 - x;
    }
    w = 0.916461398268964 - log(z);
    u = sqrt(w);
    s = (log(u) + 0.488826640273108) / w;
    t = 1.0 / (u + 0.231729200323405);
    y = u * (1.0 - s * (s * 0.124610454613712 + 0.5)) - 
    ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
       0.150689047360223) * t + 0.116065025341614) * t + 
     0.499999303439796) * t;
    t = 3.97886080735226 / (y + 3.97886080735226);
    u = t - 0.5;
    s = (((((((((0.00112648096188977922 * u + 
		 1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
	       7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
	     0.00339721910367775861) * u - 0.011274916933250487) * u - 
	   0.0118598117047771104) * u + 0.0142961988697898018) * u + 
	 0.0346494207789099922) * u + 0.00220995927012179067;
    s = ((((((((((((s * u - 0.0743424357241784861) * u - 
		   0.105872177941595488) * u + 0.0147297938331485121) * u + 
		 0.316847638520135944) * u + 0.713657635868730364) * u + 
	       1.05375024970847138) * u + 1.21448730779995237) * u + 
	     1.16374581931560831) * u + 0.956464974744799006) * u + 
	   0.686265948274097816) * u + 0.434397492331430115) * u + 
	 0.244044510593190935) * t - 
    z * exp(y * y - 0.120782237635245222);
    y += s * (y * s + 1.0);
    if (x > 1.0) 
    {
      y = -y;
    }
    return y;
  }
  double gamma(double x)
  {
    double y;
    lgama_(1, x, y);
    return y;
  }
  double gammaln(double x)
  {
    double y;
    lgama_(0, x, y);
    return y;
  }
  double digamma(double x)
  {
    double y;
    psi_(x, y);
    return y;
  }
  
  double xlogy(double x, double y)
  {
    if(x==0.0)
      return 0.0;
    else
      return x*log(y);
  }
  double rand()
  {
    return genrand_real3();
  }
  double randn()
  {
    // uses the polar form of the Box-Mueller transform. 
    double x1, x2, w;
    static bool stored=false;
    static double storeVal;
    if(stored)
    {	
      stored = false;
      return storeVal;
    }
    else
    {
      do {
	x1 = 2.0 * genrand_real1() - 1.0;
	x2 = 2.0 * genrand_real1() - 1.0;
	w = x1 * x1 + x2 * x2;
      } while ( w >= 1.0 );
      
      w = sqrt( (-2.0 * log( w ) ) / w );
      storeVal = x1 * w;
      stored=true;
      return x2 * w;
    }
  }
  
  // Give a random perumation of numbers.
  vector<unsigned long> randpermTrunc(unsigned long maxVal, unsigned long length)
  {
    BOUNDCHECK(length<=maxVal);
    vector<unsigned long> perm;
    vector<unsigned long> indices;
    for(unsigned long i=0; i < maxVal; i++)
    {
      indices.push_back(i);
    }
    for(unsigned long i=0; i<length; i++)
    {
      unsigned long ind = (unsigned long)(ndlutil::rand()*indices.size());
      perm.push_back(indices[ind]);
      indices.erase(indices.begin()+ind);
    }
    return perm;
  }
  // give a truncation of a random permutation of numbers
  vector<unsigned long> randperm(unsigned long maxVal)
  {
    return randpermTrunc(maxVal, maxVal);
  }
  /******* Mersenne Twister Random number generator. **************/
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


  /* initializes mt[N] with a seed */
  void init_genrand(unsigned long s)
  {
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
      mt[mti] = 
      (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti] &= 0xffffffffUL;
      /* for >32 bit machines */
    }
  }
  
  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  /* slight change for C++, 2004/2/26 */
  void init_by_array(unsigned long init_key[], int key_length)
  {
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
      + init_key[j] + j; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++; j++;
      if (i>=N) { mt[0] = mt[N-1]; i=1; }
      if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
      - i; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }
    
    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
  }
  
  /* generates a random number on [0,0xffffffff]-interval */
  unsigned long genrand_int32(void)
  {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
      int kk;
      
      if (mti == N+1)   /* if init_genrand() has not been called, */
	init_genrand(5489UL); /* a default initial seed is used */
      
      for (kk=0;kk<N-M;kk++) {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for (;kk<N-1;kk++) {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
      mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
      
      mti = 0;
    }
    
    y = mt[mti++];
    
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    
    return y;
  }
  
  /* generates a random number on [0,0x7fffffff]-interval */
  long genrand_int31(void)
  {
    return (long)(genrand_int32()>>1);
  }
  
  /* generates a random number on [0,1]-real-interval */
  double genrand_real1(void)
  {
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
  }
  
  /* generates a random number on [0,1)-real-interval */
  double genrand_real2(void)
  {
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
  }
  
  /* generates a random number on (0,1)-real-interval */
  double genrand_real3(void)
  {
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
  }
  
  /* generates a random number on [0,1) with 53-bit resolution*/
  double genrand_res53(void) 
  { 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
  } 
  /* These real versions are due to Isaku Wada, 2002/01/09 added */
  
  /*int main(void)
    {
    int i;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
    printf("1000 outputs of genrand_int32()\n");
    for (i=0; i<1000; i++) {
    printf("%10lu ", genrand_int32());
    if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand_real2()\n");
    for (i=0; i<1000; i++) {
    printf("%10.8f ", genrand_real2());
    if (i%5==4) printf("\n");
    }
    return 0;
    }*/

}

//#ifdef _MSC_VER
extern "C" double rand_(const int &flag)
{
  // Generate random numbers using MersenneTwister code from 
  //  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
  if (flag == 1) {
    ndlutil::init_genrand(0);
  }
  if (flag != 0) {
    ndlutil::init_genrand(flag);
  }
  return ndlutil::genrand_real2();
}
//#endif
