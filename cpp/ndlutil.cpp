#include "ndlutil.h"

namespace ndlutil {

  const double HALFSQRTTWO = 0.7071067811865476;

  double ngaussian(double x)
  {
    x *= x;
    x = exp(-.5*x);
    x = x/sqrt(2.0*M_PI);
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
      x = 1/(sqrt(2.0*M_PI)*0.5*derfcx_(-HALFSQRTTWO*x));
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

  double erfcinv(double x)
  {
    double s, t, u, w, y, z;
    
    z = x;
    if (x > 1) {
        z = 2 - x;
    }
    w = 0.916461398268964 - log(z);
    u = sqrt(w);
    s = (log(u) + 0.488826640273108) / w;
    t = 1 / (u + 0.231729200323405);
    y = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
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
    y += s * (y * s + 1);
    if (x > 1) {
        y = -y;
    }
    return y;
}


}
