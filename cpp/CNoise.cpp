#include "CNoise.h"

using namespace std;
void CNoise::getGradTransParams(CMatrix& g, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y)
{
  CMatrix params(1, nParams);
  getParams(params);
  getGradParams(g, mu, varsigma, y);
  double param = 0.0;
  double gval=0.0;
  for(int i=0; i<transIndex.size(); i++)
    {
      param = params.getVals(transIndex[i]);
      gval = g.getVals(transIndex[i]);
      param = transforms[i]->gradfact(param);
      g.setVals(param*gval, transIndex[i]);
    }
  
}
void CNoise::getTransParams(CMatrix& params) const
{
  getParams(params);
  double param = 0.0;
  for(int i=0; i<transIndex.size(); i++)
    {
      param = params.getVals(transIndex[i]);
      param = transforms[i]->xtoa(param);
      params.setVals(param, transIndex[i]);
    }
}
void CNoise::setTransParams(const CMatrix& params)
{
  double param = 0.0;
  CMatrix paramCopy(params.getRows(), params.getCols());
  paramCopy.deepCopy(params);
  for(int i=0; i<transIndex.size(); i++)
    {
      param = paramCopy.getVals(transIndex[i]);
      param = transforms[i]->atox(param);
      paramCopy.setVals(param, transIndex[i]);
    }
  setParams(paramCopy);
}

CGaussianNoise::~CGaussianNoise()
{
}
void CGaussianNoise::setInitParam()
{
  setNParams(getNProcess()+1);
  sigma2 = 1e-6;
  // transform sigma2 (the last parameter).
  addTransform(new CNegLogLogitTransform(), getNParams()-1);
  setLogConcave(true);
  setSpherical(true);
  setMissing(false);
}
void CGaussianNoise::setInitParam(const int numProcess)
{
  setNProcess(numProcess);
  CMatrix b(1, getNProcess(), 0.0);
  bias.deepCopy(b);
  setInitParam();
}
void CGaussianNoise::setInitParam(const CMatrix& y)
{
  setNProcess(y.getCols());
  bias.deepCopy(meanRow(y));
  setInitParam();
}
ostream& CGaussianNoise::display(ostream& os)
{
  double b = 0.0;
  for(int j=0; j<bias.getCols(); j++)
    {
      b = bias.getVals(j);
      os << "Gaussian bias on process " << j << ": " << b << endl; 
    }
  os << "Gaussian noise: " << sigma2 << endl;
  return os;
}
void CGaussianNoise::setParams(const CMatrix& params)
{
  for(int j=0; j<bias.getCols(); j++)
    {
      bias.setVals(params.getVals(j), j);
    }
  sigma2 = params.getVals(nParams-1);
}
void CGaussianNoise::getParams(CMatrix& params) const
{
  for(int j=0; j<bias.getCols(); j++)
    params.setVals(bias.getVals(j), j);
  params.setVals(sigma2, nParams-1);
}
 
void CGaussianNoise::getGradParams(CMatrix& g, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const
{
  double nu=0.0;
  double u=0.0;
  double gsigma2=0.0;
  for(int j=0; j<y.getCols(); j++)
    {
      double gbias = 0.0;
      double b = bias.getVals(j);
      for(int i=0; i<y.getRows(); i++)
	{
	  nu = 1/(varsigma.getVals(i, j)+sigma2);
	  u=y.getVals(i, j) - mu.getVals(i, j)-b;
	  u*=nu;
	  gbias+=u;
	  gsigma2+=nu-u*u;
	}
      g.setVals(gbias, 0, j);
    }
  g.setVals(-0.5*gsigma2, 0, nParams-1);
}
void CGaussianNoise::getGradInputs(CMatrix& dlnZ_dmu, CMatrix& dlnZ_dvs, 
				   const CMatrix& mu, const CMatrix& varsigma, 
				   const CMatrix& y) const
{
  // gradient with respect to mu and varsigma of log likelihood.  
  double muval;
  double vsval;
  for(int j=0; j<y.getCols(); j++)
    {
      for(int i=0; i<y.getRows(); i++)
	{
	  muval = -bias.getVals(j);
	  vsval = 1/(sigma2+varsigma.getVals(i, j));
	  muval += y.getVals(i, j)-mu.getVals(i, j);
	  muval *= vsval;
	  vsval = 0.5*(muval*muval - vsval);
	  dlnZ_dmu.setVals(muval, i, j);
	  dlnZ_dvs.setVals(vsval, i, j);
	}
    }
}
  
void CGaussianNoise::getNuG(CMatrix& g, CMatrix& nu, 
			    const CMatrix& mu, const CMatrix& varsigma, 
			    const CMatrix& y, const int index) const
{
  double nuval=0.0;
  double gval=0.0;
  for(int j=0; j<y.getCols(); j++)
    {
      nuval=1./(sigma2+varsigma.getVals(index, j));
      nu.setVals(nuval, index, j);
      gval=y.getVals(index, j)-mu.getVals(index, j)-bias.getVals(j);
      g.setVals(gval*nuval, index, j);
    }
}
void CGaussianNoise::updateSites(CMatrix& m, CMatrix& beta, const int actIndex, 
				 const CMatrix& g, const CMatrix& nu, 
				 const CMatrix& mu, const CMatrix& varsigma, 
				 const CMatrix& y,
				 const int index) const
{
  for(int j=0; j<y.getCols(); j++)
    {
      m.setVals(y.getVals(index, j)-bias.getVals(j), actIndex, j);
      beta.setVals(1/sigma2, actIndex, j);
    }
}
CMatrix CGaussianNoise::out(const CMatrix& mu, const CMatrix& varsigma) const
{
  CMatrix y(mu.getRows(), mu.getCols());
  y.deepCopy(mu);
  y+=bias;
  return y;
}
CMatrix CGaussianNoise::likelihood(const CMatrix& mu, const CMatrix& varsigma, 
				   const CMatrix& y) const
{
  CMatrix L(mu.getRows(), mu.getCols());
  double arg=0.0;
  double var=0.0;
  for(int i=0; i<mu.getRows(); i++)
    {
      for(int j=0; j<mu.getCols(); j++)
	{
	  arg = y.getVals(i, j) - mu.getVals(i, j) - bias.getVals(j);
	  arg *= arg;
	  var = varsigma.getVals(i, j) + sigma2;
	  arg = 1/sqrt(2*M_PI*var)*exp(-.5*arg*arg/var);
	  L.setVals(arg, i, j);
	}
    }
  return L;
}

double CGaussianNoise::logLikelihood(const CMatrix& mu, const CMatrix& varsigma, 
				     const CMatrix& y) const
{
  double arg=0.0;
  double var=0.0;
  double L=0.0;
  for(int i=0; i<mu.getRows(); i++)
    {
      for(int j=0; j<mu.getCols(); j++)
	{
	  arg = y.getVals(i, j) - mu.getVals(i, j) - bias.getVals(j);
	  arg *= arg;
	  var = varsigma.getVals(i, j) + sigma2;
	  arg = arg/var;
	  L += log(var)+arg;
	}
    }  
  L += mu.getRows()*mu.getCols()*log(2*M_PI);
  L *= -0.5;
  return L;
}
  
  
