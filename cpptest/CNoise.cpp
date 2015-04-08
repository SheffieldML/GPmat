#include "CNoise.h"

using namespace std;
void CNoise::getNuG(CMatrix& g, CMatrix& nu, const int index) const
{
  assert(g.dimensionsMatch(nu));
  assert(g.getRows()==nData);
  assert(g.getCols()==nProcesses);
  assert(index>=0);
  assert(index<nData);
  double nuval=0.0;
  double gval=0.0;
  
  for(int j=0; j<g.getCols(); j++)
    {
      getGradInputs(gval, nuval, index, j);
      nuval = gval*gval - 2*nuval;
      if(nuval<0)
	{
	  if(isLogConcave())
	    {
	      throw "Log concave noise model has value of nu < 0.";
	    }
	  else 
	    {
	      nuval=ndlutil::SMALLVAL;
	    }
	}
      if(isnan(nuval)) 
	throw "Nu is NaN";
      if(abs(nuval)<ndlutil::SMALLVAL)
	nuval=ndlutil::EPS;
      nu.setVal(nuval, index, j);
      
      g.setVal(gval, index, j);
    }
  
  
}

void CNoise::updateSites(CMatrix& m, CMatrix& beta, const int actIndex, 
				 const CMatrix& g, const CMatrix& nu, 
				 const int index) const
{
  assert(index>=0);
  assert(index<nData);
  assert(m.dimensionsMatch(beta));
  assert(actIndex>=0);
  assert(actIndex<m.getRows());
  assert(g.dimensionsMatch(nu));
  assert(m.getCols()==g.getCols());
  double nuVal=0.0;
  double gVal=0.0;
  for(int j=0; j<m.getCols(); j++)
    {
      nuVal = nu.getVal(index, j);
      m.setVal(getMu(index, j) + g.getVal(index, j)/nuVal, actIndex, j);
      beta.setVal(nuVal/(1-nuVal*getVarSigma(index, j)), actIndex, j);
    }
}
mxArray* CNoise::toMxArray() const
{
  int dims[1];
  dims[0] = 1;
  const char *fieldNames[] = {"type", "transforms", "numProcess", "spherical", "nParams", "missing", "logconcave"};

  mxArray* matlabArray = mxCreateStructArray(1, dims, 7, fieldNames);
    
  // type field.
  const char *typeName[1];
  string ty=getType();
  typeName[0] = ty.c_str();
  mxSetField(matlabArray, 0, "type", 
	     mxCreateCharMatrixFromStrings(1, typeName));
  
  // transforms field.
  mxSetField(matlabArray, 0, "transforms", transformsToMxArray());

  // inputDimension field.
  mxSetField(matlabArray, 0, "numProcess", convertMxArray((double)nProcesses));

  // nParams field.
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)nParams));
    
  // spherical field
  mxSetField(matlabArray, 0, "spherical", convertMxArray((double)spherical));

  // missing field
  mxSetField(matlabArray, 0, "missing", convertMxArray((double)missing));
  
  // logConcave field
  mxSetField(matlabArray, 0, "logconcave", convertMxArray((double)logConcave));
// priors field.
  /// if priors exist need to add the field.

  // Noise specific code.
  addParamToMxArray(matlabArray);
  return matlabArray;

}
void CNoise::fromMxArray(const mxArray* matlabArray) 
{
  string mxType = mxArrayExtractStringField(matlabArray, "type");
  if(mxType!=type)
    cerr << "Error mismatch between saved type, " << mxType << ", and Class type, " << type << "." << endl;
  
  mxArray* transformArray = mxArrayExtractMxArrayField(matlabArray, "transforms");
  if(transformArray!=NULL)
    transformsFromMxArray(transformArray);
  nProcesses = mxArrayExtractIntField(matlabArray, "numProcess");
  spherical = mxArrayExtractBoolField(matlabArray, "spherical");
  logConcave = mxArrayExtractBoolField(matlabArray, "logconcave");
  missing = mxArrayExtractBoolField(matlabArray, "missing");
  
  // TODO priors ... need to deal with priors
  extractParamFromMxArray(matlabArray);
}
  
  
void CNoise::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "nParams");
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)nParams));
  string pName;
  for(int i=0; i<nParams; i++)
    {
      pName = getParamName(i);
      mxAddField(matlabArray, pName.c_str());      
      mxSetField(matlabArray, 0, pName.c_str(), convertMxArray(getParam(i))); 
    } 
}

void CNoise::extractParamFromMxArray(const mxArray* matlabArray) 
{
  nParams = mxArrayExtractIntField(matlabArray, "nParams");
  string pName;
  for(int i=0; i<nParams; i++)
    {
      pName=getParamName(i);
      setParam(mxArrayExtractDoubleField(matlabArray, pName), i);
    }

}
bool CNoise::equals(const CNoise& noise, const double tol) const
{
  if(getType()!=noise.getType())
    return false;
  if(getNumParams()!=noise.getNumParams())
    return false;
  CMatrix params(1, getNumParams());
  getParams(params);
  CMatrix noiseParams(1, getNumParams());
  noise.getParams(noiseParams);
  if(!params.equals(noiseParams, tol))
    return false;
  return true;
}

CGaussianNoise::~CGaussianNoise()
{
}
void CGaussianNoise::setInitParam()
{
  setType("gaussian");
  setNoiseName("Gaussian");
  setNumParams(getNumProcesses()+1);
  mu.resize(y.getRows(), y.getCols());
  varSigma.resize(y.getRows(), y.getCols());
  mu.zeros();
  varSigma.zeros();
  bias.deepCopy(meanRow(y));
  for(int j=0; j<y.getCols(); j++)
    setParamName("bias" + j, j);
  sigma2 = 1e-6;
  setParamName("sigma2", getNumProcesses());
  clearTransforms();
  // transform sigma2 (the last parameter).
  addTransform(new CNegLogLogitTransform(), getNumParams()-1);
  setLogConcave(true);
  setSpherical(true);
  setMissing(false);
}
ostream& CGaussianNoise::display(ostream& os)
{
  double b = 0.0;
  os << "Gaussian Noise: " << endl; 
  for(int j=0; j<bias.getCols(); j++)
    {
      b = bias.getVal(j);
      os << "Bias on process " << j << ": " << b << endl; 
    }
  os << "Variance: " << sigma2 << endl;
  return os;
}
void CGaussianNoise::setParam(const double val, const int index)
{
  assert(index>=0);
  assert(index<getNumParams());
  if(index<getNumProcesses())
    bias.setVal(val, index);
  else
    sigma2=val;
}  
void CGaussianNoise::setParams(const CMatrix& params)
{
  assert(getNumParams()==getNumProcesses()+1);
  assert(params.getCols()==getNumParams());
  assert(params.getRows()==1);
  assert(getNumProcesses()==bias.getCols());
  for(int j=0; j<bias.getCols(); j++)
    {
      bias.setVal(params.getVal(j), j);
    }
  sigma2 = params.getVal(getNumParams()-1);
}
double CGaussianNoise::getParam(const int index) const
{
  assert(index>=0);
  assert(index<getNumParams());
  if(index<getNumProcesses())
    return bias.getVal(index);
  else
    return sigma2;

}
void CGaussianNoise::getParams(CMatrix& params) const
{
  assert(getNumParams()==getNumProcesses()+1);
  assert(params.getCols()==getNumParams());
  assert(params.getRows()==1);
  assert(getNumProcesses()==bias.getCols());
  for(int j=0; j<getNumProcesses(); j++)
    params.setVal(bias.getVal(j), j);
  params.setVal(sigma2, getNumParams()-1);
}
 
void CGaussianNoise::getGradParams(CMatrix& g) const
{
  assert(g.getCols()==getNumParams());
  assert(g.getRows()==1);
  double nu=0.0;
  double u=0.0;
  double gsigma2=0.0;
  double b=0.0;
  for(int j=0; j<y.getCols(); j++)
    {
      double gbias = 0.0;
      b = bias.getVal(j);
      for(int i=0; i<y.getRows(); i++)
	{
	  nu = 1/(getVarSigma(i, j)+sigma2);
	  u=getTarget(i, j) - getMu(i, j)-b;
	  u*=nu;
	  gbias+=u;
	  gsigma2+=nu-u*u;
	}
      g.setVal(gbias, 0, j);
    }
  g.setVal(-0.5*gsigma2, 0, getNumParams()-1);
}

void CGaussianNoise::getGradInputs(double& gmu, double& gvs, const int i, const int j) const
{
  gmu = -bias.getVal(j);
  gvs = 1/(sigma2+getVarSigma(i, j));
  gmu += getTarget(i, j)-getMu(i, j);
  gmu *= gvs;
  gvs = 0.5*(gmu*gmu - gvs);
}
  
void CGaussianNoise::getNuG(CMatrix& g, CMatrix& nu, const int index) const
{
  double nuval=0.0;
  double gval=0.0;
  for(int j=0; j<y.getCols(); j++)
    {
      nuval=1./(sigma2+getVarSigma(index, j));
      if(isnan(nuval))
	{
	  cout << "Sigma2 " << sigma2 << endl;
	  cout << "varSigma " << getVarSigma(index, j) << endl;
	}
      assert(!isnan(nuval));
      nu.setVal(nuval, index, j);
      gval=getTarget(index, j)-getMu(index, j)-bias.getVal(j);
      g.setVal(gval*nuval, index, j);
    }
}
void CGaussianNoise::updateSites(CMatrix& m, CMatrix& beta, const int actIndex, 
				 const CMatrix& g, const CMatrix& nu, 
				 const int index) const
{
  for(int j=0; j<y.getCols(); j++)
    {
      m.setVal(getTarget(index, j)-bias.getVal(j), actIndex, j);
      beta.setVal(1/sigma2, actIndex, j);
    }
}
void CGaussianNoise::out(CMatrix& yTest, const CMatrix& muTest, const CMatrix& varSigmaTest) const
{
  assert(yTest.dimensionsMatch(muTest));
  assert(muTest.dimensionsMatch(varSigmaTest));
  yTest.deepCopy(muTest);
  yTest+=bias;
}
void CGaussianNoise::likelihood(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, 
				   const CMatrix& yTest) const
{
  assert(yTest.getCols()==getNumProcesses());
  assert(L.dimensionsMatch(muTest));
  assert(yTest.dimensionsMatch(muTest));
  assert(muTest.dimensionsMatch(varSigmaTest));
  double arg=0.0;
  double var=0.0;
  for(int i=0; i<muTest.getRows(); i++)
    {
      for(int j=0; j<muTest.getCols(); j++)
	{
	  arg = yTest.getVal(i, j) - muTest.getVal(i, j) - bias.getVal(j);
	  arg *= arg;
	  var = varSigmaTest.getVal(i, j) + sigma2;
	  arg = 1/sqrt(2*M_PI*var)*exp(-.5*arg*arg/var);
	}
    }
}

double CGaussianNoise::logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, 
				     const CMatrix& yTest) const
{
  assert(yTest.getCols()==getNumProcesses());
  assert(yTest.dimensionsMatch(muTest));
  assert(yTest.dimensionsMatch(varSigmaTest));
  double arg=0.0;
  double var=0.0;
  double L=0.0;
  for(int i=0; i<muTest.getRows(); i++)
    {
      for(int j=0; j<muTest.getCols(); j++)
	{
	  arg = yTest.getVal(i, j) - muTest.getVal(i, j) - bias.getVal(j);
	  arg *= arg;
	  var = varSigmaTest.getVal(i, j) + sigma2;
	  arg = arg/var;
	  L += log(var)+arg;
	}
    }  
  L += muTest.getRows()*muTest.getCols()*log(2*M_PI);
  L *= -0.5;
  return L;
}

void CGaussianNoise::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "nParams");
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)getNumParams()));
  mxAddField(matlabArray, "bias");
  mxSetField(matlabArray, 0, "bias", bias.toMxArray());
  mxAddField(matlabArray, "sigma2");
  mxSetField(matlabArray, 0, "sigma2", convertMxArray(sigma2));
}

void CGaussianNoise::extractParamFromMxArray(const mxArray* matlabArray) 
{
  setNumParams(mxArrayExtractIntField(matlabArray, "nParams"));
  mxArray* biasField = mxArrayExtractMxArrayField(matlabArray, "bias");
  bias.fromMxArray(biasField);
  sigma2 = mxArrayExtractDoubleField(matlabArray, "sigma2");

}




CProbitNoise::~CProbitNoise()
{
}
void CProbitNoise::setInitParam()
{
  setType("probit");
  setNoiseName("Probit");
  setNumParams(getNumProcesses());
  mu.resize(y.getRows(), y.getCols());
  varSigma.resize(y.getRows(), y.getCols());
  mu.zeros();
  varSigma.zeros();
  double nClass1=0.0;
  bias.resize(1, y.getCols());
  for(int j=0; j<y.getCols(); j++)
    {
      for(int i=0; i<y.getRows(); i++)
	{
	  if(y.getVal(i, j)==1)
	    nClass1++;
	}
      bias.setVal(ndlutil::invCumGaussian(nClass1/(double)y.getRows()), j);
      setParamName("bias" + j, j);
    }
      
  // sigma2 isn't treated as a parameter.
  sigma2 = 1e-6;
  clearTransforms();
  setLogConcave(true);
  setSpherical(false);
  setMissing(false);
}
ostream& CProbitNoise::display(ostream& os)
{
  double b = 0.0;
  os << "Probit noise: " << endl;
  for(int j=0; j<bias.getCols(); j++)
    {
      b = bias.getVal(j);
      os << "Bias on process " << j << ": " << b << endl; 
    }
  return os;
}
void CProbitNoise::setParam(const double val, const int index)
{
  assert(index>=0);
  assert(index<getNumParams());
  bias.setVal(val, index);
}  
void CProbitNoise::setParams(const CMatrix& params)
{
  assert(getNumParams()==getNumProcesses());
  assert(params.getCols()==getNumParams());
  assert(params.getRows()==1);
  assert(getNumProcesses()==bias.getCols());
  for(int j=0; j<bias.getCols(); j++)
    {
      bias.setVal(params.getVal(j), j);
    }
}
double CProbitNoise::getParam(const int index) const
{
  assert(index>=0);
  assert(index<getNumParams());
  if(index<getNumProcesses())
    return bias.getVal(index);

}
void CProbitNoise::getParams(CMatrix& params) const
{
  assert(getNumParams()==getNumProcesses());
  assert(params.getCols()==getNumParams());
  assert(params.getRows()==1);
  assert(getNumProcesses()==bias.getCols());
  for(int j=0; j<getNumProcesses(); j++)
    params.setVal(bias.getVal(j), j);
}
 
void CProbitNoise::getGradParams(CMatrix& g) const
{
  assert(g.getCols()==getNumParams());
  assert(g.getRows()==1);
  double c=0.0;
  double u=0.0;
  double b=0.0;
  for(int j=0; j<y.getCols(); j++)
    {
      double gbias = 0.0;
      b = bias.getVal(j);
      for(int i=0; i<y.getRows(); i++)
	{
	  c = getTarget(i, j)/sqrt(getVarSigma(i, j)+sigma2);
	  u=c*(getMu(i, j)+b);
	  u=ndlutil::gradLnCumGaussian(u);
	  gbias+=u*c;
	}
      g.setVal(gbias, 0, j);
    }
}
void CProbitNoise::getGradInputs(double& gmu, double& gvs, const int i, const int j) const
{
  double b = bias.getVal(j);
  double c = getTarget(i, j)/sqrt(sigma2+getVarSigma(i, j));
  double u = c*(getMu(i, j) + b);
  gmu = ndlutil::gradLnCumGaussian(u)*c;
  gvs = -0.5*c*u*gmu;
}
  
void CProbitNoise::out(CMatrix& yTest, const CMatrix& muTest, const CMatrix& varSigmaTest) const
{
  assert(yTest.dimensionsMatch(muTest));
  assert(muTest.dimensionsMatch(varSigmaTest));
  for(int j=0; j<yTest.getCols(); j++)
    {
      double b = bias.getVal(j);
      for(int i=0; i<yTest.getRows(); i++)	
	{
	  if(muTest.getVal(i, j)>-b)
	    yTest.setVal(1.0, i, j);
	  else
	    yTest.setVal(0.0, i, j);
	}
    }
}
void CProbitNoise::likelihood(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, 
				   const CMatrix& yTest) const
{
  assert(yTest.getCols()==getNumProcesses());
  assert(L.dimensionsMatch(muTest));
  assert(yTest.dimensionsMatch(muTest));
  assert(muTest.dimensionsMatch(varSigmaTest));
  double arg=0.0;
  double var=0.0;
  for(int i=0; i<muTest.getRows(); i++)
    {
      for(int j=0; j<muTest.getCols(); j++)
	{
	  arg = yTest.getVal(i, j) - muTest.getVal(i, j) - bias.getVal(j);
	  arg *= arg;
	  var = varSigmaTest.getVal(i, j) + sigma2;
	  arg = 1/sqrt(2*M_PI*var)*exp(-.5*arg*arg/var);
	}
    }
}

double CProbitNoise::logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, 
				     const CMatrix& yTest) const
{
  assert(yTest.getCols()==getNumProcesses());
  assert(yTest.dimensionsMatch(muTest));
  assert(yTest.dimensionsMatch(varSigmaTest));
  double arg=0.0;
  double var=0.0;
  double L=0.0;
  for(int i=0; i<muTest.getRows(); i++)
    {
      for(int j=0; j<muTest.getCols(); j++)
	{
	  arg = muTest.getVal(i, j) + bias.getVal(j);
	  arg *= getTarget(i, j);
	  var = sqrt(varSigmaTest.getVal(i, j) + sigma2);
	  L += ndlutil::lnCumGaussian(arg/var);
	}
    }  
  return L;
}

void CProbitNoise::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "nParams");
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)getNumParams()));
  mxAddField(matlabArray, "bias");
  mxSetField(matlabArray, 0, "bias", bias.toMxArray());
  mxAddField(matlabArray, "sigma2");
  mxSetField(matlabArray, 0, "sigma2", convertMxArray(sigma2));
}

void CProbitNoise::extractParamFromMxArray(const mxArray* matlabArray) 
{
  setNumParams(mxArrayExtractIntField(matlabArray, "nParams"));
  mxArray* biasField = mxArrayExtractMxArrayField(matlabArray, "bias");
  bias.fromMxArray(biasField);
  sigma2 = mxArrayExtractDoubleField(matlabArray, "sigma2");

}

CNcnmNoise::~CNcnmNoise()
{
}
void CNcnmNoise::setInitParam()
{
  setType("ncnm");
  setNoiseName("Null Category");
  setNumParams(getNumProcesses()+2);
  mu.resize(y.getRows(), y.getCols());
  varSigma.resize(y.getRows(), y.getCols());
  mu.zeros();
  varSigma.zeros();
  double nClass1=0.0;
  double nClass2=0.0;
  double nMissing=0.0;
  bias.resize(1, y.getCols());
  for(int j=0; j<y.getCols(); j++)
    {
      for(int i=0; i<y.getRows(); i++)
	{
	  if(y.getVal(i, j)==1.0)	    
	    nClass1++;
	  else if(y.getVal(i, j)==-1.0)
	    nClass2++;
	  else
	    nMissing++;
	}
      bias.setVal(ndlutil::invCumGaussian(nClass1/(nClass1+nClass2)), j);
      setParamName("bias" + j, j);
    }
  gamman=nMissing/(double)y.getRows();
  gammap = gamman;

  sigma2 = ndlutil::EPS;
  width = 1.0;
  clearTransforms();

  // sigmoid transforms on gamman and gammap.
  addTransform(new CSigmoidTransform(), getNumParams()-1);
  addTransform(new CSigmoidTransform(), getNumParams()-2);

  // sigma2 isn't treated as a parameter.
  setLogConcave(false);
  setSpherical(false);
  setMissing(true);
}
ostream& CNcnmNoise::display(ostream& os)
{
  double b = 0.0;
  os << "Ncnm noise: " << endl;
  for(int j=0; j<bias.getCols(); j++)
    {
      b = bias.getVal(j);
      os << "Bias on process " << j << ": " << b << endl; 
    }
  os << "Missing label probability for -ve class: " << gamman << endl;
  os << "Missing label probability for +ve class: " << gammap << endl;
  return os;
}
void CNcnmNoise::setParam(const double val, const int index)
{
  assert(index>=0);
  assert(index<getNumParams());
  if(index<getNumProcesses())
    {
      bias.setVal(val, index);
      return;
    }
  if(index==getNumProcesses())
    {
      gamman=val;
      return;
    }
  if(index==getNumProcesses()+1)
    {
      gammap=val;
      return;
    }
}  
void CNcnmNoise::setParams(const CMatrix& params)
{
  assert(params.getCols()==getNumParams());
  assert(params.getRows()==1);
  int nProc = getNumProcesses();
  for(int j=0; j<nProc; j++)
    {
      bias.setVal(params.getVal(j), j);
    }
  gamman=params.getVal(nProc);
  gammap=params.getVal(nProc+1);
}
double CNcnmNoise::getParam(const int index) const
{
  assert(index>=0);
  assert(index<getNumParams());
  if(index<getNumProcesses())
    return bias.getVal(index);
  if(index==getNumProcesses())
    return gamman;
  if(index==getNumProcesses()+1)
    return gammap;

}
void CNcnmNoise::getParams(CMatrix& params) const
{
  assert(params.getCols()==getNumParams());
  assert(params.getRows()==1);
  int nProc=getNumProcesses();
  for(int j=0; j<nProc; j++)
    params.setVal(bias.getVal(j), j);
  params.setVal(gamman, nProc);
  params.setVal(gammap, nProc+1);
}
 
void CNcnmNoise::getGradParams(CMatrix& g) const
{
  assert(g.getCols()==getNumParams());
  assert(g.getRows()==1);
  double ggamman=0.0;
  double ggammap=0.0;
  double halfWidth = width/2.0;
  for(int j=0; j<y.getCols(); j++)
    {
      double gbias = 0.0;
      double b = bias.getVal(j);
      for(int i=0; i<y.getRows(); i++)
	{
	  double muAdj = getMu(i, j)+b;
	  double c = 1/sqrt(sigma2+getVarSigma(i, j));
	  double targVal = getTarget(i, j);
	  if(targVal==-1.0)
	    {
	      muAdj+=halfWidth;
	      muAdj=muAdj*c;
	      gbias-= c*ndlutil::gradLnCumGaussian(-muAdj);
	      ggamman-=1.0/(1.0-gamman);
	    }
	  else if(targVal==1.0)
	    {
	      muAdj-=halfWidth;
	      muAdj*=c;
	      gbias+=c*ndlutil::gradLnCumGaussian(muAdj);
	      ggammap-=1.0/(1.0-gammap);
	    }
	  else
	    {
	      muAdj+=halfWidth;
	      double u=muAdj*c;
	      double uprime=(muAdj-width)*c;
	      double lndenom = ndlutil::lnCumGaussSum(-u, uprime, gamman, gammap);
	      double lnNumer1 = log(gamman) - ndlutil::HALFLOGTWOPI -.5*u*u;
	      double lnNumer2 = log(gammap) - ndlutil::HALFLOGTWOPI -.5*uprime*uprime;
	      double B1 = exp(lnNumer1-lndenom);
	      double B2 = exp(lnNumer2-lndenom);
	      gbias+=c*(B2-B1);
	      ggammap+=exp(ndlutil::lnCumGaussian(uprime)-lndenom);
	      ggamman+=exp(ndlutil::lnCumGaussian(-u)-lndenom);
	    }
	      
	}
      g.setVal(gbias, 0, j);
    }
  g.setVal(ggamman, 0, getNumProcesses());
  g.setVal(ggammap, 0, getNumProcesses()+1);
}
void CNcnmNoise::getGradInputs(double& gmu, double& gvs, const int i, const int j) const
{
  double b = bias.getVal(j);
  double c = 1.0/sqrt(sigma2+getVarSigma(i, j));
  double muAdj = getMu(i, j)+b;
  double targ = getTarget(i, j);
  double halfWidth = width/2.0;
  if(targ==-1.0)
    {
      muAdj+=halfWidth;
      muAdj*=c;
      gmu=-ndlutil::gradLnCumGaussian(-muAdj)*c;
      gvs=-.5*c*muAdj*gmu;
    }
  else if(targ==1.0)
    {
      muAdj-=halfWidth;
      muAdj*=c;
      gmu=ndlutil::gradLnCumGaussian(muAdj)*c;
      gvs=-.5*c*muAdj*gmu;
    }
  else // missing data
    {
      muAdj+=halfWidth;
      double u=c*muAdj;
      double uprime=(muAdj-width)*c;
      double lndenom=ndlutil::lnCumGaussSum(-u, uprime, gamman, gammap);
      double lnNumer1 = log(gamman) - ndlutil::HALFLOGTWOPI -.5*(u*u);
      double lnNumer2 = log(gammap) - ndlutil::HALFLOGTWOPI -.5*(uprime*uprime);
      double B1 = exp(lnNumer1 - lndenom);
      double B2 = exp(lnNumer2 - lndenom);
      gmu = c*(B2-B1);
      gvs = -.5*c*c*(uprime*B2-u*B1);
      
    }
  if(isnan(gmu))
    throw "gmu is NaN";
  if(isnan(gvs))
    throw "gvs is NaN";
}
  
void CNcnmNoise::out(CMatrix& yTest, const CMatrix& muTest, const CMatrix& varSigmaTest) const
{
  assert(yTest.dimensionsMatch(muTest));
  assert(muTest.dimensionsMatch(varSigmaTest));
  for(int j=0; j<yTest.getCols(); j++)
    {
      double b = bias.getVal(j);
      for(int i=0; i<yTest.getRows(); i++)	
	{
	  double muVal = muTest.getVal(i, j);
	  if(muVal>-b)
	    yTest.setVal(1.0, i, j);
	  else
	    yTest.setVal(0.0, i, j);
	}
    }
}
void CNcnmNoise::likelihood(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, 
				   const CMatrix& yTest) const
{
  assert(yTest.getCols()==getNumProcesses());
  assert(L.dimensionsMatch(muTest));
  assert(yTest.dimensionsMatch(muTest));
  assert(muTest.dimensionsMatch(varSigmaTest));
  double halfWidth = width/2.0;
  for(int j=0; j<muTest.getCols(); j++)
    {
      double b = bias.getVal(j);
      for(int i=0; i<muTest.getRows(); i++)
	{
	  double muAdj=muTest.getVal(i, j) + b;
	  double c=1/sqrt(sigma2+varSigmaTest.getVal(i, j));	  
	  double targVal=yTest.getVal(i, j);
	  if(targVal==1.0)
	    {
	      muAdj-=halfWidth;
	      L.setVal(ndlutil::cumGaussian(muAdj*c)*(1-gammap), i, j);
	    }
	  else if(targVal==-1.0)
	    {
	      muAdj+=halfWidth;
	      L.setVal(ndlutil::cumGaussian(-muAdj*c)*(1-gamman), i, j);
	    }
	  else // missing data
	    {
	      muAdj+=halfWidth;
	      L.setVal(gamman*ndlutil::cumGaussian(-muAdj*c)+gammap*(ndlutil::cumGaussian((muAdj-width)*c)), i, j);
	    }
	}
    }
}

double CNcnmNoise::logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, 
				     const CMatrix& yTest) const
{
  assert(yTest.getCols()==getNumProcesses());
  assert(yTest.dimensionsMatch(muTest));
  assert(yTest.dimensionsMatch(varSigmaTest));
  double arg=0.0;
  double var=0.0;
  double L=0.0;
  double halfWidth = width/2.0;
  double logPosGamma = log(1.0-gammap);
  double logNegGamma = log(1.0-gamman);
  for(int j=0; j<muTest.getCols(); j++)    
    {
      double b=bias.getVal(j);
      for(int i=0; i<muTest.getRows(); i++)
	{
	  double muAdj = muTest.getVal(i, j) + b;
	  double c=1/sqrt(sigma2+varSigmaTest.getVal(i, j));
	  double targVal=yTest.getVal(i, j);
	  if(targVal==1.0)
	    {
	      muAdj-=halfWidth;
	      L+=ndlutil::lnCumGaussian(muAdj*c);
	      L+=logPosGamma;
	      
	    }
	  else if(targVal==-1.0)
	    {
	      muAdj+=halfWidth;
	      L+=ndlutil::lnCumGaussian(-muAdj*c);
	      L+=logNegGamma;
	    }
	  else // missing data.
	    {
	      muAdj+=halfWidth;
	      double u=muAdj*c;
	      double uprime=(muAdj-width)*c;
	      L+=ndlutil::lnCumGaussSum(-u, uprime, gamman, gammap);	      
	    }
	}
    }  
  return L;
}

void CNcnmNoise::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "nParams");
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)getNumParams()));
  mxAddField(matlabArray, "bias");
  mxSetField(matlabArray, 0, "bias", bias.toMxArray());
  mxAddField(matlabArray, "sigma2");
  mxSetField(matlabArray, 0, "sigma2", convertMxArray(sigma2));
  mxAddField(matlabArray, "width");
  mxSetField(matlabArray, 0, "width", convertMxArray(width));
  mxAddField(matlabArray, "gamman");
  mxSetField(matlabArray, 0, "gamman", convertMxArray(gamman));
  mxAddField(matlabArray, "gammap");
  mxSetField(matlabArray, 0, "gammap", convertMxArray(gammap));
}

void CNcnmNoise::extractParamFromMxArray(const mxArray* matlabArray) 
{
  setNumParams(mxArrayExtractIntField(matlabArray, "nParams"));
  mxArray* biasField = mxArrayExtractMxArrayField(matlabArray, "bias");
  bias.fromMxArray(biasField);
  sigma2 = mxArrayExtractDoubleField(matlabArray, "sigma2");
  width = mxArrayExtractDoubleField(matlabArray, "width");
  gamman = mxArrayExtractDoubleField(matlabArray, "gamman");
  gammap = mxArrayExtractDoubleField(matlabArray, "gammap");

}
