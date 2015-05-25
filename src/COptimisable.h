#ifndef COPTIMISABLE_H
#define COPTIMISABLE_H
#include "CMatrix.h"
#include "ndlutil.h"
#include "ndlfortran.h"

// abstract base class for making a class optimisable.
class COptimisable {
  
 public:
  enum 
  {
    CG, 
    SCG, 
    GD, 
    BFGS,
    LBFGS
  };
  COptimisable()
  {
    // set default optimisation parameters
    setVerbosity(2);
    setDefaultOptimiser(SCG);
    setFuncEvalTerminate(false);
    setIterTerminate(true);
    setMaxFuncEvals(1000);
    setMaxIters(1000);
    setObjectiveTol(1e-6);
    setParamTol(1e-6);
    setLearnRate(0.01);
    setMomentum(0.9);
    iter = 0;
    funcEval = 0;
  }
  virtual ~COptimisable() {}
  virtual inline void setVerbosity(int val) const 
  {
    verbosity = val;
  }  
  virtual inline int getVerbosity() const 
  {
    return verbosity;
  }
  virtual unsigned int getOptNumParams() const=0;
  virtual void getOptParams(CMatrix& param) const=0;
  virtual void setOptParams(const CMatrix& param)=0;
  virtual double computeObjectiveGradParams(CMatrix& g) const=0;
  virtual double computeObjectiveVal() const=0;

  inline void setDirection(const CMatrix& vals)
  {
    DIMENSIONMATCH(vals.getCols()==getOptNumParams());
    DIMENSIONMATCH(vals.getRows()==1);
    direction.deepCopy(vals);
  }
  inline void getDirection(CMatrix& vals) const
  {
    DIMENSIONMATCH(vals.getCols()==getOptNumParams());
    DIMENSIONMATCH(vals.getRows()==1);
    DIMENSIONMATCH(direction.dimensionsMatch(vals));
    vals.deepCopy(direction);
  }
  void checkGradients();  
  void gdOptimise();
  void gdPullbackOptimise();
  //void netlabScgOptimise();
  void lbfgsOptimise();
  void scgOptimise();
  void cgOptimise();
  //void lineMinimisation(const CMatrix& direction);
  double oneDObjectiveVal(double val);
  //void lineMinimisation();
  //void bracketMinimum(double& a, double& b, double& c, double& fa, unsigned int maxStep);
  void setLearnRate(double val)
  {
    learnRate = val;
  }
  double getLearnRate() const
  {
    return learnRate;
  }
  void setMomentum(double val)
  {
    momentum = val;
  }
  double getMomentum() const
  { 
    return momentum;
  }
  void setFuncEvalTerminate(bool val)
  {
    funcEvalTerminate = val;
  }
  bool isFuncEvalTerminate() const
  {
    return funcEvalTerminate;
  }
  void setIterTerminate(bool val)
  {
    iterTerminate = val;
  }
  bool isIterTerminate() const
  {
    return iterTerminate;
  }
  void setMaxFuncEvals(unsigned int val)
  {
    maxFuncVal = val;
  }
  unsigned int getMaxFuncEvals() const
  {
    return maxFuncVal;
  }
  void setMaxIters(unsigned int val)
  {
    maxIters = val;
  }
  unsigned int getMaxIters() const
  {
    return maxIters;
  }
  void setObjectiveTol(double val)
  {
    objectiveTol = val;
  }
  double getObjectiveTol() const
  {
    return objectiveTol;
  }
  void setParamTol(double val)
  {
    parameterTol = val;
  }
  double getParamTol() const
  {
    return parameterTol;
  }
  void setDefaultOptimiser(int val) const
  {
    defaultOptimiser = val;
  }
  int getDefaultOptimiser() const
  {
    return defaultOptimiser;
  }
  void setDefaultOptimiserStr(string val)
  {
    if(val == "conjgrad")
      defaultOptimiser = CG;
    else if(val == "scg")
      defaultOptimiser = SCG;
    else if(val == "quasinew")
      defaultOptimiser = BFGS;
    else if(val == "graddesc")
      defaultOptimiser = GD;
    else 
      throw ndlexceptions::NotImplementedError("Unknown optimisation");
  }
  string getDefaultOptimiserStr() const
  {
    switch(defaultOptimiser)
    {
    case CG:
      return "conjgrad";
    case SCG:
      return "scg";
    case BFGS:
      return "quasinew";
    case GD:
      return "graddesc";
    default:
      throw ndlexceptions::NotImplementedError("Unknown optimisation.");

    }
  }
  void runDefaultOptimiser()
  {
    switch(defaultOptimiser)
    {
    case CG:
      cgOptimise();
      break;
    case SCG:
      scgOptimise();
      break;
    case BFGS:
      lbfgsOptimise();
      break;
    case GD:
      gdOptimise();
      break;
    default:
      throw ndlexceptions::NotImplementedError("Unknown optimisation.");
      
    }
  }
 private:
  
  double objectiveTol;
  double parameterTol;
  
  #ifndef SWIG
    const static bool evalFunc=true;
    const static double phi;
    const static double cphi;
    const static double smallNum;
  #else
    static bool evalFunc=true;
    static double phi;
    static double cphi;
    static double smallNum;
  #endif
  double learnRate;
  double momentum;

  unsigned int iter;
  unsigned int funcEval;

  unsigned int maxIters;
  unsigned int maxFuncVal;

  mutable int verbosity; 

  mutable int defaultOptimiser;

  bool funcEvalTerminate;
  bool iterTerminate;
  
  CMatrix direction; // direction for 1-D optimisation.
  CMatrix paramStoreOne;
  CMatrix paramStoreTwo;
}; 

class CProbabilisticOptimisable : public COptimisable
{
 public:
  CProbabilisticOptimisable() : COptimisable() {}
  virtual double logLikelihood() const=0;
  virtual double logLikelihoodGradient(CMatrix& g) const=0;
  virtual double computeObjectiveGradParams(CMatrix& g) const 
  {
    double L = logLikelihoodGradient(g);
    g.negate();
    return -L;
  }
  virtual double computeObjectiveVal() const 
  {
    return -logLikelihood();
  }
};
#endif
