#ifndef CIVM_H
#define CIVM_H
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "COptimisable.h"
#include "CMatlab.h"
#include "CMatrix.h"
#include "CKern.h"
#include "CNoise.h"
using namespace std;

const double NULOW=1e-16;

class CIvm : public CMatinterface, public COptimisable {
 public:
  CIvm(const CMatrix& inData, const CMatrix& targetData, 
       CKern& kernel, CNoise& noiseModel, const int selectCrit,
       const int dVal, const int verbos=2);
  void init();
  void selectPoints(); // select active set points.
  void addPoint(const int index); // add a point to the model.
  void updateSite(const int index); // update the site parameters at index.
  void updateM(const int index); // update M at index.

  int selectPointAdd(); // select a point to add to active set.
  int entropyPointAdd(); // add a point selected by entropy change.
  int randomPointAdd();  // add a point selected randomly.
  double entropyChangeAdd(const int) const; // entropy change associated with adding a point

  int selectPointRemove();  // select a point to remove from the active set.
  int entropyPointRemove(); // remove a point selected by entropy change.
  int randomPointRemove(); // remove a point selected randomly.
  double entropyChangeRemove(const int) const; // entropy change associated with removing a point

  void out(CMatrix& yPred, const CMatrix& inData) const;
  void posteriorMeanVar(CMatrix& mu, CMatrix& varSigma, const CMatrix& X) const;
  inline void setVerbosity(const int val)
    {
      verbosity = val;
    }
  inline int getVerbosity() const
    {
      return verbosity;
    }
  inline int changeEntropy(const double val)
    {
      cumEntropy += val;
      lastEntropyChange = val;
    }

  // Gradient routines
  void updateCovGradient(int index) const;
  

  inline void setTerminate(const bool val)
    {
      terminate = val;
    }
  inline void setEpUpdate(const bool val)
    {
      epUpdate = val;
    }
  inline bool isTerminate()
    {
      return terminate;
    }
  inline bool isEpUpdate()
    {
      return epUpdate;
    }
  void updateNuG();
  // update K with the kernel computed from the active points.
  void updateK() const;
  // update invK with the inverse of the kernel plus beta terms computed from the active points.
  void updateInvK(int index=0) const;
  // compute the approximation to the log likelihood.
  double approxLogLikelihood() const;
  // compute the gradients of the approximation wrt parameters.
  void approxLogLikelihoodGradient(CMatrix& g) const;
  
  void optimise(const int maxIters=15, const int kernIters=100, const int noiseIters=100);
  void display(ostream& os) const;

  inline int getOptNumParams() const
    {
      return kern.getNumParams();
    }    
  void getOptParams(CMatrix& param) const
    {
      kern.getTransParams(param);
    }
  void setOptParams(const CMatrix& param)
    {
      kern.setTransParams(param);
    }
  void computeObjectiveGradParams(CMatrix& g) const
    {
      approxLogLikelihoodGradient(g);
      g.negate();
    }
  double computeObjectiveVal() const
    {
      return -approxLogLikelihood();
    }
  mxArray* toMxArray() const
    {
      // TODO Write this function
    }
  void fromMxArray(const mxArray* matlabArray) 
    {
      // TODO Write this function
    }
  const CMatrix& X;

  // arguably these are noise model associated.
  const CMatrix& y;
  
  CMatrix nu;
  CMatrix g;
  
  CMatrix Kstore;

  
  // these are IVM associated.
  CMatrix m;
  CMatrix beta;

  CMatrix s;
  CMatrix a;
  CMatrix ainv;

  // these really just provide local storage
  mutable CMatrix covGrad;
  mutable CMatrix invK;
  mutable double logDetK;
  mutable CMatrix K;

  CMatrix activeX;

  CMatrix* M;
  CMatrix* L;
  CMatrix* Linv;

  //  COptions options;
  
  vector<int> activeSet;
  vector<int> inactiveSet;

  CKern& kern;
  CNoise& noise;
  

 private:
  bool terminate;
  bool epUpdate;

  int numCovStruct;
  int activeSetSize;

  int numTarget;
  int numData;
  
  double lastEntropyChange;
  double cumEntropy;
  enum{ENTROPY, RENTROPY, RANDOM};
  const int selectionCriterion; // need to set this up with enum
  
  int verbosity;
};

#endif
