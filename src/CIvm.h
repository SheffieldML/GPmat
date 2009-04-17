#ifndef CIVM_H
#define CIVM_H
#include "CMltools.h"
using namespace std;

const double NULOW=1e-16;
const string IVMVERSION="0.2";

class CIvm : public CMapModel, public CProbabilisticOptimisable, public CStreamInterface, public CMatInterface
{
 public:
  CIvm();
  // Constructor given a filename.
  // CIvm(const string modelFileName, int verbos=2);
  // Constructor given a kernel and a noise model.
  CIvm(CMatrix* inData, CMatrix* targetData, 
       CKern* kernel, CNoise* noiseModel, int selectCrit,
       unsigned int dVal, int verbos=2);
  CIvm(CMatrix& actX, CMatrix& actY, 
       CMatrix& mmat, CMatrix& betamat, 
       vector<unsigned int> actSet, CKern* kernel, 
       CNoise* noiseModel, int selectCrit=ENTROPY, 
       int verbos=2);

#ifdef _NDLMATLAB
  // Constructor using file containing ivmInfo.
  CIvm(CMatrix* inData, 
       CMatrix* targetData, 
       CKern* kernel, 
       CNoise* noiseModel, 
       const string ivmInfoFile, 
       const string ivmInfoVariable, 
       int verbos=2);
#endif


  void writeParamsToStream(ostream& os) const;
  void readParamsFromStream(istream& is);
  // initialise the model.
  void init();
  // Initialise the storeage for the model.
  void initStoreage();
  // Set the initial values for the model.
  void initVals();
  void selectPoints(); // select active set points.
  void addPoint(unsigned int index); // add a point to the model.
  void updateSite(unsigned int index); // update the site parameters at index.
  void updateM(unsigned int index); // update M at index.

  unsigned int selectPointAdd(); // select a point to add to active set.
  unsigned int entropyPointAdd(); // add a point selected by entropy change.
  unsigned int randomPointAdd();  // add a point selected randomly.
  double entropyChangeAdd(unsigned int) const; // entropy change associated with adding a point

  unsigned int selectPointRemove();  // select a point to remove from the active set.
  unsigned int entropyPointRemove(); // remove a point selected by entropy change.
  unsigned int randomPointRemove(); // remove a point selected randomly.
  double entropyChangeRemove(unsigned int) const; // entropy change associated with removing a point
  void test(const CMatrix& ytest, const CMatrix& Xin) const;

  void likelihoods(CMatrix& pout, CMatrix& yTest, const CMatrix& Xin) const;
  double logLikelihood(const CMatrix& yTest, const CMatrix& Xin) const;

  // For MapModel interface.
  void out(CMatrix& yPred, const CMatrix& inData) const;
  void out(CMatrix& yPred, CMatrix& probPred, const CMatrix& inData) const;
  double outGradParams(CMatrix& g, const CMatrix& Xin, unsigned int pointNo, unsigned int outputNo) const;
  double outGradX(CMatrix& g, const CMatrix& Xin, unsigned int pointNo, unsigned int outputNo) const;

  void posteriorMeanVar(CMatrix& mu, CMatrix& varSigma, const CMatrix& X) const;
  string getNoiseName() const
  {
    return pnoise->getName();
  }
  inline void changeEntropy(double val)
  {
    cumEntropy += val;
    lastEntropyChange = val;
  }

  // Gradient routines
  void updateCovGradient(unsigned int index) const;
  

  inline void setTerminate(const bool val)
  {
    terminate = val;
  }
  inline void setEpUpdate(const bool val)
  {
    epUpdate = val;
  }
  inline bool isTerminate() const
  {
    return terminate;
  }
  inline bool isEpUpdate() const 
  {
    return epUpdate;
  }
  void updateNuG();
  // update K with the kernel computed from the active points.
  void updateK() const;
  // update invK with the inverse of the kernel plus beta terms computed from the active points.
  void updateInvK(unsigned int index=0) const;
  // compute the approximation to the log likelihood.
  double logLikelihood() const;
  // compute the gradients of the approximation wrt parameters.
  double logLikelihoodGradient(CMatrix& g) const;
  
  void optimise(unsigned int maxIters=15, unsigned int kernIters=100, unsigned int noiseIters=100);
  bool equals(const CIvm& model, double tol=ndlutil::MATCHTOL) const;
  void display(ostream& os) const;

 
  inline unsigned int getOptNumParams() const
  {
    return pkern->getNumParams();
  }    
  void getOptParams(CMatrix& param) const
  {
    pkern->getTransParams(param);
  }
  void setOptParams(const CMatrix& param)
  {
    pkern->setTransParams(param);
  }
  
  string getTypeSelection() const
  {
    switch(selectionCriterion)
    {
    case ENTROPY:
      return "entropy";
    case RENTROPY:
      return "rentropy";
    case RANDOM:
      return "random";
    default:
      cerr << "Unrecognised selection criterion." << endl;
    }
  }
  void setTypeSelection(const string val) 
  {
    if(val=="entropy")
      selectionCriterion=ENTROPY;
    else if(val=="rentropy")
      selectionCriterion=RENTROPY;
    else if(val=="random")
      selectionCriterion=RANDOM;
    else
      cerr << "Unrecognised selection criterion " << val << "." << endl;
  }
  void setTypeSelection(unsigned int val)
  {
    BOUNDCHECK(val>=ENTROPY && val<=RANDOM);
    selectionCriterion=val;
  }
  
#ifdef _NDLMATLAB
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* matlabArray);
#endif

  unsigned int getActiveSetSize() const
  {
    return activeSetSize;
  }
  
  double getActiveX(unsigned int i, unsigned int j) const
  {
    return activeX.getVal(i, j);
  }
  unsigned int getActivePoint(unsigned int i) const
  {
    return activeSet[i];
  }
  // arguably these are noise model associated.
  CMatrix* pX;
  CMatrix* py;
  
  CMatrix nu;
  CMatrix g;
  
  CMatrix Kstore;

  
  // these are IVM associated.
  CMatrix m;
  CMatrix beta;

  // these really just provide local storage
  mutable CMatrix covGrad;
  mutable CMatrix invK;
  mutable double logDetK;
  mutable CMatrix K;

  mutable CMatrix s;
  mutable CMatrix a;
  mutable CMatrix ainv;


  CMatrix activeX;
  CMatrix activeY;

  CMatrix* M;
  CMatrix* L;
  CMatrix* Linv;

  //  COptions options;
  
  vector<unsigned int> activeSet;
  vector<unsigned int> inactiveSet;

  CKern* pkern;
  CNoise* pnoise;
  enum{ENTROPY, RENTROPY, RANDOM};


 private:

  void _init();

  bool terminate;
  bool epUpdate;
  bool loadedModel;

  unsigned int numCovStruct;
  unsigned int activeSetSize;

  unsigned int numTarget;
  unsigned int numData;
  
  double lastEntropyChange;
  double cumEntropy;
  int selectionCriterion; 
  
  string type;
};

// Functions which operate on the object
void writeIvmToStream(const CIvm& model, ostream& out);
void writeIvmToFile(const CIvm& model, const string modelFileName, const string comment="");
CIvm* readIvmFromStream(istream& in);
CIvm* readIvmFromFile(const string modelfileName, int verbosity=2);

#endif
