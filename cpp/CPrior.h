#include <cmath>
#include "CMatrix.h"
#include "mex.h"


class CPrior : public CMatinterface {

 public:
  CPrior();
  ~CPrior();
  CPrior(CPrior& prior);
  virtual void setInitParam();
  
};    
