#ifndef CDIST_H
#define CDIST_H
#include <cmath>
#include <vector>
#include "CMatrix.h"
#include "CTransform.h"
#include "mex.h"


class CDist {

 public:
  CDist();
  ~CDist();
  //CDist(CDist& dist);
  void setInitParam();
  CMatrix getParams();
  void setParams(CMatrix& params);
  double logProb(CMatrix& data);
  CMatrix gradParam(CMatrix& data);
  void addTransform(CTransform* transform, const int index)
    {
      transforms.push_back(transform);
      transIndex.push_back(index);
    }
 protected:
  vector<CTransform*> transforms;
  vector<int> transIndex;
  int nParams;
};

class CGaussian : public CDist {

 public:
  CGaussian()
    {
      //CDist::CDist();
      setInitParam();
    }
  /*  CGaussian(CGaussian& gauss)
  {
      //CDist::CDist(gauss);
      precision = gauss.precision;
      }*/
  ~CGaussian()
    {
    }
  void setInitParam();
  void setParams(CMatrix& params);
  CMatrix getParams();
  double logProb(CMatrix& data);
  CMatrix CGaussian::gradParam(CMatrix& X);

 private:
  double precision;
  int nParams;
  vector<CTransform> transforms;
};

#endif

    
