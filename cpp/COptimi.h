#ifndef COPTIMI_H
#define COPTIMI_H
#include <iostream>
#include <cmath>
#include <vector>
#include "CMatrix.h"

class COptimi {
   
  virtual void gradientParams(CMatrix& g)=0;
  virtual void gradientX(CMatrix& g)=0;
  virtual void getParams(CMatrix& params);
  virtual void setParams(const CMatrix& params)=0;
  virtual void 
};

#endif
