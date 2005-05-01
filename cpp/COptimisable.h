#ifndef COPTIMISABLE_H
#define COPTIMISABLE_H
#include "CMatrix.h"
#include "ndlutil.h"
// abstract class for making a class optimisable.
class COptimisable {
  
 public:
  COptimisable()
    {
    }
  virtual int getVerbosity() const=0;
  virtual int getOptNumParams() const=0;
  virtual void getOptParams(CMatrix& param) const=0;
  virtual void setOptParams(const CMatrix& param)=0;
  virtual void computeObjectiveGradParams(CMatrix& g) const=0;
  virtual double computeObjectiveVal() const=0;

  inline void setDirection(const CMatrix& vals)
    {
      assert(vals.getCols()==getOptNumParams());
      assert(vals.getRows()==1);
      direction.deepCopy(vals);
    }
  inline void getDirection(CMatrix& vals) const
    {
      assert(vals.getCols()==getOptNumParams());
      assert(vals.getRows()==1);
      assert(direction.dimensionsMatch(vals));
      vals.deepCopy(direction);
    }
  void checkGradients();  
  void gdOptimise(double learnRate=0.01, double momentum=0.9, int display=0, int maxIters=1000, const double objectiveTol=1e-6, const double paramTol=1e-6);
  void gdPullbackOptimise(double learnRate=0.01, int display=0, int maxIters=1000, const double objectiveTol=1e-6, const double paramTol=1e-6);
  void netlabScgOptimise(int maxIters=100, const double objectiveTol=1e-6, const double paramTol=1e-6);
  void scgOptimise(int maxIters=100, const double objectiveTol=1e-6, const double paramTol=1e-6);
  void lineMinimisation(const CMatrix& direction);
  double oneDObjectiveVal(double val);
  void lineMinimisation();
  void bracketMinimum(double& a, double& b, double& c, double& fa, const int maxStep);
  
 private:
  
  double objectiveTol;
  double parameterTol;
  //int maxIters;
  //double learnRate;
  //double momentum;

  const static bool evalFunc=true;
  //int display;
  const static double phi=1.618033988749895;
  const static double cphi=-0.6180339887498949;
  const static double smallNum=1e-11;

  CMatrix direction; // direction for 1-D optimisation.
  CMatrix paramStoreOne;
  CMatrix paramStoreTwo;
}; 

#endif
