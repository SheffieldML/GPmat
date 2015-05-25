#include "ndlassert.h"

void ndlCheckBounds(bool cond)
{
  if(cond)
    return;
  else
  {
    throw std::out_of_range("Bound check failed.");
  }
}
void ndlCheckDimension(bool cond)
{
  if(cond)
    return;
  else
  {
    throw ndlexceptions::Error("Dimension match check failed.");
  }
}
void ndlCheckCharArguments(bool cond)
{
  if(cond)
    return;
  else
  {
    throw ndlexceptions::Error("Char argument check failed.");
  }
}
void ndlCheckMatrixProperties(bool cond)
{
  if(cond)
    return;
  else
  {
    throw ndlexceptions::Error("Matrix properties check failed");
  }
}
void ndlCheckSanity(bool cond)
{
  if(cond)
    return;
  else
  {
    throw ndlexceptions::Error("Sanity check failed.");
  }
}
void ndlCheckZeroOrPositive(bool cond)
{
  if(cond)
    return;
  else
  {
    throw ndlexceptions::Error("Zero or positive check failed.");
  }
}
