#ifndef NDLASSERT_H
#define NDLASSERT_H
#ifdef _NDLCASSERT
#include <cassert>
#define BOUNDCHECK(x) assert(x)
#define DIMENSIONMATCH(x) assert(x)
#define CHARARGUMENTS(x) assert(x)
#define MATRIXPROPERTIES(x) assert(x)
#define SANITYCHECK(x) assert(x)
#define CHECKZEROORPOSITIVE assert(x)
#else
#include "ndlexceptions.h"
#define BOUNDCHECK(x) ndlCheckBounds(x)
#define DIMENSIONMATCH(x) ndlCheckDimension(x)
#define CHARARGUMENTS(x) ndlCheckCharArguments(x)
#define MATRIXPROPERTIES(x) ndlCheckMatrixProperties(x)
#define SANITYCHECK(x) ndlCheckSanity(x)
#define CHECKZEROORPOSITIVE(x) ndlCheckZeroOrPositive(x)
using namespace std;
const double test = 1.0;
void ndlCheckBounds(bool cond);
void ndlCheckDimension(bool cond);
void ndlCheckCharArguments(bool cond);
void ndlCheckMatrixProperties(bool cond);
void ndlCheckSanity(bool cond);
void ndlCheckZeroOrPositive(bool cond);
#endif /* _NDLPYTHON */
#endif /* NDLASSERT_H */
