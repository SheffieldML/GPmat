// CGp.i -- SWIG interface
%module ndlgp
%{
#define NDLSWIG
#include "CGp.h"
%}

%rename(gp) CGp;
%rename(ndlmatrix) CMatrix;

%rename(ndlkern) CKern;
%rename(whiteKern) CWhiteKern;
%rename(biasKern) CBiasKern;
%rename(rbfKern) CRbfKern;
%rename(ratquadKern) CRatQuadKern;
%rename(matern32Kern) CMatern32Kern;
%rename(matern52Kern) CMatern52Kern;
%rename(linKern) CLinKern;
%rename(mlpKern) CMlpKern;
%rename(polyKern) CPolyKern;
%rename(linardKern) CLinardKern;
%rename(rbfardKern) CRbfardKern;
%rename(mlpardKern) CMlpardKern;
%rename(polyardKern) CPolyardKern;

// parse the original header file
%import "CNdlInterfaces.h"
%import "CTransform.h"
%import "CDataModel.h"
%import "CDist.h"
%import "CKern.h"
%import "CNoise.h"
%import "COptimisable.h"
%import "CMltools.h"
%include "CGp.h"
