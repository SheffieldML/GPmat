// CKern.i -- SWIG interface
%module ndlkern
%{
#define NDLSWIG
#include "CKern.h"
%}

%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectorui) vector<unsigned int>;
   %template(vectors) vector<string>;
   %template(vectord) vector<double>;
   %template(vectorm) vector<CMatrix*>;
   %template(vectork) vector<CKern*>;
};


%rename(ndlmatrix) CMatrix;
%rename(ndlkern) CKern;
%rename(cmpndKern) CCmpndKern;
%rename(tensorKern) CTensorKern;
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
%import "ndlassert.h"
%import "CTransform.h"
%import "CDataModel.h"
%import "CMatrix.h"
%import "CDist.h"
%import "ndlstrutil.h"
%import "ndlutil.h"
%include "CKern.h"
