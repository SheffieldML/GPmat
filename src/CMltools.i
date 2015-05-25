// CMltools.i -- SWIG interface
%module ndlmltools
%{
#define NDLSWIG
#include "CMltools.h"
%}

%rename(ndlmatrix) CMatrix;
%rename(mlpMapping) CMlpMapping;
%rename(linearMapping) CLinearMapping;

// parse the original header file
%import "CNdlInterfaces.h"
%import "CMatrix.h"
%import "CTransform.h"
%import "CDataModel.h"
%import "CDist.h"
%import "CKern.h"
%import "CNoise.h"
%import "COptimisable.h"
%include "CMltools.h"
