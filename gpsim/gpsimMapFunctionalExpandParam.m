function model = gpsimMapFunctionalExpandParam(model, f)

% GPSIMMAPFUNCTIONALEXPANDPARAM Expand the given function values into a GPSIMMAP structure.
% FORMAT
% DESC takes the given vector of function values and places them in
% the model structure. It then updates any stored representations that
% are dependent on those function values, for example the Hessian, etc..
% ARG model : the model structure for which parameters are to be
% updated.
% ARG vals : a vector of function values for placing in the model
% structure.
% RETURN model : a returned model structure containing the updated
% function values.
% 
% SEEALSO : gpsimMapCreate, gpsimMapFunctionalExtractParam, modelExtractParam, gpsimMapUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

model.f = f';
model = gpsimMapUpdateG(model);
model = gpsimMapUpdateYpred(model);
model = gpsimMapFunctionalUpdateW(model);
model = gpsimMapUpdatePosteriorCovariance(model);
