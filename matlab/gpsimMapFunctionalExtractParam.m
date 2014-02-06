function f = gpsimMapFunctionalExtractParam(model)

% GPSIMMAPFUNCTIONALEXTRACTPARAM Extract the function values from a GPSIMMAP model.
% FORMAT
% DESC extracts the function values from a structure containing
% the information about a Gaussian process for single input motif
% modelling.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a row vector of the function parameters from the model.
%
% SEEALSO : gpsimMapCreate, gpsimMapFunctionalExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

f = model.f';
