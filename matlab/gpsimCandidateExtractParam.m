function [param, names] = gpsimCandidateExtractParam(model)

% GPSIMCANDIDATEEXTRACTPARAM Extract the parameters of a GPSIM model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process for single input motif
% modelling.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpsimCreate, gpsimAddCandidate, gpsimCandidateExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2007

% GPSIM

if nargout>1
  [param, names] = gpsimExtractParam(model.candidate);
else
  param = gpsimExtractParam(model.candidate);
end
