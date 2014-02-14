function [params, names] = expvarMeanExtractParam(model)

% EXPVARMEANEXTRACTPARAM Extract weights and biases from an EXPVARMEAN mapping.
% FORMAT
% DESC returns a vector of all the parameters from the mapping
% associated with the mean function of a exponential variational approximation.
% ARG model : the mappings from which we wish to extract the parameters.
% and biases.
% RETURN params : vector of all the parameters.
% RETURN names : optional additional returned cell array of the
% names of the parameters.
%
% SEEALSO : expvarMeanCreate, expvarMeanExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

if nargout > 1
  [params, names] = kernExtractParam(model.kern);
else
  params = kernExtractParam(model.kern);
end
