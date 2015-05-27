function [params, names] = gibbsKernExtractParam(kern)

% GIBBSKERNEXTRACTPARAM Extract parameters from the GIBBS kernel structure.
% FORMAT
% DESC extracts parameters from the Mark Gibbs's non-stationary
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the Mark Gibbs's non-stationary
% kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO gibbsKernParamInit, gibbsKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% KERN

params = zeros(1, kern.nParams);
if nargout < 2
  params(1:end-1) = modelExtractParam(kern.lengthScaleFunc);
  params(end) = kern.variance;
else
  names = cell(1, kern.nParams);
  [params(1:end-1), names(1:end-1)] = modelExtractParam(kern.lengthScaleFunc);
  params(end) = kern.variance;
  names{end} = 'variance';
end
