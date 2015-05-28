function [params, names] = polyKernExtractParam(kern)

% POLYKERNEXTRACTPARAM Extract parameters from the POLY kernel structure.
% FORMAT
% DESC Extract parameters from the polynomial kernel structure into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and their names from the polynomial
% kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO polyKernParamInit, polyKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% KERN


params = [kern.weightVariance kern.biasVariance kern.variance];
if nargout > 1
  names = {'weight variance', 'bias variance', 'variance'};
end
%/~
if any(isinf(params))
  warning('params are infinite')
end
%~/
