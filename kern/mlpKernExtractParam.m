function [params, names] = mlpKernExtractParam(kern)


% MLPKERNEXTRACTPARAM Extract parameters from the MLP kernel structure.
% FORMAT
% DESC Extract parameters from the multi-layer perceptron kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and parameter names from the multi-layer
% perceptron kernel matrix into a vector of parameters for
% optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings giving names to the parameters.
%
% SEEALSO mlpKernParamInit, mlpKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
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
