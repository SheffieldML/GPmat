function [params, names] = polyardKernExtractParam(kern)

% POLYARDKERNEXTRACTPARAM Extract parameters from the POLYARD kernel structure.
% FORMAT
% DESC Extract parameters from the automatic relevance
% determination polynomial kernel structure into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% DESC extracts parameters and parameter names from the automatic
% relevance determination polynomial kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing the parameter names.
% SEEALSO polyardKernParamInit, polyardKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% KERN

params = [kern.weightVariance kern.biasVariance kern.variance kern.inputScales];
if nargout > 1
  names = {'weight variance', 'bias variance', 'variance'};
  for i = 1:length(kern.inputScales)
    names{3+i} = ['input scale ' num2str(i)];
  end
end
