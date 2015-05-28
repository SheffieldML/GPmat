function [params, names] = linardKernExtractParam(kern)


% LINARDKERNEXTRACTPARAM Extract parameters from the LINARD kernel structure.
% FORMAT
% DESC Extract parameters from the automatic relevance
% determination linear kernel structure into a vector of parameters
% for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and parameter names from the automatic relevance
% determination linear kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO linardKernParamInit, linardKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% KERN


params = [kern.variance kern.inputScales];
if nargout > 1
  names = {'variance'};
  for i = 1:length(kern.inputScales)
    names{1+i} = ['input scale ' num2str(i)];
  end
end
