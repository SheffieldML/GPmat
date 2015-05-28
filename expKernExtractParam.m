function [params, names] = expKernExtractParam(kern)

% EXPKERNEXTRACTPARAM Extract parameters from the EXP kernel structure.
% FORMAT
% DESC extracts parameters from the exponentiated
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the exponentiated
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
% SEEALSO expKernParamInit, expKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% KERN


if nargout > 1
  [params, names] = kernExtractParam(kern.argument);
  for i = 1:length(names)
    names{i} = ['Arg kern ' names{i}];
  end
  names = {'variance', names{:}};
else
  params = kernExtractParam(kern.argument);
end

params = [kern.variance params];
