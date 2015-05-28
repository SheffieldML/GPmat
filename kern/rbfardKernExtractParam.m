function [params, names] = rbfardKernExtractParam(kern)


% RBFARDKERNEXTRACTPARAM Extract parameters from the RBFARD kernel structure.
% FORMAT
% DESC Extract parameters from the automatic relevance determination
% radial basis function kernel structure into a vector of parameters for
% optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% DESC Extract parameters and parameter names from the automatic
% relevance determination radial basis function kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containg the parameter names.
%
% SEEALSO rbfardKernParamInit, rbfardKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% KERN


params = [kern.inverseWidth kern.variance kern.inputScales];
if nargout > 1
  names = {'inverse width', 'variance'};
  for i = 1:length(kern.inputScales)
    names{2+i} = ['input scale ' num2str(i)];
  end
end
