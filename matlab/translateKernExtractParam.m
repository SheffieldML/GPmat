function [params, names] = translateKernExtractParam(kern)

% TRANSLATEKERNEXTRACTPARAM Extract parameters from the TRANSLATE kernel structure.
% FORMAT
% DESC extracts parameters from the input space translation
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the input space translation
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
% SEEALSO translateKernParamInit, translateKernExpandParam,
% kernExtractParam, cmpndKernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% KERN

kern.nParams = kern.nParams - kern.inputDimension;
if nargout == 2
  [params, names] = cmpndKernExtractParam(kern);
  namLength = length(names);
  for i = 1:kern.inputDimension
    names{namLength+i} = ['Centre ' num2str(i)];
  end
else
  params = cmpndKernExtractParam(kern);
end
params = [params kern.centre];
