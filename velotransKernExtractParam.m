function [params, names] = velotransKernExtractParam(kern)

% VELOTRANSKERNEXTRACTPARAM Extract parameters from the VELOTRANS kernel structure.
% FORMAT
% DESC extracts parameters from the velocity translate
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the velocity translate
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
% SEEALSO velotransKernParamInit, velotransKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2011
%
% KERN

kern.nParams = kern.nParams - kern.inputDimension + 1;
if nargout == 2
  [params, names] = cmpndKernExtractParam(kern);
  namLength = length(names);
  for i = 1:(kern.inputDimension - 1)
    names{namLength+i} = ['Velocity ' num2str(i)];
  end
else
  params = cmpndKernExtractParam(kern);
end
params = [params kern.velocity];
