function [params, names] = tensorKernExtractParam(kern)


% TENSORKERNEXTRACTPARAM Extract parameters from the TENSOR kernel structure.
% FORMAT
% DESC Extract parameters from the tensor product kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and parameter names from the tensor
% product kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO tensorKernParamInit, tensorKernExpandParam, kernExtractParam, scg, conjgrad, cmpndKernExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% KERN

if nargout > 1
  [params, names] = cmpndKernExtractParam(kern);
else
  params = cmpndKernExtractParam(kern);
end
