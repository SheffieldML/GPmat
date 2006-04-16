function params = tensorKernExtractParam(kern)


% TENSORKERNEXTRACTPARAM Extract parameters from the TENSOR kernel structure.
% FORMAT
% DESC Extract parameters from the tensor product kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO tensorKernParamInit, tensorKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% KERN


params = cmpndKernExtractParam(kern);
