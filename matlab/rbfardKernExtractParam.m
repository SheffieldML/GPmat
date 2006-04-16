function params = rbfardKernExtractParam(kern)


% RBFARDKERNEXTRACTPARAM Extract parameters from the RBFARD kernel structure.
% FORMAT
% DESC Extract parameters from the automatic relevance determination radial basis function kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO rbfardKernParamInit, rbfardKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% KERN


params = [kern.inverseWidth kern.variance kern.inputScales];
