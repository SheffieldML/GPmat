function params = ardKernExtractParam(kern)


% ARDKERNEXTRACTPARAM Extract parameters from the ARD kernel structure.
% FORMAT
% DESC Extract parameters from the pre-built RBF and linear ARD kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO ardKernParamInit, ardKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004
%
% KERN


params = [kern.inverseWidth kern.rbfVariance ...
          kern.biasVariance kern.whiteVariance ...
          kern.linearVariance kern.inputScales];
