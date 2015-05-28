function [params, names] = biasKernExtractParam(kern)


% BIASKERNEXTRACTPARAM Extract parameters from the BIAS kernel structure.
% FORMAT
% DESC extracts parameters from the bias kernel structure into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the bias kernel
% structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO biasKernParamInit, biasKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% KERN


params = kern.variance;
if nargout > 1
  names = {'variance'};
end
