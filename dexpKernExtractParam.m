function [params, names] = dexpKernExtractParam(kern)

% DEXPKERNEXTRACTPARAM Extract parameters from the double exponential's
% kernel structure.
%
% FORMAT
% DESC extracts parameters from the double exponential kernel's structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN params : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the double exponential
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
% SEEALSO dexpKernParamInit, dexpKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : David Luengo, 2009

% KERN


params = [kern.decay kern.variance];
if nargout > 1
  names={'decay', 'variance'};
end
