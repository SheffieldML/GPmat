function [params, names] = ratquadKernExtractParam(kern)

% RATQUADKERNEXTRACTPARAM Extract parameters from the RATQUAD kernel structure.
% FORMAT
% DESC extracts parameters from the rational quadratic
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the rational quadratic
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
% SEEALSO ratquadKernParamInit, ratquadKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% KERN

params = [kern.alpha kern.lengthScale kern.variance];
if nargout > 1
  names={'alpha', 'length scale', 'variance'};
end
