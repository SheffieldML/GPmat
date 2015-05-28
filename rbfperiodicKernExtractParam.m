function [params, names] = rbfperiodicKernExtractParam(kern)

% RBFPERIODICKERNEXTRACTPARAM Extract parameters from the RBFPERIODIC kernel structure.
% FORMAT
% DESC extracts parameters from the RBF derived periodic
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the RBF derived periodic
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
% SEEALSO rbfperiodicKernParamInit, rbfperiodicKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% KERN


params = [kern.inverseWidth kern.variance];
if nargout > 1
  names={'inverse width', 'variance'};
end
