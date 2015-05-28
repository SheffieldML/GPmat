function [params, names] = ouKernExtractParam(kern)

% OUKERNEXTRACTPARAM Extract parameters from the OU kernel structure (see
% ouKernCompute or ouKernParamInit for a more detailed description of the
% OU kernel).
% FORMAT
% DESC extracts parameters from the Ornstein-Uhlenbeck kernel
% structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the Ornstein-Uhlenbeck
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
% SEEALSO ouKernParamInit, ouKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : David Luengo, 2009
%
% KERN


params = [kern.decay kern.variance];
if nargout > 1
  names={'decay', 'variance'};
end
