function [params, names] = simwhiteKernExtractParam(kern)

% SIMWHITEKERNEXTRACTPARAM Extract parameters from the SIM-WHITE kernel
% structure.
% FORMAT
% DESC extracts parameters from the SIM-White (Single Input Motif - White)
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the SIM-White (Single
% Input Motif - White) kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If the
% field 'transforms' is not empty in the kernel structure, the parameters
% will be transformed before optimisation (for example positive only
% parameters could be logged before being returned).
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO simwhiteKernParamInit, simwhiteKernExpandParam, kernExtractParam,
% scg, conjgrad
%
% COPYRIGHT : David Luengo, 2009
%
% KERN


params = [kern.decay kern.variance kern.sensitivity];
if nargout > 1
  names = {'decay', 'variance', 'sensitivity'};
end
