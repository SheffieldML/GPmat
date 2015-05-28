function [params, names] = rbfinfwhiteKernExtractParam(kern)

% RBFINFWHITEKERNEXTRACTPARAM Extract parameters from the RBF-WHITE kernel
% (with integration limits between minus infinity and infinity) structure.
% FORMAT
% DESC extracts parameters from the RBF-WHITE kernel structure into a vector
% of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the RBF-WHITE kernel
% structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If the
% field 'transforms' is not empty in the kernel structure, the parameters
% will be transformed before optimisation (for example positive only
% parameters could be logged before being returned).
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO rbfinfwhiteKernParamInit, rbfinfwhiteKernExpandParam, kernExtractParam,
% scg, conjgrad
%
% COPYRIGHT : David Luengo, 2009
%
% KERN


params = [kern.inverseWidth kern.variance];
if nargout > 1
  names = {'inverse width', 'variance'};
end
