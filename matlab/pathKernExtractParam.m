function [params, names] = pathKernExtractParam(kern)

% PATHKERNEXTRACTPARAM Extract parameters from the RBF kernel structure.
% FORMAT
% DESC Extract parameters from the path kernel
% structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and parameter names from the path
% kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings giving names to the parameters.
%
% SEEALSO pathKernParamInit, pathKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


[gparams,gnames]=kernExtractParam(kern.gkern);

params=[sqrt(kern.cd),kern.chv,gparams];
names=['diag step cost','hor/vert step cost',gnames];
