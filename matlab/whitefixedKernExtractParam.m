function [params, names] = whitefixedKernExtractParam(kern)


% WHITEFIXEDKERNEXTRACTPARAM Extract parameters from the WHITEFIXED kernel structure.
% FORMAT
% DESC Extract parameters from the fixed parameter white noise
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and parameter names from the fixed
% parameter white noise kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO whitefixedKernParamInit, whitefixedKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Nathaniel J. King, 2006
%
% KERN

params = [];
names = {};
