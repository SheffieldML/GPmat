function [params, names] = rbfhKernExtractParam(kern)

% RBFHKERNEXTRACTPARAM Extract parameters from the RBFH kernel structure.
% FORMAT
% DESC Extract parameters from the radial basis function heat kernel
% structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and parameter names from the radial basis
% function heat kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings giving names to the parameters.
%
% SEEALSO rbfhKernParamInit, rbfhKernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

params = [kern.inverseWidthTime kern.inverseWidthSpace];
if nargout > 1
  names={'inverse width time', 'inverse width space.'};
end
