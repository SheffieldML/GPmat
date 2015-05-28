function [params, names] = rbfperiodic2KernExtractParam(kern)

% RBFPERIODIC2KERNEXTRACTPARAM Extract parameters from the RBFPERIODIC2 kernel structure.
% FORMAT
% DESC extracts parameters from the RBF periodic covariance with variying period
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the RBF periodic covariance with variying period
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
% SEEALSO rbfperiodic2KernParamInit, rbfperiodic2KernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009
%
%
% MODIFICATIONS : Andreas C. Damianou, 2011
%
% MODIFICATIONS : Michalis K. Titsias, 2011

% KERN

params = [kern.inverseWidth kern.variance kern.factor];
if nargout > 1
  names={'inverse width', 'variance', 'factor'};
end
