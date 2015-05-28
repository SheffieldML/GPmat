function [params, transform] = sqexpKernExtractParam(kern)


% SQEXPKERNEXTRACTPARAM Extract parameters from the SQEXP kernel structure.
% FORMAT
% DESC Extract parameters from the pre-built compound squared
% exponential kernel structure into a vector of parameters for
% optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% DESC Extract parameters and parameter names from the pre-built
% compound squared exponential kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : celly array of strings containing parameter names.
%
% SEEALSO sqexpKernParamInit, sqexpKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004
%
% KERN


params = [kern.inverseWidth kern.rbfVariance kern.biasVariance ...
          kern.whiteVariance];
if nargout > 1
  names{1} = 'sqexp inverse width';
  names{2} = 'sqexp rbf variance';
  names{3} = 'sqexp bias variance';
  names{4} = 'sqexp white variance';
end
