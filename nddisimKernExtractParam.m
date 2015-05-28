function [params, names] = nddisimKernExtractParam(kern)

% NDDISIMKERNEXTRACTPARAM Extract parameters from the NDDISIM kernel structure.
% FORMAT
% DESC Extract parameters from the single input motif kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and their names from the single input
% motif kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO disimKernParamInit, disimKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2009
%
% COPYRIGHT : Jaakko Peltonen, 2011
%
% KERN


params = [kern.inverseWidth, kern.di_variance, kern.decay, kern.variance, kern.delay];

if nargout > 1
  names = {'inverse width', 'di_variance', 'decay', 'variance', 'delay'};
end


%fprintf(1, 'disimKern parameters physical values:\n');
%params
%names
