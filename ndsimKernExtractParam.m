function [params, names] = ndsimKernExtractParam(kern)

% NDSIMKERNEXTRACTPARAM Extract parameters from the SIM kernel structure.
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
% SEEALSO simKernParamInit, simKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009
%
% COPYRIGHT : Jaakko Peltonen, 2011
%
% KERN

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  if isfield(kern, 'isNegativeS') && kern.isNegativeS
    params = [kern.inverseWidth kern.sensitivity kern.initialVariance];
    if nargout > 1
      names = {'inverse width', 'sensitivity', 'initial variance'};
    end
  else
    params = [kern.inverseWidth kern.variance kern.initialVariance];
    if nargout > 1
      names = {'inverse width', 'variance', 'initial variance'};
    end
  end
else
  if isfield(kern, 'isNegativeS') && kern.isNegativeS
    params = [kern.inverseWidth kern.sensitivity];
    if nargout > 1
      names = {'inverse width', 'sensitivity'};
    end
  else
    params = [kern.inverseWidth kern.variance];
    if nargout > 1
      names = {'inverse width', 'variance'};
    end
  end
end

%fprintf(1, 'simKern parameters physical values:\n');
%params
%names
