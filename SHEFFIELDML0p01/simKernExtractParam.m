function [params, names] = simKernExtractParam(kern)

% SIMKERNEXTRACTPARAM Extract parameters from the SIM kernel structure.
%
%	Description:
%
%	PARAM = SIMKERNEXTRACTPARAM(KERN) Extract parameters from the single
%	input motif kernel structure into a vector of parameters for
%	optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = SIMKERNEXTRACTPARAM(KERN) Extract parameters and
%	their names from the single input motif kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO SIMKERNPARAMINIT, SIMKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006, 2009 Neil D. Lawrence

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  if isfield(kern, 'isNegativeS') && kern.isNegativeS
    params = [kern.decay kern.inverseWidth kern.sensitivity kern.initialVariance];
    if nargout > 1
      names = {'decay', 'inverse width', 'sensitivity', 'initial variance'};
    end
  else
    params = [kern.decay kern.inverseWidth kern.variance kern.initialVariance];
    if nargout > 1
      names = {'decay', 'inverse width', 'variance', 'initial variance'};
    end
  end
else
  if isfield(kern, 'isNegativeS') && kern.isNegativeS
    params = [kern.decay kern.inverseWidth kern.sensitivity];
    if nargout > 1
      names = {'decay', 'inverse width', 'sensitivity'};
    end
  else
    params = [kern.decay kern.inverseWidth kern.variance];
    if nargout > 1
      names = {'decay', 'inverse width', 'variance'};
    end
  end
end
