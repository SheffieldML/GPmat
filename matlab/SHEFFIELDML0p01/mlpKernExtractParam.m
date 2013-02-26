function [params, names] = mlpKernExtractParam(kern)

% MLPKERNEXTRACTPARAM Extract parameters from the MLP kernel structure.
%
%	Description:
%
%	PARAM = MLPKERNEXTRACTPARAM(KERN) Extract parameters from the
%	multi-layer perceptron kernel matrix into a vector of parameters for
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
%	[PARAM, NAMES] = MLPKERNEXTRACTPARAM(KERN) Extract parameters and
%	parameter names from the multi-layer perceptron kernel matrix into a
%	vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings giving names to the parameters.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO MLPKERNPARAMINIT, MLPKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


params = [kern.weightVariance kern.biasVariance kern.variance];
if nargout > 1
  names = {'weight variance', 'bias variance', 'variance'};
end
