function [params, names] = polyKernExtractParam(kern)

% POLYKERNEXTRACTPARAM Extract parameters from the POLY kernel structure.
%
%	Description:
%
%	PARAM = POLYKERNEXTRACTPARAM(KERN) Extract parameters from the
%	polynomial kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = POLYKERNEXTRACTPARAM(KERN) Extract parameters and
%	their names from the polynomial kernel structure.
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
%	% SEEALSO POLYKERNPARAMINIT, POLYKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2005, 2006 Neil D. Lawrence


params = [kern.weightVariance kern.biasVariance kern.variance];
if nargout > 1
  names = {'weight variance', 'bias variance', 'variance'};
end
