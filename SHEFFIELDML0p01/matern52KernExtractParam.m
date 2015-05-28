function [params, names] = matern52KernExtractParam(kern)

% MATERN52KERNEXTRACTPARAM Extract parameters from the MATERN52 kernel structure.
%
%	Description:
%
%	PARAM = MATERN52KERNEXTRACTPARAM(KERN) extracts parameters from the
%	matern kernel with nu=5/2 kernel structure into a vector of
%	parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = MATERN52KERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the matern kernel with nu=5/2 kernel
%	structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing names for each parameter.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO MATERN52KERNPARAMINIT, MATERN52KERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006 Neil D. Lawrence


params = [kern.lengthScale kern.variance];
if nargout > 1
  names={'length scale', 'variance'};
end