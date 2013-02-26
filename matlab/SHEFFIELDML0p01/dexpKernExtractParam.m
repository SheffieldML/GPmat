function [params, names] = dexpKernExtractParam(kern)

% DEXPKERNEXTRACTPARAM Extract parameters from the double exponential's
%
%	Description:
%	kernel structure.
%	
%
%	PARAMS = DEXPKERNEXTRACTPARAM(KERN) extracts parameters from the
%	double exponential kernel's structure into a vector of parameters
%	for optimisation.
%	 Returns:
%	  PARAMS - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = DEXPKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the double exponential kernel structure.
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
%	See also
%	% SEEALSO DEXPKERNPARAMINIT, DEXPKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2009 David Luengo



params = [kern.decay kern.variance];
if nargout > 1
  names={'decay', 'variance'};
end
