function [params, names] = whitehKernExtractParam(kern)

% WHITEHKERNEXTRACTPARAM Extract parameters from the WHITEH kernel structure.
%
%	Description:
%
%	PARAM = WHITEHKERNEXTRACTPARAM(KERN) Extract parameters from the
%	whiteh noise kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = WHITEHKERNEXTRACTPARAM(KERN) Extract parameters and
%	parameter names from the whiteh noise kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings giving paramter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO WHITEHKERNPARAMINIT, WHITEHKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence

params = kern.variance;
if nargout > 1
  names = {'variance'};
end

