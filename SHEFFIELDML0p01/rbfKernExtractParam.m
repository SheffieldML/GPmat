function [params, names] = rbfKernExtractParam(kern)

% RBFKERNEXTRACTPARAM Extract parameters from the RBF kernel structure.
%
%	Description:
%
%	PARAM = RBFKERNEXTRACTPARAM(KERN) Extract parameters from the radial
%	basis function kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = RBFKERNEXTRACTPARAM(KERN) Extract parameters and
%	parameter names from the radial basis function kernel structure.
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
%	See also
%	% SEEALSO RBFKERNPARAMINIT, RBFKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


params = [kern.inverseWidth kern.variance];
if nargout > 1
  names={'inverse width', 'variance'};
end
