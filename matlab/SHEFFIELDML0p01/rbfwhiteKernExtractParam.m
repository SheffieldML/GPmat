function [params, names] = rbfwhiteKernExtractParam(kern)

% RBFWHITEKERNEXTRACTPARAM Extract parameters from the RBF-WHITE kernel
%
%	Description:
%	structure.
%
%	PARAM = RBFWHITEKERNEXTRACTPARAM(KERN) extracts parameters from the
%	RBF-WHITE kernel structure into a vector of parameters for
%	optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = RBFWHITEKERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the RBF-WHITE kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing names for each parameter.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	scg, conjgrad
%	
%	
%
%	See also
%	% SEEALSO RBFWHITEKERNPARAMINIT, RBFWHITEKERNEXPANDPARAM, KERNEXTRACTPARAM, 


%	Copyright (c) 2009 David Luengo


params = [kern.inverseWidth kern.variance];
if nargout > 1
  names = {'inverse width', 'variance'};
end
