function [params, names] = rbfperiodicKernExtractParam(kern)

% RBFPERIODICKERNEXTRACTPARAM Extract parameters from the RBFPERIODIC kernel structure.
%
%	Description:
%
%	PARAM = RBFPERIODICKERNEXTRACTPARAM(KERN) extracts parameters from
%	the RBF derived periodic kernel structure into a vector of
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
%	[PARAM, NAMES] = RBFPERIODICKERNEXTRACTPARAM(KERN) extracts
%	parameters and parameter names from the RBF derived periodic kernel
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
%	% SEEALSO RBFPERIODICKERNPARAMINIT, RBFPERIODICKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2007 Neil D. Lawrence


params = [kern.inverseWidth kern.variance];
if nargout > 1
  names={'inverse width', 'variance'};
end