function [params, names] = noneKernExtractParam(kern)

% NONEKERNEXTRACTPARAM Extract parameters from the NONE kernel structure.
%
%	Description:
%
%	PARAM = NONEKERNEXTRACTPARAM(KERN) extracts parameters from the
%	dummy kernel function kernel structure into a vector of parameters
%	for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = NONEKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the dummy kernel function kernel structure.
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
%	% SEEALSO NONEKERNPARAMINIT, NONEKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2008 Neil D. Lawrence

params = [];
names = {};
