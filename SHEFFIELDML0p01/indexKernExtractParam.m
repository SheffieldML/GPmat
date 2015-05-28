function [params, names] = indexKernExtractParam(kern)

% INDEXKERNEXTRACTPARAM Extract parameters from the INDEX kernel structure.
%
%	Description:
%
%	PARAM = INDEXKERNEXTRACTPARAM(KERN) extracts parameters from the
%	index based covariance function kernel structure into a vector of
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
%	[PARAM, NAMES] = INDEXKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the index based covariance function kernel
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
%	% SEEALSO INDEXKERNPARAMINIT, INDEXKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2011 Neil D. Lawrence

  params = [kern.variance];
  names{1} = ['index variance'];
end