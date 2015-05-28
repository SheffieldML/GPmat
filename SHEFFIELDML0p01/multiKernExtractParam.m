function [params, names] = multiKernExtractParam(kern)

% MULTIKERNEXTRACTPARAM Extract parameters from the MULTI kernel structure.
%
%	Description:
%
%	PARAM = MULTIKERNEXTRACTPARAM(KERN) Extract parameters from the
%	multiple output block kernel structure into a vector of parameters
%	for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. The vector
%	   of transforms is assumed to be empty here, any transormation of
%	   parameters is assumed to be done in the component kernels.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = MULTIKERNEXTRACTPARAM(KERN) Extract parameters and
%	parameter names from the multiple output block kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. The vector
%	   of transforms is assumed to be empty here, any transormation of
%	   parameters is assumed to be done in the component kernels.
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%
%	See also
%	% SEEALSO MULTIKERNPARAMINIT, MULTIKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006 Neil D. Lawrence


if nargout > 1
  [params, names] = cmpndKernExtractParam(kern);
else
  params = cmpndKernExtractParam(kern);
end