function [params, names] = invcmpndKernExtractParam(kern)

% INVCMPNDKERNEXTRACTPARAM Extract parameters from the INVCMPND kernel structure.
%
%	Description:
%
%	PARAM = INVCMPNDKERNEXTRACTPARAM(KERN) Extract parameters from the
%	inv. precision compound kernel matrix into a vector of parameters
%	for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. The vector
%	   of 'transforms' is assumed to be empty here. Any transformations
%	   of parameters should be done in component kernels.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%
%	See also
%	% SEEALSO CMPNDKERNEXTRACTPARAM, INVCMPNDKERNPARAMINIT, INVCMPNDKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2012 Andreas C. Damianou



if nargout > 1
    [params, names] = cmpndKernExtractParam(kern);
else
    params = cmpndKernExtractParam(kern);
end