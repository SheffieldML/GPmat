function [params, names] = sdlfmvKernExtractParam(kern)

% SDLFMVKERNEXTRACTPARAM Extract parameters from the SDLFMV kernel structure.
%
%	Description:
%
%	PARAM = SDLFMVKERNEXTRACTPARAM(KERN) Extract parameters from the
%	switching dynamical LFM kernel structure for the velocities into a
%	vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = SDLFMVKERNEXTRACTPARAM(KERN) Extract parameters and
%	their names from the switching dynamical LFM kernel structure for
%	the velocities.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel.
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%
%	See also
%	% SEEALSO SDLFMKERNEXTRACTPARAM, SDLFMKERNEXPANDPARAM, KERNEXTRACTPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


% This is just a wrapper function that calls the extractParam function for
% the SDLFM kernel.

if nargout > 1
    [params, names] = sdlfmKernExtractParam(kern);
else
    params = sdlfmKernExtractParam(kern);
end