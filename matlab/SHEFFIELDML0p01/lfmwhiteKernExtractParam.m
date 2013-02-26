function [params, names] = lfmwhiteKernExtractParam(kern)

% LFMWHITEKERNEXTRACTPARAM Extract parameters from the LFM-WHITE kernel
%
%	Description:
%	structure.
%
%	PARAM = LFMWHITEKERNEXTRACTPARAM(KERN) extracts parameters from the
%	LFM-White (Latent Force Model - White) kernel structure into a
%	vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = LFMWHITEKERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the LFM-White (Latent Force Model - White)
%	kernel structure.
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
%	% SEEALSO LFMWHITEKERNPARAMINIT, LFMWHITEKERNEXPANDPARAM, KERNEXTRACTPARAM, 


%	Copyright (c) 2009 David Luengo


params = [kern.mass kern.spring kern.damper kern.variance kern.sensitivity];
if nargout > 1
  names = {'mass', 'spring', 'damper', 'variance', 'sensitivity'};
end
