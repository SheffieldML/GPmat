function [params, names] = lfmKernExtractParam(kern)

% LFMKERNEXTRACTPARAM Extract parameters from the LFM kernel structure.
%
%	Description:
%
%	PARAM = LFMKERNEXTRACTPARAM(KERN) Extract parameters from the single
%	input motif kernel structure into a vector of parameters for
%	optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = LFMKERNEXTRACTPARAM(KERN) Extract parameters and
%	their names from the single input motif kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel.
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%
%	See also
%	% SEEALSO LFMKERNPARAMINIT, LFMKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2007 Neil D. Lawrence


params = [kern.mass, kern.spring, kern.damper,  kern.inverseWidth, kern.sensitivity];
if nargout > 1
  names = {'mass', 'spring', 'damper', 'inverse width', 'sensitivity'};
end