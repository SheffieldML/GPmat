function kern = lfmvKernExpandParam(kern, params)

% LFMVKERNEXPANDPARAM Create kernel structure from LFMV kernel's parameters.
%
%	Description:
%
%	KERN = LFMVKERNEXPANDPARAM(KERN, PARAM) returns a LFMV kernel
%	structure filled with the parameters in the given vector. This is
%	used as a helper function to enable parameters to be optimised in,
%	for example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%
%	See also
%	LFMKERNPARAMINIT, LFMKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


kern = lfmKernExpandParam(kern, params);