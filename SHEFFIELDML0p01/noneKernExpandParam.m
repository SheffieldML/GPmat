function kern = noneKernExpandParam(kern, params)

% NONEKERNEXPANDPARAM Create kernel structure from NONE kernel's parameters.
%
%	Description:
%
%	KERN = NONEKERNEXPANDPARAM(KERN, PARAM) returns a dummy kernel
%	function kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
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
%	NONEKERNPARAMINIT, NONEKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2008 Neil D. Lawrence


kern = kern;